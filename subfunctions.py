"""subfunctions.py"""

# Imports
import numpy as np
import math

# Define constants and dictionaries (given from appendix A and C)
wheel = {
    'radius': 0.30, # in m
    'mass': 1.0 # in kg
    }

speed_reducer = {
    'type': "reverted",
    'diam_pinion': 0.04, # in m
    'diam_gear': 0.07, # in m
    'mass': 1.50 # in kg
    }

motor = {
    'torque_stall': 170.0, # in Nm
    'torque_noload': 0.0, #in Nm
    'speed_noload': 3.80, # in rad/s
    'mass': 5.0 # in kg
        }

chassis = {
    'mass': 659.0 # in kg
   }

science_payload = {
    'mass': 75.0 # in kg
    }

power_subsys = {
    'mass': 90.0 # in kg
    }

planet = {
    'g': 3.72 # in m/s^2
    }

wheel_assembly = {
    'wheel': wheel,
    'speed_reducer': speed_reducer,
    'motor': motor
    }

rover = {
    'wheel_assembly': wheel_assembly,
    'chassis': chassis,
    'science_payload': science_payload,
    'power_subsys': power_subsys
    }

#####

def get_mass(rover):

    if not isinstance (rover, dict):
        raise Exception ("Input must be a dictionary.")


    chassis_mass = rover["chassis"]["mass"]
    power_mass = rover["power_subsys"]["mass"]
    science_mass = rover["science_payload"]["mass"]
    wheel_mass = rover["wheel_assembly"]["wheel"]["mass"]
    speed_mass = rover["wheel_assembly"]["speed_reducer"]["mass"]
    motor_mass = rover["wheel_assembly"]["motor"]["mass"]

    total_mass = chassis_mass + power_mass + science_mass + 6*(wheel_mass + speed_mass + motor_mass)

    return total_mass

#####

def get_gear_ratio(speed_reducer):
    """
    Returns the speed reduction ratio for the speed
    reducer based on speed_reducer dict.
    """
    if isinstance(speed_reducer, dict) == False:
        raise Exception("speed_reducer must be a dictionary.")

    if speed_reducer['type'].lower() != "reverted": # include .lower() so that input letter case doesn't matter
        raise Exception("Only reverted type speed reducers are used for this phase of the project.")

    # gear ration Ng = ( d2 / d1 )^2
    Ng = ( speed_reducer['diam_gear'] / speed_reducer['diam_pinion'] ) **2

    return Ng

#####

def tau_dcmotor(omega, motor):
    
    # omega must be scalar or 1D array
    if not (np.isscalar(omega) or (isinstance(omega, np.ndarray) and omega.ndim == 1)):
        raise Exception("Omega must be a scalar or vector")
        
    if not isinstance(motor, dict):
        raise Exception("Motor must be a dictionary.")

    torque_stall = motor["torque_stall"]
    torque_noload = motor["torque_noload"]
    speed_noload = motor["speed_noload"]

    if np.isscalar(omega): # evaluating tau if omega is scalar
        if omega <0:
            torque = torque_stall

        elif omega > speed_noload:
            torque = 0

        else:
            torque = torque_stall - ((torque_stall - torque_noload)/speed_noload)*omega

    else: # evaluating tau if omega is vector
        torque = np.where( omega < 0, torque_stall, np.where(omega > speed_noload, 0, torque_stall - ((torque_stall - torque_noload)/speed_noload)*omega ) )

    return torque

#####

def F_drive(omega, rover):
    """
    Returns the force applied to the rover by the drive system given
    information about the drive system (wheel_assembly) and the motor shaft speed.
    """
    # True if input variable omega is a 1d numpy array or scalar value
    omega_type_check = False
    if isinstance(omega, np.ndarray) and omega.ndim == 1:
          omega_type_check = True
    if np.isscalar(omega):
        omega_type_check = True

    if omega_type_check == False:
        raise Exception("Omega must be a scalar or vector value")

    if isinstance(rover, dict) == False:
        raise Exception("Rover must be a dictionary.")

    motor = rover['wheel_assembly']['motor']
    sr = rover['wheel_assembly']['speed_reducer']
    wheel = rover['wheel_assembly']['wheel']

    tau = tau_dcmotor(omega, motor) # call to tau_dcmotor

    Ng = get_gear_ratio(sr)

    # Base formula: Fd = torque / wheel radius
    # Ng * tau is output torque accounting for gear ratio
    # multiply by 6 to account for all 6 wheels
    Fd = 6 * Ng * tau / wheel['radius']

    return Fd

#####

def F_gravity(terrain_angle, rover, planet):

    if not (np.isscalar(terrain_angle) or isinstance(terrain_angle, np.ndarray)):
        raise Exception ("Terrain angle(s) must be a scalar or vector.")
        
    if np.any(terrain_angle < -75.0) or np.any(terrain_angle > 75.0):
        raise Exception ("Terrain angle(s) must be between -75 and +75 degrees.")
        
    if not isinstance(rover, dict):
        raise Exception ("Rover must be a dictionary.")
        
    if not isinstance(planet, dict):
        raise Exception ("Planet must be a dictionary.")

    g_mars = planet["g"]
    mass_rover = get_mass(rover)
    alpha = terrain_angle * math.pi / 180 
    
    Fgt = - mass_rover * g_mars * np.sin(alpha) # negative sign to account for rover rolling DOWN the incline

    return Fgt

#####

def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Returns the magnitude of the force acting on the rover in the direction
    of its translational motion due to rolling resistances given the
    terrain inclination angle, rover properties, and a rolling resistance coefficient.
    """

    # .asarray() turns variable into type array if not already
    omega = np.asarray(omega)
    terrain_angle = np.asarray(terrain_angle)

    if omega.ndim > 1 or terrain_angle.ndim > 1:
        raise Exception('Omega and terrain_angle must be scalars or 1D arrays.')

    if omega.size != terrain_angle.size:
        raise Exception('Omega and terrain_angle must be the same size.')

    if np.any(terrain_angle < -75.0) or np.any(terrain_angle > 75.0):
        raise Exception ("Terrain angle(s) must be between -75 and +75 degrees.")

    if isinstance(rover, dict) == False:
        raise Exception("Rover must be a dictionary.")

    if isinstance(planet, dict) == False:
        raise Exception("Planet must be a dictionary.")

    if np.isscalar(Crr) == False or Crr <= 0:
        raise Exception("Crr must be a positive scalar.")

    # Calculate normal force
    alpha = terrain_angle * math.pi / 180 # convert degrees to radians
    m = get_mass(rover)
    Fn = m * planet['g'] * np.cos(alpha)


    # Calculate constant rolling resistance coefficient
    simple_Frr = Crr * Fn

    # Take rover velocity into account
    Ng = get_gear_ratio( rover['wheel_assembly']['speed_reducer'] )
    r = rover['wheel_assembly']['wheel']['radius']
    v_rover = r * omega / Ng # in m/s


    # Calculate erf for scalar or vector
    if v_rover.ndim == 0: # means that v_rover is a scalar value
       erf_value = math.erf(40 * float(v_rover))

    else: # means that v_rover is a vector
       erf_value = np.zeros(v_rover.size) # creates array of zeros that is same size as v_rover

       for i in range(v_rover.size):
           erf_value[i] = math.erf(40 * v_rover[i])

    Frr = -erf_value * simple_Frr # negative sign to account for retarding
    # already accounting for 6 wheels since get_mass(rover) accounts for all 6 wheels
    
    return Frr

#####

def F_net(omega, terrain_angle, rover, planet, Crr):

    omega_type_check = False
    if isinstance(omega, np.ndarray):
          omega_type_check = True
          
    if np.isscalar(omega):
        omega_type_check = True

    if omega_type_check == False:
        raise Exception("omega must be a scalar or vector value")

    angle_type_check = False
    if isinstance(terrain_angle, np.ndarray):
        angle_type_check = True
    if np.isscalar(terrain_angle):
        angle_type_check = True

    if angle_type_check == False:
        raise Exception("Array of terrain angles must be a scalar or vector")


    # .asarray() turns variable into type array if not already
    # turning omega and terrain_angle into arrays so we can compare size
    omega = np.asarray(omega)
    terrain_angle = np.asarray(terrain_angle)

    if omega.ndim > 1 or terrain_angle.ndim > 1:
        raise Exception('omega and terrain_angle must be scalars or 1D arrays')

    if omega.size != terrain_angle.size:
        raise Exception('omega and terrain_angle must be the same size')

    if np.any(terrain_angle < -75.0) or np.any(terrain_angle > 75.0):
        raise Exception ("Terrain angle must be between -75 and +75 degrees")
        
    if not isinstance(rover, dict):
        raise Exception ("Rover must be a dictionary.")
        
    if not isinstance(planet, dict):
        raise Exception ("Llanet must be a dictionary.")

    Fd = F_drive(omega, rover) 
    Fg = F_gravity(terrain_angle, rover, planet)
    Fr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fnet = Fd + Fg + Fr

    return Fnet