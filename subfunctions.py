"""subfunctions.py"""

# Imports
import numpy as np
import math
from define_experiment import experiment1
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp, trapezoid

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
    'torque_stall': 170.0,
    'torque_noload': 0.0,
    'speed_noload': 3.80,
    'mass': 5.0,
    'effcy_tau': np.array([0, 10, 20, 40, 70, 165]),
    'effcy': np.array([0, 0.55, 0.75, 0.71, 0.50, 0.05])
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

#####

# Computes the rotational speed of the motor shaft
def motorW(v, rover):
    
    v_type_check = False
    
    # Validate that v is a scalar or vector, i'm missing one check
    if v_type_check == False:
        raise Exception("omega must be a scalar or vector value")
    
    # Validate that w is a dictionary
    if not isinstance(rover, dict):
        raise Exception ("rover must be a dictionary.")
        
        
    # rotational speed of motor = (rover velocity/wheel radius)*gear ratio
    
    r = rover["wheel_assembly"]["radius"]
    G = get_gear_ratio(speed_reducer)
    
    # calculate the motor angular velocity
    w = (v/r)*G
    
    return w

#####

# rover_dynamics to compute the derivative of the state vector
def rover_dynamics(t, y, rover, planet, experiment):
    
    # Validate that t is a scalar
    if not np.isscalar(t):
        raise Exception("t must be a scalar")
    
    # Validate that y is a 1D numpy array
    if y > 1:
        raise Exception('Y must be a scalar or a 1D array.')
    
    # Validate that rover is a dictionary
    if not isinstance(rover, dict):
        raise Exception ("Rover must be a dictionary.")
    
    # Validate that planet is a dictionary
    if not isinstance(planet, dict):
        raise Exception ("Planet must be a dictionary.")
        
    # Validate that experiment is a dictionary
    if not isinstance(experiment, dict):
        raise Exception ("Experiment must be a dictionary.")
    
 
    # Caclualate terrain angle after creating interpolation function        
    alpha_dist = experiment['alpha_dist']
    alpha_deg = experiment['alpha_deg']
        
    alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate')
    terrain_angle = alpha_fun(y[1])
    
    # Rover velocity
    velocity = y[0]
    
    # Calculate motor angular velocity
    omega = motorW(velocity, rover)
    
    # Calculate net force
    Fnet = F_net(omega, terrain_angle, rover, planet, experiment['Crr'])
    
    # Calculate total rover mass
    mass = get_mass(rover)
    
    # Caclulate acceleration
    acceleration = Fnet/mass
    
    #Return derivetives
    dydt = np.array([acceleration, velocity])
    
    return dydt

#####

# Computes the instantaneous mechanical power output
def mechpower(v, rover):
    
    # Validate that v is a scalar or a 1D numpy array
    if v > 1:
        raise Exception('v must be a scalar or a 1D array.')
    
    # Validate that v contains numerical values
    if not np.issubdtype(np.asarray(v).dtype, np.number):
        raise Exception("v must contain numerical values")
    
    # Validate that rover is a dictionary
    if not isinstance(rover, dict):
        raise Exception ("Rover must be a dictionary.")
    
    # Calculate motor angular velocity
    omega = motorW(v, rover)
    
    # Calculate the motor torque
    tau = tau_dcmotor(omega, rover)
    
    # Calculate the mechanical power
    Power = tau * omega
    
    return Power

#####

# Computes the total electrical energy consumed by all six rover motors during a simulation
def battenergy(t, v, rover):

    if not isinstance(t, np.ndarray) or t.ndim != 1:
        raise Exception("t must be a 1D numpy array.")

    if not isinstance(v, np.ndarray) or v.ndim != 1:
        raise Exception("v must be a 1D numpy array.")

    if len(t) != len(v) or len(t) < 2:
        raise Exception("t and v must have equal lengths of at least 2.")

    if not np.issubdtype(t.dtype, np.number) or \
       not np.issubdtype(v.dtype, np.number):
        raise Exception("t and v must contain numerical values.")

    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(v)):
        raise Exception("t and v must contain finite values.")

    if np.any(np.diff(t) <= 0):
        raise Exception("Time must be strictly increasing.")

    if not isinstance(rover, dict):
        raise Exception("rover must be a dictionary.")

    # Get motor dictionary
    motor = rover['wheel_assembly']['motor']

    # Calculate motor rotational speed
    omega = motorW(v, rover)

    # Calculate motor torque
    tau = tau_dcmotor(omega, motor)

    # Calculate mechanical power of one motor
    P_motor = mechpower(v, rover)

    # Get motor efficiency data
    effcy_tau = motor['effcy_tau']
    effcy = motor['effcy']

    # Interpolate motor efficiency
    effcy_fun = interp1d(effcy_tau, effcy, kind='cubic')

    # Keep torque within measured efficiency range
    tau_eval = np.clip(tau, effcy_tau[0], effcy_tau[-1])

    # Calculate efficiency at each torque
    eta = effcy_fun(tau_eval)

    # Calculate total battery power for six motors
    P_batt = np.zeros_like(P_motor, dtype=float)

    active = np.abs(P_motor) > 1e-12

    if np.any(eta[active] <= 0):
        raise Exception("Motor efficiency must be positive.")

    P_batt[active] = 6 * P_motor[active] / eta[active]

    # Integrate battery power over time
    E = trapezoid(P_batt, t)

    return float(E)

####

# Simulates the rover's trajectory and calculates its telemetry data
def simulate_rover(rover, planet, experiment, end_event):

    # Validate dictionaries
    if not isinstance(rover, dict):
        raise Exception("rover must be a dictionary.")

    if not isinstance(planet, dict):
        raise Exception("planet must be a dictionary.")

    if not isinstance(experiment, dict):
        raise Exception("experiment must be a dictionary.")

    if not isinstance(end_event, dict):
        raise Exception("end_event must be a dictionary.")

    # Get initial conditions and time range
    t_range = experiment['time_range']
    y0 = experiment['initial_conditions']

    # Create events to terminate the simulation
    events = end_of_mission_event(end_event)

    # Solve the rover differential equations
    sol = solve_ivp(
        lambda t, y: rover_dynamics(t, y, rover, planet, experiment),
        t_range,
        y0,
        events=events
    )

    if not sol.success:
        raise Exception("Simulation failed: " + sol.message)

    # Get time, velocity, and position
    t = sol.t
    v = sol.y[0]
    x = sol.y[1]

    # Calculate mechanical power
    P = mechpower(v, rover)

    # Calculate battery energy
    E = battenergy(t, v, rover)

    # Calculate total distance traveled
    distance = x[-1] - x[0]

    # Calculate completion time
    completion_time = t[-1] - t[0]

    # telemetry dictionary
    rover['telemetry'] = {
        'Time': t,
        'completion_time': completion_time,
        'velocity': v,
        'position': x,
        'distance_traveled': distance,
        'max_velocity': np.max(v),
        'average_velocity': distance / completion_time,
        'power': P,
        'battery_energy': E,
        'energy_per_distance': E / distance
    }

    return rover

####

# given from Canvas
def end_of_mission_event(end_event):
    """
    Defines an event that terminates the mission simulation. Mission is over
    when rover reaches a certain distance, has moved for a maximum simulation 
    time or has reached a minimum velocity.            
    """
    
    mission_distance = end_event['max_distance']
    mission_max_time = end_event['max_time']
    mission_min_velocity = end_event['min_velocity']
    
    # Assume that y[1] is the distance traveled
    distance_left = lambda t,y: mission_distance - y[1]
    distance_left.terminal = True
    
    time_left = lambda t,y: mission_max_time - t
    time_left.terminal = True
    
    velocity_threshold = lambda t,y: y[0] - mission_min_velocity;
    velocity_threshold.terminal = True
    velocity_threshold.direction = -1
    
    # terminal indicates whether any of the conditions can lead to the
    # termination of the ODE solver. In this case all conditions can terminate
    # the simulation independently.
    
    # direction indicates whether the direction along which the different
    # conditions is reached matters or does not matter. In this case, only
    # the direction in which the velocity treshold is arrived at matters
    # (negative)
    
    events = [distance_left, time_left, velocity_threshold]
    
    return events



