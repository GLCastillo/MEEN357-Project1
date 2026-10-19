import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from define_rover import define_rover

# Get the rover dictionary
rover = define_rover()

# Get the motor efficiency and torque data
effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
effcy = rover['wheel_assembly']['motor']['effcy']

# Create the cubic interpolation function
effcy_fun = interp1d(effcy_tau, effcy, kind='cubic')

# Generate 100 evenly spaced torque values
tau_values = np.linspace(min(effcy_tau), max(effcy_tau), 100)
