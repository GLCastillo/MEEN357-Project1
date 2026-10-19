"""
efficiency_visualization.py
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from subfunctions import rover

# Get the motor efficiency and torque data
effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
effcy = rover['wheel_assembly']['motor']['effcy']

# Create the cubic interpolation function
effcy_fun = interp1d(effcy_tau, effcy, kind='cubic')

# Generate 100 evenly spaced torque values
tau_values = np.linspace(min(effcy_tau), max(effcy_tau), 100)

# Evaluate efficiency at the 100 torque alues
effcy_values = effcy_fun(tau_values)

#display settings
plt.figure(figsize=(9,6), dpi=300)

# Plot original data points as stars
plt.plot(effcy_tau, effcy, '*', label='known data points')

# Plot interpolated efficiencies as a line
plt.plot(tau_values, effcy_values, '-', label='cubic spline interpolation')

#label axes
plt.xlabel('Motor torque [N*m]')
plt.ylabel('Motor efficiency')
plt.title('Motor efficiency vs Torque (N*m)')

# Add legend
plt.legend()
plt.grid(True)
plt.show()
