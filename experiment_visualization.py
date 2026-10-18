### To plot motor torque vs. efficiency

import numpy as np
import matplotlib.pyplot as plt
from define_experiment import experiment1
from scipy.interpolate import interp1d

experiment, end_event = experiment1()
alpha_dist = experiment['alpha_dist']
alpha_deg = experiment['alpha_deg']

# generating 100 evenly space distances from 0 to 1000 m
values = np.linspace(min(alpha_dist), max(alpha_dist), 100)

# interpolation function
alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate')

alpha_values = alpha_fun(values)

plt.plot(alpha_dist, alpha_deg, label = 'original data')
plt.plot(values, alpha_values, label ='cubic interpolation')
plt.xlabel("Distance along path [m]")
plt.ylabel("Terrain angle [degrees]")
plt.show()
