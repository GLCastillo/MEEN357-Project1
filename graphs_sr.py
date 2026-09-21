"""
graphs_sr.py
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from subfunctions import tau_dcmotor, get_gear_ratio, motor, rover


# Generate input motor shaft speeds from 0 to no-load speed
omega_in = np.linspace(0, motor['speed_noload'], 100)


# Calculate input motor shaft torque
tau_in = tau_dcmotor(omega_in, motor)


# Get gear ratio
Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
# Calculate speed reducer values
omega_out = omega_in / Ng

tau_out = Ng * tau_in

power_out = tau_out * omega_out

# display settings
plt.figure(figsize=(16, 20))

# top plot: motor shaft speed [rad/s] vs. motor shaft torque [Nm] (use torque on the x-axis)

plt.subplot(3,1,1)
plt.plot(tau_out, omega_out)
plt.xlabel('Speed Reducer Output Torque [Nm]')
plt.ylabel('Speed Reducer Output Speed [rad/s]')

# middle plot: motor power [W] vs. motor shaft torque [Nm] (use torque on the x-axis)

plt.subplot(3,1,2)
plt.plot(tau_out, power_out)
plt.xlabel('Speed Reducer Output Torque [Nm]')
plt.ylabel('Speed Reducer Output Power [W]')

# bottom plot: motor power [W] vs. motor shaft speed [rad/s] (use speed on the x-axis)

plt.subplot(3,1,3)
plt.plot(omega_out, power_out)
plt.xlabel('Speed Reducer Output Speed [rad/s]')
plt.ylabel('Speed Reducer Output Power [W]')
