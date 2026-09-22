"""
graph_motor.py
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from subfunctions import tau_dcmotor, motor


# Generate motor shaft speeds from 0 to no-load speed
omega = np.linspace(0, motor['speed_noload'], 100)


# Calculate motor shaft torque
tau = tau_dcmotor(omega, motor)


# Calculate motor power
power = tau * omega

# display settings
plt.figure(figsize=(8, 10))

# top plot: motor shaft speed [rad/s] vs. motor shaft torque [Nm] (use torque on the x-axis)

plt.subplot(3,1,1)
plt.plot(tau, omega)
plt.xlabel('Motor Shaft Torque [Nm]')
plt.ylabel('Motor Shaft Speed [rad/s]')

# middle plot: motor power [W] vs. motor shaft torque [Nm] (use torque on the x-axis)

plt.subplot(3,1,2)
plt.plot(tau, power)
plt.xlabel('Motor Shaft Torque [Nm]')
plt.ylabel('Motor Power [W]')

# bottom plot: motor power [W] vs. motor shaft speed [rad/s] (use speed on the x-axis)

plt.subplot(3,1,3)
plt.plot(omega, power)
plt.xlabel('Motor Shaft Speed [rad/s]')
plt.ylabel('Motor Power [W]')
