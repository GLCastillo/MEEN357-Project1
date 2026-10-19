"""
rover_experiment1.py
"""

# Imports
import matplotlib.pyplot as plt
from subfunctions import rover, planet, simulate_rover
from define_experiment import experiment1

# Load experiment and end_event dictionaries
experiment, end_event = experiment1()

# Stopping conditions
end_event['max_distance'] = 1000
end_event['max_time'] = 10000
end_event['min_velocity'] = 0.01

# Simulate rover trajectory
rover = simulate_rover(rover, planet, experiment, end_event)

# Get telemetry data
telemetry = rover['telemetry']

t = telemetry['Time']
x = telemetry['position']
v = telemetry['velocity']
P = telemetry['power']

# display settings
plt.figure(figsize=(8, 10))

# top plot: rover position vs. time
plt.subplot(3,1,1)
plt.plot(t, x)
plt.xlabel('Time [s]')
plt.ylabel('Rover Position [m]')
plt.title('Rover Position vs. Time')

# middle plot: Rover velocity vs. time
plt.subplot(3,1,2)
plt.plot(t, v)
plt.xlabel('Time [s]')
plt.ylabel('Rover Velocity [m/s]')
plt.title('Rover Velocity vs. Time')

# bottom plot: Motor power vs. time
plt.subplot(3,1,3)
plt.plot(t, P)
plt.xlabel('Time [s]')
plt.ylabel('Motor Power [W]')
plt.title('Motor Power vs. Time')

plt.show()

# Print telemetry results
print('Simulation Results:')
print('Completion Time [s]:', telemetry['completion_time'])
print('Distance Traveled [m]:', telemetry['distance_traveled'])
print('Maximum Velocity [m/s]:', telemetry['max_velocity'])
print('Average Velocity [m/s]:', telemetry['average_velocity'])
print('Battery Energy [J]:', telemetry['battery_energy'])
print('Energy per Distance [J/m]:', telemetry['energy_per_distance'])
