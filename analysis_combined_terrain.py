"""
analysis_combined_terrain.py
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from subfunctions import F_net, get_gear_ratio, rover, planet


Crr_array = np.linspace(0.01,0.5,25)

slope_array_deg = np.linspace(-15,35,25)

CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)

VMAX = np.zeros(np.shape(CRR), dtype = float)


Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

r = rover['wheel_assembly']['wheel']['radius']

vmax_possible = 3.8 * r / Ng

# Size of matrix
N = np.shape(CRR)[0]


# Loop through all combinations of Crr and slope
for i in range(N):

    for j in range(N):

        Crr_sample = float(CRR[i, j])
        slope_sample = float(SLOPE[i, j])

        # Initial velocity bounds
        vlow = 0
        vhigh = vmax_possible

        # Convert velocity bounds to motor shaft speed
        omega_low = np.array([(vlow * Ng) / r])
        omega_high = np.array([(vhigh * Ng) / r])

        # Calculate net force at each bound
        F_low = F_net(
            omega_low,
            slope_sample,
            rover,
            planet,
            Crr_sample
        )

        F_high = F_net(
            omega_high,
            slope_sample,
            rover,
            planet,
            Crr_sample
        )


        # Check if a root exists between vlow and vhigh
        if F_low[0] * F_high[0] <= 0:

            # Bisection method
            for k in range(50):

                vmid = (vhigh + vlow) / 2

                # Convert rover velocity to motor shaft speed
                omega = np.array([(vmid * Ng) / r])

                net_force = F_net(
                    omega,
                    slope_sample,
                    rover,
                    planet,
                    Crr_sample
                )

                # If net force is positive, rover is still accelerating
                if net_force[0] > 0:

                    vlow = vmid

                else:

                    vhigh = vmid


            # Store maximum velocity
            VMAX[i, j] = vmid


        else:

            # No terminal velocity found
            VMAX[i, j] = np.nan


# Create 3D surface plot
figure = plt.figure(figsize=(8, 10))

ax = figure.add_subplot(
    111,
    projection='3d'
)

ax.view_init(
    elev=25,
    azim=-135
)

ax.plot_surface(
    CRR,
    SLOPE,
    VMAX
)

ax.set_xlabel('Coefficient of Rolling Resistance')
ax.set_ylabel('Terrain Slope [deg]')
ax.set_zlabel('Maximum Velocity [m/s]')

ax.set_title(
    'Maximum Rover Velocity vs. Rolling Resistance and Terrain Slope'
)

plt.show()
