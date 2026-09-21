import numpy as np
import matplotlib.pyplot as plt
from subfunctions import F_net, rover, planet

Crr = 0.15

slope_array_deg = np.linspace(-15, 35, 25)
v_max = np.zeros(len(slope_array_deg))

Ng = (0.07 / 0.04)**2
r = 0.3

for i in range(len(slope_array_deg)):
    slope = slope_array_deg[i]
    
    vlow = 0
    vhigh = 0.372
    
    for j in range(50):
        vmid = (vhigh + vlow) / 2
      
        omega = np.array([((vmid * Ng) / r)])

        net_force = F_net(omega, slope, rover, planet, Crr)
        
        if net_force[0] > 0:
            vlow = vmid
        else:
            vhigh = vmid
    
    v_max[i]=vmid
    
plt.plot(slope_array_deg, v_max)
plt.xlabel("Terrain slope [degrees]")
plt.ylabel("Maximum velocity [m/s]")
plt.show()
