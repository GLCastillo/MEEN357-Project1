#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:05:52 2026

@author: lauraalvarez
"""
import numpy as np
import matplotlib.pyplot as plt
from subfunctions import F_net, rover, planet

slope = 0

Crr_array = np.linspace(0.01, 0.5, 25)
vmax = np.zeros(len(Crr_array))

Ng = (0.07/0.04)**2
r = 0.3

for i in range(len(Crr_array)):
    Crr = float(Crr_array[i])
    
    vlow = 0
    vhigh = 3.8
    
    for j in range(50):
        vmid = (vhigh + vlow)/2
        
        omega = np.array([(vmid * Ng)/r])
        net_force = F_net(omega, slope, rover, planet, Crr)
        
        if net_force[0]>0:
            vlow = vmid
        else:
            vhigh = vmid
    
    vmax[i] = vmid

plt.plot(Crr_array, vmax)
plt.xlabel("Coefficient of rolling resistance")
plt.ylabel("Maximum velocity [m/s]")
plt.show()
