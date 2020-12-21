"""
Power v3.1 with sun

Calculates power based on satellite position with respect to sun
"""
import math
import numpy as np

def flux(long, sunAngle, i0, j0, k0):
    """Calculates solar power received by satellite
    Input: satellite longitude, sun angle and current satellite orientations in i0, j0, k0 vectors
    Output: power (Watt)
    k0 points upwards from top of satellite (solar panel exists there)"""
    
    efficiency2 = 0.307 * 0.88 * (1 - (75-28) * 0.0022) # efficiency of Pumpkin solar panel
    efficiency1 = 0.295 * 0.88 * (1 - (75-28) * 0.0022) # efficiency of Endurosat solar panel
    
    Area2U = 0.01076664    # Area of one 2U Pumpkin panel in m^2
    Area1U = 0.00603       # Area of one 1U Endurosat panel in m^2
    # using arrays for inertial i and k unit vectors
    
    phi = 1373*(math.cos(23)*np.array([1,0,0]) - math.sin(23)*np.array([0,0,1])) 
    # solar flux vector
    
    powerTop = phi.dot(k0)*Area1U
    if powerTop < 0:   # top is sunlit
        powerTop = efficiency1*abs(powerTop)
    else:
        powerTop = 0

    powerSidei0 = efficiency2*abs(phi.dot(i0)*Area2U)  # don't double count
    powerSidej0 = efficiency2*abs(phi.dot(j0)*Area2U)*0.75    #reduced to 1U on a side
    
    powerAll = powerTop + powerSidei0 + powerSidej0    # before considering dark side
    
    if abs(sunAngle - long) < np.pi/2:
        powerAll = powerAll               # lighted up
    else:
        powerAll = 0                      # in shadow of Earth
    
    return powerAll
