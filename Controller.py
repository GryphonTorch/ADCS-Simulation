"""
Controller
Calculates deviation from target attitude and outputs control signal
Proposed Two Axis (pitch and roll) stabilization code assuming circular orbit
@author: Yu Jun
"""
import numpy as np
import math

def RateFieldControl(omegaX,omegaY,omegaZ,Bx,By,Bz,K,turns,Area):
    """Generates currents in coils based on both magnetic field and torque
    desired. Assume proportional control in setting torque using angular 
    velocity (error)"""
    
    TorqueX = -K*omegaX   # Proportional control
    TorqueY = -K*omegaY
    TorqueZ = -K*omegaZ

    Ix = K/(turns*Area)*(By*TorqueZ - Bz*TorqueY)  # from paper
    Iy = K/(turns*Area)*(Bz*TorqueX - Bx*TorqueZ)
    Iz = K/(turns*Area)*(Bx*TorqueY - By*TorqueX)

    return [Ix,Iy,Iz]
