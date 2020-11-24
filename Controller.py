"""
Controller File (24 Nov 2020)
Calculates output torque based on input error
Proposed Two Axis (pitch and roll) stabilization code assuming circular orbit
"""
import numpy as np
import math

def nominalTorqueBFPA(OmegaX, OmegaY, OmegaZ, BfieldBFPA):
    """
    Input: current errors and Bfield to find corrective 
    torques (all in BPFA). Omegas are floats, BfieldBFPA is array
    Output: Torque array in BFPA, to use directly in Solver module.
    
    0.25 Amp*m2 is the nominal working dipole strength. This method is 
    BANG-BANG CONTROL : magnetorquer switches on or off only, at 0.25 Ampm^2
    dipole strength. This value from Nanoavionics spec sheet, projeted to 
    consume 140 mW nominally. 
    """
    
    TorqueX = np.array([0,0,0])   # initialize local variables
    TorqueY = np.array([0,0,0])   # default bang bang buffer zone
    TorqueZ = np.array([0,0,0])
    dipoleX = np.array([0,0,0])   # vectors
    dipoleY = np.array([0,0,0])
    dipoleZ = np.array([0,0,0])
    omega = np.array([OmegaX, OmegaY, OmegaZ])

    dipoleX = np.cross(np.array([1,0,0]),BfieldBFPA)
    if dipoleX.dot(omega) < -10E-7:
        TorqueX = 0.25*dipoleX        # run in +ive i0 direction
#        print("dipoleX dot:", dipoleX.dot(omega), "+i0")
    elif dipoleX.dot(omega) > 10E-7:
        TorqueX = -0.25*dipoleX       # run in -ve i0 direction
#        print("dipoleX dot:", dipoleX.dot(omega), "-i0")
    
    dipoleY = np.cross(np.array([0,1,0]),BfieldBFPA)
    if dipoleY.dot(omega) < -10E-7:
        TorqueY = 0.25*dipoleY        # run in +ive i0 direction
#        print("dipoleY dot:", dipoleY.dot(omega), "+j0")
    elif dipoleY.dot(omega) > 10E-7:
        TorqueY = -0.25*dipoleY       # run in -ve i0 direction
#        print("dipoleY dot:", dipoleY.dot(omega), "-j0")

    dipoleZ = np.cross(np.array([0,0,1]),BfieldBFPA)
    if dipoleZ.dot(omega) < -10E-7:
        TorqueZ = 0.25*dipoleZ        # run in +ive i0 direction
#        print("dipoleZ dot:", dipoleZ.dot(omega), "+k0")
    elif dipoleZ.dot(omega) > 10E-7:
        TorqueZ = -0.25*dipoleZ       # run in -ve i0 direction
#        print("dipoleZ dot:", dipoleZ.dot(omega), "-k0")

    return TorqueX + TorqueY + TorqueZ  # do all at once

Bfield = np.array([-8.11481912e-06,  5.47955177e-06,  2.39722085e-05])
#nominalTorqueBFPA(1, 0, 0, Bfield)

def RateFieldControl(omegaX,omegaY,omegaZ,Bx,By,Bz,K,turns,Area):
    """Generates currents in coils based on both magnetic field and torque
    desired. Assume proportional control in setting torque using angular 
    velocity (error)
    Do not use for now. """
    
    TorqueX = -K*omegaX   # Proportional control
    TorqueY = -K*omegaY
    TorqueZ = -K*omegaZ

    Ix = K/(turns*Area)*(By*TorqueZ - Bz*TorqueY)  # from paper
    Iy = K/(turns*Area)*(Bz*TorqueX - Bx*TorqueZ)
    Iz = K/(turns*Area)*(Bx*TorqueY - By*TorqueX)

    return [Ix,Iy,Iz]
