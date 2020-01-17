"""
Controller
Calculates deviation from target attitude and outputs control signal
Proposed Two Axis (pitch and roll) stabilization code assuming circular orbit
@author: Yu Jun
"""
import numpy as np
import math

def RateBasedControl(omegaX,omegaY,omegaZ, Kp):
    """Simple Control Method for Two Axis Stablisation (Roll and Pitch)
    Inputs: Spacecraft Attitude and State Vector, Kp Gain 
    Outputs: Roll and Pitch Control Currents
    
    Assume current omegas known, either from B-dot or accelerometer IRL"""
    Roll = - omegaX*Kp      # proportional gain. counter error
    Pitch = - omegaY*Kp
    return Roll,Pitch

def RollPitchControl(i0,j0,k0,x,y,z,vx,vy,vz, Kp):
    """Simple Control Method for Two Axis Stablisation (Roll and Pitch)
    Inputs: Spacecraft Attitude and State Vector, Kp Gain 
    Outputs: Roll and Pitch Control Currents 
    
    Note in this simulation, we *know* where the satellite is supposed to point
    at every moment in time in its orbit. So calculating an attitude "error" 
    makes sense. In real life, we must pre-load the satellite attitude history,
    then time-sync since it probably *cannot calculate* where it should point.
    
    Error defined as: Real - Target, with counterclockwise as +ve """    
    
    """Calculate Target Attitude"""                                                             
    #sub 1's: orbit frame unit vectors
    i1 = [vx,vy,vz]/np.linalg.norm([vx,vy,vz])   # spacecraft forward x direction in GCI coord
    k1 = [x,y,z]/np.linalg.norm([x,y,z])         # spacecraft out vector in GCI    
    #circular orbit assumed?
    j1 = np.cross(k1,i1)/np.linalg.norm(np.cross(k1,i1))    # orthogonal bases
    
    """Roll Error - along satellite i1 axis; align in j1k1 plane"""
    proj_jk = np.dot([i0,j0,k0],k1) + np.dot([i0,j0,k0],j1)    # no i1 component
    #print(np.dot(proj_jk,i0))     # should be zero, diagnostic
    Roll_cos = np.dot(proj_jk,k1)/np.linalg.norm(proj_jk)
    if np.dot(np.cross(k1, proj_jk),i1) > 0 :    # projection clockwise from k1
        Roll = math.acos(Roll_cos)               # Roll error > 0
    else:                                        # below, so -ve error
        Roll = -math.acos(Roll_cos)
  
    """PITCH ERROR: Torque vector along satellite j1 axis; align in i1k1 plane"""
    proj_ik = np.dot([i0,j0,k0],i1) + np.dot([i0,j0,k0],k1)    # no j1 component
    Pitch_cos = np.dot(proj_ik,i1)/np.linalg.norm(proj_ik)
    if np.dot(np.cross(i1,proj_ik),j1) > 0 :    # projection clockwise from i1
        Pitch = math.acos(Pitch_cos)            # Pitch error > 0 "nose up"
    else:                                       # below, so -ve error
        Pitch = - math.acos(Pitch_cos)

    #print("Roll angle error:", Roll, "Pitch angle error:", Pitch)
    MagTorque1 = -Roll*Kp      # proportional gain. counter error
    MagTorque2 = -Pitch*Kp
    
    return [MagTorque1,MagTorque2]

'''
#sub 0's: actual spacecraft orientations
i0 = np.array([0,1,0])   # forward. test values, to call in function! 
j0 = np.array([-1,0,0]) 
k0 = np.cross(i0,j0)     # radially out
x = 7000000
y = 400000
z = 1000
vx = 400
vy = 7000
vz = 200
Kp = 2 
print(PID(i0,j0,k0,x,y,z,vx,vy,vz,Kp))
'''
