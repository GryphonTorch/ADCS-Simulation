"""
CubeSat Attitude Determination and Control System Simulation
@author: Yu Jun
version 1.2: fixed missing B-field dependence; torque computed in spacecraft frame
"""

import math
import numpy as np
import csv

import Propagator
import MagneticField
import Controller
import Solver

iniPos = np.array([2315480.356,6240325.846,1359050.958,\
                   -7259.977,2459.492,1284.609])   
# based on ISS orbit
Pos = iniPos

Ix = 1   # moment of inertia along *principal* axes
Iy = 1
Iz = 1
turns = 10
area = 0.001
Kp = 0.2
dt = 0.001 # unified throughout
q = 1000   # data record rate (every 100th frame)
Duration = 10 # seconds
totalSteps = int(Duration/dt)

i = np.array([1,0,0])     # unit vectors in Earth non-rotating inertial frame
j = np.array([0,1,0]) 
k = np.array([0,0,1]) 

i0 = np.array([0,1,0])    # initial attitude of spacecraft, in inertial coord
j0 = np.array([-1,0,0])   # these '0' vectors must be orthogonal and unit mag.
k0 = np.cross(i0,j0)
omegaX = 1    # Starting test values, in spacecraft frame. rad/s
omegaY = 2    # so correspond to roll/pitch/yaw
omegaZ = 0.1
TestData = np.array([omegaX,omegaY,omegaZ,Kp])

History = []   # initialize records
i0data = []
Time = 0       # seconds

# in the loop
for n in range(totalSteps):
    """I. 2 Body Forward Propagation"""
    newPos = Propagator.RK4(Pos[0],Pos[1],Pos[2],Pos[3],Pos[4],Pos[5],dt)
    # using classical 2-body only (good enough for now)
    x  = newPos[0]
    y  = newPos[1]
    z  = newPos[2]  
    vx = newPos[3]
    vy = newPos[4]
    vz = newPos[5]  
    Pos = newPos                            # close loop
    
    r_vectorMag = (x**2 + y**2 + z**2)**0.5 # magnitude of radius vector
    lat = np.arcsin(z/r_vectorMag)          # latitude, radian
    long = np.arctan2(y,x)                  # longitude; y=0 is Greenwich?  
    
    """II. Magnetic Field Calculation"""
    BfieldGCI = MagneticField.TiltedDipole(lat, long, r_vectorMag)  # in Earth non-rotating frame
    Bfield = MagneticField.GCItoBFPAtransform(i,j,k,i0,j0,k0,BfieldGCI[0],\
                                              BfieldGCI[1],BfieldGCI[2])
    # transforms to satellite principal axes frame
    # sub 0's: actual spacecraft orientation, unit vectors
    
    """III. Magnetorquer Output"""
    '''The B-field information is used to calculate torque by working out
    the current that gives the required torque i.e. we choose a control
    architecture for the torque, and then there is a formula to get currents
    that takes into account B. We must check if power usage is okay.'''
    
    Currents = Controller.RateFieldControl(omegaX,omegaY,omegaZ,Bfield[0],\
                                           Bfield[1],Bfield[2],Kp,turns,area)
    
    NetTorque = - Kp*np.array([omegaX, omegaY, omegaZ])

    """IV. Numerical Integration for omegas"""
    nextOmega = Solver.EulerEqnSolver(omegaX, omegaY, omegaZ, NetTorque[0],\
                               NetTorque[1],  NetTorque[2], Ix,Iy,Iz,dt)
    # set yaw torque (Mz) = 0 in two axis control
    omegaX = nextOmega[0]
    omegaY = nextOmega[1]
    omegaZ = nextOmega[2]
    
    """V. Effecting omega"""  # ensure unit vectors remain unit magnitude, orthogonal
    # omegaX. rotate the principal axes accordingly. omegaX = Roll. x0 supposed to be forward facing
    d0X = omegaX*dt   # radian, small angle 
    i0new = i0
    j0new = j0*math.cos(d0X) + k0*math.sin(d0X)
    k0new = k0*math.cos(d0X) - j0*math.sin(d0X)
    
    i0 = i0new
    j0 = j0new
    k0 = k0new   # put in the new values

    # omegaY = Pitch
    d0Y = omegaY*dt   # radian, small angle 
    j0new = j0
    i0new = i0*math.cos(d0Y) - k0*math.sin(d0Y)
    k0new = k0*math.cos(d0Y) + i0*math.sin(d0Y)
    
    i0 = i0new
    j0 = j0new
    k0 = k0new   # put in the new values

    # omegaZ = Yaw
    d0Z = omegaZ*dt   # radian, small angle 
    k0new = k0
    i0new = i0*math.cos(d0Z) + j0*math.sin(d0Z)
    j0new = j0*math.cos(d0Z) - i0*math.sin(d0Z)
    
    i0 = i0new
    j0 = j0new
    k0 = k0new   # put in the new values

    if (n//q)*q == n:
        Time = Time + q*dt
        Data = [Time,Pos[0],Pos[1],Pos[2],Pos[3],Pos[4],Pos[5],Bfield[0],\
            Bfield[1],Bfield[2],NetTorque[0],NetTorque[1],omegaX, omegaY, omegaZ]
        History.append(Data)     # B vectors experienced

print(omegaX)
print("Simulation done. Timestep", dt, "sec. Data recorded every", q, "frames.")
print("Total time:", dt*totalSteps, "sec")
print("[omegaX, omegaY, omegaZ, Kp]:", TestData)

with open('SimulationData.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Time (s)","x (m)", "y (m)", "z (m)", "vx (m/s)" ,"vy (m/s)",\
                     "vz (m/s)","B_x (T)","B_y (T)","B_z (T)",\
                     "RollTorque (Nm)","PitchTorque (Nm)","OmegaX (rad/s)",\
                     "OmegaY (rad/s)", "OmegaZ (rad/s)",]) 
    for row in History:
        writer.writerow(row)
f.close()    
