"""
CubeSat Attitude Determination and Control System Simulation
@author: Yu Jun
version 2.2: Bang Bang Control with Sun added
Next up: Solar power function (block sun in shadow region), 
documentation (diagram of coordinate systems)
Orbit checked with STK
Implement power calculations (flux, Earth horizon w/ Michael)
"""

import math
import numpy as np
from numpy import linalg as LA
import csv
import Propagator
import MagneticField
import Controller
import Solver
import Power

'''=========Test conditions (change this part only========='''

iniPos = np.array([-6719.4, 385.319, 2.669, -0.272368, -4.77507, 6.03443])
#iniPos = np.array([-5801232.89, 3520987.41, 6.751, -2136.353, -3519.9, 6120.767])
#ISS from STK

dt = 1 # unified throughout
q = 10   # data record rate (every q frames)
Duration = 4#160*60 # seconds
i0 = np.array([1,0,0])    # initial attitude of spacecraft, in inertial coord
j0 = np.array([0,1,0])   # these '0' vectors must be orthogonal and unit mag.
k0 = np.cross(i0,j0)
omegaX = 0.2    # Starting test values, in spacecraft frame. rad/s
omegaY = 0.1    # so correspond to roll/pitch/yaw
omegaZ = 0.3
'''========================================================'''


'''Satellite constants (input once the design is finalised)'''
Ix = 1   # moment of inertia along *principal* axes
Iy = 1
Iz = 1
turns = 10
area = 0.001
Kp = 0.02

'''======Computation constants (don't need to change)======'''
totalSteps = int(Duration/dt)
Pos = iniPos              # initialize state vector
i = np.array([1,0,0])     # unit vectors in Earth non-rotating inertial frame
j = np.array([0,1,0]) 
k = np.array([0,0,1]) 
TestData = np.array([omegaX,omegaY,omegaZ,Kp])   # save initial test data
Jx = 0 # initialise current in x torque coil
Jy = 0
Jz = 0
M_old = [0.0, 0.0, 0.0] # initialize torque history
sunAngle = 0   # assume in ecliptic 
History = []   # initialize records
i0data = []
Time = 0       # seconds

orbitDebug = []  # initialize empty list for testing

'''*****======Start loop======*****'''

for n in range(totalSteps):
    """I. 2 Body Forward Propagation"""
    print(Pos)
    newPos = Propagator.RK4(Pos[0],Pos[1],Pos[2],Pos[3],Pos[4],Pos[5],dt)
    # using classical 2-body only (good enough for now)
    x  = newPos[0]
    y  = newPos[1]
    z  = newPos[2]  
    vx = newPos[3]
    vy = newPos[4]
    vz = newPos[5]  
    Pos = newPos                            # close loop
    print(newPos)   # debug
    
    r_vectorMag = (x**2 + y**2 + z**2)**0.5 # magnitude of radius vector
    lat = np.arcsin(z/r_vectorMag)          # latitude, radian
    long = np.arctan2(y,x)                  # longitude; y=0 is Greenwich?  
    
    """II. Magnetic Field Calculation"""
    BfieldGCI = MagneticField.TiltedDipoleB(lat, long, r_vectorMag)  # in Earth non-rotating frame
    BfieldBFPA = MagneticField.GCItoBFPAtransform(i,j,k,i0,j0,k0,BfieldGCI[0],\
                                              BfieldGCI[1],BfieldGCI[2])
    BfieldNED = MagneticField.TiltedDipoleNED(lat, long, r_vectorMag)
    # transforms to satellite principal axes frame
    # sub 0's: actual spacecraft orientation, unit vectors
    
    """III. Magnetorquer Output"""
    '''The B-field information is used to calculate torque by working out,
    under nominal Bang Bang control. Then convert to BFPA frame. '''

    NetTorque = Controller.nominalTorqueBFPA(omegaX,omegaY,omegaZ,BfieldBFPA)

    """IV. Numerical Integration for omegas"""
    nextOmega = Solver.EulerEqnSolver(omegaX,omegaY,omegaZ, NetTorque[0], 
                                      NetTorque[1], NetTorque[0],Ix,Iy,Iz,dt)
      
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


    """ VI. Power calculation with dark side"""
    sunAngle = sunAngle + 2*np.pi/(24*60*60)*dt  # sun moves
    if sunAngle >= 2*np.pi:
        sunAngle = sunAngle - 2*np.pi            # keep to within 0, 2pi range
    # because at ISS inclinations the sun's 23 deg tilt won't affect coverage     
    power = Power.flux(long, sunAngle, i0, j0, k0) # compute power

    if (n//q)*q == n:
        Time = Time + q*dt
        Data = [Time,Pos[0]/1000,Pos[1]/1000,Pos[2]/1000,Pos[3]/1000,Pos[4]/1000,Pos[5]/1000,BfieldNED[0],\
            BfieldNED[1],BfieldNED[2],NetTorque[0],NetTorque[1],NetTorque[2],\
            omegaX, omegaY, omegaZ, power]
        History.append(Data)     # B vectors experienced
    
    orbitDebug.append([x,y,z,vx,vy,vz])    
    
'''*****======End loop======*****'''
    
print("Simulation done. Timestep used: ", dt, "sec. Data recorded every", q, "frames.")
print("Total time:", dt*totalSteps, "sec")
print("[omegaX, omegaY, omegaZ, Kp]:", TestData)

with open('SimulationData.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Time (s)","x (m)", "y (m)", "z (m)", "vx (m/s)" ,\
                     "vy (m/s)","vz (m/s)","B_North (nT)","B_East (nT)",\
                     "B_Down (nT)","RollTorque (Nm)","PitchTorque (Nm)", \
                     "YawTorque (Nm)","OmegaX (rad/s)",\
                     "OmegaY (rad/s)", "OmegaZ (rad/s)","Power (Watt)"]) 
    for row in History:
        writer.writerow(row)
f.close()    

#debug
with open('AttSim_orbitDebug.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Position"]) 
    for row in orbitDebug:
        writer.writerow(row)
f.close()    