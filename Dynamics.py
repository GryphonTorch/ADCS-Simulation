# -*- coding: utf-8 -*-
"""
3-D Rotational Dynamics
4 Dec 2019
@author: user
"""
import numpy as np

#i0, j0, k0 are unit vectors of satellite in inertial frame
i0 = np.array([0,1.00,0])   # test values. actually 2 out of 3 are needed
j0 = np.array([-1.00,0,0]) 
k0 = np.cross(i0,j0)
i = np.array([1,0,0])   # "correct" axis values
j = np.array([0,1,0])
k = np.array([0,0,1])
Torque = 0  # net torque
AtHistory = []

mass = 2  # kilogram
Ixx = 1/12*mass*(0.1**2+0.2**2)
Iyy = Ixx
Izz = 1/12*mass*(0.1**2+0.2**2)
MomentBFPA = np.array([[Ixx, 0, 0], [0, Iyy, 0], [0, 0, Izz]])
# moment of inertia tensor in fixed axis frame

TransformToInertial = np.array([[np.dot(i0,i), np.dot(j0,i), np.dot(k0,i)],\
                      [np.dot(i0,j), np.dot(j0,j), np.dot(k0,j)],\
                      [np.dot(i0,k), np.dot(j0,k), np.dot(k0,k)]])
# transformation matrix from old system vector to new system

Moment1 = np.matmul(TransformToInertial, MomentBFPA)
MomentIn= np.matmul(Moment1,TransformToInertial)

#print(MomentBFPA, "Along principal axes")
print("===Moment of Inertia in Inertial Axes===")
print(MomentIn)
#print(i0)

# keep track of two histories: MomentIn and Omega

dt = 1
Omega = [0,0.2,0.3]    # angular velocity disturbance vector
for n in range(5):   # check letter
    Torque = 0
    Omega0 = Omega    # record past omega
    print(Omega0 , "omega0 start")
    MomentIn0 = MomentIn
    Omega = (Torque + MomentIn[0]*Omega0 + MomentIn[1]*Omega0 + MomentIn[2]*Omega0) * \
            np.linalg.inv(2*MomentIn+MomentIn0)
    Omega = np.array([Omega.item(0),Omega.item(4),Omega.item(8)])
    # IDK if this is correct. seems only diagonal present?

    i0 = i0 + np.cross(Omega,i0)*dt
    j0 = j0 + np.cross(Omega,j0)*dt 
    k0 = k0 + np.cross(Omega,k0)*dt
    AtHistory.append([i0,j0,k0]) 
    #print(i0)
    
    
    n = n + 1
    #print(i0)