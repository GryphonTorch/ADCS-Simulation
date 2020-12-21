"""
Propagator File (21 Dec 2020)
J6 Orbit propagator using [x,y,z,vx,vy,vz] state vector, Runge Kutta 4 integration
For Earth orbits. Distances all in km, not meters. 

More info on Runge Kutta method: 
https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods 
Simple check at the bottom (#commented)
"""
from numpy import *
import numpy as np
import csv

GM = 398600.      # Standard Gravitational Parameter, in units of km^3 s^-2
rE = 6378.137     # Radius of Earth, in km
# Earth zonal harmonics, dimensionless
J2 = 0.0010826266835531513
J3 = -0.0000025
J4 = -0.0000016
J5 = -0.00000015
J6 = 0.00000057
dt=1.     # timestep

def grad(p0,p1,p2,p3,p4,p5):     # RK4 gradient function
    r = sqrt(p0**2 + p1**2 + p2**2)   # Earth radius for J term calculations

    Jx = 1 - J2*(3./2.)*(rE/r)**2*(5*p2**2/r**2-1) + \
        J3*(5./2.)*(rE/r)**3*(3*p2/r-7*p2**3/r**3) - \
        J4*(5./8.)*(rE/r)**4*(3-42*p2**2/r**2+63*p2**4/r**4) - \
        J5*(3./8.)*(rE/r)**5*(35*p2/r-210*p2**3/r**3+231*p2**5/r**5) + \
        J6*(1./16.)*(rE/r)**6*(35-945*p2**2/r**2+3465*p2**4/r**4-3003*p2**6/r**6)
    Jz = 1 + J2*(3./2.)*(rE/r)**2*(3-5*p2**2/r**2) + \
        J3*(3./2.)*(rE/r)**3*(10*p2/r-(35./3.)*p2**3/r**3-r/p2) - \
        J4*(5./8.)*(rE/r)**4*(15-70*p2**2/r**2+63*p2**4/r**4) - \
        J5*(1./8.)*(rE/r)**5*(315*p2/r-945*p2**3/r**3+693*p2**5/r**5-15*p2/r) + \
        J6*(1./16.)*(rE/r)**6*(315-2205*p2**2/r**2+4851*p2**4/r**4-3003*p2**6/r**6)
    thetaP = 0.00007292115
    v = sqrt((p3+thetaP*p1)**2+(p4-thetaP*p0)**2+p5**2)
    return [p3,p4,p5,-GM*(p0)/r**3*Jx,-GM*(p1)/r**3*Jx,-GM*(p2)/r**3*Jz]

def RK4(u0,u1,u2,u3,u4,u5,dt):   # standard RK4 implementation
    k1 = grad(u0,u1,u2,u3,u4,u5)
    k2 = grad(u0+k1[0]*dt/2, u1+k1[1]*dt/2, u2+k1[2]*dt/2, \
              u3+k1[3]*dt/2, u4+k1[4]*dt/2, u5+k1[5]*dt/2)
    k3 = grad(u0+k2[0]*dt/2, u1+k2[1]*dt/2, u2+k2[2]*dt/2,\
              u3+k2[3]*dt/2, u4+k2[4]*dt/2, u5+k2[5]*dt/2)
    k4 = grad(u0+k3[0]*dt, u1+k3[1]*dt, u2+k3[2]*dt, \
              u3+k3[3]*dt, u4+k3[4]*dt,u5+k3[5]*dt)
    res = [u0 + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), \
           u1 + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]), \
           u2 + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]), \
           u3 + dt/6*(k1[3]+2*k2[3]+2*k3[3]+k4[3]), \
           u4 + dt/6*(k1[4]+2*k2[4]+2*k3[4]+k4[4]), \
           u5 + dt/6*(k1[5]+2*k2[5]+2*k3[5]+k4[5])]
    return res

#Simple check
#print(RK4(-6719.4, 385.319, 2.669, -0.272368, -4.77507, 6.03443,1))
