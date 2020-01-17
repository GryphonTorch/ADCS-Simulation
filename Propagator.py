"""
Propagator
Orbit propagator using two body physics and standard Runge Kutta 4 integration
Use "import Propagator.oneStepAhead(x,y,z,vx,vy,vz,dt)" to calculate next step
state vector. 
"""

from numpy import *
import csv

GM = 3.986*(10**14)   # SI units
rE = 6378137          # meters
u = array([2315480.356,6240325.846,1359050.958,-7259.977,2459.492,1284.609],float)
#    test case^
J2 = 0.0010826266835531513
J3 = -0.0000025
J4 = -0.0000016
J5 = -0.00000015
J6 = 0.00000057
totalTime=860
dt=1.

def grad(p0,p1,p2,p3,p4,p5):   
    # faster version
    r = sqrt(p0**2+p1**2+p2**2)
    thetaP = 0.00007292115
    v = sqrt((p3+thetaP*p1)**2+(p4-thetaP*p0)**2+p5**2)
    return [p3,p4,p5,-GM*(p0)/r**3,-GM*(p1)/r**3,-GM*(p2)/r**3]

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

"""def grad(p0,p1,p2,p3,p4,p5):   
    r = sqrt(p0**2+p1**2+p2**2)
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
    return [p3,p4,p5,-GM*(p0)/r**3*Jx,-GM*(p1)/r**3*Jx,-GM*(p2)/r**3*Jz]"""


'''
outputList = [u]         
for i in range(0,int(totalTime/dt)):
    u = RK4(u[0],u[1],u[2],u[3],u[4],u[5],dt)
    outputList.append(u)

with open('SimulationData2.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["x", "y", "z", "vx" ,"vy", "vz"]) 
    for row in outputList:
        writer.writerow(row)
f.close()
'''   