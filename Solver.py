# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 15:37:18 2019
Solver module
Goal: output accurate, fast omegax/y/z values after one timestep
input: this step's omega, dt
Design notes: At first we used single timestep forward march, now upgraded to
Runge Kutta 4 -> much faster when timestep is 1 sec instead of 1 ms before 
IN SATELLITE BFPA
@author: user
"""

Ixx = 1   # moment of inertia along *principal* axes
Iyy = 1
Izz = 1  

def EulerEqnSolver(omegaX,omegaY,omegaZ, MomentX, MomentY, MomentZ, Ixx, Iyy, Izz, dt):
    """Numerically solves Euler equations of motion
    Inputs: current ang. velocity, torques, & moments of inertia in Principal Axes
    Outputs: next timestep's angular velocities in Principal Axes frame
    Direct forward march numerical integration
    Linearize omega dot across one timestep; expect some error over time"""
    
    omegaX0 = omegaX   # old variable, nth step
    omegaY0 = omegaY
    omegaZ0 = omegaZ
    
    omegaX = omegaX0 + (MomentX/Ixx + (Iyy-Izz)/Ixx*omegaY0*omegaZ0)*dt #(n+1)th step
    omegaY = omegaY0 + (MomentY/Iyy + (Izz-Ixx)/Iyy*omegaZ0*omegaX0)*dt
    omegaZ = omegaZ0 + (MomentZ/Izz + (Ixx-Iyy)/Izz*omegaX0*omegaY0)*dt
    return [omegaX,omegaY,omegaZ]


# def grad(p0,p1,p2,p3,p4,p5, M_old, M_new, dt):   
#     """Input: omegaX, omegaY, omegaZ, omegaXdot, omegaYdot, omegaZdot, previous
#     torque component (array), current torque component (array), dt"""
#     p0dot = p3   # by definition of omega dot
#     p1dot = p4 
#     p2dot = p5   
#     p3dot = 1/Ixx*(Iyy - Izz)*(p1*p5 + p4*p2) + (M_new[0] - M_old[0])/dt
#     p4dot = 1/Iyy*(Izz - Ixx)*(p2*p3 + p5*p0) + (M_new[1] - M_old[1])/dt
#     p5dot = 1/Izz*(Ixx - Iyy)*(p0*p4 + p3*p1) + (M_new[2] - M_old[2])/dt
    
#     """FOROGOT ABOUT CROSS TERMS !! 16 Feb 2020. Fix this"""
    
#     return [p0dot,p1dot,p2dot,p3dot,p4dot,p5dot]

# def EulerRK4(u0,u1,u2,u3,u4,u5, M_old, M_new, dt):   # standard RK4 implementation
#     """Outputs omegaX, omegaY and omegaZ in BFPA frame (required by Euler eqns) 
#     Code copied from propagator module.
#     Is M_old, M_new, dt placed in the correct position? a "Constant" in integration"""
    
#     k1 = grad(u0,u1,u2,u3,u4,u5, M_old, M_new, dt)
    
#     k2 = grad(u0+k1[0]*dt/2, u1+k1[1]*dt/2, u2+k1[2]*dt/2, \
#               u3+k1[3]*dt/2, u4+k1[4]*dt/2, u5+k1[5]*dt/2, M_old, M_new, dt)
    
#     k3 = grad(u0+k2[0]*dt/2, u1+k2[1]*dt/2, u2+k2[2]*dt/2,\
#               u3+k2[3]*dt/2, u4+k2[4]*dt/2, u5+k2[5]*dt/2, M_old, M_new, dt)
    
#     k4 = grad(u0+k3[0]*dt, u1+k3[1]*dt, u2+k3[2]*dt, \
#               u3+k3[3]*dt, u4+k3[4]*dt,u5+k3[5]*dt, M_old, M_new, dt)
    
#     res = [u0 + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), \
#            u1 + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]), \
#            u2 + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]), \
#            u3 + dt/6*(k1[3]+2*k2[3]+2*k3[3]+k4[3]), \
#            u4 + dt/6*(k1[4]+2*k2[4]+2*k3[4]+k4[4]), \
#            u5 + dt/6*(k1[5]+2*k2[5]+2*k3[5]+k4[5])]
    
#     return res