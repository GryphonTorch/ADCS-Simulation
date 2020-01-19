# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 15:37:18 2019
Solver
Goal: output accurate, fast omegax/y/z values after one timestep
input: this step's omega, dt

@author: user
"""

def EulerEqnSolver(omegaX,omegaY,omegaZ, MomentX, MomentY, MomentZ, Ix, Iy, Iz, dt):
    """Numerically solves Euler equations of motion
    Inputs: current ang. velocity, torques, & moments of inertia in Principal Axes
    Outputs: next timestep's angular velocities in Principal Axes frame
    Direct forward march numerical integration
    Linearize omega dot across one timestep; expect some error over time"""
    
    omegaX0 = omegaX   # old variable, nth step
    omegaY0 = omegaY
    omegaZ0 = omegaZ
    
    omegaX = omegaX0 + (MomentX/Ix + (Iy-Iz)/Ix*omegaY0*omegaZ0)*dt #(n+1)th step
    omegaY = omegaY0 + (MomentY/Iy + (Iz-Ix)/Iy*omegaZ0*omegaX0)*dt
    omegaZ = omegaZ0 + (MomentZ/Iz + (Ix-Iy)/Iz*omegaX0*omegaY0)*dt
    return omegaX,omegaY,omegaZ



def functionX(omegaX,omegaY,omegaZ,omegaXprev,MomentX,Ix, Iy, Iz, dt):
    """For Newton method. needs previous step omegaX0 to calculate omega dot"""
    fx = Ix*(omegaX-omegaXprev)/dt - (Iy-Iz)*omegaY*omegaZ - MomentX
    return fx

def functionY(omegaX,omegaY,omegaZ,omegaYprev,MomentY,Ix, Iy, Iz, dt):
    """For Newton method. needs previous step omegaX0 to calculate omega dot"""
    fy = Iy*(omegaY-omegaYprev)/dt - (Iz-Ix)*omegaZ*omegaX - MomentY
    return fy

def functionZ(omegaX,omegaY,omegaZ,omegaZprev,MomentZ,Ix, Iy, Iz, dt):
    """For Newton method. needs previous step omegaX0 to calculate omega dot"""
    fz = Iz*(omegaZ-omegaZprev)/dt - (Ix-Iy)*omegaX*omegaY - MomentZ
    return fz

"""
    fX = Solver.functionX(omegaX,omegaY,omegaZ,omegaXprev,\
                          MomentX,Ix, Iy, Iz, dt)
    fprimeX = (fX - fXprev)/dt
    omegaXprev = omegaX            # save n step 
    omegaX = omegaX - fX/fprimeX   # Newton's method for n+1 step
    
    fY = Solver.functionX(omegaX,omegaY,omegaZ,omegaYprev,\
                          MomentY,Ix, Iy, Iz, dt)
    fprimeY = (fY - fYprev)/dt
    omegaYprev = omegaY
    omegaY = omegaY - fY/fprimeY

    fZ = Solver.functionZ(omegaX,omegaY,omegaZ,omegaZprev,\
                          MomentZ,Ix, Iy, Iz, dt)
    fprimeZ = (fZ - fZprev)/dt
    omegaZprev = omegaZ
    omegaZ = omegaZ - fZ/fprimeZ

    # Update variables    
    fXprev = fX    # for Newton method
    fYprev = fY
    fZprev = fZ
    
    #INITIALIZE AT THE START
    fXprev = 0   # initialize values
    fYprev = 0
    fZprev = 0
"""