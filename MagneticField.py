# -*- coding: utf-8 -*-
"""
Magnetic Field File (24 Nov 2020)
Contains two B field calculations: one in rotating Earth frame (North, East, Down)
and another in the x,y,z frame (with a vector transformation)

See MIT notes for tilted dipole mode (Slide 34):
https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-851-satellite-engineering-fall-2003/lecture-notes/l9_acs.pdf
"""
import numpy as np
    
def TiltedDipoleXYZ(lat,long, r_vectorMag):
    """outputs magnetic field in ECI -> ""x,y,z" frame"""
    matrix = np.array([[-np.cos(lat), np.sin(lat)*np.cos(long), np.sin(lat)*np.sin(long)],\
                       [0,np.sin(long), -np.cos(long)],\
                       [-2*np.sin(lat), -2*np.cos(lat)*np.cos(long), -2*np.cos(lat)*np.sin(long)]])
    vector = np.array([-29900, -1900, 5530])

    Bfield = (6378/r_vectorMag)**3*matrix.dot(vector)
    #In north, east and down currently. Use 6378 which is Earth radius in km.
    
    Bfieldx = Bfield[0]*(-np.sin(lat)*np.cos(long)) + \
              Bfield[1]*np.sin(long) + Bfield[2]*(-np.cos(lat)*np.cos(long))
    Bfieldy = Bfield[0]*(-np.sin(lat)*np.sin(long)) + Bfield[1]*np.cos(long) +\
              Bfield[2]*(-np.cos(lat)*np.sin(long))
    Bfieldz = Bfield[0]*np.cos(lat) - Bfield[2]*np.sin(lat)
    # nanoTesla
    return np.array([Bfieldx*10**(-9) , Bfieldy*10**(-9) , Bfieldz*10**(-9) ])   
    
def GCItoBFPAtransform(i,j,k,i0,j0,k0,x,y,z):
    """Transforms vector x,y,z in coordinate frame with unit vectors i,j,k
    into vector x0,y0,z0 in coordinate frame with unit vectors i0,j0,k0. 
    Use: transform B field from Geocentric Inertial Frame to spacecraft
    Body-Fixed Principal axes frame. See MIT Dynamics lecture for math.  
    Outputs: new vector x0,y0,z0."""
    
    x0 = np.dot(i0,i)*x + np.dot(i0,j)*y + np.dot(i0,k)*z
    y0 = np.dot(j0,i)*x + np.dot(j0,j)*y + np.dot(j0,k)*z
    z0 = np.dot(k0,i)*x + np.dot(k0,j)*y + np.dot(k0,k)*z

    return np.array([x0,y0,z0])

def TiltedDipoleNED(lat,long, r_vectorMag):
    """Debug function to compare against STK. 
    Returns array of B field, in North, East, Down (NED) components
    in nanoTesla"""
    
    B_row1 = np.array([-np.cos(lat),np.sin(lat)*np.cos(long),
                       np.sin(lat)*np.sin(long)])
    B_row2 = np.array([0,np.sin(long),-np.cos(long)])
    B_row3 = np.array([-2*np.sin(lat),-2*np.cos(lat)*np.cos(long),
                       -2*np.cos(lat)*np.sin(long)])
    B_column = np.array([-29900,-1900,5530])   # from physics
    
    Mat1 = np.multiply(B_row1, B_column)
    B_north = Mat1[0] + Mat1[1] + Mat1[2]   # matrix multiplication for B_north
    Mat2 = np.multiply(B_row2, B_column)
    B_east = Mat2[0] + Mat2[1] + Mat2[2]    # B_east
    Mat3 = np.multiply(B_row3, B_column)
    B_down = Mat3[0] + Mat3[1] + Mat3[2]    # B_down
    BfieldRot = np.array([B_north,B_east,B_down])*(6378/r_vectorMag)**3*10**(-9)    # in nanoTesla       
    
    return BfieldRot

#Check for difference
#print(TiltedDipoleXYZ(1,2, 7000))
#print(TiltedDipoleNED(1,2, 7000))