# -*- coding: utf-8 -*-
"""
Magnetic field calculato in rotating Earth frame
    See MIT notes for tilted dipole model 
@author: user
"""
import numpy as np

def TiltedDipole(lat,long, r_vectorMag):
    """Analytical Model of Earth's Magnetic Field as a Tilted Dipole.
    Input: latitude and longitude radians, radial distance from Earth core
    Output: Magnetic field in inertial Earth x, y, z component
    (pick z is geographic North pole, y=0 at Greenwich)"""
    
    B_row1 = np.array([-np.cos(lat),np.sin(lat)*np.cos(long),np.sin(lat)*np.sin(long)])
    B_row2 = np.array([0,np.sin(long),-np.cos(long)])
    B_row3 = np.array([-2*np.sin(lat),-2*np.cos(lat)*np.cos(long),-2*np.cos(lat)*np.sin(long)])
    B_column = np.array([-29900,-1900,5530])   # from physics
    
    Mat1 = np.multiply(B_row1, B_column)
    B_north = Mat1[0] + Mat1[1] + Mat1[2]   # matrix multiplication for B_north
    Mat2 = np.multiply(B_row2, B_column)
    B_east = Mat2[0] + Mat2[1] + Mat2[2]    # B_east
    Mat3 = np.multiply(B_row3, B_column)
    B_down = Mat3[0] + Mat3[1] + Mat3[2]    # B_down
    BfieldRot = 10**(-9)*np.array([B_north,B_east,B_down])*(6378000/r_vectorMag)**3   
    # this outputs North, East, Down vectors. in Teslas
    
    # Now, in inertial frame
    Bfieldx = BfieldRot[0]*(-np.sin(lat)*np.cos(long)) + \
              BfieldRot[1]*np.cos(long) + BfieldRot[2]*(-np.cos(lat)*np.cos(long))
    Bfieldy = BfieldRot[0]*(-np.sin(lat)*np.sin(long)) + BfieldRot[1]*np.sin(long) +\
              BfieldRot[2]*(-np.cos(lat)*np.sin(long))
    Bfieldz = BfieldRot[0]*np.cos(lat) - BfieldRot[2]*np.sin(lat)
        
    return [Bfieldx,Bfieldy,Bfieldz]

