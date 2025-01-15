#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 09:44:42 2021

@author: ngp
"""

'''
Set of routines to create rotation matrices between various coordinate systems 
with the added functionality of extranting the Euler angles, as well as rotating
gauss coefficients.

The coordinate sysmtems included here come from spacepy functionality and are 
defined as follows:

Earth-centered Inertial Systems
-------------------------------
    * **ECI2000** Earth-centered Inertial, J2000 epoch
    * **ECIMOD** Earth-centered Inertial, mean-of-date
    * **ECITOD** Earth-centered Inertial, true-of-date
    * **GEI** Geocentric Equatorial Inertial (IRBEM approximation of TOD)

Magnetospheric Systems
----------------------
    * **GSM** Geocentric Solar Magnetospheric
    * **GSE** Geocentric Solar Ecliptic
    * **SM** Solar Magnetic
    * **MAG** Geomagnetic Coordinate System (aka CDMAG)

Earth-fixed Systems
-------------------
    * **GEO** Geocentric geographic, aka Earth-centered Earth-fixed
    * **GDZ** Geodetic coordinates

By convention *all* systems are treated as natively Cartesian except
geodetic (GDZ), which is defined in [altitude, latitude, longitude]
where altitude is relative to a reference ellipsoid. Similarly, distance
units are assumed to be Earth radii (Re) in all systems except GDZ, where
altitude is given in km. Conversions to GDZ will output altitude in km
regardless of the input distance units and conversions from GDZ will
output in Re regardless of input units. In all other cases, the distance
units will be preserved.
'''


import pyshtools
import numpy as np
# from datetime import datetime as dt
import math
from spacepy.coordinates import Coords      
from spacepy.time import Ticktock

def rotationMatrixToEulerAngles(R) :
    assert(isRotationMatrix(R))
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
    return np.array([x, y, z])

def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

def RM(t,from_str='GEO',to_str='SM'): 
    
    y = Coords([[1,0,0],[0,1,0],[0,0,1]], from_str, 'car',
               units=['km', 'km', 'km'], 
               use_irbem=False)
    
    if np.size(t)==1:    
        y.ticks= Ticktock([t,t,t], 'ISO')
        x = y.convert(to_str,'car')
        return np.vstack([x.x,x.y,x.z]).T
    elif np.size(t)>1:
        RMs=[]
        for ti in t:
            y.ticks= Ticktock([ti,ti,ti], 'ISO')
            x = y.convert(to_str,'car')            
            RMs.append(np.vstack([x.x,x.y,x.z]).T)
        return RMs
    
def EulerAngles(t,from_str='GEO',to_str='SM',Ms=[]):   
    if len(Ms)==0:
        Ms=RM(t,from_str,to_str)
    
    if np.size(t)==1:    
        return rotationMatrixToEulerAngles(Ms)

    elif np.size(t)>1:
        EA=[]
        for m in Ms:         
            EA.append(rotationMatrixToEulerAngles(m))
            
        return EA  
    
def rotate_gauss(clm,t,params,Ms=[],from_str='GEO',to_str='SM'):
  
    # dj=pyshtools.shtools.djpi2(params.n_max)
    dj=pyshtools.rotate.djpi2(params.n_max)
    try: 
        len(t)
        time_steps=len(t)
        if time_steps!=len(clm):
            print('The time array and the gauss coeffs are not the same length')
            return 
    except:
        time_steps=1
    
    if len(Ms)==0:
        Ms=RM(t,from_str=from_str,to_str=to_str)
    

    
    if time_steps>1: 
        clm_prime=[]
        for i in range(time_steps):
            Euler=EulerAngles(t[i],Ms=Ms[i])
            # print(t[i])
            clm_ppyf=pyshtools.rotate.SHRotateRealCoef(
                clm_format(clm[i].data,params), Euler, dj)
            clm_prime.append(clm_format(clm_ppyf,params))#, params.n_max))
        return clm_prime
       
    else:
        Euler=EulerAngles(t)
        # print(t)
        # return pyshtools.shtools.SHRotateRealCoef(
        clm_ppyf=pyshtools.rotate.SHRotateRealCoef(    
            clm_format(clm,params),Euler, dj)#, params.n_max)
        return clm_format(clm_ppyf,params)
     



def clm_format(clm,params):
    
    lmax=params.n_max
    ems=params.ms
    elles=params.ns    
    
    #Transform a list of gauss coefficient to the matrix style expected by pyshtools 
    if len(clm.shape)==1:   
        clim_out=np.zeros([2,lmax+1,lmax+1])
        
        for i in range(len(elles)):
            l=elles[i]
            m=ems[i]
            
            if m>=0:
                clim_out[0,l,m]=clm[i]
            else:
                clim_out[1,l,abs(m)]=clm[i]
                
        return clim_out
    else:
    # Transform the matrix type coefficients from matrix shape of pyshtools
    # to 1D arrays.     
        clim_out=np.zeros(len(elles))
        for i in range(len(elles)):
            l=elles[i]
            m=ems[i]
            if m>=0:
                clim_out[i]=clm[0,l,m]
            else:
                clim_out[i]=clm[1,l,abs(m)]
                
        return clim_out    

def get_MagLat(lat_geo,lon_geo,t,R=[]):
    
    deg2rad=np.pi/180
    R=RM(np.mean(t),'GEO','MAG')
    x=np.cos(lat_geo.values*deg2rad)*np.cos(lon_geo.values*deg2rad)
    y=np.cos(lat_geo.values*deg2rad)*np.sin(lon_geo.values*deg2rad)
    z=np.sin(lat_geo.values*deg2rad)
    A=np.vstack((x,y,z))
    
    A_mag=R @ A
    hmag=A_mag[0,:]**2+ A_mag[1,:]**2
    
    MagLat=[math.degrees(math.atan2(A_mag[2,i],hmag[i])) for i in range(len(hmag))]


    return MagLat

