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

    # time_strings = np.array([ti.isoformat() for ti in t])
    # spacepy_time = Ticktock(time_strings, "UTC")
    # positions = np.stack([np.ones(np.shape(lat_geo)),lat_geo, lon_geo]).T
    # geo_coords = Coords(positions, 'GEO', 'sph', 
    #                                     units=['Re','deg', 'deg'], 
    #                                     ticks=spacepy_time, 
    #                                     use_irbem=False)
    # mag_coords = geo_coords.convert("MAG", "sph")

    # return mag_coords.lati

    return MagLat


    
#    
# #            

#     cilm_SM = pyshtools.shtools.SHRotateRealCoef (gh_clim, Euler, dj, lmax)





# def SM2GEOC(rad_sm,theta_sm,phi_sm,time):
#     if len(time)==0:
#         time = spacepy_Ticktock(dt.now(), 'UTC')
#         time_strings = np.array([time.isoformat()])
#     else:
#         time_strings = np.array([_t.isoformat() for _t in time])
    
#     spacepy_time = spacepy_Ticktock(time_strings, "UTC",)
#     positions = np.stack([rad_sm, theta_sm, phi_sm]).T
#     spacepy_sm_coords = spacepy_Coords(positions, 'SM', 'sph', 
#                                         units=['deg', 'deg', 'm'], 
#                                         ticks=spacepy_time, 
#                                         use_irbem=False)
#     spacepy_geoc_coords = spacepy_sm_coords.convert("GEO", "sph")
#     # lat = spacepy_mag_coords.lati
#     # mlt = convert_longitude_to_local_time(spacepy_mag_coords.long)
#     return spacepy_geoc_coords.long, spacepy_geoc_coords.long



# def GEO_to_MAG_spacepy(GEO_lat, GEO_lon, GEO_rad, time):
#     time_strings = np.array([_t.isoformat() for _t in time])
#     spacepy_time = spacepy_Ticktock(time_strings, "UTC",)
#     positions = np.stack([GEO_rad, GEO_lat, GEO_lon]).T
#     spacepy_geo_coords = spacepy_Coords(positions, 'GEO', 'sph', units=['deg', 'deg', 'm'], ticks=spacepy_time, use_irbem=False)
#     spacepy_mag_coords = spacepy_geo_coords.convert("CDMAG", "sph")
#     mlat = spacepy_mag_coords.lati
#     mlt = convert_longitude_to_local_time(spacepy_mag_coords.long)
#     return mlat, mlt

# def convert_longitude_to_local_time(longitude_deg):
#     longitude_deg = (longitude_deg + 180)%360
#     return longitude_deg/15
















# #def rotation(params,t,gh)
 
# def GEO2SM(gh,qs,t,params):
    
        
#     lmax=params.n_max
# #    mmax=params.m_max
# #
# #    ems=params.ms
# #    elles=params.ns    
# #    
# #    wigner = spherical.Wigner(lmax)
# #    
# #    #The rotation matrix:
#     dj=pyshtools.shtools.djpi2(lmax)

#     gh_SM=np.zeros(np.size(gh[0]))
#     qs_SM=np.zeros(np.size(qs[0]))    
    
    
#     try: get_MagLat
#         len(t)
#         time_steps=len(t)
#     except:
#         time_steps=1
    
    
#     if time_steps>1: 
#         for i in range(time_steps):
#             t0=t[i]
#             gh_i=gh[i]
#             qs_i=qs[i]

#             mjd=datetoday2k(t0)
#             print(mjd)
#             R=GEO_SM_matrix(mjd)
#             Euler=rotationMatrixToEulerAngles(R)
# #        gh_SM_i=From_GEO_to_SM(gh_i,Euler,lmax,dj)
#             gh_SM=np.vstack([gh_SM,From_GEO_to_SM(gh_i,Euler,params,dj)])
# #                           ignore_index = True)
#             qs_SM=np.vstack([qs_SM,From_GEO_to_SM(qs_i,Euler,params,dj)])
       
#     else:
#         t0=t
#         gh_i=gh
#         qs_i=qs

#         mjd=datetoday2k(t0)
#         print(mjd)
#         R=GEO_SM_matrix(mjd)
#         Euler=rotationMatrixToEulerAngles(R)
# #        gh_SM_i=From_GEO_to_SM(gh_i,Euler,lmax,dj)
#         gh_SM=np.vstack([gh_SM,From_GEO_to_SM(gh_i,Euler,params,dj)])
# #                           ignore_index = True)
#         qs_SM=np.vstack([qs_SM,From_GEO_to_SM(qs_i,Euler,params,dj)])
# #       
    
#     return gh_SM[1:],qs_SM[1:]
    

# def SM2GEO(gh,qs,t,params):
    
        
#     lmax=params.n_max
# #    mmax=params.m_max
# #
# #    ems=params.ms
# #    elles=params.ns    
# #    
# #    wigner = spherical.Wigner(lmax)
# #    
# #    #The rotation matrix:
#     dj=pyshtools.shtools.djpi2(lmax)

#     gh_SM=np.zeros(np.size(gh[0]))
#     qs_SM=np.zeros(np.size(qs[0])) 
    
#     try: 
#         len(t)
#         time_steps=len(t)
#     except:
#         time_steps=1
    
    
#     if time_steps>1: 
#         for i in range(time_steps):
            
#             t0= t[i]
#             gh_i=gh[i]
#             qs_i=qs[i]
#             mjd=datetoday2k(t0)
#             print(mjd)
#             R=GEO_SM_matrix(mjd)
#             Euler=rotationMatrixToEulerAngles(np.transpose(R))
#     #        gh_SM_i=From_GEO_to_SM(gh_i,Euler,lmax,dj)
#             gh_SM=np.vstack([gh_SM,From_GEO_to_SM(gh_i,Euler,params,dj)])
#     #                           ignore_index = True)
#             qs_SM=np.vstack([qs_SM,From_GEO_to_SM(qs_i,Euler,params,dj)])
#     else:
#         t0=t
#         gh_i=gh
#         qs_i=qs
#         mjd=datetoday2k(t)
#         print(mjd)
#         R=GEO_SM_matrix(mjd)
#         Euler=rotationMatrixToEulerAngles(np.transpose(R))
# #        gh_SM_i=From_GEO_to_SM(gh_i,Euler,lmax,dj)
#         gh_SM=From_GEO_to_SM(gh_i,Euler,params,dj)
# #                           ignore_index = True)
#         qs_SM=From_GEO_to_SM(qs_i,Euler,params,dj)
        
# #       
#     return gh_SM[1:],qs_SM[1:]
        

# def From_GEO_to_SM(gh_i,Euler,params,dj):
    
#     lmax=params.n_max
# #    mmax=params.m_max

#     ems=params.ms
#     elles=params.ns   
# #    N=len(elles)
# #    
# #    A_lm=np.zeros(np.size(elles))
# #    
#     gh_clim=np.zeros([2,lmax+1,lmax+1])
    
#     for i in range(len(elles)):
#         l=elles[i]
#         m=ems[i]
        
#         if m>=0:
#             gh_clim[0,l,m]=gh_i[i]
#         else:
#             gh_clim[1,l,abs(m)]=gh_i[i]
# #    
# #            

#     cilm_SM = pyshtools.shtools.SHRotateRealCoef (gh_clim, Euler, dj, lmax)
    
#     gh_SM=np.zeros(np.size(gh_i))

#     for i in range(len(elles)):
#         if ems[i]>=0:
#             gh_SM[i]=cilm_SM[0,elles[i],ems[i]]

#         else:
#             gh_SM[i]=cilm_SM[1,elles[i],abs(ems[i])]
            
#     return  gh_SM
    


# def GEO_SM_matrix(mjd):            
#     M_GEO_SM= makeRotMatrix('GSM_SM',mjd) @ \
#         makeRotMatrix('GSE_GSM',mjd) @ \
#         makeRotMatrix('GEI_GSE',mjd) @ \
#         makeRotMatrix('GEI_GEO',mjd).transpose()
        
#     return M_GEO_SM

    
# def makeRotMatrix(direction,mjd):
    
#     i=[1,0,0]
#     j=[0,1,0]
#     k=[0,0,1] 
    
#     if direction=='GSM_SM':
#         igrf=get_IGRF_coef(mjd).values
# #        31 and 13 have the wrong sign:
#         R =np.array([coordt.GSM_SM(mjd,i,igrf),\
#             coordt.GSM_SM(mjd,j,igrf),\
#             coordt.GSM_SM(mjd,k,igrf)])
#     elif direction == 'GSE_GSM':
#         igrf=get_IGRF_coef(mjd).values
# #        31 and 13 have the wrong sign:
#         R =np.array([coordt.GSE_GSM(mjd,i,igrf),\
#             coordt.GSE_GSM(mjd,j,igrf),\
#             coordt.GSE_GSM(mjd,k,igrf)]    )
#     elif direction == 'GEI_GSE':
#         R =np.array([coordt.GEI_GSE(mjd,i),\
#             coordt.GEI_GSE(mjd,j),\
#             coordt.GEI_GSE(mjd,k)]    )
#     elif direction == 'GEI_GEO':
#         R =np.array([coordt.GEI_GEO(mjd,i),\
#             coordt.GEI_GEO(mjd,j),\
#             coordt.GEI_GEO(mjd,k)]    )
    
#     return R.transpose()
    
# def get_IGRF_coef(mjd):
    
#     import glob
#     import os    
#     import gmpylib.get_data.IGRF_reader as read
    
#     list_of_files = glob.glob('/users/swarm/L2PSInput/Auxiliary/*IGR*DBL') # * means all if need specific format then *.csv
#     latest_file = max(list_of_files, key=os.path.getctime)
#     fpath=latest_file.split('/')
#     IGRF=read.read_IGRF(fpath[-1])
    
#     '''
#     Find the igrf coefficients before and after the desired date
#     '''
#     [n,r]=divmod(day2ktodate(mjd).year,5)
#     y0=5*n
#     '''
#     number of days in 5 years
#     '''
#     Dmjd5y=datetoday2k(dt(y0+5,1,1))-datetoday2k(dt(y0,1,1))
#     '''
#     earlier value of the igrf used
#     '''
#     igrf0=IGRF[str(y0)+'-01-01'][0:3]
#     '''
#     Number of days since the earlier value of the igrf values
#     '''
#     Dmjd=mjd-datetoday2k(dt(y0,1,1))
    
#     '''
#     linear interpolation fo the desired date
#     '''
    
#     return ((IGRF[str(y0+5)+'-01-01'][0:3]-igrf0)*Dmjd/Dmjd5y) + igrf0
    


# def GEOC2SM_R(time):
#     B0=1
    
# def GEOC2SM_gauss(clms,time):
#     from scipy.spatial.transform import Rotation   
#     mlat,mlon=GEO2SM_rtp(1000,0,0,time)
    
    
#     ### first transform the matrix to euler angles
#     r =  Rotation.from_matrix(rotation_matrix)
#     angles = r.as_euler("zyx",degrees=True)
    

# def GEO2MAG_rtp(rad_gg,theta_gg,phi_gg,time,i_trans=1):
#     if len(time)==0:
#         time = spacepy_Ticktock(dt.now(), 'UTC')
#         time_strings = np.array([time.isoformat()])
#     else:
#         time_strings = np.array([_t.isoformat() for _t in time])
    
#     spacepy_time = spacepy_Ticktock(time_strings, "UTC",)
#     positions = np.stack([rad_gg, theta_gg, phi_gg]).T
#     spacepy_geo_coords = spacepy_Coords(positions, 'GEO', 'sph', 
#                                         units=['deg', 'deg', 'm'], 
#                                         ticks=spacepy_time, 
#                                         use_irbem=False)
#     spacepy_mag_coords = spacepy_geo_coords.convert("MAG", "sph")
#     mlat = spacepy_mag_coords.lati
#     mlon = spacepy_mag_coords.long
#     return mlat, mlon
    
    
# def GEO2SM_rtp(rad_gg,theta_gg,phi_gg,time,i_trans=1):
#     if len(time)==0:
#         time = spacepy_Ticktock(dt.now(), 'UTC')
#         time_strings = np.array([time.isoformat()])
#     else:
#         time_strings = np.array([_t.isoformat() for _t in time])
    
#     spacepy_time = spacepy_Ticktock(time_strings, "UTC",)
#     positions = np.stack([rad_gg, theta_gg, phi_gg]).T
#     spacepy_geo_coords = spacepy_Coords(positions, 'GEO', 'sph', 
#                                         units=['deg', 'deg', 'm'], 
#                                         ticks=spacepy_time, 
#                                         use_irbem=False)
#     spacepy_mag_coords = spacepy_geo_coords.convert("SM", "sph")
#     mlat = spacepy_mag_coords.lati
#     mlt = convert_longitude_to_local_time(spacepy_mag_coords.long)
#     return mlat, mlt

# def SM2GEOC_rtp(rad_sm,theta_sm,phi_sm,time,i_trans=1):
#     if len(time)==0:
#         time = spacepy_Ticktock(dt.now(), 'UTC')
#         time_strings = np.array([time.isoformat()])
#     else:
#         time_strings = np.array([_t.isoformat() for _t in time])
    
#     spacepy_time = spacepy_Ticktock(time_strings, "UTC",)
#     positions = np.stack([rad_sm, theta_sm, phi_sm]).T
#     spacepy_sm_coords = spacepy_Coords(positions, 'SM', 'sph', 
#                                         units=['deg', 'deg', 'm'], 
#                                         ticks=spacepy_time, 
#                                         use_irbem=False)
#     spacepy_geoc_coords = spacepy_sm_coords.convert("GEO", "sph")
#     # lat = spacepy_mag_coords.lati
#     # mlt = convert_longitude_to_local_time(spacepy_mag_coords.long)
#     return spacepy_geoc_coords.long, spacepy_geoc_coords.long