#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:22:45 2024

@author: ngp
"""
from viresclient import SwarmRequest
import numpy as np
import chaosmagpy as cp
import pandas as pd
import os
from datetime import datetime as dt


from utils.Coord_Trans import get_MagLat


ResidualsMC4='/users/swarm/L2PSProcessing/MMA_SHA_2E/development_bham/'
resname='_magneto_15sec_geo_only_alldata_nonans_2024.txt'

def LoadData(params,source='Vires'):
    

    dsa = get_Data(params.tini,params.tfin,'a')
    dsb = get_Data(params.tini,params.tfin,'b')
    # dsm = get_Data(params.tini,params.tfin,'m')

    data= SelectData(pd.concat([dsa,dsb]),params)
    
    return data

        
        
        
def get_Data(ts,te,sat,source='Vires',iono=True):

    
    if source == 'Vires':
        if iono==False:
            m= [
                    """
                    'Model' = 'CHAOS-Core'
                            + 'CHAOS-Static'
                    """
                 ]
        else:
            m= [
                    """
                    'Model' = 'CHAOS-Core'
                            + 'CHAOS-Static'
                            + 'MIO_SHA_2C'
                    """
                 ]
    
        col='SW_OPER_MAG'+sat.upper()+'_LR_1B'
        # for col in colectionlist:
        request = SwarmRequest()
        request.set_collection(col)
        request.set_products(
            measurements=["B_NEC"],
            models=m,
            residuals=True,
            auxiliaries=["MLT", "QDLat", "Dst", "QDBasis",
            'DipoleAxisVector'],
            sampling_step="PT25S" 
            )


        print('Getting data from Swarm'+sat.upper())
        d = request.get_between(ts,te)
        data = d.as_xarray()
        
        df=pd.DataFrame({  't'      : cp.data_utils.mjd2000(data.Timestamp),
                            'r'      : data.Radius/1000,                 #Rad in km
                            'phi'    : data.Longitude,
                            'theta'  : 90-data.Latitude,      #Colatitude in degrees
                            'B_rtp_1': -data.B_NEC_res_Model.data[:,2],  # Rad is -Center
                            'B_rtp_2': -data.B_NEC_res_Model.data[:,0],  # Theta is -North
                            'B_rtp_3': data.B_NEC_res_Model.data[:,1],
                            'sat'    : data.Spacecraft,
                            # 'MLT'    : data.MLT,
                          })
        
    else:
        
        filename = get_filename_2Fres(sat)
        jd_s=cp.data_utils.mjd2000(ts)
        jd_e=cp.data_utils.mjd2000(te)

        """------------------------------
        Load data from satellite 
        """
     
        df = pd.read_csv(filename,
                    skiprows=14,
                    header=None,
                    sep='\s+',
                    names=['t', 'r', 'phi', 'theta',
                            'B_rtp_1', 'B_rtp_2', 'B_rtp_3',
                            'r_noniono_residual',
                            'theta_noniono_residual',
                            'phi_noniono_residual', 'sat'
                            ])
        df.reset_index()
        """
        filter to the desired dates
        """
        df = df.drop(df.columns[[7,8,9,10]], axis=1)
        df = df[np.logical_and(df.t>=jd_s,df.t<=jd_e)]

        ''' Change angles to degrees'''
        r2d=180/np.pi
        df.phi=df.phi*r2d
        df.theta=df.theta*r2d
        df['sat']=sat.upper()
    
    df['time']=cp.data_utils.timestamp(df.t).astype(dt)
    return df.reset_index(drop=True)

# def DataSelection(params,data):if source=='MMA2F':
    
def get_filename_2Fres(sat):
    
    if sat.upper()=='M':
        name='MSS1A'
    else:
        name='swarm'+sat.upper()
        
    return os.path.join(ResidualsMC4,name,name[0].upper()+name[1:]+resname)


def SelectData(data,params):
    
    # data['time']=cp.data_utils.timestamp(data.t).astype(dt)

    ''' Build the mask for locat time '''
    LT=convert_longitude_to_local_time(data.phi,data.t)
    mask_lt= np.logical_or(LT < (params.LT_limit) , \
                                LT > ((24-params.LT_limit)))
    
    lat_mag=get_MagLat(90-data.theta,data.phi,data.time)
    mask_maglat= np.logical_and(np.abs(lat_mag) >=params.min_gm_lat, 
                                np.abs(lat_mag)<= params.max_gm_lat)

    return data[np.logical_and(mask_maglat,mask_lt)]
        

def convert_longitude_to_local_time(longitude_deg,t):
    longitude_deg = (longitude_deg/360) + (t%1)
    return (longitude_deg%1)*24