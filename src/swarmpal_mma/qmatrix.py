#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:57:46 2024

@author: ngp
"""

import numpy as np
# import pandas as pd
import xarray as xr
import h5py
import statsmodels.api as sm


from swarmpal_mma.utils.SHA_utils import SHA_coeff_index_to_nm
from swarmpal_mma.utils.SHA_utils import design_SHA

def estimate_SH_coefficients_1D(data,params):
    '''    
    # % Estimate external SH coefficients using robust LS and 1-D impulse responses
    # % to model induction effect of the mantle
    # %
    # % Details are published here https://doi.org/10.1029/2020JA028672
    # %
    # % June 2021, Alexander Grayver, ETH
    Rewritten between July 2021 to Jan 2022 to be part of the production code (ngp)
    '''
    
    
    # % General constants
    a   = params['R_earth'] #6371.2
    rad = np.pi/180
    paramsdt=params['delt']/24 # dt in hours, convert to days
    '''
    Initialisation
    '''
    # model = pd.DataFrame(index=[0])

    
    n_e_terms = range(1,1+(params.n_max*(params.n_max+2)))
    nm_ext =np.array( SHA_coeff_index_to_nm(n_e_terms))
    terms_e = (abs(nm_ext[1]) <= params.m_max)
    


# % Exclude invalid data
    data=data.dropna(subset=['B_rtp_1','B_rtp_2','B_rtp_3'])
    

    # Assume you have data in NEC system:
    theta_tr = data.theta
    phi_tr = data.phi
    B_r_tr = data.B_rtp_1
    B_theta_tr = data.B_rtp_2
    B_phi_tr = data.B_rtp_3

        
        

    if len(data.t)>0 :
        print('Data range and integration time',max(data.t),min(data.t),paramsdt)
        n_bins = int(np.round((max(data.t) - min(data.t))/paramsdt))
    else:
        n_bins=0

    
    t1=np.array([i*paramsdt for i in range(n_bins)])+min(data.t)
    t2=t1+paramsdt
    
# 
    #load the Q matrix:
        
    with h5py.File(params.q_file, "r") as f:
        # Get the data
        Q_kernels = np.asarray(np.array(list(f['Q_kernels'])))
        
    
    kernel_length = np.round(params.n_lag_days / paramsdt)
    

    ds_t=t1 #+(paramsdt/2)
    
    Lm=len(terms_e)
    nan_array=np.empty(n_bins)
    nan_array[:]=np.nan
    
    da=xr.Dataset({"qs":([n_bins,Lm],np.zeros([n_bins, Lm])),
                  "gh":([n_bins,Lm],np.zeros([n_bins, Lm])),
                  "R2_e":np.zeros(n_bins),
                  "R2_i":np.zeros(n_bins),
                  "data_per_bin":nan_array.copy(),
                  "Cm_cond_e":nan_array.copy(),
                  "Cm_cond_i":nan_array.copy(),
                  "MSER_i":nan_array.copy(),
                  "MSER_e":nan_array.copy()
                  },coords={"time": ds_t})
    
    for i in range(n_bins):
        if(np.mod(i, 100) == 0):
            print('Done '+str(i)+' out of '+str(n_bins)+'\n')
    
        index = np.logical_and(data.t >= t1[i] , data.t < t2[i])
        if len(index[index])>15:

            t_bin = data.t[index]
            r_bin = data.r[index]
            theta_bin = theta_tr[index]
            phi_bin = phi_tr[index]
            dBt_bin = B_theta_tr[index]
            dBp_bin = B_phi_tr[index]
            dBr_bin = B_r_tr[index]
        
            if(len(t_bin) < sum(terms_e)):
                print('Time bin #'+str(i)+ ' has length of '+str(len(t_bin))+'!')
                continue        
            
            A_r_e, A_theta_e, A_phi_e = design_SHA(r_bin/a, \
                    theta_bin*rad, phi_bin*rad, params.n_max, 'ext')
            
            A_r_i, A_theta_i, A_phi_i = design_SHA(r_bin/a,\
                    theta_bin*rad, phi_bin*rad, params.n_max, 'int')
            
            B_theta_int = np.zeros(len(t_bin))
            B_phi_int = np.zeros(len(t_bin))
            
            if kernel_length<=len(Q_kernels.T):
                t_past = np.arange(i,max(i-kernel_length+2,0),-1,dtype='int')
            else:
                t_past = np.arange(len(Q_kernels.T),0,-1,dtype='int') \
                        + (i-len(Q_kernels.T))
            
            G =np.array(np.zeros(len(A_r_e)*2))
     
              
            idx = 0
            for n in range(1,params.n_max+1):
                n_indices = np.logical_and(nm_ext[0,:] == n, \
                            abs(nm_ext[1,:]) <= params.m_max)
                
                
                G_n = np.vstack([A_theta_e[:, n_indices] + \
                    Q_kernels[0,n-1]*A_theta_i[:, n_indices], \
                    A_phi_e[:, n_indices] + \
                    Q_kernels[0,n-1]*A_phi_i[:, n_indices]])
                
            
                
                
                Qt = Q_kernels[0:len(t_past),n-1]
                for m in range(0,min(params.m_max, n)+1):
       
                    q_nm = da.qs.data[t_past,idx]
                    B_theta_int = B_theta_int + np.multiply(Qt,q_nm).sum(axis=0) * \
                         A_theta_i[:, idx]
                    B_phi_int = B_phi_int + np.multiply(Qt,q_nm).sum(axis=0) * \
                        A_phi_i[:, idx]
                    
                    idx = idx + 1;
                    if m > 0:
                        s_nm = da.qs.data[t_past,idx]
                        B_theta_int = B_theta_int + \
                            np.multiply(Qt,s_nm).sum(axis=0) * \
                            A_theta_i[:, idx]
                        B_phi_int = B_phi_int +\
                            np.multiply(Qt,s_nm).sum(axis=0) * \
                            A_phi_i[:, idx]
                        idx = idx + 1
                
                G=np.column_stack((G, G_n))
                
            G=G[:,1:]
            
            yt = dBt_bin - B_theta_int
            yp = dBp_bin - B_phi_int
            
            da.data_per_bin.data[i]=len(yt)+len(yp)
            
            mod = sm.RLM(np.hstack([yt,yp]), G)
            da.qs.data[i] = mod.fit().params
            
            yr = dBr_bin - np.multiply(A_r_e,mod.fit().params).sum(axis=1)
            
            mod2= sm.RLM(yr,A_r_i)
            da.gh.data[i] = mod2.fit().params
    
            
            wls=sm.WLS(mod.endog, mod.exog, weights=mod.fit().weights).fit()
            da.R2_e.data[i]=wls.rsquared_adj
            da.Cm_cond_e.data[i]=wls.condition_number
            da.MSER_e.data[i]=wls.mse_resid
            '''
            check RMS
            '''
            wls=sm.WLS(mod2.endog, mod2.exog, weights=mod2.fit().weights).fit()
            da.R2_i[i].data=wls.rsquared_adj
            da.Cm_cond_i.data[i]=wls.condition_number
            da.MSER_i.data[i]=wls.mse_resid
        
    return da

