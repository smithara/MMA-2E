#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 09:59:20 2024

@author: ngp
"""
import numpy as np
import pyshtools


def SHA_coeff_index_to_nm(indices):
    # % function [n,m] = SHA_coeff_index_to_nm(indices)
    # % Compute degrees (n) and orders (m - +/-) for the given
    # % coefficient indices corresponding to the output of design_SHA
        
    n = np.floor(np.sqrt(indices))
      
    #% Compute orders using three step approach
    i = indices - (n-1)*(n+1)
    j = ((i > 1) & (np.mod(i,2) == 1))
    m = np.floor(i/2)
    m[j] = -m[j]
      
    return n, m


def design_SHA(r, theta, phi, N, i_e_flag='int'):
    '''
    % [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N)
    %
    % Calculates design matrices A_i that connects the vector 
    % of (Schmidt-normalized) spherical harmonic expansion coefficients, 
    % x = (g_1^0; g_1^1; h_1^1; g_2^0; g_2^1; h_2^1; ... g_N^N; h_N^N) 
    % and the magnetic component B_i, where "i" is "r", "theta" or "phi":
    %        B_i = A_i*x
    % Input: r(:)      radius vector (in units of the reference radius a)
    %        theta(:)  colatitude    (in radians)
    %        phi(:)    longitude     (in radians)
    %        N         maximum degree/order
    %
    % [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N, i_e_flag)
    % with i_e_flag = 'int' for internal sources (g_n^m and h_n^m)
    %                 'ext' for external sources (q_n^m and s_n^m) 
    %
    % Uses MEX file design_SHA_m if available, and Matlab program design_SHA_matlab.m else
    
    % January 2003, Nils Olsen, DSRI
    Modified Jan 2022, ngp
    '''
    
    '''% determine size of input arrays'''
    max_size = max(np.size(r, 0),np.size(theta, 0),np.size(phi, 0))
    max_length = max(r.size, theta.size, phi.size)

# % convert to matrix if input parameter is scalar
    if len(r)     == 1:
        r = r*np.ones(max_size)
    if len(theta) == 1:
        theta = theta*np.ones(max_size)
    if len(phi)   == 1: 
        phi = phi*np.ones(max_size)
        
# % check for equal length of all input parameters
    if np.size(r) != np.size(theta) or np.size(r) != np.size(phi):
        print('Variables must be of equal size (or scalars)')
    

# % convert to row vector
    r = np.reshape(r, max_length)
    theta = np.reshape(theta, max_length)
    phi = np.reshape(phi, max_length)

    A_r, A_theta, A_phi = design_SHA_matlab(r, theta, phi, N, i_e_flag)

    return A_r, A_theta, A_phi

# % -----------------------------------------------------------------------

def design_SHA_matlab(r, theta, phi, N, i_e_flag='int'):
    '''
    % [A_r, A_theta, A_phi] = design_SHA_matlab(r, theta, phi, N)
    %
    % Calculates design matrices A_i that connects the vector 
    % of (Schmidt-normalized) spherical harmonic expansion coefficients, 
    % x = (g_1^0; g_1^1; h_1^1; g_2^0; g_2^1; h_2^1; ... g_N^N; h_N^N) 
    % and the magnetic component B_i, where "i" is "r", "theta" or "phi":
    %        B_i = A_i*x
    % Input: r(:)      radius vector (in units of the reference radius a)
    %        theta(:)  colatitude    (in radians)
    %        phi(:)    longitude     (in radians)
    %        N         maximum degree/order
    %
    % [A_r, A_theta, A_phi] = design_SHA(r, theta, phi, N, i_e_flag)
    % with i_e_flag = 'int' for internal sources (g_n^m and h_n^m)
    %                 'ext' for external sources (q_n^m and s_n^m) 
    %
    
    % March 2001, Nils Olsen, DSRI
    
    Modified for python in January 2022 (ngp)
    '''

    N_koeff=((N+1)**2)-1
    
    k=0
    nm=np.zeros([N+1,N+1])-100
    mat_gnm=nm.copy()
    mat_hnm=nm.copy()
    
    for n in range(0,N+1):
        for m in range(0,n+1):
            if m==0:
                nm[n,0]=k
                k=k+1
            else:
                nm[n,m]=k
                k=k+1
    k=0
    p=1

    for n in range(1,N+1):
   
      for m in range(n+1):
         
         if m == 0:
            k = k+1
            mat_gnm[n,0]=k
#            i_pnm[p]=1
            p=p+1
            # % terms corresponding to g_n^0
         else:
            k = k+1
            mat_gnm[n,m]=k
#            i_pnm[p]=(-1)**m
            p=p+1
            
            k = k+1
            mat_hnm[n,m]=k

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    N_data    = len(theta)

    A_r       = np.asarray(np.zeros([N_data, N_koeff]))
    A_theta   = np.asarray(np.zeros([N_data, N_koeff]))
    A_phi     = np.asarray(np.zeros([N_data, N_koeff]))

    Pnm = []
    dPnm= []

    for z in cos_theta:
        ar, dar = pyshtools.legendre.PlmSchmidt_d1(N,z,csphase=1)     # % P_n^m and derivatives vrt. theta
        Pnm.append(ar)
        dPnm.append(-(dar)*np.sqrt(1-z*z))
    Pnm=np.array(Pnm)
    dPnm=np.array(dPnm)
    
    for n in range(1,N+1):
       if i_e_flag =='int':
          r_n       = r**(-(n+2))
       elif i_e_flag == 'ext':
          r_n       = r**(n-1)
       else:
          print( 'i_e_flag neither "int" nor "ext". Assumed "int".')
          r_n       = r**(-(n+2))
       
       if i_e_flag == 'int':  # % internal sources by default
          for m in range(n+1):
             cos_phi   = np.cos(m*phi)
             sin_phi   = np.sin(m*phi)
             nm_i      = int(nm[n,m])
             gnm       = int(mat_gnm[n,m])-1
             hnm       = int(mat_hnm[n,m])-1
             
             if m == 0:
                
                A_r[:,gnm]       =  (n+1)*r_n*Pnm[:,nm_i]
                A_theta[:,gnm]   = -r_n *dPnm[:,nm_i]
                A_phi[:,gnm]     =  r_n*0
                
             else:
                
                A_r[:,gnm]       =  (n+1)*r_n*cos_phi*Pnm[:,nm_i]
                A_theta[:,gnm]   = -r_n*cos_phi*dPnm[:,nm_i]
                A_phi[:,gnm]     =  r_n*m*sin_phi*Pnm[:,nm_i]/sin_theta
                
                A_r[:,hnm]       =  (n+1)*r_n*sin_phi*Pnm[:,nm_i]
                A_theta[:,hnm]   = -r_n*sin_phi*dPnm[:,nm_i]
                A_phi[:,hnm]     = -r_n*m*cos_phi*Pnm[:,nm_i]/sin_theta 
             
                
       else:  # % external sources, i_e_flag ='ext'
          for m in range(n+1):
             cos_phi   = np.cos(m*phi)
             sin_phi   = np.sin(m*phi)
             nm_i      = int(nm[n,m])
             gnm       = int(mat_gnm[n,m])-1
             hnm       = int(mat_hnm[n,m])-1
             
             if m == 0:
                
                A_r[:,gnm]       = -n*r_n*Pnm[:,nm_i]
                A_theta[:,gnm]   = -r_n*dPnm[:,nm_i]
                A_phi[:,gnm]     =  r_n*0
             else:
                 
                A_r[:,gnm]       = -n*r_n*cos_phi*Pnm[:,nm_i]
                A_theta[:,gnm]   = -r_n*cos_phi*dPnm[:,nm_i]
                A_phi[:,gnm]     =  r_n*m*sin_phi*Pnm[:,nm_i]/sin_theta
                
                A_r[:,hnm]       = -n*r_n*sin_phi*Pnm[:,nm_i]
                A_theta[:,hnm]   = -r_n*sin_phi*dPnm[:,nm_i]
                A_phi[:,hnm]     = -r_n*m*cos_phi*Pnm[:,nm_i]/sin_theta
             
    return A_r, A_theta, A_phi
