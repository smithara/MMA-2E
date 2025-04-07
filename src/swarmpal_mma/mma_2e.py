#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:07:02 2024

@author: ngp
"""
from datetime import datetime as dt
from scipy import stats
import pandas as pd
import numpy as np


import utils.configuration as config
from utils.GetResiduals import LoadData
from qmatrix import estimate_SH_coefficients_1D

def mma2e(ts=[],te=[]):
    
    params=config.BasicConfig()
    params.fullreset()
    
    if isinstance(ts,dt): params.tini=ts
    if isinstance(te,dt): params.tfin=te
    
    
    data= LoadData(params)
    
    threshold = 15
    outliers = pd.Series(data=False,index=data.index)
    for vi in range(1,4):
        vec_c='B_rtp_'+str(vi)
        z=np.abs(stats.zscore(data[vec_c]))
        outliers= (z > threshold) | outliers
    data=data[~outliers]

    mma = estimate_SH_coefficients_1D(data, params)
    
    return mma