#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:07:02 2024

@author: ngp
"""

# from utils.Config import BasicConfig as config
import utils.configuration as config
from datetime import datetime as dt

from utils.GetResiduals import LoadData

def mma2e(ts,te):
    
    params=config.BasicConfig()
    params.fullreset()
    
    params.tini=ts
    params.tfin=te
    
    
    data= LoadData(params)
    
    mma = estimate_SH_coefficients_1D(data, params)