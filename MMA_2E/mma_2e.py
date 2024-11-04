#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:07:02 2024

@author: ngp
"""
from datetime import datetime as dt

import utils.configuration as config
from qmatrix import estimate_SH_coefficients_1D

from utils.GetResiduals import LoadData

def mma2e(ts=[],te=[]):
    
    params=config.BasicConfig()
    params.fullreset()
    
    if isinstance(ts,dt): params.tini=ts
    if isinstance(te,dt): params.tfin=te
    
    
    data= LoadData(params)
    
    mma = estimate_SH_coefficients_1D(data, params)