
from datetime import datetime as dt
from scipy import stats
import pandas as pd
import numpy as np
import chaosmagpy as cp
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import calendar


import sys
ROOT = os.path.join(os.path.dirname(__file__),'..','MMA_2E')
sys.path.append(ROOT)
import utils.configuration as config
from utils.GetResiduals import LoadData
from qmatrix import estimate_SH_coefficients_1D




def cleanSat(data):
    
    threshold = 15
    outliers = pd.Series(data=False,index=data.index)
    for vi in range(1,4):
        vec_c='B_rtp_'+str(vi)
        z=np.abs(stats.zscore(data[vec_c]))
        outliers= (z > threshold) | outliers
    return data[~outliers]


def getData(ts=dt(2024,10,1),te=dt(2024,11,1)):
    params=config.BasicConfig()
    params.fullreset()
    
    

    
    if isinstance(ts,dt): params.tini=ts
    if isinstance(te,dt): params.tfin=te
    
    
    data_2F= LoadData(params,source='MMA2F',mss=True)
    data_2F=cleanSat(data_2F)
    
    data_Vi= LoadData(params,source='Vires')
    data_Vi=cleanSat(data_Vi)
    return data_Vi,data_2F,params

ts=dt(2024,10,1)
te=dt(2024,11,1)
data_Vi,data_2F,params =getData(ts,te)
