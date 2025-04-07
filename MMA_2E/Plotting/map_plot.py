        
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:51:27 2024

@author: ngp
"""
'''
Usage:
im=maps(clm,fig,ax,ll='SM',lmax=lmax,vi=vi,cm=cm,lim=lim,cobar=False

set up a map axis as:
fig, axes = plt.subplots(subplot_kw={'projection':ccrs.EqualEarth(180)},
                      figsize=(5,3))


INPUT:

    clm     array of the gauss coefficients
    fig     figure handle
    ax      axis to use
    ll      coordinate system, use 'SM' in solar magnetic
    lmax    is the max degree in the array clm
    vi      is the component you want to plot vi could be 'rad', 'theta','phi', 'total','horz'
    cm      is the colormap cm='viridis'
    lim     is defined if you want to set the limits in the colormap, do not defined to select limits automatically
    cobar   boolean variable to plot a color bar

'''


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


import pyshtools as pysh
import cartopy.crs as ccrs
# from cartopy.feature.nightshade import Nightshade
RE=6.371e6

def list2matrix(clm,lmax):
    coeffs = np.zeros((2, lmax+1, lmax+1))
    j=0
    for l in range(1,lmax+1):
         coeffs[0,l,0]=clm[j]
         j=j+1
         for m in range(1,l+1):
             coeffs[0,l,m]=clm[j]
             coeffs[1,l,m]=clm[j+1]
             
             j=j+2
    return coeffs

def myround(x, base=10,up=1):
    return base * (round(x/base)+ up)


def map_surface_rtp(ext,fig=False,ax=False,ll='',field='',cm='viridis',lmax=3,vi='theta',lim=0,cobar=True):

    coeffs= list2matrix(ext,lmax)
    cilm=pysh.SHMagCoeffs.from_array(coeffs,r0=RE,csphase=-1) 
    grid = cilm.expand(lmax_calc=lmax,lmax=15)
    
    if vi=='theta':
        data=grid.theta.to_xarray()
        field='{\\theta}'
        cl=1
    elif vi=='rad':
        data=grid.rad.to_xarray()
        field='r'
        cl=-1
    elif vi=='phi':
        data=grid.phi.to_xarray()
        field='{\\phi}'
        cl=-1
    elif vi=='total':
        data=grid.total.to_xarray()
        field='T'
        cl=1        
    elif vi=='horz':
        data=np.sqrt((grid.phi.to_xarray())**2+(grid.theta.to_xarray())**2)
        field='H'
        cl=1
    
    

    if lim==0:
        lim=max(myround(data.max().values),myround(abs(data.min().values)))
        lev=np.arange(-lim, lim+0.1, lim/10)
    else:
        if cl==-1:
            lev=np.arange(-lim, lim+0.1, lim/5)
        elif cl==1:
#            lim=myround(data.max().values)
            lev=np.arange(0, lim+0.1, lim/10)

    if not(fig) or not(ax):
        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.EqualEarth(180))



    img=ax.contourf(data.lon, data.lat, data.values,levels=lev,
                transform=ccrs.PlateCarree(),
                cmap=cm)
    if len(lev)<10:
        ax.contour(data.lon, data.lat, data.values,levels=lev,
                transform=ccrs.PlateCarree())#,

    ax.add_patch(patches.Rectangle(xy=[270, -90], width=90, height=180,
                                   linewidth=0,
                                   fill=True,
                                   fc='w',
                                   alpha=0.5,
                                   zorder=10,
                                   hatch='///',
                                   transform=ccrs.PlateCarree())
                 )
    ax.add_patch(patches.Rectangle(xy=[0, -90], width=90, height=180,
                                   linewidth=0,
                                   fill=True,
                                   fc='w',
                                   alpha=0.5,
                                   zorder=10,
                                   hatch='///',
                                   transform=ccrs.PlateCarree())
                 )

    if cobar:
        cbar=fig.colorbar(img,ax=ax, ticks=[lev[0],0,lev[-1]])
        cbar.set_label(r'$B_'+field+'$ (nT)')

    ax.set_global()

    return img  
    