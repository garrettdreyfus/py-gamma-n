import xarray as xr
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from profile import Profile
import gsw 
import pickle
from lib.bottle_to_cast import bottle_to_cast

refdata = xr.open_dataset('util/refprofiles.nc')

def replaceNAN(m):
    col_mean = np.nanmean(m, axis=0)
    inds = np.where(np.isnan(m))
    m[inds] = np.take(col_mean, inds[1])
    return m


def gamma_n(refdata,s,t,p,lon,lat):
    if isinstance(lat,float) or isinstance(lat,float)  :
        lat = np.asarray([lat]*len(s))
        lon = np.asarray([lon]*len(s))
    else:
        lat = np.asarray(lat)
        lon = np.asarray(lon)
    lats = np.asarray(refdata.coords["lat"]).flatten() 
    lons = np.asarray(refdata.coords["lon"]).flatten() 
    lati = np.floor(1 + (len(lats)-1)*(lat-lats[0])/(lats[-1]-lats[0])).astype(int)
    lati = [lati,lati-1]
    loni = ((np.floor(1+ (len(lons) - 1)*(lon-lons[0])/((lons[-1]-lons[0]))))-1).astype(int)
    loni = [loni,loni+1]
    ref_s = []
    ref_t = []
    ref_p = []
    ref_gamma = []
    #hard coding combinations in to enforce order
    refpres = np.empty((len(lat),len(refdata.coords["pres"])))
    for coords in [[lati[1],loni[0]],[lati[0],loni[0]],[lati[0],loni[1]],[lati[1],loni[1]]]:
        l1 = xr.DataArray(coords[1])
        l2 = xr.DataArray(coords[0])
        reftemp = np.asarray(refdata["t"][l1,l2])
        refsal = np.asarray(refdata["s"][l1,l2])
        refgamma = np.asarray(refdata["gamma"][l1,l2])
        refpres[:] = refdata.coords["pres"]
        sref,tref,pref,gammaref = bottle_to_cast(s,t,p,refsal,reftemp,refpres,refgamma)
        ref_s.append(sref)
        ref_t.append(tref)
        ref_p.append(pref)
        ref_gamma.append(gammaref)
            
    ref_s,ref_p,ref_t,ref_gamma  = replaceNAN(np.asarray(ref_s)),\
            replaceNAN(np.asarray(ref_p)),replaceNAN(np.asarray(ref_t)),replaceNAN(np.asarray(ref_gamma))

    rx = (lon-lons[loni[0]])/(lons[lati[0]+1]-lons[loni[0]])
    ry = (lat-lats[lati[0]])/(lats[lati[0]+1]-lats[lati[0]])
    
    gamma_n = (1-ry)*(ref_gamma[0] + rx*(ref_gamma[1] - ref_gamma[0])) + ry*(ref_gamma[3] + rx*(ref_gamma[2] - ref_gamma[3]));

    return gamma_n, [0,0,0]
    

gamma_n = partial(gamma_n,refdata)

def neutralsurfaces(s,t,p,gamma_n,surfaces):
    sals=[]
    temps=[]
    pres = []
    for surf in surfaces:
        solp = np.interp(surf,np.asarray(gamma_n).flatten(),np.asarray(p).flatten())
        pres.append(solp)
        sals.append(np.interp(solp,p,s))
        temps.append(np.interp(solp,p,t))
    return sals,temps,pres

   
