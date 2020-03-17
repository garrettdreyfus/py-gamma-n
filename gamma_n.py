import xarray as xr
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from profile import Profile
import gsw 
import pickle
from lib.bottle_to_cast import bottle_to_cast

refdata = xr.open_dataset('util/refprofiles.nc')
print(refdata.coords["lat"])

def loadProfiles(refdata):
    ## gamma_n goes lon, lat, pres
    profiles =[]
    lats = np.asarray(refdata.coords["lat"]).flatten() 
    lons = np.asarray(refdata.coords["lon"]).flatten() 
    for lon in range(len(lons)):
        profiles.append([])
        for lat in range(len(lats)):
            sals = np.asarray(refdata["SP_ref"][lon][lat][:]).flatten()
            if not np.isnan(sals).all():
                temps = np.asarray(refdata["t_ref"][lon][lat][:]).flatten()
                gamma = np.asarray(refdata["gamma_n_ref"][lon][lat][:]).flatten()
                pres = np.asarray(refdata["p_ref"][:]).flatten()
                #print(sals,temps,gamma,pres)
                profiles[-1].append(Profile(sals,temps,pres,gamma,lats[lat],lons[lon]))
            else:
                profiles[-1].append(np.nan)
    return profiles
            


def singlegamma_n(refdata,s,t,p,lon,lat):
    lats = np.asarray(refdata.coords["lat"]).flatten() 
    lons = np.asarray(refdata.coords["lon"]).flatten() 
    lati = int(np.floor(1 + (len(lats)-1)*(lat-lats[0])/(float(lats[-1]-lats[0]))))
    lati = [lati,lati-1]
    loni = int(np.floor(1+ (len(lons) - 1)*float(lon-lons[0])/(float(lons[-1]-lons[0]))))-1
    loni = [loni,loni+1]
    ref_s = []
    ref_t = []
    ref_p = []
    ref_gamma = []
    ds = []
    #hard coding combinations in to enforce order
    for coords in [[lati[1],loni[0]],[lati[0],loni[0]],[lati[0],loni[1]],[lati[1],loni[1]]]:
        coords = (int(coords[1]),int(coords[0]))
        reftemp = np.asarray(refdata.variables["t"][coords[0]][coords[1]])
        refsal = np.asarray(refdata.variables["s"][coords[0]][coords[1]])
        refgamma = np.asarray(refdata.variables["gamma"][coords[0]][coords[1]])
        refpres = np.asarray(refdata.coords["pres"])
        #refsal=refsal[~np.isnan(reftemp)]
        #refpres=refpres[~np.isnan(reftemp)]
        #refgamma=refgamma[~np.isnan(reftemp)]
        #reftemp=reftemp[~np.isnan(reftemp)]
        if ~np.isnan(reftemp).all():
            sref,tref,pref,gammaref = bottle_to_cast(s,t,p,refsal,reftemp,refpres,refgamma)
            #sref,tref,pref,gammaref = prof.neutralDepth(s,t,p) 
            ref_s.append(sref)
            ref_t.append(tref)
            ref_p.append(pref)
            ref_gamma.append(gammaref)
        else:
            ref_s.append(np.nan)
            ref_t.append(np.nan)
            ref_p.append(np.nan)
            ref_gamma.append(np.nan)
            
    ref_s = np.asarray(ref_s)
    ref_p = np.asarray(ref_p)
    ref_t = np.asarray(ref_t)
    ref_gamma = np.asarray(ref_gamma)

    #ref_s[np.isnan(ref_s)] = np.nanmean(ref_s)
    #ref_p[np.isnan(ref_p)] = np.nanmean(ref_p)
    #ref_t[np.isnan(ref_t)] = np.nanmean(ref_t)
    #ref_gamma[np.isnan(ref_gamma)] = np.nanmean(ref_gamma)

    rx = (lon-lons[loni[0]])/(lons[lati[0]+1]-lons[loni[0]])
    ry = (lat-lats[lati[0]])/(lats[lati[0]+1]-lats[lati[0]])
    
    gamma_n = (1-ry)*(ref_gamma[0] + rx*(ref_gamma[1] - ref_gamma[0])) + ry*(ref_gamma[3] + rx*(ref_gamma[2] - ref_gamma[3]));



    #scaling = (np.asarray(ds))/np.linalg.norm(np.asarray(ds))

    #if len(scaling) != 0:
        #return np.dot(scaling,ref_gamma)/np.sum(scaling), [np.dot(scaling,ref_s)/np.sum(scaling),np.dot(scaling,ref_t)/np.sum(scaling),np.dot(scaling,ref_p)/np.sum(scaling)]
    #else:
        #return np.nan,[np.nan]*3
    return gamma_n, [0,0,0]
    

singlegamma_n = partial(singlegamma_n,refdata)

def gamma_n(s,t,p,lon,lat):
    return singlegamma_n(s,t,p,lon,lat)

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

   
