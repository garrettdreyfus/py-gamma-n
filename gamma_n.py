import h5py
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from profile import Profile
import gsw 
import pickle

refdata = h5py.File('refdata.mat', 'r')

def loadProfiles(refdata):
    ## gamma_n goes lon, lat, pres
    profiles =[]
    lats = np.asarray(refdata["lat_ref"]).flatten() 
    lons = np.asarray(refdata["long_ref"]).flatten() 
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
            
            
lon = 187.317;
lat = -41.6667;

SP =[35.066,35.086,35.089,35.078,35.025,34.851,34.696,34.572,34.531,34.509,34.496,34.452,34.458,34.456,34.488,34.536,34.579,34.612,34.642,34.657,34.685,34.707,34.72,34.729]

t = [12.25,12.21,12.09,11.99,11.69,10.54,9.35,8.36,7.86,7.43,6.87,6.04,5.5,4.9,4.04,3.29,2.78,2.45,2.211,2.011,1.894,1.788,1.554,1.38]

p = [1.0,48.0,97.0,145.0,194.0,291.0,388.0,485.0,581.0,678.0,775.0,872.0,969.0,1066.0,1260.0,1454.0,1647.0,1841.0,2020.0,2216.0,2413.0,2611.0,2878.0,3000.0]

gamma_n_known = [26.65,26.682830469203406,26.710963096615604,26.723242299110460,26.741488538021695,26.825445912051336,26.918689217252997,26.989761790054338,27.039067923101946,27.089143151019517,27.166567035269665,27.260376554533835,27.343619695291586,27.421578895148251,27.557338511940429,27.698188932980081,27.798443363873236,27.866285802482334,27.920185440895871,27.959264296723756,27.997866000490600,28.031679411184577,28.079958980601589,28.117372360538731]

SP_ns_known =[34.90641,34.630793019073700,34.724828398226208]

t_ns_known = [10.90640,2.300296198425469,1.460659594687408]

p_ns_known = [260.109,1953.131680473056,2943.451620399693]

gamma_n_surfaces = [26.8, 27.9, 28.1]



def singlegamma_n(refprofiles,refdata,s,t,p,lon,lat):
    lats = np.asarray(refdata["lat_ref"]).flatten() 
    lons = np.asarray(refdata["long_ref"]).flatten() 
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
        prof = refprofiles[coords[0]][coords[1]]
        if prof:
            sref,tref,pref,gammaref = prof.neutralDepth(s,t,p) 
            ref_s.append(sref)
            ref_t.append(tref)
            ref_p.append(pref)
            ref_gamma.append(gammaref)
            ds.append(np.linalg.norm(((lon-lons[coords[0]])/(lons[loni[1]]-lons[loni[0]]),(lat-lats[coords[1]])/(lats[lati[1]]-lats[lati[0]]))))
        else:
            ref_s.append(np.nan)
            ref_t.append(np.nan)
            ref_p.append(np.nan)
            ref_gamma.append(np.nan)
    ref_s = np.asarray(ref_s)
    ref_p = np.asarray(ref_p)
    ref_t = np.asarray(ref_t)
    ref_gamma = np.asarray(ref_gamma)
    ref_s[np.isnan(ref_s)] = np.nanmean(ref_s)
    ref_p[np.isnan(ref_p)] = np.nanmean(ref_p)
    ref_t[np.isnan(ref_t)] = np.nanmean(ref_t)
    ref_gamma[np.isnan(ref_gamma)] = np.nanmean(ref_gamma)

    rx = (lon-lons[loni[0]])/(lons[lati[0]+1]-lons[loni[0]])
    ry = (lat-lats[lati[0]])/(lats[lati[0]+1]-lats[lati[0]])

    gamma_n = (1-ry)*(ref_gamma[0] + rx*(ref_gamma[1] - ref_gamma[0])) + ry*(ref_gamma[3] + rx*(ref_gamma[2] - ref_gamma[3]));




    scaling = (np.asarray(ds))/np.linalg.norm(np.asarray(ds))

    return np.dot(scaling,ref_gamma)/np.sum(scaling), [np.dot(scaling,ref_s)/np.sum(scaling),np.dot(scaling,ref_t)/np.sum(scaling),np.dot(scaling,ref_p)/np.sum(scaling)]
    #return gamma_n, [0,0,0]
    


with open('profiles.pickle', 'wb') as outfile:
    pickle.dump(loadProfiles(refdata),outfile)
with open('profiles.pickle', 'rb') as outfile:
    profiles = pickle.load(outfile)

singlegamma_n = partial(singlegamma_n,profiles,refdata)

def gamma_n(s,t,p,lon,lat):
    solutions=[]
    solutiont = []
    solutionp = []
    solutiongamma =[]
    for l in range(len(p)):
        outg,[outs,outt,outp]= singlegamma_n(s[l],t[l],p[l],lon,lat)
        solutions.append(outs)
        solutiont.append(outt)
        solutionp.append(outp)
        solutiongamma.append(outg)
    return solutiongamma, [solutions,solutiont,solutionp]

def neutralsurfaces(s,t,p,gamma_n,surfaces):
    sals=[]
    temps=[]
    pres = []
    for surf in surfaces:
        solp = np.interp(surf,gamma,p)
        pres.append(solp)
        sals.append(np.interp(solp,p,s))
        temps.append(np.interp(solp,p,t))
    return sals,temps,pres

    
SP = gsw.SA_from_SP(SP,p,lon,lat)
t = gsw.CT_from_t(SP,t,p)
gamma,debug = gamma_n(SP,t,p,lon,lat)

knowns,knownt,knownp = neutralsurfaces(SP,t,p,gamma,gamma_n_surfaces)
 
fig, axs = plt.subplots(1,3)
fig.suptitle('comparisons')
axs[0].scatter(range(3),knowns,label="python")
axs[0].scatter(range(3),gsw.SA_from_SP(SP_ns_known,p_ns_known,lon,lat),label="matlab")
axs[1].scatter(range(3),knownt,label="python")
axs[1].scatter(range(3),gsw.CT_from_t(gsw.SA_from_SP(SP_ns_known,p_ns_known,lon,lat),t_ns_known,p_ns_known),label="matlab")
axs[2].scatter(range(3),knownp,label="python")
axs[2].scatter(range(3),p_ns_known,label="matlab")
axs[0].legend()
axs[1].legend()
axs[2].legend()
axs[0].set_ylabel("Salinity")
axs[1].set_ylabel("Temperature")
axs[2].set_ylabel("Pressure")
plt.show()


