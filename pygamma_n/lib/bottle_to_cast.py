import numpy as np
import gsw
#import matplotlib.pyplot as plt

def diagMult(a,b):
    return np.nansum((a * b),axis=-1)

def zeroCross(Es,Esvalues,startEs,values):
    dP = diagMult(Es,values)
    dE = diagMult(Es,Esvalues)
    firstE =  diagMult(startEs,Esvalues)
    firstP =  diagMult(values,-startEs)
    psol = firstE * (dP/dE) + firstP
    return psol

def bottle_to_cast(s,t,p,s_ref,t_ref,p_ref,gamma_ref):
    #assuming s,t,p are vectors
    #we'd like s_ref,t_ref, and p_ref to be len(p)xlen(p_ref)
    if s.shape != t.shape and t.shape != p.shape:
        raise ValueError("s,t and p must all have the same dimensions")
    if isinstance(s,float):
        s=np.asarray([s])
        t=np.asarray([t])
        p=np.asarray([p])
    Sref = s_ref
    Tref = t_ref
    Gammaref = gamma_ref
    Pref = p_ref

    Sdata = np.column_stack([s]*p_ref.shape[1])
    Tdata = np.column_stack([t]*p_ref.shape[1])
    Pdata = np.column_stack([p]*p_ref.shape[1])

    refdens = gsw.rho(Sref,Tref,(Pref+Pdata)/2.0)
    datadens = gsw.rho(Sdata,Tdata,(Pref+Pdata)/2.0)
    ##Differences in potential densities denotated as E
    Esvalues = refdens-datadens
    ###Find points before zero crossing
    Es = (np.diff(np.sign(Esvalues)))/2.0
    ##Add the points after the zero crossing
    shiftedEs = np.insert(Es,0,0,axis=1)
    Es = np.append(-Es,np.zeros([Es.shape[0],1]),1)
    startEs = Es
    Es = Es+shiftedEs

    if np.max(np.nansum(np.abs(Es),axis=1)/2) >2:
        print("Multiple Solutions")
    
    #instead of just bisecting lets take a weighted average
    Emask = np.zeros_like(Es)
    Pmask = np.zeros_like(Es)

    valid = np.full_like(Es,np.nan)
    valid[np.where(Es)] = Es[np.where(Es)]
    mask = np.nanmax(valid,axis=1)

    Es[np.where(Es == np.inf)] = 0
    Es[mask!=1] = 0
 
    psol = zeroCross(Es,Esvalues,startEs,Pref)
    ssol = zeroCross(Es,Esvalues,startEs,Sref)
    tsol = zeroCross(Es,Esvalues,startEs,Tref)
    gsol = zeroCross(Es,Esvalues,startEs,Gammaref)
    Es[np.where(Es == np.inf)] = 0
    #set to nan points without a solution or that have multiple zerocrossings seperated by over 100m
    psol[np.logical_or(psol==0,psol==np.inf)] = np.nan
    ssol[np.logical_or(ssol==0,ssol==np.inf)] = np.nan
    tsol[np.logical_or(tsol==0,tsol==np.inf)] = np.nan
    gsol[np.logical_or(gsol==0,gsol==np.inf)] = np.nan
    #print(psol)
    return ssol,tsol,psol,gsol

