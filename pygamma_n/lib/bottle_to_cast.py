import numpy as np
import gsw

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
    Es = (np.diff(np.sign(Esvalues)))
    ##Add the points after the zero crossing
    shiftedEs = np.insert(Es,0,0,axis=1)
    Es = np.append(Es,np.zeros([Es.shape[0],1]),1)
    Es = Es+shiftedEs

    if np.max(np.nansum(Es,axis=1)/2) >2:
        print("Multiple Solutions")
    
    #instead of just bisecting lets take a weighted average
    z = np.zeros_like(Es)
    valid = np.full_like(Es,np.nan)
    ## fill the indexs with their z values
    valid[np.where(Es)] = np.abs(Pref)[np.where(Es)]
    z[np.where(Es)] = np.abs(Esvalues)[np.where(Es)]
    Es = np.reciprocal(z)
    mask = np.ptp(valid,axis=1)
    Es[np.where(Es == np.inf)] = 0

 
    ##Make Sure each row sums to one
    Es = Es/np.nansum(Es,axis=1)[:,None]
    #This multiplies but only cares for the diagonals
    psol = np.nansum((Es * Pref),axis=-1)
    ssol = np.nansum((Es * Sref),axis=-1)
    tsol = np.nansum((Es * Tref),axis=-1)
    gsol = np.nansum((Es * Gammaref),axis=-1)
    #set to nan points without a solution or that have multiple zerocrossings seperated by over 100m
    psol[np.where(mask>100)] = np.nan
    psol[psol==0] = np.nan
    ssol[np.where(mask>100)] = np.nan
    ssol[ssol==0] = np.nan
    tsol[np.where(mask>100)] = np.nan
    tsol[tsol==0] = np.nan
    gsol[np.where(mask>100)] = np.nan
    gsol[gsol==0] = np.nan
    return ssol,tsol,psol,gsol

