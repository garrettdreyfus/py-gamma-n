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
    dims = (len(p),len(p_ref))
    Sref = np.empty(dims)
    Sref[:] = s_ref
    Tref = np.empty(dims)
    Tref[:] = t_ref
    Pref = np.empty(dims)
    Pref[:] = p_ref
    Gammaref = np.empty(dims)
    Gammaref[:] = gamma_ref

    Sdata = np.column_stack([s]*len(p_ref))
    Tdata = np.column_stack([t]*len(p_ref))
    Pdata = np.column_stack([p]*len(p_ref))

    #print(Sref)
    #print(Sdata)
    refdens = gsw.rho(Sref,Tref,(Pref+Pdata)/2.0)
    datadens = gsw.rho(Sdata,Tdata,(Pref+Pdata)/2.0)
    Es = refdens-datadens
    cross = np.where(np.diff(np.sign(Es)))
    psol = Pref[cross]
    ssol = Sref[cross]
    tsol = Tref[cross]
    gsol = Gammaref[cross]
    return ssol,tsol,psol,gsol

