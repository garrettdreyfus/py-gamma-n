from gamma_n import *
import scipy.io as sio
from progress.bar import Bar
def nepbCTDExtract(fname):
    ctddata = sio.loadmat(fname)
    lats = np.asarray(ctddata["lat"])[0][:10]
    lons = np.asarray(ctddata["lon"])[0][:10]
    pres = np.asarray(ctddata["Pint"]).T[0]
    sals = np.asarray(ctddata["Sint"]).T[:10]
    temps = np.asarray(ctddata["Tint"]).T[:10]
    CT = np.asarray(ctddata["CT"]).T[:10]
    SA = np.asarray(ctddata["Sint_abs"]).T[:10]
    nspres = np.asarray(ctddata["P_gamma"]).T[:10]
    PV = np.asarray(ctddata["PV"]).T[:10]
    ns = np.asarray(ctddata["P_gref"])
    profiles = []
    gammas = []
    print(pres.shape)
    print(sals.shape)
    newpres = np.empty((len(sals),6000))
    newpres[:] = pres
    pres = newpres
    lons = np.asarray([lons]*6000).T
    lats = np.asarray([lats]*6000).T
    abssal = gsw.SA_from_SP(sals,pres,lons,lats)
    ct = gsw.CT_from_t(abssal,temps,pres)
    gamma = np.empty(sals.shape)
    indexs = np.asarray(list(np.ndindex(gamma.shape)))
    print(sals[indexs])
    print(indexs.shape)
    print(pres.shape)
    print(lats.shape)
    print(lons.shape)
    print(indexs)
    gamma,debug = gamma_n(abssal[indexs],ct[indexs],pres[indexs],lons[indexs],lats[indexs])
    gamma[indexs] = gamma
    plt.plot(gamma,data["pres"])
    plt.show()
    print("boop")
    gammas.append(np.asarray(gamma))
    return np.asarray(gammas)

gammas = nepbCTDExtract("testing/data/newnepbdata.mat")



