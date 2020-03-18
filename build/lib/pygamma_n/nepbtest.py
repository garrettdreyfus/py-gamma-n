from gamma_n import *
import scipy.io as sio
from progress.bar import Bar
def nepbCTDExtract(fname):
    ctddata = sio.loadmat(fname)
    lats = np.asarray(ctddata["lat"])[0][:1]
    lons = np.asarray(ctddata["lon"])[0][:1]
    pres = np.asarray(ctddata["Pint"]).T[0]
    sals = np.asarray(ctddata["Sint"]).T[:1]
    temps = np.asarray(ctddata["Tint"]).T[:1]
    CT = np.asarray(ctddata["CT"]).T[:1]
    SA = np.asarray(ctddata["Sint_abs"]).T[:1]
    nspres = np.asarray(ctddata["P_gamma"]).T[:1]
    PV = np.asarray(ctddata["PV"]).T[:1]
    ns = np.asarray(ctddata["P_gref"])
    profiles = []
    gammas = []
    newpres = np.empty((len(sals),6000))
    newpres[:] = pres
    pres = newpres
    lons = np.asarray([lons]*6000).T
    lats = np.asarray([lats]*6000).T
    abssal = gsw.SA_from_SP(sals,pres,lons,lats)
    ct = gsw.CT_from_t(abssal,temps,pres)
    gamma = np.empty(sals.shape)
    print(gamma.shape)
    indexs = np.asarray(list(np.ndindex(gamma.shape)))
    print(indexs.shape)
    xindexs = indexs.T[0][500:2000]
    yindexs = indexs.T[1][500:2000]
    gamma,debug = gamma_n(abssal[xindexs][yindexs],ct[xindexs][yindexs],\
            pres[xindexs][yindexs],lons[xindexs][yindexs],lats[xindexs][yindexs])
    gamma[xindexs][yindexs] = gamma
    plt.plot(gamma,data["pres"])
    plt.show()
    print("boop")
    gammas.append(np.asarray(gamma))
    return np.asarray(gammas)

gammas = nepbCTDExtract("testing/data/newnepbdata.mat")



