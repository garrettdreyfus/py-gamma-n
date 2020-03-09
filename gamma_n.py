import h5py
from functools import partial
import numpy as np
import matplotlib.pyplot as plt

refdata = h5py.File('refdata.mat', 'r')

def gamma_n(refdata,s,t,p,lon,lat):
    print(len(refdata["lat_ref"][:]))
    print(len(refdata["long_ref"][:]))
    for l in range(len(refdata["long_ref"][:])):
        for j in range(len(refdata["lat_ref"][:])):
            plt.plot(np.asarray(refdata["gamma_n_ref"])[l][j][:])
    plt.show()

gamma_n = partial(gamma_n,refdata)
gamma_n(0,0,0,0,0)
