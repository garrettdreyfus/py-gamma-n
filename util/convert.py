import gsw
import h5py
import xarray as xr
import numpy as np

def convertToDictionary(refdata):
    ## gamma_n goes lon, lat, pres
    data = {}
    data["lats"] = np.asarray(refdata["lat_ref"])
    data["lons"] = np.asarray(refdata["long_ref"])
    p_ref = np.empty((91,45,33))
    for i in range(p_ref.shape[0]):
        for j in range(p_ref.shape[1]):
            p_ref[i][j] = np.asarray(refdata["p_ref"]).flatten()
    data["s"] = []
    data["t"] = []
    data["p"] = []
    for j in range(refdata["SP_ref"].shape[0]):
        data["s"].append(np.asarray(gsw.SA_from_SP(np.asarray(refdata["SP_ref"][j]),p_ref[j],np.full_like(data["lats"],data["lons"][j]),data["lats"])))
        data["t"].append(np.asarray(gsw.CT_from_t(data["s"][-1],np.asarray(refdata["t_ref"][j]),p_ref[j])))
        data["p"].append(p_ref[j])
    data["t"] =np.asarray(data["t"])
    print(data["t"].shape)
    data["p"] =np.asarray(refdata["p_ref"]).flatten()
    data["lats"] = np.asarray(refdata["lat_ref"]).flatten()
    data["lons"] = np.asarray(refdata["long_ref"]).flatten()
    data["gamma"] = np.asarray(refdata["gamma_n_ref"])
    return data
refdata = h5py.File('refdata.mat', 'r')
d = convertToDictionary(refdata)    
print(d["t"].shape)
ds = xr.DataSet(data_vars={"t":(["x","y","z"],d["t"]), 
                           "s":(["x","y","z"],d["s"]),
                           "gamma":(["x","y","z"],d["gamma"])}, 
                coords={"lat": (["y"], d["lats"]), 
                        "lon": (["x"], d["lons"]), 
                        "pres": (["z"],d["p"])})

ds.to_netcdf('refprofiles.nc')
