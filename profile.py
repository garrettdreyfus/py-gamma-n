import numpy as np
import operator
import gsw 
from mpl_toolkits.mplot3d import Axes3D
import datetime
import matplotlib.pyplot as plt
from scipy import interpolate
import mygsw
from scipy.interpolate import UnivariateSpline
#import gswmatlab.pyinterface as matgsw


class Profile:
    def __init__(self,eyed, data,tempunit,salunit):
        ##id of profiles plus info
        if tempunit not in ["insitu","conservative","potential"]:
            raise ValueError("This temperature unit is not supported")
        if salunit not in ["practical","absolute","insitu"]:
            raise ValueError("This salinity unit is not supported")
        if not {"sal","temp","pres","lat","lon"}.issubset(data.keys()):
            raise ValueError("This does not contain the required information")
        if abs(max(data["pres"])-min(data["pres"])) <50:
            print(data["pres"])
            raise ValueError("This does not contain enough pressure information ")

        self.eyed = eyed
        self.lat = data["lat"]
        self.lon = data["lon"]
	#Temerature Salinity and Pressure
        self.temps = np.asarray(data["temp"])
        self.sals = np.asarray(data["sal"])
        self.pres = np.abs(np.asarray(data["pres"]))

        if salunit == "practical":
            self.sals = gsw.SA_from_SP(self.sals,self.pres,self.lon,self.lat)

        if tempunit == "potential":
            self.temps = gsw.CT_from_pt(self.sals,self.temps)

        if tempunit == "insitu":
            self.temps = gsw.CT_from_t(self.sals,self.temps,np.abs(self.pres))

        s = np.argsort(self.pres)
        self.temps = self.temps[s]
        self.sals = self.sals[s]
        self.pres = self.pres[s]
        ##Interpolated Temperature, Salinity, and Pressure
        self.itemps = []
        self.isals = []
        self.ipres = []
       	self.idensities = []
        self.neutraldepth = {}
        self.interpolate()

   ##apply a salinty offset 
    def interpolate(self):
        self.ipres = range(int(min(self.pres)),int(max(self.pres)))
        if len(self.pres)>4:

            self.isals = np.interp(self.ipres,self.pres,self.sals)
            self.itemps = np.interp(self.ipres,self.pres,self.temps)
            self.ialpha = gsw.alpha(self.isals,self.itemps,self.ipres)
            self.ibeta = gsw.beta(self.isals,self.itemps,self.ipres)
            self.idalphadtheta = gsw.cabbeling(self.isals,self.itemps,self.ipres)
            self.idalphadp = gsw.thermobaric(self.isals,self.itemps,self.ipres)

            ###using gsw
            self.n2 = gsw.Nsquared(self.isals,self.itemps,self.ipres,self.lat)[0]

    def neutralDepth(self,s,t,p,debug=False,searchrange=100,depthname=None):
        depth = int(depth)
        plowerbound = min(self.ipres[0])
        pupperbound = min(self.ipres[-1])
        at = np.where(np.asarray(self.ipres) == depth)[0][0]

        depths =  (np.asarray(self.ipres[:]) +p)/2.0
        
        selfdensities = gsw.rho(self.isals[:],\
                self.itemps[:],\
                depths)

        refdensities = gsw.rho([s]*(len(depths)),\
                [t]*(len(depths)),\
                depths)

        Es = selfdensities-refdensities

        zero_crossings = np.where(np.diff(np.sign(Es)))[0]
        smallest = np.argmin(np.abs(Es))
        if len(zero_crossings)>=1 :
            if abs(p2.ipres[zero_crossings[0]] - p2.ipres[zero_crossings[-1]])>100:
                return None
            a  =np.asarray(zero_crossings)
            p2.neutraldepth[depthname] = np.mean(np.asarray(self.ipres)[a])
            #print("More than one crossing")
            return np.mean(self.itemps[a]),np.mean(self.isals[a]),self.neutraldepth[depthname]
        else:
            return None
