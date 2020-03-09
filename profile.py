import numpy as np
import gsw 
import datetime
import matplotlib.pyplot as plt
from scipy import interpolate


class Profile:
    def __init__(self,sals,temps,pres,gamma,lat,lon):
        ##id of profiles plus info
        tempunit = "insitu"
        salunit = "practical"
        if tempunit not in ["insitu","conservative","potential"]:
            raise ValueError("This temperature unit is not supported")
        if salunit not in ["practical","absolute","insitu"]:
            raise ValueError("This salinity unit is not supported")
        self.lat = lat
        self.lon = lon
        nanmask = ~np.isnan(temps) 
        nanmask = np.logical_and(~np.isnan(sals),nanmask)
        nanmask = np.logical_and(~np.isnan(pres),nanmask)
        nanmask = np.logical_and(~np.isnan(gamma),nanmask)

	#Temerature Salinity and Pressure
        self.temps = np.asarray(temps)[nanmask]
        self.sals = np.asarray(sals)[nanmask]
        self.pres = np.asarray(np.abs(pres))[nanmask]
        self.gamma = np.asarray(gamma)[nanmask]

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
        self.gamma = self.gamma[s]
        ##Interpolated Temperature, Salinity, and Pressure
        self.itemps = []
        self.isals = []
        self.ipres = []
        self.igamma = []
       	self.idensities = []
        self.neutraldepth = {}
        self.interpolate()

   ##apply a salinty offset 
    def interpolate(self):
        self.ipres = range(int(min(self.pres)),int(max(self.pres)))
        if len(self.pres)>4:

            self.isals = np.interp(self.ipres,self.pres,self.sals)
            self.itemps = np.interp(self.ipres,self.pres,self.temps)
            self.igamma = np.interp(self.ipres,self.pres,self.gamma)
            self.ialpha = gsw.alpha(self.isals,self.itemps,self.ipres)
            self.ibeta = gsw.beta(self.isals,self.itemps,self.ipres)
            self.idalphadtheta = gsw.cabbeling(self.isals,self.itemps,self.ipres)
            self.idalphadp = gsw.thermobaric(self.isals,self.itemps,self.ipres)

            ###using gsw
            self.n2 = gsw.Nsquared(self.isals,self.itemps,self.ipres,self.lat)[0]

    def neutralDepth(self,s,t,p,debug=False,searchrange=100,depthname=None):
        at = np.where(np.asarray(self.ipres) == p)[0][0]

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
            if abs(self.ipres[zero_crossings[0]] - self.ipres[zero_crossings[-1]])>100:
                return None
            a  =np.asarray(zero_crossings)
            #print("More than one crossing")
            return np.mean(self.isals[a]),np.mean(self.itemps[a]),np.mean(np.asarray(self.ipres)[a]),np.mean(self.igamma[a])
        else:
            return None
