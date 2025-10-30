#!/usr/bin/env python
from glob import glob
from matplotlib import pyplot
import spiceypy as spice
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import pickle
import time
import sys
import re
from coordinates import *

#Hard coded resolution
latres = 3 #degrees
lonres = 3 #degrees
altres = 10 #km
print("Warning! Resolution is assumed to be {} lon x {} lat x {} alt".format(lonres,latres,altres))
time.sleep(1)

dpr = spice.dpr()

file = sys.argv[1]
#Load Necessary spice Kernels
spice.furnsh('mhdMetaK.txt')

f = open(file,'r')

showkernels = False
if showkernels:
    n = spice.ktotal("ALL")
    print("Kernels loaded:")
    for i in range(n):
        file, filtyp, srcfil, handle = spice.kdata(i, "ALL", 512, 32, 512)
        print(f"{i:02d} {filtyp:>4}  {file}")

started = False
while not started:
    temp = f.readline()
    if "timestamp" in temp:
        ctime = re.search(r'(\d{4}/\d{2}/\d{2})\s+(\d{2}:\d{2}:\d{2}\.\d{3})', temp).group(0)
        etime = spice.str2et(ctime)
        print(f"Converting MHD file {file} date: {spice.et2utc(etime,'ISOC',3)}")

    if "VARIABLES" in temp:
        started = True

#Get relevant information from kernels based on et
#Open file for writing
file2 = 'MHDGEO_'+ctime[0:4]+ctime[5:7]+ctime[8:10]+'_'+ctime[11:13]+ctime[14:16]\
    +ctime[17:19]+'.dat'

f2 = open(file2,'w')

radius = spice.bodvrd("MARS","RADII",3)[1][0]
method = "Intercept/ellipsoid"
#Get subsolar point for current et
output = spice.subslr(method,"mars",etime,"iau_mars","None","mars")
spoint = output[0]
sppc = np.array(spice.reclat(spoint)) #r, lon, lat
ls = spice.lspcn('MARS',etime,'None')
sslongitude = sppc[1]
sslatitude = sppc[2]
# tilt =  sppc[2]
# sslat = sppc[2]
# sslong =  sppc[1]
inclination = 25.91 / dpr
longitude = np.arange(0,360,lonres)
latitude = np.arange(-90,91,latres)
altitude = np.arange(100,301,altres)
bTotalEast = np.zeros((len(longitude),len(latitude),len(altitude)))
bTotalNorth = np.zeros((len(longitude),len(latitude),len(altitude)))
bTotalUp= np.zeros((len(longitude),len(latitude),len(altitude)))
# bMagnitude = np.zeros((len(longitude),len(latitude),len(altitude)))
# bElevation = np.zeros((len(longitude),len(latitude),len(altitude)))
bType = np.zeros((len(longitude),len(latitude),len(altitude)),int)


## Read the MHD file
endoffile = False
line = 0
maxlon = -1e9
minlon = 1e9
maxlong = -10
maxlat = -1e9
minlat = 1e9
bmax = 0
B = []
loc = []
msox = []
msoy = []
msoz = []
bmag = []
btyp = []

while not endoffile:
    temp1 = f.readline()
    if len(temp1) == 0:
        endoffile = True
    else:
        temp2 = [float(t) for t in temp1.split()]
        z = temp2[2]
        x = temp2[0]
        y = temp2[1]
        alt = temp2[3]
        type = int(temp2[6])
        bx = temp2[7]
        by = temp2[8]
        bz = temp2[9]

        # msox.append(bx)
        # msoy.append(by)
        # msoz.append(bz)


        Xgcm = convertMSO2GCM(x,y,z,sslongitude,sslatitude,inclination,ls)
        #
        # breakpoint()
        alts = np.linalg.norm(Xgcm)
        lon = np.arctan2(Xgcm[1],Xgcm[0])#-sslong
        if lon < 0:
            lon += 2*np.pi
        lat = np.arcsin(Xgcm[2]/alts)
        alts = (alt*radius)-radius

        Bgcm = convertMSO2GCM(bx,by,bz,sslongitude,sslatitude,inclination,ls)
        totalFieldUp,totalFieldNorth,totalFieldEast = convertVector(Bgcm,lat,lon)
        totalFieldNorth = -totalFieldNorth

        if alt > 100 and alt <= 300:
            ilon = (np.abs(longitude-lon*dpr)).argmin()
            ilat = (np.abs(latitude-lat*dpr)).argmin()
            ialt = (np.abs(altitude-alt)).argmin()
            if lon > 358.5/dpr and longitude[0] == 0:
                ilon = 0
                lon = 0
            if (np.abs(longitude[ilon]-lon*dpr) > 2 or np.abs(latitude[ilat]-lat*dpr >2) or \
                np.abs(altitude[ialt]-alt) > 5):
                print("issue with grid?")
                breakpoint()
            bTotalUp[ilon,ilat,ialt] = totalFieldUp
            bTotalEast[ilon,ilat,ialt] = totalFieldEast
            bTotalNorth[ilon,ilat,ialt] = totalFieldNorth
            bType[ilon,ilat,ialt] = type


        if line % 1000 == 0:
            print("Processing line {}".format(line))

        line += 1
f.close()

f2.write("Lon(rad)  Lat(rad)  Alt(km) Bup  Bnorth  Beast  Btype\n")
f2.write("#START\n")
i = 0
izero = 0
for ilon in range(len(longitude)):
    for ilat in range(len(latitude)):
        for ialt in range(len(altitude)):


            if ilon > 0 and ilat > 0 and ilon < len(longitude)-1 and ilat < len(latitude)-1:
                bMagnitude =np.linalg.norm([bTotalUp[ilon,ilat,ialt],\
                    bTotalEast[ilon,ilat,ialt],bTotalNorth[ilon,ilat,ialt]])

                if bMagnitude == 0:
                    izero += 1

                    # East
                    x = [longitude[ilon-1], longitude[ilon+1]]
                    y = [latitude[ilat-1], latitude[ilat+1]]
                    Z = [
                        [bTotalEast[ilon-1, ilat-1, ialt], bTotalEast[ilon-1, ilat+1, ialt]],
                        [bTotalEast[ilon+1, ilat-1, ialt], bTotalEast[ilon+1, ilat+1, ialt]]
                    ]

                    f = RegularGridInterpolator((x, y), Z, method='linear', bounds_error=False, fill_value=None)
                    bTotalEast[ilon,ilat,ialt] = float(f((longitude[ilon], latitude[ilat])))

                    # North
                    Z = [
                        [bTotalNorth[ilon-1, ilat-1, ialt], bTotalNorth[ilon-1, ilat+1, ialt]],
                        [bTotalNorth[ilon+1, ilat-1, ialt], bTotalNorth[ilon+1, ilat+1, ialt]]
                    ]

                    f = RegularGridInterpolator((x, y), Z, method='linear', bounds_error=False, fill_value=None)
                    bTotalNorth[ilon,ilat,ialt] = float(f((longitude[ilon], latitude[ilat])))

                    # Up
                    Z = [
                        [bTotalUp[ilon-1, ilat-1, ialt], bTotalUp[ilon-1, ilat+1, ialt]],
                        [bTotalUp[ilon+1, ilat-1, ialt], bTotalUp[ilon+1, ilat+1, ialt]]
                    ]

                    f = RegularGridInterpolator((x, y), Z, method='linear', bounds_error=False, fill_value=None)
                    bTotalUp[ilon,ilat,ialt] = float(f((longitude[ilon], latitude[ilat])))


                    if abs(longitude[ilon] - longitude[ilon-1]) < abs(longitude[ilon+1] - longitude[ilon]):
                        thislon = ilon-1
                    else:
                        thislon = ilon+1
                    if abs(latitude[ilat] - latitude[ilat-1]) < abs(latitude[ilat+1] - latitude[ilat]):
                        thislat = ilat-1
                    else:
                        thislat = ilat+1
                    bType[ilon,ilat,ialt] = bType[thislon,thislat,ialt]

            if bType[ilon,ilat,ialt] == 0:
                bType[ilon,ilat,ialt] = 4
            bTotalEast[ilon,ilat,ialt] = 1e-3 if  bTotalEast[ilon,ilat,ialt] == 0 else  bTotalEast[ilon,ilat,ialt]
            bTotalNorth[ilon,ilat,ialt] = 1e-3 if bTotalNorth[ilon,ilat,ialt] == 0 else bTotalNorth[ilon,ilat,ialt]
            bTotalUp[ilon,ilat,ialt] = 1e-3 if    bTotalUp[ilon,ilat,ialt] == 0 else    bTotalUp[ilon,ilat,ialt]

            f2.write("{:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:d}\n".format(\
            longitude[ilon]/dpr,latitude[ilat]/dpr,\
            altitude[ialt],bTotalEast[ilon,ilat,ialt],bTotalNorth[ilon,ilat,ialt],\
            bTotalUp[ilon,ilat,ialt],\
            bType[ilon,ilat,ialt]))
            i+=1
print("{} zeros out of {}".format(izero,i))

f2.close()
