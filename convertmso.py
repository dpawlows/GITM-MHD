#!/usr/bin/env python
import spiceypy as spice
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay
import numpy as np
import time
import sys
import re
from coordinates import *

#Hard coded resolution
latres = 3  # degrees
lonres = 3  # degrees
altres = 10  # km
print("Warning! Resolution is assumed to be {} lon x {} lat x {} alt".format(lonres,latres,altres))
time.sleep(1)


def nearest_index(value, start, step, count):
    """Return the index of the closest value on a regularly spaced grid."""
    idx = int(np.floor((value - start) / step + 0.5))
    return min(max(idx, 0), count - 1)


def fill_with_interpolation(values, mask, known_points, known_values, missing_coords, tri):
    """Fill missing grid cells in a single altitude slice."""
    if known_points.size == 0 or missing_coords.size == 0:
        return values

    filled = values.copy()

    linear_values = np.full(missing_coords.shape[0], np.nan)
    if tri is not None:
        linear_interp = LinearNDInterpolator(tri, known_values)
        linear_values = linear_interp(missing_coords)

    need_nearest = np.isnan(linear_values)
    if np.any(need_nearest):
        nearest_interp = NearestNDInterpolator(known_points, known_values)
        linear_values[need_nearest] = nearest_interp(missing_coords[need_nearest])

    filled[~mask] = linear_values
    return filled

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
longitude = np.arange(0, 360, lonres)
latitude = np.arange(-90, 91, latres)
altitude = np.arange(100, 301, altres)
grid_shape = (len(longitude), len(latitude), len(altitude))

bTotalEast = np.zeros(grid_shape)
bTotalNorth = np.zeros(grid_shape)
bTotalUp = np.zeros(grid_shape)
bType = np.zeros(grid_shape, int)
sample_counts = np.zeros(grid_shape, dtype=np.int32)

lon_grid, lat_grid = np.meshgrid(longitude, latitude, indexing='ij')


## Read the MHD file
endoffile = False
line = 0


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
        lon = np.arctan2(Xgcm[1], Xgcm[0])
        if lon < 0:
            lon += 2 * np.pi
        lat = np.arcsin(Xgcm[2] / np.linalg.norm(Xgcm))

        Bgcm = convertMSO2GCM(bx,by,bz,sslongitude,sslatitude,inclination,ls)
        totalFieldUp,totalFieldNorth,totalFieldEast = convertVector(Bgcm,lat,lon)
        totalFieldNorth = -totalFieldNorth

        if 100 <= alt <= 300:
            lon_deg = (lon * dpr) % 360.0
            lat_deg = lat * dpr
            ilon = int(np.floor((lon_deg + lonres / 2.0) / lonres)) % len(longitude)
            ilat = nearest_index(lat_deg, latitude[0], latres, len(latitude))
            ialt = nearest_index(alt, altitude[0], altres, len(altitude))

            sample_counts[ilon, ilat, ialt] += 1
            bTotalUp[ilon, ilat, ialt] += totalFieldUp
            bTotalEast[ilon, ilat, ialt] += totalFieldEast
            bTotalNorth[ilon, ilat, ialt] += totalFieldNorth
            bType[ilon, ilat, ialt] = type


        if line % 1000 == 0:
            print("Processing line {}".format(line))

        line += 1
f.close()

valid_mask = sample_counts > 0

with np.errstate(invalid='ignore'):
    bTotalEast[valid_mask] /= sample_counts[valid_mask]
    bTotalNorth[valid_mask] /= sample_counts[valid_mask]
    bTotalUp[valid_mask] /= sample_counts[valid_mask]

for ialt in range(len(altitude)):
    mask = valid_mask[:, :, ialt]
    missing = ~mask
    if not np.any(missing):
        continue

    known_points = np.column_stack((lon_grid[mask], lat_grid[mask]))
    if known_points.size == 0:
        continue

    missing_coords = np.column_stack((lon_grid[missing], lat_grid[missing]))
    tri = Delaunay(known_points) if known_points.shape[0] >= 3 else None

    bTotalEast[:, :, ialt] = fill_with_interpolation(
        bTotalEast[:, :, ialt],
        mask,
        known_points,
        bTotalEast[:, :, ialt][mask],
        missing_coords,
        tri,
    )
    bTotalNorth[:, :, ialt] = fill_with_interpolation(
        bTotalNorth[:, :, ialt],
        mask,
        known_points,
        bTotalNorth[:, :, ialt][mask],
        missing_coords,
        tri,
    )
    bTotalUp[:, :, ialt] = fill_with_interpolation(
        bTotalUp[:, :, ialt],
        mask,
        known_points,
        bTotalUp[:, :, ialt][mask],
        missing_coords,
        tri,
    )

    type_interp = NearestNDInterpolator(known_points, bType[:, :, ialt][mask])
    filled_types = bType[:, :, ialt].astype(float)
    filled_types[missing] = type_interp(missing_coords)
    bType[:, :, ialt] = np.rint(filled_types).astype(int)

magnitude = np.sqrt(bTotalEast**2 + bTotalNorth**2 + bTotalUp**2)
izero = int(np.count_nonzero(magnitude == 0))

f2.write("Lon(rad)  Lat(rad)  Alt(km) Bup  Bnorth  Beast  Btype\n")
f2.write("#START\n")
total_cells = np.prod(grid_shape)
for ilon in range(len(longitude)):
    for ilat in range(len(latitude)):
        for ialt in range(len(altitude)):
            if bType[ilon, ilat, ialt] == 0:
                bType[ilon, ilat, ialt] = 4

            if bTotalEast[ilon, ilat, ialt] == 0:
                bTotalEast[ilon, ilat, ialt] = 1e-3
            if bTotalNorth[ilon, ilat, ialt] == 0:
                bTotalNorth[ilon, ilat, ialt] = 1e-3
            if bTotalUp[ilon, ilat, ialt] == 0:
                bTotalUp[ilon, ilat, ialt] = 1e-3

            f2.write(
                "{:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:d}\n".format(
                    longitude[ilon] / dpr,
                    latitude[ilat] / dpr,
                    altitude[ialt],
                    bTotalEast[ilon, ilat, ialt],
                    bTotalNorth[ilon, ilat, ialt],
                    bTotalUp[ilon, ilat, ialt],
                    bType[ilon, ilat, ialt],
                )
            )

print("{} zeros out of {}".format(izero, total_cells))

f2.close()
