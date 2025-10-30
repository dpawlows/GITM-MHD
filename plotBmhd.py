#!/usr/bin/env python

import sys
import numpy as np
from mhd import readmhd
from matplotlib import pyplot as pp

def usagemhd():
    print('Usage: plotBmhd.py filename altitude')
    print('   Plot the magnetic field given an ascii MHD file')
    exit(1)

try:
    args = sys.argv[2]
except:
    usagemhd()

try:
    sys.argv.index("-h")
    usagemhd()
except:
    pass

args = sys.argv
file = args[1]
pos = file.rfind('.dat')-13
outfile = 'Bup_'+file[pos:pos+13]+'.png'
print('Working on {}...'.format(outfile))
plotaltitude = int(args[2])

data = readmhd(file)
total = data["B"]

palt = min(range(len(data['Altitude'])), key=lambda i: abs(data['Altitude'][i]-plotaltitude))

totalmag = np.linalg.norm(total[:,:,palt,:],axis=2)
# indmag = np.linalg.norm(induced[:,:,palt,:],axis=2)
fig,ax = pp.subplots(nrows=2,ncols=1,figsize=(8,4))
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE =16

pp.rc('font', size=SMALL_SIZE)          # controls default text sizes
#pp.rc('axes', titlesize=20)     # fontsize of the axes title
#pp.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# ax[1].xaxis.label.set_size(MEDIUM_SIZE)
# ax[0].yaxis.label.set_size(MEDIUM_SIZE)
# ax[1].yaxis.label.set_size(MEDIUM_SIZE)
pp.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pp.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pp.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
pp.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# Plot only B up:
totalmag = total[:,:,palt,2]
# indmag = induced[:,:,palt,2]
# crustalmag = totalmag - indmag

cmap = 'bwr'
cbar_num_format = "%d"
# crustalmax = np.max(abs(crustalmag))
totalmax = np.max(abs(totalmag))
# crustallevels =  np.linspace(-crustalmax,crustalmax,90)
totallevels =  np.linspace(-totalmax,totalmax,90)

# lon1 = np.where(data['Longitude'] < 0, data['Longitude']+360,data['Longitude'])
lon1 = data['Longitude']
# imin = np.argmin(abs(lon1-0))
# lon2 = np.roll(lon1,imin)
data0 = np.transpose(data['bType'][:,:,palt])
#data1 = np.roll(data0,imin,1)
# breakpoint()
cont1=ax[0].pcolor(lon1,data['Latitude'],\
data0,cmap=cmap,shading='nearest')
cb1 = pp.colorbar(cont1,ax=ax[0])
# cbticks = np.linspace(-crustalmax,crustalmax,9)
# breakpoint()
#
# cb1.set_label('Crustal B$_{up}$ (nT)')
# cb1.set_ticks(cbticks)
# cb1.set_ticklabels(['{:.0f}'.format(x) for x in cbticks])
# ax[0].set_ylabel("Latitude")
# data2 = np.roll(np.transpose(totalmag),imin,1)
data2 = np.transpose(totalmag)
cont2=ax[1].contourf(lon1,data['Latitude'],\
    data2,levels=totallevels,cmap=cmap)
cbticks = np.linspace(-totalmax,totalmax,9)
cb2 = pp.colorbar(cont2,ax=ax[1])
cb2.set_ticks(cbticks)
cb2.set_ticklabels(['{:.0f}'.format(x) for x in np.linspace(-totalmax,totalmax,9)])

cb2.set_label('B$_{up}$ (nT)')

pp.ylabel("Latitude")
pp.xlabel("Longitude")

pp.savefig(outfile)
# pp.contour(lon,lat,total[:,])
