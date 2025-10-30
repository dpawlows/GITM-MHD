import numpy as np


def readmhd(file):

    f = open(file,'r')
    started = False

    while not started:
        temp = f.readline()
        if  temp[0:6] == "#START":
            started = True


    longitude = []
    latitude = []
    altitude = []

    longitude = np.arange(0,360,3)
    latitude = np.arange(-90,91,3)
    altitude = np.arange(100,301,10)
    data = {}
    B = np.zeros((len(longitude),len(latitude),len(altitude),3))
    bType = np.zeros((len(longitude),len(latitude),len(altitude)),dtype='int')
    for ilon in range(len(longitude)):
        for ilat in range(len(latitude)):
            for ialt in range(len(altitude)):
                temp = f.readline()
                temp1 = temp.split()

                lon = float(temp1[0])*180/np.pi
                # lon = lon - 360 if lon > 180 else lon
                lat = float(temp1[1])*180/np.pi
                alt = float(temp1[2])
                if np.round(lon) != longitude[ilon]:
                    print("Lon mismatch?")
                    breakpoint()
                if round(lat) != latitude[ilat]:
                    print("Lat mismatch?")
                    breakpoint()
                if round(alt) != altitude[ialt]:
                    breakpoint()
                    print("Alt mismatch?")
                B[ilon,ilat,ialt,:] = np.array([float(temp1[3]),float(temp1[4]),float(temp1[5])])
                bType[ilon,ilat,ialt] = int(temp1[6])

    data["Longitude"] = longitude
    data["Latitude"] = latitude
    data["Altitude"] = altitude
    data["B"] = B
    data["bType"] = bType

    return data
