'''
This script is to compute sgs TKE in form of k from high resolved inflow data for coarser grids

k_sgs = 1/2 (<u'^2_x> + <u'^2_y> + <u'^2_z>)

author: mahansinger

'''

import numpy as np
import pandas as pd
import os, sys
from scipy.interpolate import griddata
from numba import jit, njit
#import matplotlib.pyplot as plt

# VARIABLES

data_points_org = 28980
data_points_new = 3300


idEnd = 36956 	#24956       #ANPASSEN an precursor data 20001steps=0.0005s
idStart = 0     #python starts with 0


#OK
def readPointsVector(filename='points_precursor',header=0,data_points=data_points_org):
    # read in the points data and returns an array with y and z coordinates, x is not needed as it is y-z plane
    data = pd.read_csv(filename, header=header , names=['raw'])
    data['raw'] = data['raw'].map(lambda x: x.strip('()'))

    try:
        assert len(data.raw) == data_points
    except AssertionError:
        data = data.drop(data.index[-1])

    # empty list
    x = []
    y = []
    z = []

    for row, val in data['raw'].iteritems():
        text = data['raw'].iloc[row]
        text = text.split()
        x.append(float(text[0]))
        y.append(float(text[1]))
        z.append(float(text[2]))

    # numpy points array
    #np_y = np.asarray(y)
    np_yz = np.asarray((y,z))
    return np_yz

def distance(dns_y,dns_z,les_y,les_z):
    # helper function
    return np.sqrt((dns_y-les_y)**2 + (dns_z - les_z)**2)


def matchCenters(DNS_points,LES_points,factor=9):
    # this routine finds the nearest DNS cell centers for each LES cell
    # factor determines how many DNS centers are considered for each LES cell center
    # e.g. factor = 9: assign the 9 closest DNS cell/face centers to an LES cell

    LES_match_list = []

    for l in range(0,LES_points.shape[1]):
        les_y = LES_points[0,l]
        les_z = LES_points[1,l]

        # generate a list to store the nearest points
        # filled with tuples: 1. entry is distance, 2. the point number of DNS
        nearest_points = [(1e5,i) for i in range(1,factor+1)]
        #print(nearest_points)

        for d in range(0,DNS_points.shape[1]):
            dns_y = DNS_points[0,d]
            dns_z = DNS_points[1,d]

            dist = distance(dns_y,dns_z,les_y,les_z)
            #print(dist)

            # now check if distance is closer then the farest point (last one in list of tuples)?
            # if so overwrite that one and sort it in descending order (nearest is first entry in list)
            if dist<nearest_points[-1][0]:
                nearest_points[-1] = (dist,d)
                nearest_points.sort(key = lambda tup : tup[0])

        LES_match_list.append(nearest_points)

        print('Done with LES point nr: ',l)

    return LES_match_list

