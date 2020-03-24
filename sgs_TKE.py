'''
This script is to compute sgs TKE in form of k from high resolved inflow data for coarser grids

k_sgs = 1/2 (<u'_x^2> + <u'_y^2> + <u'_z^2>)

author: mahansinger

'''

import numpy as np
import pandas as pd
import os, sys
from scipy.interpolate import griddata
from numba import jit, njit
#import matplotlib.pyplot as plt
from numba import jit

###################################
# VARIABLES

data_points_DNS = 28980
data_points_LES = 3300


idEnd = 36956 	#24956       #ANPASSEN an precursor data 20001steps=0.0005s
idStart = 0     #python starts with 0
###################################

#OK
def readPointsVector(filename='points_precursor',header=0,data_points=data_points_DNS,factor=1.0):
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
    # numpy points array
    np_x = np.asarray(x) * factor
    np_y = np.asarray(y) * factor
    np_z = np.asarray(z) * factor
    return np_x, np_y, np_z


def readPointsScalar(filename='points_precursor',header=0,data_points=data_points_DNS):
    # read in the points_precursor data
    data = pd.read_csv(filename, header=header , names=['raw'])
    data['raw'] = data['raw'].map(lambda x: x.strip('()'))

    try:
        assert len(data.raw) == data_points
    except AssertionError:
        data = data.drop(data.index[-1])

    # empty list
    x = []

    for row, val in data['raw'].iteritems():
        text = data['raw'].iloc[row]
        text = text.split()
        x.append(float(text[0]))

    # numpy points array
    np_x = np.asarray(x)

    return np_x


def writePointsScalar(header_name='headerPoints.txt', out_path='BC_final/FUEL/points', Value=None):
    size_x = len(Value)
    # file to read
    fout = open(out_path, 'w')
    try:
        with open(header_name, encoding='utf-8') as f:
            for line in f.readlines():
                fout.write(line)
        #fout.write('\n')
        fout.write('\n')
        fout.write(str(int(size_x)) + '\n')
        fout.write('(\n')

        for i in range(0, len(Value)):
            #mystring = '( ' + str(Value[i]) + ')\n'
            mystring = str(Value[i]) + '\n'
            # print(mystring)
            fout.write(mystring)
        fout.write(')')
    finally:
        fout.close()


def writePointsVector(header_name='headerPoints.txt',out_path='BC_final/FUEL/points', x=None,y=None,z=None):

    size_x=len(x)
    # file to read
    fout = open(out_path,'w')
    try:
        with open(header_name, encoding='utf-8') as f:
            for line in f.readlines():
                fout.write(line)
        #fout.write('\n')
        fout.write('\n')
        fout.write(str(int(size_x))+'\n')
        fout.write('(\n')

        for i in range(0,len(x)):
            mystring = '( '+str(x[i])+' '+str(y[i])+' '+str(z[i])+' )\n'
            #print(mystring)
            fout.write(mystring)
        fout.write(')')

    finally:
        fout.close()


def writeU(data_points=data_points_LES,pfy=None,pfz=None,y=None,z=None):
    # do grid interpolation and write the data to new U files
    names = os.listdir('BC_precursor/FUEL/')
    names = [f for f in names if os.path.isdir('BC_precursor/FUEL/')]
    names.sort()

    for id in range(idStart,idEnd):
        U_path ='BC_precursor/FUEL/'+str(names[id])+'/U'

        try:
            Ux, Uy, Uz = readPointsVector(filename=U_path,header=1,data_points=data_points_DNS)
        except IndexError:
            print('Check the header lines in readPoints of writeU')

        # grid interpolation
        # nur z und y weil x eh x=0..
        pfUx = griddata((y, z), Ux, (pfy, pfz), method='nearest')
        pfUy = griddata((y, z), Uy, (pfy, pfz), method='nearest')
        pfUz = griddata((y, z), Uz, (pfy, pfz), method='nearest')

        # create directory if not exists
        out_path = 'BC_final/FUEL/'+str(names[id])
        try:
            if not os.path.exists(out_path):
                os.makedirs(out_path)
        except:
            print('Check the output directory for the new U fields')

        # write the new U points

        writePointsVector(header_name='headerU.txt',out_path=out_path+'/U',x=pfUx,y=pfUy,z=pfUz)
        # print progress
        print('Percent U: ', round((id+1)/idEnd,3)*100)


def writeScalar(data_points=data_points_LES,scalar = 'CH4',pfy=None,pfz=None,y=None,z=None):
    # do grid interpolation and write the data to new U files
    names = os.listdir('BC_precursor/FUEL/')
    names = [f for f in names if os.path.isdir('BC_precursor/FUEL/')]
    names.sort()

    for id in range(idStart,idEnd):
        Scalar_path ='BC_precursor/FUEL/'+str(names[id])+'/'+scalar

        try:
            Scalar = readPointsScalar(filename=Scalar_path,header=1,data_points=data_points_DNS)
            # grid interpolation
            # nur z und y weil x eh x=0..
            pfScalar = griddata((y, z), Scalar, (pfy, pfz), method='nearest')
        except IndexError:
            print('Check the header lines in readPoints of writeScalar')

        # create directory if not exists
        out_path = 'BC_final/FUEL/'+str(names[id])
        try:
            if not os.path.exists(out_path):
                os.makedirs(out_path)
        except:
            print('Check the output directory for the new U fields')

        # write the new U points

        writePointsScalar(header_name='headerScalar.txt',out_path=out_path+'/'+scalar,Value=pfScalar)
        # print progress
        print('Percent '+scalar, round((id+1)/idEnd,3)*100)


@jit
def distance(dns_y,dns_z,les_y,les_z):
    # helper function
    return np.sqrt((dns_y-les_y)**2 + (dns_z - les_z)**2)


def matchCenters(DNS_points,LES_points,cells=9):
    # this routine finds the nearest DNS cell centers for each LES cell
    # factor determines how many DNS centers are considered for each LES cell center
    # e.g. factor = 9: assign the 9 closest DNS cell/face centers to an LES cell

    LES_match_list = []

    for l in range(0,LES_points.shape[1]):
        les_y = LES_points[0,l]
        les_z = LES_points[1,l]

        # generate a list to store the nearest points
        # filled with tuples: 1. entry is distance, 2. the point number of DNS
        nearest_points = [(1e5,i) for i in range(1,cells+1)]
        #print(nearest_points)

        for d in range(0,DNS_points.shape[1]):
            dns_y = DNS_points[0,d]
            dns_z = DNS_points[1,d]

            dist = distance(dns_y,dns_z,les_y,les_z)
            #print(dist)

            # now check if distance is closer then the most far point (last one in list of tuples)?
            # if so overwrite that one and sort it in descending order (nearest is first entry in list)
            if dist < nearest_points[-1][0]:
                nearest_points[-1] = (dist,d)
                nearest_points.sort(key = lambda tup : tup[0])

        LES_match_list.append(nearest_points)

        print('Done with LES point nr: ',l)

    return LES_match_list


@jit
def computeTKE_LES(LES_match_list,U_DNS_list,max_TKE,cells=9,nrLES = data_points_LES):

    TKE_list = np.zeros(nrLES)

    for i, pair in enumerate(LES_match_list):
        DNS_points = np.zeros(cells)
        U_points = np.zeros(cells)
        V_points = np.zeros(cells)
        W_points = np.zeros(cells)

        for p in range(0,cells):
            dns_point = pair[p][1]
            DNS_points[p]= pair[p][1]
            # CHECKEN!
            U_points[p] = U_DNS_list[dns_point,0]
            V_points[p] = U_DNS_list[dns_point,1]
            W_points[p] = U_DNS_list[dns_point,2]

        uprime = vprime = wprime = 0

        for p in range(0,cells):
            uprime += (U_points[p] - U_points.mean()) ** 2
            vprime += (V_points[p] - V_points.mean()) ** 2
            wprime += (W_points[p] - W_points.mean()) ** 2

        TKE = (uprime/cells + vprime/cells + wprime/cells) / 2

        if TKE > max_TKE:
            TKE = max_TKE

        TKE_list[i] = TKE

    return TKE_list



def writeTKE(LES_match_list,max_TKE=1):
    # write the TKE for each time step
    names = os.listdir('BC_precursor/FUEL/')
    names = [f for f in names if os.path.isdir('BC_precursor/FUEL/')]
    names.sort()

    this_maxTKE = max_TKE

    for id in range(idStart,idEnd):
        U_path ='BC_precursor/FUEL/'+str(names[id])+'/U'

        try:
            Ux, Uy, Uz = readPointsVector(filename=U_path,header=1,data_points=data_points_DNS)
            # compute TKE

            U_vector =np.array([Ux,Uy,Uz]).T
            TKE_list = computeTKE_LES(LES_match_list=LES_match_list, U_DNS_list=U_vector, max_TKE=this_maxTKE, cells=9,nrLES=data_points_LES)
        except IndexError:
            print('Check the header lines in readPoints of writeU')

        # create directory if not exists
        out_path = 'BC_final/FUEL/'+str(names[id])
        try:
            if not os.path.exists(out_path):
                os.makedirs(out_path)
        except:
            print('Check the output directory for the new U fields')

        # write the new U points

        writePointsScalar(header_name='headerScalar.txt',out_path=out_path+'/k',Value = TKE_list)
        # print progress
        print('Percent k', round((id+1)/idEnd,3)*100)


########################
# MAIN PART!
#######################

# read in the DNS points
dns_x, dns_y, dns_z = readPointsVector('points_precursor', header=0, data_points=data_points_DNS,factor=1.0)

# read in the LES points
les_x,les_y, les_z = readPointsVector('points_LES', header=0, data_points=data_points_LES)

# write the new points file
writePointsVector(header_name='headerPoints.txt',out_path='BC_final/FUEL/points', x=les_x,y=les_y,z=les_z)

# interpolate and write the new U field
writeU(data_points=data_points_LES,pfy=les_y,pfz=les_z,y=dns_y,z=dns_z)

# interpolate and write the new scalar fields
writeScalar(data_points=data_points_LES,scalar='CH4', pfy=les_y,pfz=les_z,y=dns_y,z=dns_z)
writeScalar(data_points=data_points_LES,scalar='O2', pfy=les_y,pfz=les_z,y=dns_y,z=dns_z)
writeScalar(data_points=data_points_LES,scalar='N2', pfy=les_y,pfz=les_z,y=dns_y,z=dns_z)

#LES_points = readPointsVector(filename='points_LES',header=0,data_points=data_points_LES)
#DNS_points = readPointsVector(filename='points_final',header=0,data_points=data_points_DNS,factor=0.001)

#
LES_yz = np.array([les_y,les_z])
DNS_yz  =np.array([dns_y,dns_z])
LES_match_list = matchCenters(DNS_points=DNS_yz ,LES_points=LES_yz,cells=int(data_points_DNS/data_points_LES))
writeTKE(LES_match_list,max_TKE=100)

