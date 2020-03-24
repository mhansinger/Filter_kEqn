'''

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


def writePointsScalar(header_name='headerPoints.txt', out_path='BC_final/points', Value=None):
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


def writePointsVector(header_name='headerPoints.txt',out_path='BC_final/points', x=None,y=None,z=None):

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
        out_path = 'BC_final/'+str(names[id])
        try:
            if not os.path.exists(out_path):
                os.makedirs(out_path)
        except:
            print('Check the output directory for the new U fields')

        # write the new U points

        writePointsVector(header_name='headerU.txt',out_path=out_path+'/U',x=pfUx,y=pfUy,z=pfUz)
        # print progress
        print('Percent U: ', round((id+1)/idEnd,3)*100)


def Average_Scalar(scalar = 'CH4',y=None,z=None):
    # do grid interpolation and write the data to new U files
    names = os.listdir('BC_precursor/FUEL/')
    names = [f for f in names if os.path.isdir('BC_precursor/FUEL/')]
    names.sort()

    # set up the array for the averaging
    Scalar_Av = np.zeros(data_points_DNS)

    for id in range(idStart,idEnd):
        Scalar_path ='BC_precursor/FUEL/'+str(names[id])+'/'+scalar

        try:
            Scalar = readPointsScalar(filename=Scalar_path,header=1,data_points=data_points_DNS)

            # update the values
            Scalar_Av =+ Scalar
        except IndexError:
            print('Check the header lines in readPoints of writeScalar')

        # print progress
        print('Percent '+scalar, round((id+1)/idEnd,3)*100)

    # finally divide through number of time steps
    Scalar_Av = Scalar_Av/ idEnd

    Scalar_field_2D = griddata((y, z), Scalar_Av, (y, z), method='nearest')

    np.savetxt(scalar+'_time_Av_2D.txt',Scalar_field_2D)

    writePointsScalar(header_name='headerScalar.txt', out_path=scalar+'_time_Av.txt', Value=Scalar_Av)


########################
# MAIN PART!
#######################

# read in the DNS points
dns_x, dns_y, dns_z = readPointsVector('points_precursor', header=0, data_points=data_points_DNS,factor=1.0)

# write the new points file
writePointsVector(header_name='headerPoints.txt',out_path='BC_final/points', x=les_x,y=les_y,z=les_z)

# interpolate and write the new scalar fields
Average_Scalar(scalar='CH4',y=dns_y,z=dns_z)
writeScalar(data_points=data_points_LES,scalar='O2', pfy=les_y,pfz=les_z,y=dns_y,z=dns_z)
writeScalar(data_points=data_points_LES,scalar='N2', pfy=les_y,pfz=les_z,y=dns_y,z=dns_z)


