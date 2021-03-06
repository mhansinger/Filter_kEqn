# for python3

import numpy as np
import pandas as pd
import os, sys
from scipy.interpolate import griddata
#import matplotlib.pyplot as plt

# VARIABLES

data_points_org = 9681
data_points_new = 2241


idEnd = 10530   #24956       #ANPASSEN an precursor data 20001steps=0.0005s
idStart = 0     #python starts with 0


def readPointsVector(filename='points_precursor',header=1,data_points=data_points_org):
    # read in the points_precursor data
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
    np_x = np.asarray(x)
    np_y = np.asarray(y)
    np_z = np.asarray(z)
    return np_x, np_y, np_z


def readPointsScalar(filename='points_precursor',header=1,data_points=data_points_org):
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


def writeU(data_points=data_points_new,pfy=None,pfz=None,y=None,z=None):
    # do grid interpolation and write the data to new U files
    names = os.listdir('BC_precursor/FUEL/')
    names = [f for f in names if os.path.isdir('BC_precursor/FUEL/')]
    names.sort()

    for id in range(idStart,idEnd):
        U_path ='BC_precursor/FUEL/'+str(names[id])+'/U'

        try:
            Ux, Uy, Uz = readPointsVector(filename=U_path, header=1,data_points=data_points_org)
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


def writeScalar(data_points=data_points_new,scalar = 'CH4',pfy=None,pfz=None,y=None,z=None):
    # do grid interpolation and write the data to new U files
    names = os.listdir('BC_precursor/FUEL/')
    names = [f for f in names if os.path.isdir('BC_precursor/FUEL/')]
    names.sort()

    for id in range(idStart,idEnd):
        Scalar_path ='BC_precursor/FUEL/'+str(names[id])+'/'+scalar

        try:
            Scalar = readPointsScalar(filename=Scalar_path,header=2,data_points=data_points_org)
        except IndexError:
            print('Check the header lines in readPoints of writeScalar')

        # grid interpolation
        # nur z und y weil x eh x=0..
        pfScalar = griddata((y,z), Scalar, (pfy, pfz), method='nearest')

        # create directory if not exists
        out_path = 'BC_final/'+str(names[id])
        try:
            if not os.path.exists(out_path):
                os.makedirs(out_path)
        except:
            print('Check the output directory for the new U fields')

        # write the new U points

        writePointsScalar(header_name='headerScalar.txt',out_path=out_path+'/'+scalar,Value=pfScalar)
        # print progress
        print('Percent '+scalar, round((id+1)/idEnd,3)*100)


########################
# MAIN PART!
#######################

# read in the original points
x, y, z = readPointsVector('points_precursor', header=20, data_points=data_points_org)

# read in the patch field
pfx,pfy, pfz = readPointsVector('points_final', header=20, data_points=data_points_new)

# write the new points file
writePointsVector(header_name='headerPoints.txt',out_path='BC_final/points', x=pfx,y=pfy,z=pfz)

# interpolate and write the new U field
writeU(data_points=data_points_new,pfy=pfy,pfz=pfz,y=y,z=z)

# interpolate and write the new scalar fields
#writeScalar(data_points=data_points_new,scalar='CH4', pfy=pfy,pfz=pfz,y=y,z=z)
#writeScalar(data_points=data_points_new,scalar='O2', pfy=pfy,pfz=pfz,y=y,z=z)
#writeScalar(data_points=data_points_new,scalar='N2', pfy=pfy,pfz=pfz,y=y,z=z)
