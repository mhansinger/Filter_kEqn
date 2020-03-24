# create symbolic links CH4 -> f

import os
import numpy
from os.path import join

# where are the time steps
t_dirs = '/home/hansinger/TNF_Lr75-80/case-30_80_FPV/constant/boundaryData/FUEL'

# get list of time steps
t_list = os.listdir(t_dirs)

t_list = [t for t in t_list if os.path.isdir(join(t_dirs,t))]

# go to t_dirs
os.chdir(t_dirs)

# loop over the time steps
for t in t_list:
    try:
        os.chdir(join(t_dirs,t))
        if os.path.exists('f'):
            print('f already exists')
        else:
            os.symlink('CH4', 'f')
        os.chdir('..')
    except FileExistsError:
        print('File exists already...ok')
        pass

print('Done!')
