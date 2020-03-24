# create symbolic links CH4 -> f

import os
import numpy
from os.path import join
import pandas as pd

# where are the time steps
t_dirs = './' #'/home/hansinger/TNF_Lr75-80/case-30_80_FPV/constant/boundaryData/FUEL'

# get list of time steps
t_list = os.listdir(t_dirs)

#t_list = [t for t in t_list if os.path.isdir(join(t_dirs,t))]

# go to t_dirs
# os.chdir(t_dirs)

length=1280

# loop over the time steps

for t in t_list:
    try:
        os.chdir(t)
        if os.path.exists('O2'):
            print('O2 exists already')
        else:
            os.symlink('f', 'CH4')
            this_f = pd.read_csv('f')
            # this_f.columns = ['data']
            this_O2 = this_f.copy()
            this_N2 = this_f.copy()

            f = this_f.iloc[10:-1].values.astype(float)
            air = 1 - f

            O2 = air * 0.209
            N2 = air - O2

            this_O2.iloc[10:-1] = O2
            this_N2.iloc[10:-1] = N2

            this_O2.to_csv('O2', index=False)
            this_N2.to_csv('N2', index=False)

        os.chdir('../')

        print('Done with %s ' % t)
    except:
        print('Problem occured')
        #pass

print('All Done!')
