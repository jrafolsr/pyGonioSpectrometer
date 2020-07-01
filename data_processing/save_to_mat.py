# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:24:17 2020

@author: JOANRR
"""
import numpy as np
from scipy.io import loadmat, savemat
from os.path import join as pjoin

# This is manual, so far, needs to be a function
# fcal = r'.\calibration_files'  
# EL_SIM2 = loadmat(pjoin(fcal, 'ErrLands-320nm-0.01di.mat'))
# EL_SIM2['qnormSpecRadInt_sim_2D'] = EL_SIM2.pop('ans')
# EL_SIM2['ipos_sim'] = np.arange(0.01, 1.0, 0.01)
# EL_SIM2['dAL_sim'] = np.array([380])
# EL_SIM2['wl_q'] = np.arange(380,781, 1)
# EL_SIM2['angle_q'] = np.arange(0,90, 10)


# EL_SIM2['ipos_sim'] = EL_SIM2['ipos_sim'].reshape((1, len(EL_SIM2['ipos_sim'])))
# EL_SIM2['dAL_sim'] = EL_SIM2['dAL_sim'].reshape((1, len(EL_SIM2['dAL_sim'])))
# EL_SIM2['wl_q'] = EL_SIM2['wl_q'].reshape((1, len(EL_SIM2['wl_q'])))
# EL_SIM2['angle_q'] = EL_SIM2['angle_q'].reshape((1, len(EL_SIM2['angle_q'])))

# savemat(pjoin(fcal, 'ErrLands-320nm-0.01di_restructured.mat'), EL_SIM2,
#         oned_as = 'row', do_compression=True)

fcal = r'.\calibration_files'  
EL_SIM2 = loadmat(pjoin(fcal, 'ErrLands-380nm-0.01di.mat'))

for i, row in enumerate(EL_SIM2['qnormSpecRadInt_sim_2D'][:,0]):
    EL_SIM2['qnormSpecRadInt_sim_2D'][i,0] = EL_SIM2['qnormSpecRadInt_sim_2D'][i,0][::5,:] 

EL_SIM2['wl_q'] = EL_SIM2['wl_q'][0,::5]


savemat(pjoin(fcal, 'ErrLands-320nm-0.01di_low-resolution.mat'), EL_SIM2,
        oned_as = 'row', do_compression=True)