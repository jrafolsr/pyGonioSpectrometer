# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:24:17 2020

@author: JOANRR
"""
import numpy as np
from scipy.io import loadmat, savemat
from os.path import join as pjoin

# This is manual, so far, needs to be a function
fcal = r'.\calibration_files'  
EL_SIM2 = loadmat(pjoin(fcal, 'ErrorLandscape_2D-OLED_Hres_v1'))
EL_SIM2.pop('test')
EL_SIM2['qnormSpecRadInt_sim_2D'] = EL_SIM2.pop('ans')
EL_SIM2['ipos_sim'] = np.arange(0.01, 1.0, 0.01)
EL_SIM2['dAL_sim'] = np.arange(50, 205, 5)
EL_SIM2['wl_q'] = np.arange(380,781, 1)
EL_SIM2['angle_q'] = np.arange(0,90, 10)


EL_SIM2['ipos_sim'] = EL_SIM2['ipos_sim'].reshape((1, len(EL_SIM2['ipos_sim'])))
EL_SIM2['dAL_sim'] = EL_SIM2['dAL_sim'].reshape((1, len(EL_SIM2['dAL_sim'])))
EL_SIM2['wl_q'] = EL_SIM2['wl_q'].reshape((1, len(EL_SIM2['wl_q'])))
EL_SIM2['angle_q'] = EL_SIM2['angle_q'].reshape((1, len(EL_SIM2['angle_q'])))

savemat(pjoin(fcal, 'ErrorLandscape_2D-OLED_Hres_v1_restructured.mat'), EL_SIM2,
        oned_as = 'row', do_compression=True)