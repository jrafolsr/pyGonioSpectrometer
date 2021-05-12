# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 09:26:42 2020

@author: JOANRR
"""
from os.path import join as pjoin
from os.path import dirname
# Calibration variables for the data_processing.py and data_fitting.py
fcal = pjoin(dirname(__file__), 'calibration_files')
# Instrument response function, 160504 ML
path_IRF = pjoin(fcal, 'IRF_FlameBlueFiber_wLens2945K.txt')

# Eye response function
path_eye_response = pjoin(fcal, 'CIE1988photopic.txt')

ABS_CALFACTOR = 7.10e6 # Calibration from 15/07/2020
#ABS_CALFACTOR = 1.44e6 # Calibration factor 26/02/2021 for the Raspberry Setup
# ABS_CALFACTOR = 5.4178E6 # OLD Mattias' absolute numbers calibration factor, 170224 ML

# Pixel size of the McScience substrate
PIXEL_SIZE = 4e-6 # m^2
