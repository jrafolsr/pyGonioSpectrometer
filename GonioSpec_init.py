# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 09:32:02 2020

@author: OPEGLAB
"""


#-------------------------------- GUI -----------------------------------------
# Port where the GUI will be deployed, just change if it conflicts with another app
PORT = 8051

#-------------------------------- Spectrometer --------------------------------
# Max number of counts accepted but he Spectromenter, now it corresponds to the OceanOptics Flame value
SATURATION_COUNTS = 65535
# Default integration time and number of spectra taken
INITIAL_INTEGRATION_TIME = 100.0
INITIAL_NSPECTRA = 10

#------------------------------------ Gonio -----------------------------------
# Port where arduino is connected
ARDUINO_PORT = 'ASRL7::INSTR' 
# Angular step, can't be zero or bigger than INITIAL_MAX_ANGLE
INITIAL_STEP = 10 
# Max angle gonio will do
INITIAL_MAX_ANGLE = 80 
 # Waiting time between any gonio movement in seconds. 3 is quite a safe value, could be possibly lowered for faster measruements
WAIT_TIME = 3

#--------------------------------- Saving data---------------------------------
# Default path to the directotry where the file will be saved, unless later modified by user
INITIAL_PATH = r'C:\Users\OPEGLAB\Documents\data\goniospectrometer'
# Default filename where the data will be saved, unless later modified by user
INITIAL_FILENAME = 'angular_sri'