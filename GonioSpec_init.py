# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 09:32:02 2020

@author: OPEGLAB
"""


#-------------------------------- GUI -----------------------------------------
# Port where the GUI will be deployed, just change if it conflicts with another app
PORT = 8052

#-------------------------------- Spectrometer --------------------------------
# Max number of counts accepted but he Spectromenter, now it corresponds to the OceanOptics Flame value
SATURATION_COUNTS = 65535
# Default integration time and number of spectra taken
INITIAL_INTEGRATION_TIME = 100.0
INITIAL_NSPECTRA = 10

#------------------------------------ Gonio -----------------------------------
# Port where arduino is connected
#ARDUINO_PORT = 'ASRL7::INSTR' 
# Angular step, can't be zero or bigger than INITIAL_MAX_ANGLE
INITIAL_STEP = 10.8 
# Max angle gonio will do
INITIAL_MAX_ANGLE = 86.4

 # Waiting time between any gonio movement in seconds. Deprecated, as teh RaspberryMotorController takes care of it, let's put it to 0.1 s just in case
WAIT_TIME = 0.1

#--------------------------------- Saving data---------------------------------
# Default path to the directotry where the file will be saved, unless later modified by user
INITIAL_PATH = r'/home/pi/Documents/data'
# Default filename where the data will be saved, unless later modified by user
INITIAL_FILENAME = 'angular_sri'