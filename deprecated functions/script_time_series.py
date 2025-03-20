import sys
sys.path.append('/home/pi/Documents/PythonScripts/')
from pyGonioSpectrometer import gonio_time_series
from pyGonioSpectrometer.instrumentation import list_spectrometers
from os.path import join as pjoin
import os
#%%
# --------------- DATA saving -------------------------------------------------
# Define folder and filename (without extension)
filename = 'OLED_DX' # FILL
main_folder = '/home/pi/Documents/data/joan/2025/test-updated-software' # FILL
# --------------- GONIO parameters --------------------------------------------
# Define step angle and max angle
angle_step = 5.4 # FILL
angle_max = 86.4 # FILL
# -------------- SPECTROMETER paramaters --------------------------------------
integration_time =  80 # FILL
n_spectra = 2 # FILL
max_time = 1 # Max. time allowed per total adquisition (number of spectra x integration time) in s

# -------------- TIME-SERIES paramaters ---------------------------------------
# Define variables for the time series
interval_luminance = 2 # FILL, interval to take luminance measurements in seconds
interval_gonio = 180 # FILL, interval to take gonio measurements in seconds
stop_luminance_after = 60 # FILL, time at which the program will stop measuring the luminance measurements



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# ----- DO NOT TOUCH BEYOND THIS POINTS UNLESS YOU KNOW WHAT YOU'RE DOING -----
# Assuming there is only one spectrometer, so taking the first element
name_spectrometer = list_spectrometers()[0]

folder = pjoin(main_folder, filename)
filename = filename + '_time-series'
# This code snippet only makes sure that the folder exists
if not os.path.exists(folder):
    os.mkdir(folder)
print (f'INFO: The data will be save at:\n\t{folder}')
#%%
# The function that calls the measurement
gonio_time_series(filename, folder,\
                  integration_time, n_spectra, interval_gonio,\
                  name_spectrometer,\
                  interval_luminance = interval_luminance,\
                  angle_step = angle_step, angle_max = angle_max,\
                  stop_luminance_after = stop_luminance_after,
                  max_time = max_time)


