# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:12:33 2020

@author: OPEGLAB
"""

from pyGonioSpectrometer import gonio_measurement
from pyGonioSpectrometer.instrumentation import list_ports, list_spectrometers

from time import sleep, time



# This will only be used if the script is called as it is, not if it is used as a function
# Assume the first port is the Arduino
name_motor = list_ports()[1]
# Assuming there is only one spectrometer, so taking the first element
name_spectrometer = list_spectrometers()[0]
# Define measurement variables
folder = r'C:\Users\OPEGLAB\Documents\data\xiaoying'
filename = 'LEC_time_series'
angle_step = 10
angle_max = 80
# Define variables for the spectrometer measurements
integration_time = 250
n_spectra = 10
disable_gonio = False
plot = False
# Define variables for the time series
interval = 300 # interval to wait in seconds


def sleep_tuned(seconds):
    for i in range(int(seconds*10)):
        sleep(0.1)
        
def gonio_time_series(interval, name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = 'gonio_mesurement', folder = '.',\
                      disable_gonio = False, plot = False):
    k = 0
    while True:
        try:
            start_time = time()
            print(f'INFO: Taking measurement #{k:d}')
            gonio_measurement(name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = filename, folder = folder,\
                      disable_gonio = disable_gonio, plot = plot)
            ellapsed_time = time() - start_time
            sleeping_time = interval - ellapsed_time
            print(f'INFO: The measurement lasted {ellapsed_time:.1f} s')
            print(f'      Going to sleep for {sleeping_time:.0f} s')
            sleep_tuned(sleeping_time)
            k += 1
        except KeyboardInterrupt:
            print('INFO: Measurement interrupted by the user')
            break
        except Exception as e:
            print(e)
            break
        
if __name__ == '__main__':
    gonio_time_series(interval, name_motor, angle_max, angle_step, name_spectrometer, integration_time, n_spectra,filename, folder, disable_gonio, plot)