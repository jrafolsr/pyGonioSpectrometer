# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:12:33 2020

@author: OPEGLAB
"""

from pyGonioSpectrometer import gonio_measurement
from pyGonioSpectrometer.instrumentation import list_ports, list_spectrometers

from time import sleep



# This will only be used if the script is called as it is, not if it is used as a function
# Assume the first port is the Arduino
name_motor = list_ports()[1]
# Assuming there is only one spectrometer, so taking the first element
name_spectrometer = list_spectrometers()[0]
# Define measurement variables
folder = r'.'
filename = 'CL-random_time-series'
angle_step = 10
angle_max = 80
# Define variables for the spectrometer measurements
integration_time = 100
n_spectra = 10
disable_gonio = False
plot = False
# Define variables for the time series
interval = 180 # interval to wait in seconds


def sleep_tuned(seconds):
    for i in range(int(seconds)):
        sleep(1)
        
def gonio_time_series(interval, name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = 'gonio_mesurement', folder = '.',\
                      disable_gonio = False, plot = False):
    k = 0
    while True:
        try:
            print(f'INFO: Taking measurement #{k:d}')
            gonio_measurement(name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = filename, folder = folder,\
                      disable_gonio = disable_gonio, plot = plot)
            
            print(f'INFO: Going to sleep for {interval:.0f} s')
            sleep_tuned(interval)
            k += 1
        except KeyboardInterrupt:
            print('INFO: Measurement interrupted by the user')
            break
        except Exception as e:
            print(e)
            break
        
if __name__ == '__main__':
    gonio_time_series(interval, name_motor, angle_max, angle_step, name_spectrometer, integration_time, n_spectra,filename, folder, disable_gonio, plot)