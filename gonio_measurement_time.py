# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:12:33 2020

@author: OPEGLAB
"""
from pyGonioSpectrometer import gonio_measurement
from pyGonioSpectrometer.instrumentation import list_spectrometers

from time import sleep, time

# This will only be used if the script is called as it is, not if it is used as a function
name_motor = 'ASRL7::INSTR'
# Assuming there is only one spectrometer, so taking the first element
name_spectrometer = list_spectrometers()[0]
# Define measurement variables
folder = r'C:\Users\OPEGLAB\Documents\data\goniospectrometer\joan'
filename = 'LEC_E01D2_time-series'
angle_step = 10
angle_max = 80
# Define variables for the spectrometer measurements
integration_time = 200
n_spectra = 20
disable_gonio = False
plot = False
# Define variables for the time series
interval = 600 # interval to wait in seconds

    
def gonio_time_series(interval, name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = 'gonio_mesurement', folder = '.',\
                      disable_gonio = False, plot = False):
    k = 0
    while True:
        try:
            start_time = time()
            
            print(f'\n\n<<<<<INFO: Taking measurement #{k:d}>>>>>')
                  
            gonio_measurement(name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = filename, folder = folder,\
                      disable_gonio = disable_gonio, plot = plot)
            
            ellapsed_time = time() - start_time
            
            
            print(f'INFO: The measurement lasted {ellapsed_time:.1f} s')
            
            sleeping_time = interval - ellapsed_time
            
            print(f'      Going to sleep for {sleeping_time:.0f} s')
            
            
            while ellapsed_time < interval:
                sleep(0.1)
                ellapsed_time = time() - start_time
                
            k += 1
            
        except KeyboardInterrupt:
            print('INFO: Measurement interrupted by the user')
            break
        except Exception as e:
            print(e)
            break
        
if __name__ == '__main__':
    gonio_time_series(interval, name_motor, angle_max, angle_step, name_spectrometer, integration_time, n_spectra,filename, folder, disable_gonio, plot)
    