# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:12:33 2020

@author: OPEGLAB
"""
#%%
from pyGonioSpectrometer import gonio_measurement, write_to_file
from pyGonioSpectrometer.instrumentation import list_spectrometers, SpectraMeasurement, ArduinoMotorController

from time import sleep, time
import numpy as np
from datetime import datetime
from os.path import join as pjoin

SATURATION_COUNTS = 65535 # Saturation limit of the spectrometer
UPPER_LIM = 58000 # Max. number of counts allowed before reducing the integration time
LOWER_LIM = 40000 # Min. number of counts allowed before incresing the integration time
MAX_TIME = 5000 # Max. time allowed per total adquisition (number of spectra x inegration time) in ms
#%%

# Define default vairables for gonio measurement
angle_step = 5
angle_max =80
# Define variables for the time series
interval_luminance = 10 # interval to take luminance measurements
interval_gonio = 30 # # interval to take gonio measurements


def gonio_time_series(filename, folder,\
                      integration_time, n_spectra, interval_gonio,\
                      name_motor, name_spectrometer,\
                      angle_step = angle_step, angle_max = angle_max,\
                      interval_luminance = interval_luminance,
                      ):

    # Define general variables
    disable_gonio = False # Gonio always active
    plot = False  # No plotting
    
    # Raise error if the time to sample the spectra is shorter than the time to take it
    if (n_spectra * integration_time/1000) >= (interval_luminance * 0.9):
        raise Exception('ERROR: The interval_luminance is less than 90 % of the time need to take the spectra, consider reducing n_spectra, integration_time or increasing interval_luminance')

    time_zero = time() # General inital timer
    
    # Initialize counters and flags
    shutter_open = False
    k = 0
    
    while True:
        k += 1
        
        # Adaptative interval_luminance
        total_ellapsed_time = time() - time_zero
        if total_ellapsed_time > 1200: # After 20 min, take forward luminance every min
            interval_luminance = min(60, interval_gonio)
        elif total_ellapsed_time > 300:  # After 2 min, take forward luminance every 30 seconds min
            interval_luminance = 30     
#            
        ninterval = interval_gonio // interval_luminance
        
        try:
            print(f'\n\t<<<<< INFO: Measuring L0  every {interval_luminance:.0f} s >>>>>')
            gonio = ArduinoMotorController(name_motor)
            
            if not gonio.motor.bytes_in_buffer == 0:
                print(gonio.motor.read())
                
            flame = SpectraMeasurement(name_spectrometer, integration_time = integration_time, n_spectra = n_spectra)
            
            # Timestamps for the header and filename
            itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
            timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
            start_luminance_adquisition = time()
            path = pjoin(folder, timestamp + filename + '_L0.dat')
        
            with open(path, 'a') as f:
                f.write(itimestamp + ' # Timestamp at the beginning of the measurement\n')
                f.write(f'{integration_time:.0f} # Integration time in (ms)\n')
                f.write(f'{n_spectra} # Number of spectra taken\n')
                f.write(f'# Rel.time\t Integration time\t Wavelengths\n')
                
            # Take the dark spectra at zero
            wavelengths = flame.get_wavelengths()
            write_to_file(np.nan, integration_time, wavelengths, path)
            intensities = flame.get_averaged_intensities()
            write_to_file(np.nan, integration_time, intensities, path)
            
            
            print('INFO: Opening shutter')
            gonio.move_shutter(delay = 0.5)
            shutter_open = True
            sleep(0.5)
            
            for i in range(ninterval):
                start_time = time()
                print(f'INFO: Taking spectra n.{i + 1: 2d} at forward luminance.')
                
                ellapsed_time_since_bkg = time() - start_luminance_adquisition
                intensities = flame.get_averaged_intensities()
                            
                write_to_file(ellapsed_time_since_bkg, integration_time, intensities, path)
                ellapsed_time = time() - start_time
                
                # Wait the amount of time specified by interval_luminance
                while ellapsed_time < interval_luminance:
                    ellapsed_time = time() - start_time
                    sleep(0.001)

            # Prepare for the gonio measurement.       
            # Checking if the integration time needs to be increased or decreased
            max_counts = flame.get_intensities().max()
            
            if max_counts <= LOWER_LIM:
                print('INFO: Lower limit reach, incresing integration time')
                integration_time = min(int(UPPER_LIM / max_counts * integration_time), MAX_TIME)
                if integration_time * n_spectra > MAX_TIME:
                    n_spectra = max(1, int(MAX_TIME / integration_time))
                    
            elif max_counts >= 0.9 * SATURATION_COUNTS:
                print('INFO: Saturation reach, decreasing integration time')
                while max_counts >= 0.9 * SATURATION_COUNTS:
                    integration_time *= 0.95
                    flame.config(integration_time)
                    max_counts = flame.get_intensities().max()
                    
                integration_time = max(int(integration_time), 10)
                n_spectra = min(20, int(MAX_TIME / integration_time))
            
            print(f'INFO: The adquisition is set to be  {n_spectra} x {integration_time} ms') 
            # Close the sutter as the next step will be the gonio measurement
            print('INFO: Closing shutter')
            gonio.move_shutter()
            shutter_open = False
            gonio.close()   
            flame.close()
            
            
            sleep(5.0)
            
            # The gonio measurement itself
            shutter_open = True
            print(f'\n\t<<<<< INFO: Taking measurement #{k:d} >>>>>')
            # Gonio measurements
            gonio_measurement(name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = filename, folder = folder,\
                      disable_gonio = disable_gonio, plot = plot)
            shutter_open = False

            
            
        except KeyboardInterrupt:
            print('INFO: Measurement interrupted by the user')
            break
        except Exception as e:
            print(e)
            break
        
    if shutter_open:
        gonio.move_shutter()
        
    flame.close()
    gonio.close()
    
if __name__ == '__main__':
    # Define folder and filename
    folder = pjoin(r'C:\Users\OPEGLAB\Documents\data\goniospectrometer')    
    filename = 'default_time-series'
    # Gonio measurement
    name_motor = 'ASRL7::INSTR'
    # Define variables for the spectrometer measurements
    # Assuming there is only one spectrometer, so taking the first element
    name_spectrometer = list_spectrometers()[0]
    integration_time = 100
    n_spectra = 1
    
    gonio_time_series(filename, folder,\
                      integration_time, n_spectra, interval_gonio,\
                      name_motor, name_spectrometer,\
                      interval_luminance = interval_luminance,\
                      angle_step = angle_step, angle_max = angle_max
                      )