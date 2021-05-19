# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:12:33 2020

@author: OPEGLAB
"""
#%%
from pyGonioSpectrometer import gonio_measurement, write_to_file
from pyGonioSpectrometer.instrumentation import list_spectrometers, SpectraMeasurement, RaspberryMotorController

from time import sleep, time
import numpy as np
from datetime import datetime
from os.path import join as pjoin

SATURATION_COUNTS = 65535 # Saturation limit of the spectrometer
UPPER_LIM = 58000 # Max. number of counts allowed before reducing the integration time
LOWER_LIM = 20000 # Min. number of counts allowed before incresing the integration time

#%%


def gonio_time_series(filename, folder,\
                      integration_time, n_spectra, interval_gonio,\
                      name_spectrometer,\
                      angle_step = 5, angle_max = 80,\
                      interval_luminance = 10,\
                      stop_luminance_after = 7200,\
                      max_time = 5):
    """
    Performs a time series measurement. It takes first the forward luminance every interval_luminance seconds. Every interval_gonio, it does a complete goniometer measurement (with the specified angles and step). There are a few options to control the time intervals.
    
    Parameters
    ----------
    filename : str
        A string containing the output filename. The default is 'gonio_measurement'.
    folder : str
            Path, relative or absolute, to the directory where to save the data. It should exists, the program does not check if it does.
    integration_time: int or float
        Sets the integration time in ms.
    n_spectra : int
        Number of spectra that will be averaged.
    name_spectrometer: seabreeze.spectrometers.Spectrometer class
        The spectrometer resource as the specified class. You can get it from list_spectrometers().
    angle_step : int or float, optional
        Angular step that the motor will perform in degrees. Preferably a divisor of angle_max and an integer, the program does not check for these conditions to be fullfilled. The default is 5.
    angle_max : int or float, optional
        Maximum angle to scan with the goniometer in degrees. The default is 80.
    interval_luminance : int, optional
        Interval in seconds at which teh forward luminance will be taken. The default is 10 s.
    interval_gonio : int, optional
        Interval in seconds at which teh forward luminance will be taken. The default is 120 s.
    stop_luminance_after : int, optional
        Time in seconds at which the program will stop performing the forward luminance measurements. Passed this time, only a goniometer measurement every interval_gonio will be performed. The default is 7200.
    max_time : int, optional
        Maximum time in seconds allowed per angle step (number of spectra x integration time). The default is 5000 s.
    """
    # Making sure that the angle step and the max angel are compatible, otherwise, raise an error
    _ , residu = divmod(angle_max, angle_step)
    if residu != 0.0:
            raise Exception("The angle step is not a divisor of the max. angle!")
            
    # Raise error if the time to sample the spectra is shorter than the interval_luminance
    if (n_spectra * integration_time / 1000) >= (interval_luminance * 0.9):
        raise Exception('The interval_luminance is less than 90 % of the time need to take the spectra, consider:\n 1. Reducing n_spectra or integration_time or \nor\n2. Increasing interval_luminance')
        
    # Convert max_time into milliseconds
    MAX_TIME = int(max_time * 1000)
    
    # General initial timer
    time_zero = time()
    
    # Initialize counters and flags and gonio
    gonio = RaspberryMotorController()
    shutter_open = False # Keeps track of the shutter status on the function level.
    close_resources = True # Keeps track of the shutter status of the resources
    # Loop counter
    k = 0
    
    
    
    while True:
        k += 1
        # Gonio measurements
        start_time_gonio = time()
        
        # Adjust the interval_luminance
        total_ellapsed_time = time() - time_zero
        if total_ellapsed_time > 1800: # After 30 min take spectra min
            interval_luminance = max(interval_gonio, 60)
        elif total_ellapsed_time > 300:  # After 5 min, take forward luminance every 30 seconds min
            interval_luminance = max(60, interval_luminance)   
#            
        ninterval = int(interval_gonio // interval_luminance)
        
        try:

            flame = SpectraMeasurement(name_spectrometer,\
                                       integration_time = integration_time,\
                                       n_spectra = n_spectra)
            close_resources = True
            
            # We perform only forward luminance measurements, i.e. at angle 0.
            if total_ellapsed_time <= stop_luminance_after:
                print(f'\n\t<<<<< INFO: Measuring L0  every {interval_luminance:.0f} s >>>>>')
                # Timestamps for the header and filename
                itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
                timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
                # Initialized the inner timer for the L0 interval
                start_luminance_adquisition = time()
                # Path where the data will be saved
                path = pjoin(folder, timestamp + filename + '_L0.dat')
            
                with open(path, 'a') as f:
                    f.write(itimestamp + ' # Timestamp at the beginning of the measurement\n')
                    f.write(f'{integration_time:.0f} # Integration time in (ms)\n')
                    f.write(f'{n_spectra} # Number of spectra taken\n')
                    f.write('# Rel.time\t Integration time\t Wavelengths\n')
                    
                # Take the dark spectra at zero
                wavelengths = flame.get_wavelengths()
                write_to_file(np.nan, integration_time, wavelengths, path)
                intensities = flame.get_averaged_intensities()
                write_to_file(np.nan, integration_time, intensities, path)
                
                sleep(0.25)
                gonio.move_shutter()
                shutter_open = True
                sleep(0.5)
                
                for i in range(ninterval):
                    print(f'\rINFO: Taking spectra n.{i + 1: 2d} at forward luminance...      ', end = '')
                    start_time = time()
                    
                    ellapsed_time_since_bkg = start_time - start_luminance_adquisition
                    
                    intensities = flame.get_averaged_intensities()
                    
                    print(f'\rINFO: Taking spectra n.{i + 1: 2d} at forward luminance... Done!', end = '')
                    
                    # Check for any values higher than saturation
                    if np.any(intensities > SATURATION_COUNTS):
                        print('\n! WARNING: Some values are saturating. Consider lowering the integration time.')
                    elif intensities.max() < 10000:
                        print('\n! WARNING: The max. count is less than 10000. Consider increasing the integration time')
                    
                    write_to_file(ellapsed_time_since_bkg, integration_time, intensities, path)
                    ellapsed_time = time() - start_time
                    
                    # Wait the amount of time specified by interval_luminance
                    while ellapsed_time < interval_luminance:
                        ellapsed_time = time() - start_time
                        sleep(0.001)
                   
                # Close the sutter as the next step will be the gonio measurement
                print('\nINFO: Closing shutter')
                gonio.move_shutter()
                shutter_open = False

            # If the time is >= stop_luminance_after, skip the forward luminance measurement and just do the check to adapt the int_time
            sleep(1.0)   
            print('INFO: Self-adjusting the integration time.')

            gonio.move_shutter()
            shutter_open = True
            
            integration_time, n_spectra = flame.adjust_integration_time(max_time = MAX_TIME,\
                                                                        lower_limit = LOWER_LIM,
                                                                        upper_limit = UPPER_LIM)

            print('INFO: Closing shutter')
            gonio.move_shutter()
            shutter_open = False
            close_resources = False
                
            flame.close()
            
            sleep(2.0)
            
            # The gonio measurement itself
            shutter_open = True
            
            print(f'\n\t<<<<< INFO: Taking measurement #{k:d} >>>>>')


            gonio_measurement(angle_max, angle_step,\
                      gonio,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = filename, folder = folder,\
                      disable_gonio = False, plot = False)

            
            shutter_open = False
            
            # Wait the amount of time specified by gonio
            while (time() - start_time_gonio) < interval_gonio:
                minutes, seconds =  divmod(interval_gonio - (time() - start_time_gonio), 60)  
                print(f'\rINFO: Next gonio measurement in {minutes:02.0f}:{seconds:02.0f}...', end = '')
                sleep(0.5)
            print('\n')
            
        except KeyboardInterrupt:
            print('\nINFO: Measurement interrupted by the user')
            break
        except Exception as e:
            print(e)
            break

    if shutter_open & close_resources:
        gonio.move_shutter()
    if close_resources:
        flame.close()
        gonio.close()

    
if __name__ == '__main__':
    # Define folder and filename
    folder = pjoin('/home/pi/Documents/data')    
    filename = 'default_time-series'
    # Define variables for the spectrometer measurements
    # Assuming there is only one spectrometer, so taking the first element
    name_spectrometer = list_spectrometers()[0]
    integration_time = 100
    n_spectra = 1
    interval_gonio = 120 # # interval to take gonio measurements
    
    
    gonio_time_series(filename, folder,\
                      integration_time, n_spectra, interval_gonio,\
                      name_spectrometer)