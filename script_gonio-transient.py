#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 20:45:22 2022

@author: JOANRR
"""
from pyGonioSpectrometer import GonioLogger
from time import sleep,monotonic
from datetime import datetime
from pathlib import Path
import numpy as np
from pyGonioSpectrometer import Timer
import matplotlib.pyplot as plt
import traceback

# Initilize the data saving parameters
folder = Path('/home/pi/Documents/data/joan/2025/test-updated-software')
file_id = 'example'

# Initialize the GONIO measurement parameters0
integration_time = 220             # Initial guess for the integration time
n_spectra = 2                       # I would go for one spectra to average t have a quick initial luminance
max_n_spectra = 2                   # Maximum number of spectra to take (to limit the time to adquire a spectra)
angle_step = 5.4                    # Angle step in deg, needs to be a multiple of 1.8 deg
angle_max = 86.4                   # Max angle in deg, needs t be a multple of angle_step
max_time_per_fwd_luminance = 2000 #2500   # Max time allowed to integrate the forward luminance, in ms, it basically sets the time min step for the luminance
max_time_per_angle = 440 #3000           # Maximum time in miliseconds allowed per angle step (number of spectra x integration time)

luminance_interval = 2 # Time interval at which to take forward luminance steps in seconds
stop_luminance_after = 60
gonio_interval = 300    # Time interval at which to take a full gonio scan in seconds
stop_gonio_after = 3600

max_intensity_angle = 0.0

if not folder.exists(): folder.mkdir()


# Initialize the Goniologger object with the given parameters
goniospectrometer = GonioLogger(file_id, folder,\
                                angle_step = angle_step, angle_max = angle_max,\
                                max_time_per_angle = max_time_per_fwd_luminance,\
                                integration_time = integration_time, n_spectra = n_spectra)

goniospectrometer.max_intensity_angle = max_intensity_angle

# Take a dark spectra before starting the measurement
shutter_is_closed = bool(int(input('Is the shutter closed?\n\t0 --> No\n\t1 -> Yes\n')))
if not shutter_is_closed:
    goniospectrometer.gonio.shutter_is_closed = False
    goniospectrometer.gonio.close_shutter()
    
#print('OBS! Make sure that the initial shutter position is closed!! (LEC at 180 deg with respect the spectrometer)')
goniospectrometer.take_dark_spectra()

# Set the timer object to track the luminance measurements
gonio_timer_luminance = Timer(min_time_step = luminance_interval, fix_step = True)
gonio_timer_fullscan = Timer(min_time_step=gonio_interval, max_time_step = 300, fix_step=True)
#                             intervals=(3600,3620,3630, 3640), fix_step=False)   #FOR OLED ANGE OFFSET TEST

timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")

try:
        gonio_timer_luminance.initialize() 
        gonio_timer_fullscan.initialize()
        k = 0
        print(f'\n\t<<<<< INFO: Measuring L0  every {luminance_interval:.0f} s >>>>>')
        while gonio_timer_luminance.ellapsed_time() < stop_luminance_after:
            # Take the chance to check if everythin is within the limits
            limit_hit = goniospectrometer.check_low_and_high_limits()
            if limit_hit == 1:
                # Upper limit hitted, update integration time only if we are not close to the min laready (unlikely)
                if goniospectrometer.integration_time > 1: # one mllisecond is the mi the spec can take
                    goniospectrometer.update_integration_time()
            elif limit_hit == -1:
                # Lower limit hitted, update integration time only if we are not close to the max already 
                  if goniospectrometer.integration_time < max_time_per_fwd_luminance: # one mllisecond is the mi the spec can take
                    goniospectrometer.update_integration_time()
            else:
                # don't do anything, all good
                pass

            # This block take a full gonio scan
            if gonio_timer_luminance.istime2measure():
                print(f'\rINFO: Taking spectra n.{k + 1: 2d} at forward luminance...      ', end = '')
                plt.close()
                goniospectrometer.save_spectra(gonio_timer_fullscan.ellapsed_time())
                k+=1
            sleep(0.01)
        
        print()
        
        while gonio_timer_fullscan.ellapsed_time() < stop_gonio_after:
            _temp = int(gonio_interval) - int(gonio_timer_fullscan.ellapsed_time()) % int(gonio_interval) 
            print(f'\rINFO: Next gonio measurement in... {_temp : 4.0f} s', end = '')
            if gonio_timer_fullscan.istime2measure():
                # Update the integration times but with the reference at the angle 43.2 deg, which is higher in teh half-cylinder configuraton
                goniospectrometer.update_integration_time(angle = goniospectrometer.max_intensity_angle)
                print(f'\nINFO: Moved to angle {goniospectrometer.max_intensity_angle:.2f} deg to perform the integration time adjustment.')
                
                print('\nINFO: Taking a full goniometer scan!\n')
                # Header and suffix for the gonio files
                header = f'# Ellapsed time (s): {gonio_timer_fullscan.ellapsed_time():10.2f}\n'
                suffix = f'_time-series'
                plt.close()
                goniospectrometer.take_gonio_measurement(suffix = suffix, header=header)
            sleep(0.001)
            
    
except KeyboardInterrupt:
    print('\nINFO: Terminating program')  
except Exception as e:
    print(e)
    # Print the traceback information
    traceback.print_exc()
    print('\nINFO: Terminating program')  


goniospectrometer.shutdown()

print('INFO: Program terminated.')