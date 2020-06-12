# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:25:55 2020

@author: JOANRR
"""
#%%
from os.path import join as pjoin
from time import sleep
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import winsound
from pyGonioSpectrometer.instrumentation import list_ports, ArduinoMotorController, SpectraMeasurement, list_spectrometers

# Parameters for the beeping, quite irrellevant, I'll consider to remove it
frequency = 2000  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second

# Fixed parameters
WAIT_TIME = 3
N_WAVELENGTHS = 2028

# This will only be used if the script is called as it is, not if it is used as a function
# Assume the first port is the Arduino
name_motor = list_ports()[0]
# Assuming there is only one spectrometer, so taking the first element
name_spectrometer = list_spectrometers()[0]
# Define measurement variables
folder = r'.'
filename = '_goniomeasurement'
angle_step = 20
angle_max = 80
# Define variables for the spectrometer measurements
integration_time = 100
n_spectra = 10
plot = True

def gonio_measurement(name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = 'gonio_mesurement', folder = '.',\
                      disable_gonio = False, plot = True):
    # Creates the object that will control the steppers
    gonio = ArduinoMotorController(name_motor)
    # Create the object instance taht will control the spectrometer, assuming it is the first of the list
    flame = SpectraMeasurement(name_spectrometer, integration_time=integration_time, n_spectra= n_spectra)
    #%%%
    
    n_angles = int(angle_max*100) // int(angle_step*100) + 1
    n_steps = (n_angles - 1) *2
    n_columns = 2*n_angles -1 + 4
    
    # Matrices where to save the data, based on Mattias' scheme (can be improved...)
    first_row = np.zeros((1,n_columns))  
    first_row[0,0] = integration_time
    first_row[0,1] = n_spectra
    # First columns will be the wavellengths, second the dark spectra, the rest the angles
    data = np.zeros((N_WAVELENGTHS, n_columns))
    
    # Prepraring the plot
    if plot:
        plt.ion() # Will update during the measurement
        fig, ax = plt.subplots()
        ax.set_ylabel('Counts')
        ax.set_xlabel('Wavelength (nm)')
        plot_measurement(fig, ax, [], [])
    
    # Take the dark spectra at zero
    wavelengths = flame.get_wavelengths()
    data[:,0] = wavelengths
    
    winsound.Beep(frequency, duration*2) # It will remind that the measureement starts
    
    print('INFO: Taking dark spectra....')
    temp = flame.get_averaged_intensities()
    if plot: plot_measurement(fig, ax, wavelengths, temp, 'dark')
    data[:,1] = temp
    sleep(0.25)
    # Open the shutter, assuming it is closed
    gonio.move_shutter()
    [winsound.Beep(frequency,200) for i in range (5)] # Reminder of measurement starting
    sleep(WAIT_TIME*2)
    
    # Take spectra at zero
    print('>>>>>>>>>>>>>>> STEP 0 <<<<<<<<<<<<<<<')
    print('INFO: Taking spectra at 0.00°...')
    temp = flame.get_averaged_intensities()
    if plot: plot_measurement(fig, ax, wavelengths, temp, '0.00°')
    data[:,2] = temp
    first_row[0,2] = 0.0
    
    
    # Move to the last position
    print(f'>>>>>>>>>>>>>>> STEP 1 <<<<<<<<<<<<<<<')
    out_angle = gonio.move_angle(-1.0 * angle_max)
    # Initialize the error made in each step
    error = (angle_max - out_angle)
    total = 0
    current_angle = -out_angle
    print(f'INFO: Moved {out_angle:.2f}°. Total angle moved: {total:.2f}')
    
    sleep(WAIT_TIME * 2) # Long enough time to make sure that it waits until the ed of the movement
    
    
    for k in range(n_steps):        
    
        print(f'INFO: Taking spectra at {current_angle:.2f}°...')
        
        # Save the angle
        first_row[0, k + 3] = current_angle
        # Get the whole spectra (wl and I)
        temp = flame.get_averaged_intensities()
        data[:, k + 3] = temp
    
        # Check for any values higher than saturation
        if np.any(temp > 65535): print('WARNING: Some values are saturating. Consider lowering the integration time.')
        
        # Plot the data
        if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.2f}°')
        
        # Calculating the current step to make
        current_step = angle_step + error
        # Moving the gonio
        out_angle = gonio.move_angle(current_step)
        # New error made
        error = (current_step - out_angle)
        total += abs(out_angle)
        current_angle += out_angle
        
        print(f'>>>>>>>>>>>>>>> STEP {k+2:d} <<<<<<<<<<<<<<<')
        print(f'INFO: Moved {out_angle:.2f}°. Total angle moved: {total:.2f}')
        
        sleep(WAIT_TIME)
    
    # Take last angle spectra
    print(f'INFO: Taking spectra at {current_angle:.2f}°...')
    first_row[0,k + 4] = current_angle
    temp = flame.get_averaged_intensities()
    data[:,k + 4] = temp
    
    if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.2f}°')
    sleep(WAIT_TIME)
    
    # Go back to zero
    back_angle = -1.0 * abs(current_angle)
    out_angle = gonio.move_angle(back_angle)
    current_angle -= out_angle
    print(f'>>>>>>>>>>>>>>> STEP {k+3:d} <<<<<<<<<<<<<<<')
    print(f'INFO: Moved {out_angle:.2f}°. Total angle moved: {total:.2f}')
    
    # Wait longer time, as the angle is larger and take the spectra
    sleep(WAIT_TIME*2)
    temp = flame.get_averaged_intensities()
    data[:, k + 5] = temp
    
    if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.2f}°')
    
    # Closing shutter
    gonio.move_shutter()
    
    # Saving the data code snippet
    timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
    data = np.vstack([first_row, data])
    timestamp2 = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
    path = pjoin(folder, timestamp2 + filename + '.dat')
    np.savetxt(path, data, fmt = '% 8.2f', header= timestamp)
    
    print('INFO: Measurement finished data saved at\n\t' + path)
    
    if disable_gonio: gonio.disable_gonio()
    
    gonio.close()
    flame.close()
    
def plot_measurement(fig, ax, x, y, label = None):
    ax.plot(x,y, label = label)
    ax.legend()
    plt.pause(0.25)
    
if __name__ == '__main__':
    gonio_measurement(name_motor, angle_max, angle_step, name_spectrometer, integration_time, n_spectra)