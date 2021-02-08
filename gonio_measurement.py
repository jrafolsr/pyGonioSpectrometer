# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:25:55 2020

@author: JOANRR
"""
#%%
from os.path import join as pjoin
from time import sleep, time
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import winsound
from pyGonioSpectrometer.instrumentation import list_ports, ArduinoMotorController, SpectraMeasurement, list_spectrometers

# Parameters for the beeping, quite irrellevant, I'll consider to remove it
frequency = 2000  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second

# Fixed parameters
WAIT_TIME = 1.5
N_WAVELENGTHS = 2028


def gonio_measurement(name_motor,angle_max, angle_step,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = 'gonio_mesurement', folder = '.',\
                      disable_gonio = False, plot = True):
    """
    Performs a complete measurement for the goniospectrometer setup, by taking spectra at every specified angle, for the whole forward hemisphere.
    
    Parameters
    ----------
    name_motor : str
        String containing the COM address where the motor driver is located (the Arduino). You can get the available ports by using the list_ports() function.
    angle_max : int or float
        Maximum angle to scan with the goniometer in deg.
    angle_step : int or float
        Angular step that the motor will perform in deg. Preferably a divisor of angle_max and integer, the program does not check for this conditions to be fullfilled.
    name_spectrometer : seabreeze.spectrometers.Spectrometer class
        The spectrometer resource as the specified class. You can get it from list_spectrometers().
    integration_time : int or float
        Sets the integration time in ms.
    n_spectra: int
        Number of spectra that will be averaged.
    filename : str, optional
        A string containing the output filename. The default is 'gonio_measurement'.
    folder : str, optional
        Path, relative or absolute, to the directory where to save the data. Should exists, the program does not check if it does. The default is '.'.
    disable_gonio : boolean, optional
        Whether to disable the goniometer motor after the measruement. The default is False.
    plot : boolean, optional
        Whether to plot or not the data. It is better to set it to false for a time series measurement, otherwise one will end with to many open windows. The default is True.

    """ 
    # Initializing the gonio and flame objects.
    
    gonio = None
    flame = None
    
    try:
        # Creates the object that will control the steppers
        gonio = ArduinoMotorController(name_motor)
        # Create the object instance taht will control the spectrometer, assuming it is the first of the list
        flame = SpectraMeasurement(name_spectrometer, integration_time=integration_time, n_spectra= n_spectra)
        flame.open()
       
    
        n_angles = int(angle_max*100) // int(angle_step*100) + 1
        n_steps = (n_angles - 1) *2

        
        # Matrices where to save the data, based on Mattias' scheme (can be improved...)
#        n_columns = 2*n_angles -1 + 4        
#        first_row = np.zeros((1,n_columns))  
#        first_row[0,0] = integration_time
#        first_row[0,1] = n_spectra
        # First columns will be the wavellengths, second the dark spectra, the rest the angles
#        data = np.zeros((N_WAVELENGTHS, n_columns))
        
        # Prepraring the plot
        if plot:
            plt.ion() # Will update during the measurement
            fig, ax = plt.subplots()
            ax.set_ylabel('Counts')
            ax.set_xlabel('Wavelength (nm)')
            plot_measurement(fig, ax, [], [])
        
        # Timestamps for the header and filename
        itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
        timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
        start_time = time()
        
        path = pjoin(folder, timestamp + filename + '.dat')
        
        with open(path, 'a') as f:
            f.write(itimestamp + ' # Timestamp at the beginning of the measurement\n')
            f.write(f'{integration_time:.0f} # Integration time in (ms)\n')
            f.write(f'{n_spectra:d} # Number of spectra taken\n')
              
        # Take the dark spectra at zero
        wavelengths = flame.get_wavelengths()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, np.nan, wavelengths, path)
        # Saving the data in Mattias's scheme
#        data[:,0] = wavelengths
        
        winsound.Beep(frequency, duration*2) # It will remind that the measureement starts
        
        print('INFO: Taking dark spectra....')
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, np.nan, temp, path)
      
        if plot: plot_measurement(fig, ax, wavelengths, temp, 'dark')
        # Saving the data in Mattias's scheme
#        data[:,1] = temp
        sleep(0.25)
        
        # Open the shutter, assuming it is closed
        print('INFO: Opening shutter')
        gonio.move_shutter()

        
        [winsound.Beep(frequency,100) for i in range (5)] # Reminder of measurement starting
        sleep(WAIT_TIME)
        
        # Take spectra at zero
        print('>>>>>>>>>>>>>>> STEP 0 <<<<<<<<<<<<<<<')
        print('INFO: Taking spectra at 0.00°')
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, 0.0, temp, path)
        
        if plot: plot_measurement(fig, ax, wavelengths, temp, '0.0°')
        
        # Saving the data in Mattias's scheme
#        data[:,2] = temp
#        first_row[0,2] = 0.0
        
        
        # Move to the last position
        print('>>>>>>>>>>>>>>> STEP 1 <<<<<<<<<<<<<<<')
        out_angle = gonio.move_angle(-1.0 * angle_max)
        # Initialize the error made in each step
        error = (angle_max - out_angle)
        total = 0
        current_angle = -out_angle
        print(f'INFO: Moved {out_angle:.1f}°. Total angle moved: {total:.1f}')
        
        sleep(WAIT_TIME * 2) # Long enough time to make sure that it waits until the ed of the movement
        
        k = 0
        for k in range(n_steps):        
        
            print(f'INFO: Taking spectra at {current_angle:.1f}°')
            

            # Get the whole spectra (wl and I)
            temp = flame.get_averaged_intensities()
            # Saving the data in the new scheme
            write_to_file(time()-start_time, current_angle, temp, path)
            
#            # Save the angle in Mattias' scheme
#            first_row[0, k + 3] = current_angle
#            # Saving the data in Mattias's scheme
#            data[:, k + 3] = temp
        
            # Check for any values higher than saturation
            if np.any(temp > 65535): print('WARNING: Some values are saturating. Consider lowering the integration time.')
            
            # Plot the data
            if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.1f}°')
            
            # Calculating the current step to make
            current_step = angle_step + error
            # Moving the gonio
            out_angle = gonio.move_angle(current_step)
            # New error made
            error = (current_step - out_angle)
            total += abs(out_angle)
            current_angle += out_angle
            
            print(f'>>>>>>>>>>>>>>> STEP {k+2:d} <<<<<<<<<<<<<<<')
            print(f'INFO: Moved {out_angle:.1f}°. Total angle moved: {total:.1f}')
            
            sleep(WAIT_TIME)
        
        # Take last angle spectra
        print(f'INFO: Taking spectra at {current_angle:.1f}°...')
        
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, current_angle, temp, path)
#        # Saving the data in Mattias's scheme
#        first_row[0,k + 4] = current_angle        
#        data[:,k + 4] = temp
        
        if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.1f}°')
        
        # Go back to zero
        back_angle = -1.0 * abs(current_angle)
        out_angle = gonio.move_angle(back_angle)
        current_angle -= out_angle
        
        print(f'>>>>>>>>>>>>>>> STEP {k+3:d} <<<<<<<<<<<<<<<')
        print(f'INFO: Moved {out_angle:.1f}°. Total angle moved: {total:.1f}')
        
        # Wait longer time, as the angle is larger and take the spectra
        sleep(WAIT_TIME * 2)
        print('INFO: Taking last 0° spectra....')
        
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, current_angle, temp, path)
        
        # Saving the data in Mattias's scheme 
#        first_row[0,k + 5] = current_angle 
#        data[:, k + 5] = temp
        
        if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.1f}°')
        
        # Closing shutter
        print('INFO: Closing shutter')
        gonio.move_shutter()
        sleep(0.5)
#        # Saving the data code snippet in with Mattias' scheme
#        ftimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
#        path2 = pjoin(folder, timestamp + '_'+ filename + '_Mattias.dat')
#        data = np.vstack([first_row, data])       
#        header = itimestamp + '\n' + ftimestamp
#        np.savetxt(path2, data, fmt = '% 8.2f', header= header) 
#        print('INFO: Measurement finished data saved at\n\t' + path2)
        
        if disable_gonio: gonio.disable_gonio()
               
    
    except KeyboardInterrupt:
        print('INFO: The angle scan has been cancelled by the user. Going back to the zero position.')
        if gonio != None:
            # Go back to since the spectrogoniometer movement has been cancelled.
            back_angle = -1.0 * abs(current_angle)           
    
    except Exception as e:
        print(e)
        print('INFO: Some error has ocurred during the angle scan. Going back to the zero position.')
        if gonio != None:
            # Go back to since some error has occurred during the gonio measurement
            back_angle = -1.0 * abs(current_angle)

    if gonio != None: 
        gonio.close()
    if flame != None:
        flame.close()
        
        
def plot_measurement(fig, ax, x, y, label = None):
    """
    Plots the spectra given the figure and axis handlers together with the x, and y data.
    The pause is needed to allow the live-plotting.
    """
    ax.plot(x,y, label = label)
    if label is not None:
        ax.legend()
    plt.pause(0.25)
    
def write_to_file(etime, angle, data, file, debug = False):
    """
    Writes into a file the data taken at each gonio step. It takes as inputs the etime (ellapsed time) and angle as scalers, the data (vector) and the file (as the path where to save the data). The debug is just to print it on screen or not.
    """
    t  = np.hstack((etime, angle, data))
    t = t.reshape((1, t.shape[0]))
    with open(file, 'a') as f:
       np.savetxt(f, t, fmt = '% 8.2f')
    if debug:
        print(f'INFO: Data saved at \n\t{file:s}')

    
if __name__ == '__main__':
    # This will only be used if the script is called as it is, not if it is used as a function
    # Assume the first port is the Arduino
    name_motor = list_ports()[1]
    # Assuming there is only one spectrometer, so taking the first element
    name_spectrometer = list_spectrometers()[0]
    # Define measurement variables
    folder = r'.'
    filename = 'goniomeasurement'
    angle_step = 20
    angle_max = 80
    # Define variables for the spectrometer measurements
    integration_time = 100
    n_spectra = 20
    plot = True
    
    gonio_measurement(name_motor, angle_max, angle_step, name_spectrometer, integration_time, n_spectra)