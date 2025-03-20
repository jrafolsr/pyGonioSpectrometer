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
#import winsound # Does not work in Linux
from pyGonioSpectrometer.instrumentation import RaspberryMotorController, SpectraMeasurement, list_spectrometers
from pyGonioSpectrometer.GonioSpec_init import WAIT_TIME
from pathlib import Path

# Parameters for the beeping, quite irrellevant, I'll consider to remove it
frequency = 2000  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second

# Fixed parameters
#WAIT_TIME = 1.5
N_WAVELENGTHS = 2028


def gonio_measurement(angle_max, angle_step,\
                      gonio,\
                      name_spectrometer, integration_time, n_spectra,\
                      filename = 'gonio_mesurement', folder = '.',\
                      disable_gonio = False, plot = True, header = ''):
    """
    Performs a complete measurement for the goniospectrometer setup, by taking spectra at every specified angle, for the whole forward hemisphere.
    
    Parameters
    ----------
    angle_max : int or float
        Maximum angle to scan with the goniometer in deg.
    angle_step : int or float
        Angular step that the motor will perform in deg. Preferably a divisor of angle_max and integer, the program does not check for this conditions to be fullfilled.
    gonio: RaspberryMotorController or ArduinoMotorControllerobject
        Object to controll the motor.
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
    # making sure that the angle step and the max angel are compatible
    _ , residu = divmod(angle_max, angle_step)
    if residu != 0.0:
            raise Exception("The angle step is not a divisor of the max. angle!")
    # Initializing flame objects.
    flame = None
    # Initalizing some variables
    current_angle = 0.0
        
    try:
        # Create the object instance taht will control the spectrometer, assuming it is the
#        first of the list
        flame = SpectraMeasurement(name_spectrometer, integration_time=integration_time, n_spectra= n_spectra)
        flame.open()
         
        n_angles = int(round(angle_max / angle_step, 0)) + 1
        n_steps = (n_angles - 1) * 2

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
            file_config = Path(__file__).parent / 'local-config.txt'
            if file_config.exists():
                with open(file_config) as fc:
                    for l in fc.readlines():
                        f.write('# ' + l)
                        f.write('\n') if not l.endswith('\n') else None
            f.write(header if header.startswith('#') else '# ' + header)
            f.write(f'# Timestamp at the beginning of the measurement: {itimestamp}\n')
            f.write(f'# Integration time in (ms): {integration_time:.0f}\n')
            f.write(f'# Number of spectra taken: {n_spectra:d}\n')
            f.write(f'# Time(s) Angle(deg) Wavelength(nm) -> (row 1) Counts() -> from row 2\n')
              
        # Take the dark spectra at zero
        wavelengths = flame.get_wavelengths()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, np.nan, wavelengths, path)
        # Saving the data in Mattias's scheme
#        data[:,0] = wavelengths
        
#        winsound.Beep(frequency, duration*2) # It will remind that the measureement starts
        for i in range(2):
            sleep(1)
        
        
        print('INFO: Taking dark spectra....')
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, np.nan, temp, path)
      
        if plot: plot_measurement(fig, ax, wavelengths, temp, 'dark')
        # Saving the data in Mattias's scheme
#        data[:,1] = temp
        sleep(0.25)
        
        # Open the shutter, assuming it is closed
        gonio.move_shutter()

        sleep(WAIT_TIME)
        
        # Take spectra at zero
        print(f'\rINFO Step # 0, moved  {0.0: >4.1f}°, position {0.0: >+5.1f}° |' + ' '* (n_steps + 3) + f'| {1/(n_steps + 3)*100:3.0f} % Done...', end =''  )
              
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, 0.0, temp, path)
        
        if plot: plot_measurement(fig, ax, wavelengths, temp, '0.0°')
        

        # Move to the last position
        
        out_angle = gonio.move_angle(-1.0 * angle_max)
#        # Initialize the error made in each step
#        error = (angle_max - out_angle)
        total = 0
        current_angle = -out_angle
        print(f'\rINFO: Step # 1, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |#' + ' '* (n_steps + 2) + f'| {2/(n_steps + 3)*100:3.0f} % Done...', end =''  )
              
        sleep(WAIT_TIME * 2) # Long enough time to make sure that it waits until the ed of the movement
        
        k = 0
        for k in range(n_steps):        
            # Get the whole spectra (wl and I)
            temp = flame.get_averaged_intensities()
            # Saving the data in the new scheme
            write_to_file(time()-start_time, current_angle, temp, path)
        
            # Check for any values higher than saturation
            if np.any(temp > 65535): print('WARNING: Some values are saturating. Consider lowering the integration time.')
            
            # Plot the data
            if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.1f}°')
            
            # Moving the gonio
            out_angle = gonio.move_angle(angle_step)
            
            total += abs(out_angle)

            current_angle += out_angle
            
            print(f'\rINFO: Step #{k+2:2d}, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |' + '#'* (k + 2) + ' '* (n_steps - k + 1) + f'| {(k+3)/(n_steps + 3)*100:3.0f} % Done...' , end =''  )
            
            sleep(WAIT_TIME)
        
        # Take last angle spectra
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, current_angle, temp, path)

        if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.1f}°')
        
        # Go back to zero
        back_angle = -1 * abs(current_angle)
        # Moving back the exact angle we moved to set everything to the initial position, so correct_drift = False
        out_angle = gonio.move_angle(back_angle, correct_drift = False)
        current_angle -= out_angle
        
        print(f'\rINFO: Step #{k+3:2d}, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |' + '#'* (n_steps + 3) + f'| {(k+4)/(n_steps + 3)*100:3.0f} % Done...' , end ='\n'  )
        
        # Wait longer time, as the angle is larger and take the spectra
        sleep(WAIT_TIME * 2)
        
        temp = flame.get_averaged_intensities()
        # Saving the data in the new scheme
        write_to_file(time()-start_time, current_angle, temp, path)
        
        
        if plot: plot_measurement(fig, ax, wavelengths, temp, f'{current_angle:.1f}°')
        
        # Closing shutter
        gonio.move_shutter()
        sleep(0.5)
        
        if disable_gonio: gonio.disable_gonio()
               
    
    except KeyboardInterrupt:
        print('INFO: The angle-scan has been cancelled by the user. Going back to 0°.')
        if gonio != None:
            # Go back to since the spectrogoniometer movement has been cancelled.
            back_angle = -1 * current_angle          
            out_angle = gonio.move_angle(back_angle)
            gonio.move_shutter()
            
    except Exception as e:
        print(e)
        print('INFO: Some error has ocurred during the angle-scan. Going back to 0°.')
        if gonio != None:
            # Go back to since some error has occurred during the gonio measurement
            back_angle = -1 * current_angle
            out_angle = gonio.move_angle(back_angle)
            gonio.move_shutter()
    
    finally:
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
    
    gonio_measurement(angle_max, angle_step, name_spectrometer, integration_time, n_spectra)