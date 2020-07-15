# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 08:23:32 2020

@author: JOANRR

Containts a wrapper to call the Spectrometer class from the seabreeze module and perform measurements.
"""
from seabreeze.spectrometers import list_devices, Spectrometer
from time import sleep


def list_spectrometers():
    """
    This function will list all the spectrometers availables in the system.
    
    Returns:
    --------
    ldevices: a list with all the devices available. Each element of the list will be a seabreeze.spectrometers.Spectrometer class. If none is available, it returns an empty list.
    """
    ldevices = list_devices()
    if ldevices == []:
        print('WARNING: No spectrometer found. Check the connections.')
    # else:
    #     # print('The available spectrometers are:')
    #     for device in ldevices:
    #         print(device)
    return ldevices
        

class SpectraMeasurement():
    """
    Wrapper for the seabreeze.spectrometers.Spectrometer class that helps to perform measurements, from configurating the spectrometers to return spectra averages.
    
    Parameters:
    -----------
    device: seabreeze.spectrometers.Spectrometer class
        You can get it from list_spectrometers().
    integration_time: int or float, optional
        Sets the integration time in ms. The default is 1 ms.
    n_spectra: int, optional 
        Number of spectra that will be averaged when using the method get_averaged_spectra(). The default is 1.
    """
    def __init__(self, device, integration_time = 10, n_spectra = 1):
        """Initalize all values"""
        if device is None:
            raise Exception('ERROR: USB spectrometer port not provided')
        elif type(device) is str:
            self.spec = Spectrometer.from_serial_number(device)
        else:
            self.spec = Spectrometer(device)
        self.integration_time = integration_time * 1000 # From ms in the input to microseconds to the Specrometer class 
        self.spec.integration_time_micros(self.integration_time)
        self.wavelengths = None
        self.background = None
        self.intensities = None
        self.n_spectra = n_spectra
        self.correct_nonlinearity = True
        self.correct_dark_counts  = False
        
    def config(self, integration_time, n_spectra = 1, correct_nonlinearity = True, correct_dark_counts = False):
        """
        Method to further configure the measurement once initialized. Initialy though to change only the integration time, although some other configuration parameters can be changed.
        
        Parameters:
        -----------
        integration_time: int or float
            Sets the integration time in ms.
        n_spectra: int, optional 
            Number of spectra that will be averaged when using the method get_averaged_spectra(). The default is 1.
        correct_nonlinearity: boolean, optional
            Whether to apply or no the built-in non-linearity correction when getting the intensities. The default is True
        correct_dark_counts: boolean, optional
            Whether to apply or no the built-in dark counts correction when getting the intensities. The default is False
        """
        self.integration_time = integration_time * 1000
        self.spec.integration_time_micros(self.integration_time)
        self.n_spectra = n_spectra 
        self.correct_nonlinearity = correct_nonlinearity
        self.correct_dark_counts = correct_dark_counts
        
    def get_wavelengths(self):
        """
        Method that return the an array with the wavelengths
        
        Returns:
        --------
        wavelengths: np.array with the wavelengths from the spectrometer
        
        """
        if self.wavelengths is None:
            self.wavelengths = self.spec.wavelengths()
        return self.wavelengths[20:] # The first 20 pixels are not valid
    
    def get_intensities(self):
        """
        Method that return the an array with the intensities
        
        Returns:
        --------
        intensities: np.array with the intensities from the spectrometer
            OBS! The n_spectra does not apply here
        """
        
        try:
            self.intensities = self.spec.intensities(self.correct_dark_counts, self.correct_nonlinearity)
        except Exception as e:
            print(e)
            self.spec.close()
            
        return self.intensities[20:] # The first 20 pixels are not valid
    
    def get_averaged_intensities(self):
        """
        Method that return the an array with the averaged intensities over n_spectra times, set by the initial configuration or the config() method
        
        Returns:
        --------
        intensities: np.array with the averaged intensities from the spectrometer

        """
        
        intensities = self.get_intensities()
        
        if self.n_spectra > 1:
            for i in range(self.n_spectra-1):
                intensities += self.get_intensities()
                sleep(0.010)
            intensities = intensities / self.n_spectra
        return intensities
    
   
    def set_background(self, bkg):
        """
        DEPRECATED
        """
        
        self.background = bkg
        
    def get_spectra(self):
        """
        DEPRECATED
        """
        
        if self.background is None:
            return self.get_wavelengths(), self.get_intensities()
        else:
            return self.get_wavelengths(), self.get_intensities() - self.background
        
    def close(self):
        """
        Method that closes and frees the resource.
        """
        self.spec.close()
    def open(self):
        """
        Method that opens the resource.
        """
        self.spec.open()
        self.config(self.integration_time, self.n_spectra)