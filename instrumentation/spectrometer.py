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
    Lists all the spectrometers availables in the system.
    
    Returns
    -------
    ldevices: list
        A list with all the devices available. Each element of the list will be a seabreeze.spectrometers.Spectrometer class. If none is available, it returns an empty list.
        
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
    
    Parameters
    ----------
    device : seabreeze.spectrometers.Spectrometer class
        You can get it from list_spectrometers().
    integration_time : int or float, optional
        Sets the integration time in ms. The default is 1.
    n_spectra : int, optional 
        Number of spectra that will be averaged when using the method get_averaged_spectra(). The default is 1.
    
    Attributes
    ----------
    spec : seabreeze.spectrometers.Spectrometer class
    integration_time : int
    n_spectra : int
    wavelengths : numpy.array
    background : numpy.array
    intensities : numpy.array
    correct_nonlinearity : bool
    correct_dark_counts : bool
    saturation_counts : int  
    
    Methods
    -------
    config(integration_time, n_spectra=1, correct_nonlinearity=True, correct_dark_counts=False)
        Configures the measurement once initialized. 
    get_wavelengths()
        Returns an array with the wavelengths.
    get_intensities()
        Returns an array with the intensities.
    get_averaged_intensities()
        Returns an array with the averaged intensities over n_spectra times.
    adjust_integration_time(max_time=5000,lower_limit=10000,upper_limit=58000,noise_level=2700)
    Automatically adjusts the integration time and number of spectra.
        
    """
    def __init__(self, device, integration_time = 1, n_spectra = 1):
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
        self.saturation_counts = int(self.spec.max_intensity)
        
    def config(self, integration_time, n_spectra = 1, correct_nonlinearity = True, correct_dark_counts = False):
        """
        Configures the measurement once initialized. Initialy meant to change only the integration time, although some other configuration parameters can be changed.
        
        Parameters
        ----------
        integration_time : int or float
            Sets the integration time in ms.
        n_spectra : int, optional 
            Number of spectra that will be averaged when using the method get_averaged_spectra(). The default is 1.
        correct_nonlinearity : boolean, optional
            Whether to apply or no the built-in non-linearity correction when getting the intensities. The default is True.
        correct_dark_counts : boolean, optional
            Whether to apply or no the built-in dark counts correction when getting the intensities. The default is False.
        
        """
        
        self.integration_time = integration_time * 1000
        self.spec.integration_time_micros(self.integration_time)
        self.n_spectra = n_spectra 
        self.correct_nonlinearity = correct_nonlinearity
        self.correct_dark_counts = correct_dark_counts
        
    def get_wavelengths(self):
        """
       Returns an array with the wavelengths.
        
        Returns
        -------
        wavelengths : 1D numpy.array
            Vector with the wavelengths from the spectrometer
        
        """
        if self.wavelengths is None:
            self.wavelengths = self.spec.wavelengths()
        return self.wavelengths[20:] # The first 20 pixels are not valid
    
    def get_intensities(self):
        """
        Returns an array with the intensities.
        ** OBS! The n_spectra does not apply here. **
        
        Returns
        -------
        intensities : 1D numpy.array
            Vector with the intensities from the spectrometer
            
        """
        
        try:
            self.intensities = self.spec.intensities(self.correct_dark_counts, self.correct_nonlinearity)
        except Exception as e:
            print(e)
            self.spec.close()
            
        return self.intensities[20:] # The first 20 pixels are not valid
    
    def get_averaged_intensities(self):
        """
        Returns an array with the averaged intensities over n_spectra times, set by the initial configuration or the config() method
        
        Returns
        -------
        intensities : 1D numpy.array
            Vector with the averaged intensities from the spectrometer

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
        
    def adjust_integration_time(self, max_time = 5000,\
                                lower_limit = 10000, upper_limit = 58000,
                                noise_level = 2700):
        """
        Automatically adjusts the integration time and number of spectra to an optimal value.
        At the moment, the max n_spectra is 20. In the future it will be adjustable.

        Parameters
        ----------
        max_time : int, optional
            Maximum time allowed in ms, it corresponds to the integration_time * number_of_spectra. The default is 5000.
        lower_limit : int, optional
            Lower count limit before increasing the integration time. The default is 10000.
        upper_limit : TYPE, optional
            Upper count limit before decreasing the integration time.. The default is 58000.
        noise_level : TYPE, optional
            Noise level counts. The default is 2700.

        Returns
        -------
        integration_time : int
            Updated integration time.
        n_spectra : int
            Update number of spectra taken.

        """
        # OBS! Remember I am mixing us with ms. The config function should recieved ms,
        # but the self.integration_time property is stored in us (should fix this mess...)
        
        max_counts = self.get_intensities().max()
        integration_time = self.integration_time / 1000 # To ms
        n_spectra = self.n_spectra
        
        if max_counts <= lower_limit:

            integration_time = min(int(upper_limit / (max_counts - noise_level) * integration_time), max_time)
            
            if integration_time * n_spectra > max_time:
                n_spectra = max(1, int(max_time / integration_time))
                
            print('INFO: Lower limit reached. Adquisition increased to: {:.0f} ms x {:d}'.format(integration_time, n_spectra))

        elif max_counts >= 0.9 * self.saturation_counts:

            while max_counts >= 0.9 * self.saturation_counts and integration_time > 1:
                integration_time *= 0.9
                self.config(max(integration_time, 1))
                max_counts = self.get_intensities().max()
                
            integration_time = max(int(integration_time), 10)
            
            n_spectra = min(20, int(max_time / integration_time))
            
            print('INFO: Upper limit reached. Adquisition decreased to: {:.0f} ms x {:.0f}'.format(integration_time, n_spectra))
            
        else:
            pass
            print('INFO: The aquisition parameters are ok. Aquisition: {:.0f} ms x {:.0f}'.format(integration_time, n_spectra))
            
        self.config(integration_time, n_spectra)
            
        return integration_time, n_spectra
        
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