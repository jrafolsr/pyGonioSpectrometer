# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 08:23:32 2020

@author: JOANRR

Containts a wrapper to call the Spectrometer class from the seabreeze module and perform measurements.
"""
from seabreeze.spectrometers import list_devices, Spectrometer
from time import sleep


def list_spectrometers():
    ldevices = list_devices()
    if ldevices == []:
        print('WARNING: No spectrometer found. Check the connections.')
    # else:
    #     # print('The available spectrometers are:')
    #     for device in ldevices:
    #         print(device)
    return ldevices
        

class SpectraMeasurement():
    """Wrapper to for the Spectrometer class that helps to perform measurements, from configurating the spectrometers to return spetra averages."""
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
        
    def config(self, integration_time, n_spectra = 1, correct_nonlinearity = True, correct_dark_counts = False):
        self.integration_time = integration_time * 1000
        self.spec.integration_time_micros(self.integration_time)
        self.n_spectra = n_spectra 
        self.correct_nonlinearity = correct_nonlinearity
        self.correct_dark_counts = correct_dark_counts
        
    def get_wavelengths(self):
        if self.wavelengths is None:
            self.wavelengths = self.spec.wavelengths()
        return self.wavelengths[20:] # The first 20 pixels are not valid
    
    def get_intensities(self):
        try:
            self.intensities = self.spec.intensities(self.correct_dark_counts, self.correct_nonlinearity)
        except Exception as e:
            print(e)
            self.spec.close()
            
        return self.intensities[20:] # The first 20 pixels are not valid
    
    def get_averaged_intensities(self):
        intensities = self.get_intensities()
        if self.n_spectra > 1:
            for i in range(self.n_spectra-1):
                intensities += self.get_intensities()
                sleep(0.010)
            intensities = intensities / self.n_spectra
        return intensities
    
   
    def set_background(self, bkg):
        self.background = bkg
        
    def get_spectra(self):
        if self.background is None:
            return self.get_wavelengths(), self.get_intensities()
        else:
            return self.get_wavelengths(), self.get_intensities() - self.background
        
    def close(self):
        self.spec.close()
    def open(self):
        self.spec.open()