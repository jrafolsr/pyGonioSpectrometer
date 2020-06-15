# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:13:35 2020

@author: JOANRR
"""
from os.path import join as pjoin
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.io  import loadmat

fcal = r'.\calibration_files'

# Instrument response function, 160504 ML
temp = loadmat(pjoin(fcal, 'GonioSpec_IRF_160504'))
IRF = temp['GonioSpec_IRF_160504'][:,[1]]

# Eye response function
temp = loadmat(pjoin(fcal, 'CIE1988photopic.mat'))
EyeResponse = temp['CIE1988Photopic']

ABS_CALFACTOR = 5.4178e6 # Absolute numbers calibration factor, 170224 ML
PIXEL_SIZE = 4e-6 # m^2, Size of a McScience substrate pixel. Equipment is designed to be used with these substrates only.

def read_gonio_spec(file, correct_offset = True, current = None, plot = False):
    """Read the file created by the setup Gonio spectrometer 3.0"""

    # First 3-row contains the starting time, the  integration time and the averaged times
    iTime = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
    IntTime, Nscans = np.loadtxt(file, skiprows = 1, max_rows = 2, unpack = True)

    # Get the rest, vTimes, Angles, Wavelengths, DarkSpectra and MeasSpectra
    data = np.loadtxt(file, skiprows = 3)
    vTime = data[:,0]
    Angles = data[:,1]
    Wavelengths, DarkSpectra, MeasSpectra =  data[0,2:], data[[1],2:], data[2:,2:]
    Spectra = MeasSpectra - DarkSpectra
    # Sort the Data
    isort =  Angles.argsort() # Get the indices that will sort the data based on angles
    Angles = Angles[isort]
    Spectra = Spectra[isort, :]
    
    if correct_offset:
        # Automatically find the symmetry axis assuming parabola (rought)
        angle_offset= find_symmetry(Wavelengths, Angles, Spectra, plot = plot)
        Angles -=  angle_offset
    else:
        print('INFO: No correction for the zero-angle offset is applied.')
        
    # Processing the data using all the corrections needed
    # The spectral radiant intensity [W sr-1 nm-1]
    SpecRadInt = Spectra * IRF * PIXEL_SIZE / ABS_CALFACTOR / (IntTime/1000)
    
    # The spectral luminous intensity [lm sr-1 nm-1]
    SpecLumInt = SpecRadInt * photopic_eye_response(Wavelengths)     
    
    # Luminous intensity, the integral of contributions from all wavelengths [lm sr-1 = cd]
#     print(SpecLumInt.shape, Wavelengths.shape)
    LumInt = np.trapz(SpecLumInt,Wavelengths, axis =  1)
    L0 =  np.interp(0.0, Angles, LumInt)    
    LumIntNorm = LumInt / L0   # Normalization for the 0 deg interpolated  
#     LumIntNorm = LumInt / LumInt[zeros_ff].mean()      # Normalization for the 0 deg measurements 
#     print(LumIntNorm, LumInt / LumInt[zeros_ff].mean())
    
    # Luminance, Projected area correction that gives Luminance [cd m-2]
    Luminance = LumInt / np.cos(Angles * np.pi / 180.0) / PIXEL_SIZE
    # Normalization for the 0 deg measurements
    LuminanceNorm = LumIntNorm / np.cos(Angles * np.pi / 180.0)
    
    # Forward direction luminance (mean of all "zero" degree angles) [cd m-2]
    Fw_Luminance = L0
    # Fw_Luminance = Luminance[zeros_ff].mean()
    # Forward direction current efficacy [cd/A]
    Fw_Current_eff = (Fw_Luminance * PIXEL_SIZE ) / current if current != None else np.nan
    
    # EQE calculation
    EQE = eqe_calculator(Wavelengths,Angles, SpecRadInt, current) if current != None else np.nan
    # Lambertian factor calculation
    Lamb_Corr_Factor = lambertian_correction_factor(Angles, LumIntNorm, plot  = plot)

    print(f'L0 =  {Fw_Luminance:.2f} cd m-2')
    print(f'EQE =  {EQE:.2f} %')
    print(f'CE =  {Fw_Current_eff:.2f} cd A-1')
    print(f'Lambertian Corr. Factor =  {Lamb_Corr_Factor:.4f}')
    
    if plot:
        fig, [ax,ax1] = plt.subplots(ncols=2, figsize = (8,4))
        ff_wl = (Wavelengths >= 380) &( Wavelengths <= 780)
        ax.pcolormesh(Angles,Wavelengths[ff_wl],SpecRadInt[ff_wl,:], cmap = 'inferno')
        ax1.pcolormesh(Angles,Wavelengths[ff_wl],SpecLumInt[ff_wl,:], cmap = 'inferno')
        ax.set_ylabel('Wavelength (nm)')
        ax1.set_ylabel('Wavelength (nm)')
        ax.set_title('Spectral Radiant Intensity\n($W\cdot sr^{-1}\cdot nm^{-1}$)')
        ax1.set_title('Spectral Luminous Intensity\n($lm\cdot sr^{-1}\cdot nm^{-1}$)')
        ax.set_xlabel('Angle (°)')
        ax1.set_xlabel('Angle (°)')
        ax.set_xlim(-90, 90)
        ax1.set_xlim(-90,90)
    
    return IntTime, Nscans, Angles, Wavelengths,\
            SpecRadInt,SpecLumInt, LumIntNorm, Luminance,LuminanceNorm,\
            Fw_Luminance,Fw_Current_eff, EQE

def find_symmetry(wavelengths,angles,spectra, plot = False):
    # Using the radiant intensity instead of the spectral radiant intensity (uncalibrated)
    # to calculate the symmetry point
    Ie = np.array([np.trapz(row, wavelengths) for row in spectra])
    Ie /= Ie.max()
   
#     # This ensures a symmetric parabola
    pol2_sym = lambda x,a,c, x0: a*(x-x0)**2 + c
    p0 = [-1e-6, 1, 0]
    popt, _ = curve_fit(pol2_sym,angles, Ie, p0 = p0)
       
    angle_offset = popt[2]
    
    if plot:
        fig, ax = plt.subplots()
        x = np.arange(-90,95,5)
        ax.plot(angles, Ie, 'o')
        ax.plot(x, np.polyval(popt, x))
        ax.axvline(angle_offset,  c= 'black', ls = '--')

    print('INFO: Angle offset correction applied. Offset = {:.2f}°'.format(angle_offset))
    return angle_offset



def photopic_eye_response(wl):
    return (683.002 * np.interp(wl,EyeResponse[:,0], EyeResponse[:,1])).reshape(len(wl),1)

def lambertian_correction_factor(angle, Iph, plot = False):
    """Calculates the Lambertian correction factor, see Lindh's work
    angles: vector with the angles
    Iph: normalized luminous intensity
    """
    # First, artificially add the -90 and 90° values, with value zero
    angle = (np.hstack([-90, angle, 90])) * np.pi / 180
    Iph = np.hstack([0,Iph,0])
    # It is assumeed that the angles vector is sorted
    
    # Interpolate the data
    iangle = np.linspace(-np.pi/2, np.pi/2, 90)
    iIph = np.interp(iangle, angle, Iph)
       
    # Calculate the lambertian correction factor for the positive angles and negative angle
    # As a quick workaround, take the abs of the integrand and divide by two
    # OBS! Double check that.
    Lambertian_CorrFact= 2 * np.trapz(np.abs(iIph * np.sin(iangle)), iangle) / 2

    if plot:
        fig, ax = plt.subplots()
        ax.plot(angle, Iph,'x')
        ax.plot(iangle, np.cos(iangle))
        ax.set_xlabel('Angle (°)')
        ax.set_ylabel('Norm. luminous intensity (a.u.')
        ax.set_title('Comparison with Lambertian emitter')
        
    return Lambertian_CorrFact

def eqe_calculator(wl, angles, sr_intensity, current):
    # wl: wavelength[nm]
    # angles [°]
    # sr_intensity: spectral radian intensity [W sr-1 nm-1]
    
    h = 6.626e-34  # Planck constant
    c = 2.997e8;   # Speed of light
    e = 1.602e-19  # Elementary charge
    
    # Reshape the wl vector to (N,1), for the latter integration of the wl vector
    N = len(wl)
    wl = wl.reshape((N,))
    
    # Add artificially the +/- 90°
    theta = (np.hstack([-90, angles, 90])) * np.pi / 180
    
    # Artificially adding  to a value of zero and convert  and sorting 
    # OBS! Sorting should have been already done, consider removing
    sr_intensity = np.vstack([np.zeros((1,N)), sr_intensity,np.zeros((1,N))])
    
    # Calculting the total energy in 2 steps, first, integrate over wavelengths
    integral_wl = np.trapz(sr_intensity * wl.reshape(1,N), wl * 1E-9, axis = 1)
    
    # Second, integrate over angles, theta goes form 0 to 90° only, so we integrate both -90 to 0 and 0 to 90 and take the abs value of the integrand and divide by two, as a quick work around
    total_energy = 0.5 * np.trapz(np.abs(integral_wl * np.sin(theta)), theta)
    
    # Calculate the total number of photons
    Nphotons = total_energy / h / c * 2* np.pi
    
    # Number of electrons
    Nelectrons = current / e
    
#     print(f'N_photons = {Nphotons:.2e}, N_electrons = {Nelectrons:.2e}')
    EQE = Nphotons /  Nelectrons * 100

#     print(f'EQE =  {EQE:.2f} %')
    return EQE