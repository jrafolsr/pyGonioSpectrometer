# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:13:35 2020

@author: JOANRR
"""
from os.path import join as pjoin
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import os 
import seaborn as sns


fcal = pjoin(os.path.dirname(__file__), 'calibration_files')
# Instrument response function, 160504 ML
path_IRF = pjoin(fcal, 'IRF_FlameBlueFiber_wLens2945K.txt')

# Eye response function
path_eye_response = pjoin(fcal, 'CIE1988photopic.txt')

ABS_CALFACTOR = 7.10e6 # Calibration from 15/07/2020
# ABS_CALFACTOR = 5.4178E6 # OLD Mattias' absolute numbers calibration factor, 170224 ML
PIXEL_SIZE = 4e-6 # m^2, Size of a McScience substrate pixel. Equipment is designed to be used with these substrates only.

def process_goniodata(file, correct_offset = True, current = None, plot = False, path_IRF = path_IRF, path_eye_response = path_eye_response, check_time_stability = False):
    """Read the file created by the setup Gonio spectrometer 3.0""" 
    # Load calibration files
    IRF = np.loadtxt(path_IRF, usecols = 1, unpack=True)
    
    # EyeResponse = loadmat(path_eye_response)['CIE1988Photopic']
    EyeResponse = np.loadtxt(path_eye_response)
    
    # First 3-row contains the starting time, the  integration time and the averaged times
    iTime = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
    IntTime, Nscans = np.loadtxt(file, skiprows = 1, max_rows = 2, unpack = True)

    # Get the rest, vTimes, Angles, Wavelengths, DarkSpectra and MeasSpectra
    data = np.loadtxt(file, skiprows = 3)
    vTimes = data[2:,0]
    # vTimes -= vTimes[0] # Zero the time vector to the first spectra adquisition
    Angles = data[2:,1]
    Wavelengths, DarkSpectra, MeasSpectra =  data[0,2:], data[[1],2:], data[2:,2:]
     
    Spectra = MeasSpectra - DarkSpectra
    
    
    #Check the time stability of the forward spectra, before the sorting
    if check_time_stability:
        fig, ax = plt.subplots()
        nangles = len(Angles)
        slicing = [0,int((nangles - 1)/2), -1]
        zero_spectra = Spectra[slicing, :] * IRF * PIXEL_SIZE / ABS_CALFACTOR / (IntTime/1000)
        zero_ri = np.trapz(zero_spectra,Wavelengths, axis =  1)
        ax.plot(vTimes[slicing], zero_ri/zero_ri[0] * 100, 'o-')
        ax.set_xlabel('Ellapsed time (s)')
        ax.set_ylabel('Rel. change radiant intensity (%)')
        ax.set_title('Time stability of forward angle')
    
    # Get the indices that will sort the data based on angles and sort the data
    isort =  Angles.argsort() 
    vTimes = vTimes[isort]
    Angles = Angles[isort]
    Spectra = Spectra[isort, :]
    

    
    
    if correct_offset:
        # Automatically find the symmetry axis assuming parabola (rought)
        angle_offset= find_symmetry(Wavelengths, Angles[2:-2], Spectra[2:-2,:], plot = plot)
        Angles -=  angle_offset
    else:
        print('INFO: No correction for the zero-angle offset is applied.')
        
    # Processing the data using all the corrections needed
    # The spectral radiant intensity [W sr-1 nm-1]
    SpecRadInt = Spectra * IRF * PIXEL_SIZE / ABS_CALFACTOR / (IntTime/1000)
    
    # Cut innecessary wavelengths assuming the range 450 - 800 nm to more than enough
    filt_by_wl = (Wavelengths >=449) & (Wavelengths <= 801)
    Wavelengths = Wavelengths[filt_by_wl]
    SpecRadInt = SpecRadInt[:, filt_by_wl]
    
    
    # The spectral luminous intensity [lm sr-1 nm-1]
    # Calculate first the photopic eye response for the Wavelength vector
    photopic_eye_response = (683.002 * np.interp(Wavelengths, EyeResponse[:,0], EyeResponse[:,1])).reshape(1, len(Wavelengths))
    SpecLumInt = SpecRadInt * photopic_eye_response  
    
    # Luminous intensity, the integral of contributions from all wavelengths [lm sr-1 = cd]
#     print(SpecLumInt.shape, Wavelengths.shape)
    LumInt = np.trapz(SpecLumInt,Wavelengths, axis =  1)
    L0 =  np.interp(0.0, Angles, LumInt)
    LumIntNorm = LumInt / L0   # Normalization for the 0 deg interpolated  
         
    # Luminance, Projected area correction that gives Luminance [cd m-2]
    Luminance = LumInt / np.cos(Angles * np.pi / 180.0) / PIXEL_SIZE
    # Normalization for the 0 deg measurements
    LuminanceNorm = LumIntNorm / np.cos(Angles * np.pi / 180.0)
    
    # Forward direction luminance [cd m-2]
    Fw_Luminance = np.interp(0.0, Angles, Luminance)
    # Forward direction current efficacy [cd/A]
    Fw_Current_eff = (Fw_Luminance * PIXEL_SIZE ) / current if current != None else np.nan
    
    # EQE calculation
    EQE = eqe_calculator(Wavelengths,Angles, SpecRadInt, current) if current != None else np.nan
    # Lambertian factor calculation
    Lamb_Corr_Factor = lambertian_correction_factor(Angles, LumIntNorm, plot  = plot)
    text = f'# IRF file: {path_IRF:s}\n'
    text += f'# ABS_CALFACTOR: {ABS_CALFACTOR:4.2e}\n'
    text += f'# Angle_offset = {angle_offset:.2f}°\n' if correct_offset else '# Angle_offset = nan\n'
    text += f'# I = {current*1000:.2f} mA\n' if current != None else '# Current = nan mA\n'
    text += f'# L0 =  {Fw_Luminance:.2f} cd m-2\n'
    text += f'# EQE =  {EQE:.2f} %\n'
    text += f'# CE =  {Fw_Current_eff:.2f} cd A-1\n'
    text += f'# Lambertian Corr. Factor =  {Lamb_Corr_Factor:.4f}\n'
    text += str(iTime) + '# Zero timestamp\n'
    print(text)
    
    # Saving the integrated part of the data
    fintegrated = file[:-4] + '.integrated'
    with open(fintegrated, 'w') as f:
        f.write(text)
    tdata = np.vstack((vTimes, Angles, LumInt, LumIntNorm, Luminance, LuminanceNorm)).T
                      
    with open(fintegrated, 'a') as f:
        header = 'eTime\t Angles \t LuminousIntensity \t NormLuminousIntensity \t Luminance\t NormLuminance\n'
        header += 's \t ° \t cd \t a.u. \t  cd m-2\t a.u.'
        np.savetxt(f, tdata, fmt = '% 10.6g\t', header = header)
    
    # Saving the spectral part of the data
    fradiometric = file[:-4] + '.sri'
    n = len(Angles)
    with open(fradiometric, 'w') as f:
        tdata = np.hstack((np.nan, vTimes))
        np.savetxt(f, tdata.reshape((1, n+1)), fmt = '%10.2f', header = 'ellapsedTime (s)', delimiter='\t')
        tdata = np.hstack((np.nan, Angles))
        np.savetxt(f, tdata.reshape((1, n+1)), fmt = '%10.2f', header = 'angles (°)',delimiter='\t')
        tdata = np.vstack((Wavelengths.reshape(1, len(Wavelengths)), SpecRadInt)).T
        fmt = ['%10.2f'] + ['%10.6e' for i in range(n)]
        np.savetxt(f, tdata, fmt = fmt, header = 'Wavelengths\t AngularSpectralRadiantIntensity', delimiter='\t')
        
    # Saving the spectral part of the data
    fphotometric = file[:-4] + '.sli'
    n = len(Angles)
    with open(fphotometric, 'w') as f:
        tdata = np.hstack((np.nan, vTimes))
        np.savetxt(f, tdata.reshape((1, n+1)), fmt = '%10.2f', header = 'ellapsedTime (s)', delimiter='\t')
        tdata = np.hstack((np.nan, Angles))
        np.savetxt(f, tdata.reshape((1, n+1)), fmt = '%10.2f', header = 'angles (°)',delimiter='\t')
        tdata = np.vstack((Wavelengths.reshape(1, len(Wavelengths)), SpecLumInt)).T
        fmt = ['%10.2f'] + ['%10.6e' for i in range(n)]
        np.savetxt(f, tdata, fmt = fmt, header = 'Wavelengths\t AngularSpectralLuminousIntensity', delimiter='\t')
    
    
        
    if plot:
        fig, [ax,ax1] = plt.subplots(ncols=2, figsize = (8,4))
        ff_wl = (Wavelengths >= 380) &( Wavelengths <= 780)
        ax.pcolormesh(Angles,Wavelengths[ff_wl],SpecRadInt[:, ff_wl].T, cmap = 'inferno')
        ax1.pcolormesh(Angles,Wavelengths[ff_wl],SpecLumInt[:,ff_wl].T, cmap = 'inferno')
        ax.set_ylabel('Wavelength (nm)')
        ax1.set_ylabel('Wavelength (nm)')
        ax.set_title('Spectral Radiant Intensity\n($W\cdot sr^{-1}\cdot nm^{-1}$)')
        ax1.set_title('Spectral Luminous Intensity\n($lm\cdot sr^{-1}\cdot nm^{-1}$)')
        ax.set_xlabel('Angle (°)')
        ax1.set_xlabel('Angle (°)')
        ax.set_xlim(-90, 90)
        ax1.set_xlim(-90,90)
    
    return iTime, IntTime, Nscans, Angles, Wavelengths,\
            SpecRadInt,SpecLumInt,LumInt, LumIntNorm, Luminance,LuminanceNorm,\
            Fw_Luminance,Fw_Current_eff, EQE


def process_L0(files, t0 = None, path_IRF = path_IRF, path_eye_response = path_eye_response, plot = False):
    """Read the L0-files L0 created by the setup Gonio spectrometer 3.0""" 
    # Load calibration files
    IRF = np.loadtxt(path_IRF, usecols = 1, unpack=True)
    # EyeResponse = loadmat(path_eye_response)['CIE1988Photopic']
    EyeResponse = np.loadtxt(path_eye_response)
    
    vtimes = np.array([], dtype = np.datetime64)
    vluminances = np.array([])
    
    if t0 is None:
        print('INFO: No t0 is provided, t0 taken from the first input file.')
        t0 = np.loadtxt(files[0], max_rows = 1, dtype = np.datetime64)
        
    for file in files:
        t1 = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
        # Get the rest, vTimes, Angles, Wavelengths, DarkSpectra and MeasSpectra
        data = np.loadtxt(file, skiprows = 3)
        aTimes = data[2:,0]
        integration_times = data[2:,1]
        integration_times = integration_times.reshape((len(integration_times), 1))
        Wavelengths, DarkSpectra, MeasSpectra =  data[0,2:], data[[1],2:], data[2:,2:]
        Spectra = MeasSpectra - DarkSpectra
        
        # The spectral radiant intensity [W sr-1 nm-1]
        SpecRadInt = Spectra * IRF * PIXEL_SIZE / ABS_CALFACTOR / (integration_times/1000)
        
        # Cut innecessary wavelengths assuming the range 450 - 800 nm to more than enough
        filt_by_wl = (Wavelengths >=449) & (Wavelengths <= 801)
        Wavelengths = Wavelengths[filt_by_wl]
        SpecRadInt = SpecRadInt[:, filt_by_wl]
        
        # The spectral luminous intensity [lm sr-1 nm-1]
        # Calculate first the photopic eye response for the Wavelength vector
        photopic_eye_response = (683.002 * np.interp(Wavelengths, EyeResponse[:,0], EyeResponse[:,1])).reshape(1, len(Wavelengths))
        SpecLumInt = SpecRadInt * photopic_eye_response
        
        Luminance = np.trapz(SpecLumInt, Wavelengths, axis =  1) / PIXEL_SIZE
        rtimes = [t1 + int(s*1e6) for s in aTimes]
        vtimes = np.concatenate([vtimes,rtimes])
        vluminances = np.concatenate([vluminances, Luminance])
    
    vtimes = np.float64(vtimes - t0)/1e6
    
    data = np.vstack([vtimes, vluminances]).T
    header = str(t0) + '\nRel.Time(s) \t  Luminance(cd/m2)'
    folder = os.path.dirname(files[0])
    np.savetxt(pjoin(folder, 'time-luminance.dat'), data, header = header, fmt = '%10.6f')
    return vtimes, vluminances

def plot_spectral_evolution(file, t0 = None, N = 10, skiprows = 0):
    """Read the L0-files L0 and plot its spectra radiant intensity""" 
    # Load calibration files
    IRF = np.loadtxt(path_IRF, usecols = 1, unpack=True)
    
    
    if t0 is None:
        print('INFO: No t0 is provided, t0 taken from the first input file.')
        t0 = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
        

    t1 = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
    # Get the rest, vTimes, Angles, Wavelengths, DarkSpectra and MeasSpectra
    data = np.loadtxt(file, skiprows = 3)
    aTimes = data[2:,0]
    integration_times = data[2:,1]
    integration_times = integration_times.reshape((len(integration_times), 1))
    Wavelengths, DarkSpectra, MeasSpectra =  data[0,2:], data[[1],2:], data[2:,2:]
    Spectra = MeasSpectra - DarkSpectra
    
    # The spectral radiant intensity [W sr-1 nm-1]
    SpecRadInt = Spectra * IRF * PIXEL_SIZE / ABS_CALFACTOR / (integration_times/1000)
    
    # Cut innecessary wavelengths assuming the range 450 - 800 nm to more than enough
    filt_by_wl = (Wavelengths >=450) & (Wavelengths <= 800)
    Wavelengths = Wavelengths[filt_by_wl]
    SpecRadInt = SpecRadInt[:, filt_by_wl]
    
    NormSpecRadInt = SpecRadInt / SpecRadInt.max(axis = -1, keepdims=True)
    
    rtimes = [t1 + int(s*1e6) for s in aTimes]

    rtimes = np.float64(rtimes - t0)/1e6
    
    with sns.color_palette("coolwarm", N):
        fig, ax = plt.subplots()
        
        for i in range(skiprows, N + skiprows):
            ax.plot(Wavelengths, NormSpecRadInt[i,:], label = f'{rtimes[i]:.1f}s')
        
        ax.legend()
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Norm. EL (a.u.)')
        
    return rtimes, NormSpecRadInt

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
        ax.set_title(f'Symmetry axis = {angle_offset:.2f}°')
        x = np.arange(-90,95,5)
        ax.plot(angles, Ie, 'o')
        ax.plot(x, pol2_sym(x, *popt))
        ax.axvline(angle_offset,  c= 'black', ls = '--')

    print('INFO: Angle offset correction applied. Offset = {:.2f}°'.format(angle_offset))
    
    return angle_offset


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
        ax.plot(angle*180/np.pi, Iph,'x')
        ax.plot(iangle*180/np.pi, np.cos(iangle))
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
    
    # Calculate the total number of photons per second
    Nphotons = total_energy / h / c * 2* np.pi
    
    # Number of electrons per second
    Nelectrons = current / e
    
#     print(f'N_photons = {Nphotons:.2e}, N_electrons = {Nelectrons:.2e}')
    EQE = Nphotons /  Nelectrons * 100

#     print(f'EQE =  {EQE:.2f} %')
    return EQE
