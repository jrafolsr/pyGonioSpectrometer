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
from scipy.interpolate import interp1d
from pyGonioSpectrometer.data_processing.calibration import path_IRF, path_eye_response, WL_MAX, WL_MIN, PIXEL_SIZE, ABS_CALFACTOR


def process_goniodata(file, angle_offset = 0.0, current = None, plot = False, path_IRF = path_IRF, path_eye_response = path_eye_response, correct_time_drift = False):
    """
    Reads and process the file generated by the Goniospectrometer software, version 3.0. Saves the processed dat in the files filename.integrated, filename.sri and filename.sli.
    
    Parameters:
    ----------
    file: str or path
        File path to the file generated by the Goniospectrometer software. It should have the following structure: the first three rows whould contain the timestamp at the beginning of the measurement (fmt "%Y-%m-%dT%H:%M:%S.%f"), the integration time in ms and the number of spectra averaged. After that it should a matrix data organized as:
                data[:,0]: Time ellapsed in s relative to the timestamp
                data[:,1]: Angle at which the spectra was taken
                data[0,2:]: Wavelength vector
                data[1,2:]: Intensity vector for the background spectra
                data[2:, 2:]: matrix containing the intensities for each angle (rows) and wavelength (columns)
    
    angle_offset: str or float/int, optional
        Whether to substract and offset to the angle data or not. If a 'auto' is passed, the offset is automatically derived form the symmetry axis between the left and right hemisphere of the angle resolved data, if a float or integer is passed, that value is consider to be the offset to be substracted. The default is 0.0.
    
    current : float or int, optional
        The current used to drive the device in A. Only used to calculate the efficiencies (CE and EQE). The default is None (so the CE adn EQE are not defined).

    plot: boolean, optional
        Whether to plot or not a mini-report of the processed file. The default is False.
        
    path_IRF: str or path, optional
        Path to the file containing the Instrument Response Function to correct the spectra. Should contain two columns, one with the wl and one with the counts. The wavelength vector should match the one correspinding to the experimental data. The default is defined with the path_IRF variable when loading the module.
        
    path_eye_response: str or path, optional
        Path to the file containing the Eye Response Function to correct the spectra. Should contain two columns, one with the wl and one with the counts. The program interpolates the vector to the correct values.
        
    correct_time_drift: boolean, optional
        If set to True, the radiant intensity will be corrected for the temporal drift between the first measurement and the last one, based on a linear interpolation using the 3 values taken at 0°. The default is False.
        
    Returns: (could nbe improved, a bit messy)
        iTime, IntTime, Nscans, L0, CE, EQE, Angles, Wavelengths, SpecRadInt,SpecLumInt,LumInt, LumIntNorm, Luminance,LuminanceNorm
        
    """ 
    # Load calibration files
    IRF = np.loadtxt(path_IRF, usecols = 1, unpack=True)
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
    
    # Cut innecessary (and noisy) wavelengths, taking only the range define by (WL_MIN, WL_MAX)
    filt_by_wl = (Wavelengths >= WL_MIN) & (Wavelengths <= WL_MAX)
    Wavelengths = Wavelengths[filt_by_wl]
    Spectra = Spectra[:, filt_by_wl]
    IRF = IRF[filt_by_wl]
    
    # Check the time stability of the forward spectra, before the sorting
    if correct_time_drift:
        nangles = len(Angles)
        # Zero angles 
        slicing = [0, int((nangles - 1)/2), -1]
        # Apply the IRF
        zero_spectra = Spectra[slicing, :] * IRF
        # Get the radiant intensity 
        zero_ri = np.trapz(zero_spectra, Wavelengths, axis =  1)
        # Get the interpolation function from the interp1d function
        fdrift = interp1d(vTimes[slicing], zero_ri/zero_ri[0], kind = 'linear')
        
        uncorrected = np.trapz(Spectra*IRF, Wavelengths, axis =  1) / zero_ri[0] *100
        
        # Correct the spectra using the interpolation function
        Spectra /= fdrift(vTimes).reshape((len(vTimes), 1))
        
        corrected = np.trapz(Spectra*IRF, Wavelengths, axis =  1) / zero_ri[0] *100
        
        print('INFO: Data has been corrected by the time drift derived from the 0° data')
        # Plot the report if asked
        if plot:
            fig, ax = plt.subplots()
            ax.plot(vTimes[slicing], zero_ri/zero_ri[0] * 100, 'oC0', label  = 'I0 evolution')
            ax.plot(vTimes, fdrift(vTimes) * 100, '--C0', label = 'interpolation')
            ax.plot(vTimes, uncorrected, 'xC0', label = 'uncorrected')
            ax.plot(vTimes, corrected, '.C2', label = 'corrected')
            ax.set_xlabel('Ellapsed time (s)')
            ax.set_ylabel('Rel. change radiant intensity (%)')
            ax.set_title('Time stability of forward angle')
            ax.legend()
            fig.savefig(file[:-4] + '_time-drift.png', bbox_inches = 'tight')
         
    # Get the indices that will sort the data based on angles and sort the data
    isort =  Angles.argsort() 
    vTimes = vTimes[isort]
    Angles = Angles[isort]
    Spectra = Spectra[isort, :]
    
    # Offset correction if it applies (automatic, fixed or None)
    if angle_offset == 'auto':
        # Automatically find the symmetry axis assuming parabola (rought)
        angle_offset = find_symmetry(Wavelengths, Angles[2:-2], Spectra[2:-2,:], plot = plot)
            
    elif  isinstance(angle_offset, int) or isinstance(angle_offset, float):       
        print(f'INFO: Manual angle offset correction applied. Offset = {angle_offset:.2f}°')
    else:
        raise Exception("The argument angle_offset must be either an integer, a float or 'auto'.")
        
    Angles -=  angle_offset
        
    # Processing the data using all the corrections needed
    # The spectral radiant intensity [W sr-1 nm-1]
    SpecRadInt = Spectra * IRF * PIXEL_SIZE / ABS_CALFACTOR / (IntTime/1000)
     
    # The spectral luminous intensity [lm sr-1 nm-1]
    # Calculate first the photopic eye response for the data Wavelength vector
    photopic_eye_response = (683.002 * np.interp(Wavelengths, EyeResponse[:,0], EyeResponse[:,1])).reshape(1, len(Wavelengths))
    SpecLumInt = SpecRadInt * photopic_eye_response  
    
    # Luminous intensity, the integral of contributions from all wavelengths [lm sr-1 = cd]
#     print(SpecLumInt.shape, Wavelengths.shape)
    LumInt = np.trapz(SpecLumInt, Wavelengths, axis =  1)
    LumIntNorm = LumInt / np.interp(0.0, Angles, LumInt)   # Normalization for the 0 deg interpolated  
         
    # Luminance, Projected area correction that gives Luminance [cd m-2]
    Luminance = LumInt / np.cos(Angles * np.pi / 180.0) / PIXEL_SIZE
    # Norm. luminance,pProjected area correction
    LuminanceNorm = LumIntNorm / np.cos(Angles * np.pi / 180.0)
    
    # Forward direction luminance [cd m-2]
    L0 = np.interp(0.0, Angles, Luminance) # Interpolation could be improve to a quadratic kind
    # Forward direction current efficacy [cd/A]
    CE = (L0 * PIXEL_SIZE ) / current if current != None else np.nan
    
    # EQE calculation
    EQE = eqe_calculator(Wavelengths, Angles, SpecRadInt, current) if current != None else np.nan
    
    # Lambertian factor calculation
    Lamb_Corr_Factor = lambertian_correction_factor(Angles, LumIntNorm, plot  = plot, file_id = file[:-4])
    
    text = f'# IRF_file: {path_IRF:s}\n'
    text += f'# correct_time_drift: {correct_time_drift}\n'
    text += f'# ABS_CALFACTOR: {ABS_CALFACTOR:4.2e}\n'
    text += f'# angle_offset = {angle_offset:.2f}°\n'
    text += f'# I = {current*1000:.2f} mA\n' if current != None else '# Current = nan mA\n'
    text += f'# L0 =  {L0:.2f} cd m-2\n'
    text += f'# EQE =  {EQE:.2f} %\n'
    text += f'# CE =  {CE:.2f} cd A-1\n'
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
    
    # I do not need to save the photometric spectral part of the data, do I?
    # # Saving the spectral part of the data
    # fphotometric = file[:-4] + '.sli'
    # n = len(Angles)
    # with open(fphotometric, 'w') as f:
    #     tdata = np.hstack((np.nan, vTimes))
    #     np.savetxt(f, tdata.reshape((1, n+1)), fmt = '%10.2f', header = 'ellapsedTime (s)', delimiter='\t')
    #     tdata = np.hstack((np.nan, Angles))
    #     np.savetxt(f, tdata.reshape((1, n+1)), fmt = '%10.2f', header = 'angles (°)',delimiter='\t')
    #     tdata = np.vstack((Wavelengths.reshape(1, len(Wavelengths)), SpecLumInt)).T
    #     fmt = ['%10.2f'] + ['%10.6e' for i in range(n)]
    #     np.savetxt(f, tdata, fmt = fmt, header = 'Wavelengths\t AngularSpectralLuminousIntensity', delimiter='\t')
    
    
        
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
        fig.savefig(file[:-4] + '_map.png', bbox_inches = 'tight')
    
    return iTime, IntTime, Nscans, L0,CE, EQE,\
            Angles, Wavelengths, SpecRadInt,SpecLumInt, \
            LumInt, LumIntNorm, Luminance, LuminanceNorm
            


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
   
    
    for i, file in enumerate(files):
        t1 = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
        # Get the rest, vTimes, Angles, Wavelengths, DarkSpectra and MeasSpectra
        data = np.loadtxt(file, skiprows = 3)
        aTimes = data[2:,0]
        
        # Distinguish the reading mode between a L0 (forward luminance file) and the gonio file (where I take only the 0 angle measurements (0, middle and last))
        if file[-6:] == 'L0.dat':
            integration_times = data[2:,1]
            integration_times = integration_times.reshape((len(integration_times), 1))
            Wavelengths, DarkSpectra, MeasSpectra =  data[0,2:], data[[1],2:], data[2:,2:]
        else:
            integration_times = np.loadtxt(file, skiprows = 1, max_rows = 1)
            Wavelengths, DarkSpectra, MeasSpectra =  data[0,2:], data[[1],2:], data[2:,2:]
            filter_vector = [0, int((aTimes.shape[0]-1)/2), -1]
            MeasSpectra = MeasSpectra[filter_vector, :]
            aTimes = aTimes[filter_vector]
        

        Spectra = MeasSpectra - DarkSpectra
        
        # The spectral radiant intensity [W sr-1 nm-1]
        SpecRadInt = Spectra * IRF * PIXEL_SIZE / ABS_CALFACTOR / (integration_times/1000)
        
        # Cut innecessary wavelengths assuming the range 450 - 800 nm to more than enough
        filt_by_wl = (Wavelengths >=WL_MIN) & (Wavelengths <= WL_MAX)
        Wavelengths = Wavelengths[filt_by_wl]
        SpecRadInt = SpecRadInt[:, filt_by_wl]
        
        if i == 0:
            data_to_save = Wavelengths.reshape((1, len(Wavelengths)))
        data_to_save = np.concatenate([data_to_save,SpecRadInt])
        
        # The spectral luminous intensity [lm sr-1 nm-1]
        # Calculate first the photopic eye response for the Wavelength vector
        photopic_eye_response = (683.002 * np.interp(Wavelengths, EyeResponse[:,0], EyeResponse[:,1])).reshape(1, len(Wavelengths))
        SpecLumInt = SpecRadInt * photopic_eye_response
        
        Luminance = np.trapz(SpecLumInt, Wavelengths, axis =  1) / PIXEL_SIZE
        rtimes = [t1 + int(s*1e6) for s in aTimes]
        vtimes = np.concatenate([vtimes,rtimes])
        vluminances = np.concatenate([vluminances, Luminance])
    
    vtimes = np.float64(vtimes - t0)/1e6
    vtimes = np.concatenate([[np.nan], vtimes]) 
    vluminances =  np.concatenate([[np.nan], vluminances]) 
    
    data = np.vstack([vtimes, vluminances])
       
    data_to_save = np.hstack([data.T, data_to_save])
    
    header = str(t0) 
    header += '\nFirst row correspond to the wavelength vector for the Spectral Radiant Intensity'
    header += '\nRel.Time(s) \t  Luminance(cd/m2) \t SRI(W/nm/sr)'
    
    folder = os.path.dirname(files[0])
    
    np.savetxt(pjoin(folder, 'spectral_evolution.evolution'), data_to_save, header = header, fmt = '%10.6g')
   
    return vtimes, vluminances

def plot_spectral_evolution(file, tmin = 0.0, tmax = np.inf, times = None, path = None, title = None):
    """Read the spectral evolution-file L0 and plot its spectra radiant intensity""" 

    # Get the rest, vTimes, Angles, Wavelengths, DarkSpectra and MeasSpectra
    data = np.loadtxt(file)
    rtimes = data[1:,0]
    # luminance = data[1:,1]   # I don't need it
    # The spectral radiant intensity [W sr-1 nm-1]
    Wavelengths =  data[0,2:] 
    SpecRadInt = data[1:,2:]   
    
    NormSpecRadInt = SpecRadInt / SpecRadInt.max(axis = -1, keepdims=True)
    
     # Filter out the times outside the limits
    ftime = (rtimes >= tmin) & (rtimes <=tmax)
    
    rtimes = rtimes[ftime]
    NormSpecRadInt = NormSpecRadInt[ftime, :]
    
    
    
    # Second filter for specific times
    ftime = []
    if times is not None:
        for t in times:
            ftime.append(np.abs(rtimes - t).argmin())
        rtimes = rtimes[ftime]
        NormSpecRadInt = NormSpecRadInt[ftime, :]
        
    N = len(rtimes)  
    
    with sns.color_palette("coolwarm", N):
        fig, ax = plt.subplots()
    
    lines= []
    
    for i in range(N):
        line, = ax.plot(Wavelengths, NormSpecRadInt[i,:], label = f'{rtimes[i]/60:.1f}')
        lines.append(line)
    
    if times is None:  
        lines = [lines[0], lines[int((N-1)/2)], lines[-1]]
        
    ax.legend(lines, [l.get_label() for l in lines], title = 'Time (min)')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Norm. EL (a.u.)')
    
    if title is None:
        ax.set_title('Spectral evolution')
    else:
        ax.set_title(title)
    if path is None:
        folder = os.path.dirname(file)
        path = pjoin(folder, f'spectral-evolution_times={tmin:.0f}-{tmax:.0f}s.png')
    
    fig.savefig(path, bbox_inches  = 'tight')
        
    return rtimes, NormSpecRadInt

def find_symmetry(wavelengths, angles, spectra, plot = False):
    """
    Uses the integrated radiant intensity (no need to be calibrated, from the spectra)  to calculate the symmetry point assuming a symmetric parabola as a reference (as a first approximation).
    
    Parameters
    ----------
    wavelengths: 1D np.array
        Array with the input wavelengths, [nm] should match the length of the first dimension of the sri.
        
    angles: 1D np.array
        Array with the input angles, should match the length of the second dimension of the sri.
    spectra: 2D np.array
        The array with the spectral radiant intensity in aribitray units [a.u.] (wavelength as rows) for each of the collected angles (as columns).
    
    Returns:
    --------
    angle_offset: float
        The offset by which the data should be corrected.
    """
    
    # Using the radiant intensity instead of the spectral radiant intensity (uncalibrated)
    # to calculate the symmetry point
    
    Ie = np.trapz(spectra,  wavelengths, axis = 1)
    Ie /= Ie.max()
   
#     # This ensures a symmetric parabola
    pol2_sym = lambda x,a,c, x0: a*(x-x0)**2 + c
    p0 = [-1e-6, 1, 0]
    popt, _ = curve_fit(pol2_sym,angles, Ie, p0 = p0)
       
    angle_offset = popt[2]
    
    if plot:
        fig, ax = plt.subplots()
        ax.set_title(f'Symmetry axis = {angle_offset:.2f}°')
        x = np.arange(-90, 95, 5)
        ax.plot(angles, Ie, 'o')
        ax.plot(x, pol2_sym(x, *popt))
        ax.axvline(angle_offset,  c= 'black', ls = '--')

    print(f'INFO: Automatic angle offset correction applied. Offset = {angle_offset:.2f}°')
    
    return angle_offset


def lambertian_correction_factor(angle, Iph, plot = False, file_id = 'lambertian_correction_factor'):
    """
    Calculates the Lambertian correction factor, see Lindh's work.
    
    angles: vector with the angles
    Iph: normalized luminous intensity
    """
    # First, artificially add the -90 and 90° values, with value zero and average the zeros (as there are three)
    N = int((len(angle) - 1) / 2)
    neg_a = slice(0,N-1)
    zero_a = slice(N-1, N+2)
    pos_a = slice(N+2,None)

    angle = (np.hstack([-90, angle[neg_a], angle[zero_a].mean(), angle[pos_a], 90])) * np.pi / 180

    Iph = np.hstack([0,Iph[neg_a], Iph[zero_a].mean(), Iph[pos_a],0])
    # It is assumeed that the angles vector is sorted
    
    # Interpolate the data (can be improved with quadratic interpolation) DONE
    iangle = np.linspace(-np.pi/2, np.pi/2, 90)
    f = interp1d(angle, Iph, 'quadratic')
    iIph = f(iangle)
       
    # Calculate the lambertian correction factor for the positive angles and negative angle
    # As a quick workaround, take the abs of the integrand and divide by two
    # OBS! Double check that.
    Lambertian_CorrFact = 2 * np.trapz(np.abs(iIph * np.sin(iangle)), iangle) / 2

    if plot:
        fig, ax = plt.subplots()
        ax.plot(angle * 180 / np.pi, Iph,'x')
        ax.plot(iangle * 180 / np.pi, np.cos(iangle))
        ax.set_xlabel('Angle (°)')
        ax.set_ylabel('Norm. luminous intensity (a.u.)')
        ax.set_title('Comparison with Lambertian emitter')
        fig.savefig(file_id + '_lambertianCF.png', bbox_inches = 'tight')
        
    return Lambertian_CorrFact

def eqe_calculator(wavelengths, angles, sri, current):
    """Calculates the EQE based on the input angular spectral radiant intensity, sri, without having to assume a Lambertian emission. Uses the whole -90° to 90° range to calculate it.
    
        Parameters
        ----------
        wavelengths: 1D np.array
            Array with the input wavelengths, [nm] should match the length of the first dimension of the sri.
            
        angles: 1D np.array
            Array with the input angles, should match the length of the second dimension of the sri.
        sri: 2D np.array
            The array with the spectral radiant intensity  [W sr-1 nm-1] (wavelength as rows) for each of the collected angles (as columns).

        Returns:
        --------
        EQE: float
            The EQE value.

    """ 
    
    h = 6.626e-34  # Planck constant [SI]
    c = 2.997e8   # Speed of light [SI]
    e = 1.602e-19  # Elementary charge [SI]
    
    # Reshape the wavelengths vector to (1,N), for the latter integration of the wavelengths vector
    
    N = len(wavelengths)
    wavelengths = wavelengths.reshape((1,N))
    
    # Add artificially the +/- 90°
    theta = (np.hstack([-90, angles, 90])) * np.pi / 180
    
    # Artificially adding  +/-90° sri to a value of zero and convert  and sorting 
    # OBS! Sorting should have been already done, consider removing
    sri = np.vstack([np.zeros((1,N)), sri, np.zeros((1,N))])
    
    # Calculting the total energy in 2 steps, first, integrate over wavelengths
    integral_wl = np.trapz(sri * wavelengths, wavelengths * 1E-9, axis = 1)
    
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