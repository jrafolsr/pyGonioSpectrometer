# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:21:52 2020

@author: JOANRR
"""
from pyGonioSpectrometer.data_processing import process_goniodata, process_L0, plot_spectral_evolution,\
    load_simdata, error_landscape, compare_data
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import shutil
from pyGonioSpectrometer.data_processing import __file__ as calibration_dir

plt.ioff()

calibration_dir = Path(calibration_dir).parent / 'calibration_files'
# Calibration files
calibrations ={'arduino_gonio_1': dict(path_IRF = calibration_dir / 'IRF_FlameBlueFiber_wLens2945K.txt',\
                                       abs_calfactor = 7.10e6),\
               'raspberry_gonio2':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference.txt',\
                                       abs_calfactor = 6.08E6),\
               'raspberry_gonio2_old':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens2945K.txt',\
                                       abs_calfactor = 1.44E6)}
              
# Output folders

processed_folder = 'processed-data'
fits_folder = 'fits'
plots_folder = 'plots'
raw_folder = 'raw-data'
    
def create_structure(folder):
    # Create folders with output of processed data, plots and fits
    for subfolder in [raw_folder, processed_folder, fits_folder, plots_folder]:
        t = Path(folder) /  subfolder
        t.mkdir() if not t.exists() else None  
    
    
def process_gonio_time_series(folder, current = None, angle_offset = 0.0, calibration = 'raspberry_gonio2', correct_time_drift = True):
    """
    Processes in an automated way all the time-series gonio data. For that, it restructures the folder into four subfolders: raw-data, processed-data, fits and plots.
    
    Parameters
    ----------
    folder : stri or Path
        folder that contains the data.
    current : float, optional
        Current in mA. The default is None.
    angle_offset : float, optional
        Offset angle in deg. Input 'auto' for automatic calculation. The default is 0.0.
    calibration : str, optional
        Either 'arduino_gonio_1' or 'raspberry_gonio2'. The default is 'raspberry_gonio2'.
    correct_time_drift : bool, optional
        Whether to correct or not the time drift. The default is True.

    Raises
    ------
    Exception
        If the fileiv is not found.

    Returns
    -------
    None.

    """
    folder = Path(folder)
    create_structure(folder)
    
    files = list((folder / raw_folder).glob('*_time-series.dat'))
    filesL0 = list((folder / raw_folder).glob('*_time-series*.dat'))
    
    if len(files) == 0:
        files = folder.glob('*_time-series*.dat')
        for f in files:
            shutil.move(f, folder / raw_folder / f.name)
        files = list((folder / raw_folder).glob('*_time-series.dat'))
        filesL0 = list((folder / raw_folder).glob('*_time-series*.dat'))
        
    # Assuming there is only one txt file and that it is the log
    iv_files = list(folder.glob('*.txt'))
    if len(iv_files) == 0:
        raise Exception('No iv file found.')
    fileiv = iv_files[0] 
    
    # output files
    file_log = folder  / 'performance_indicators.log'
    file_spectra = folder  / 'spectral_evolution.evolution'
    file_luminance = folder / 'time-luminance.dat'
    # Sample and measurements parameters

     # in deg, put 'auto' if you think it was misaligned
    current = current /1000 if current != None else None  # mA
    
    path_IRF = calibrations[calibration]['path_IRF']
    abs_calfactor = calibrations[calibration]['abs_calfactor']

    # Initial timestamp at the starting of the measurement
    t0 = np.loadtxt(fileiv, max_rows = 1, dtype = np.datetime64)
            
    # print(f'{len(files):d} gonio files')
    print(f'Processing {len(filesL0):d} files...')

    
    vtimes, vluminances = process_L0(filesL0, t0, folder = folder, path_IRF =  path_IRF, abs_calfactor=abs_calfactor)
    header = 'Luminance obtained from the forward spectrum\n'
    header += f'IRF_file: {path_IRF.name:s}\n'
    header += f'correct_time_drift: {correct_time_drift}\n'
    header += f'abs_calfactor: {abs_calfactor:4.2e}\n'
    header += f'angle_offset = {angle_offset:.2f} deg\n'
    header +=  f't0 = {t0}\n'
    header += 'EllapsedTime(s)\tLuminance(cd/m2)'
    np.savetxt(folder / 'time-luminance.dat', np.vstack([vtimes[1:], vluminances[1:]]).T, fmt = '%10.4f', header = header)

# ----------------- Goniometer data processing --------------------------
    for i, file in enumerate(files):
        plot = (i <= 3) or (i % 10 ==0)
      
        _ = process_goniodata(file, current = current,\
                                    angle_offset= angle_offset,\
                                    plot = plot,\
                                    correct_time_drift = correct_time_drift,\
                                    verbose =  False,\
                                    folder = str(folder / processed_folder),
                                    path_IRF = path_IRF,
                                    abs_calfactor  = abs_calfactor)
        plt.close()

    print('Done!')
 

    # Find the minimum voltage and max luminance ellapsed times and absolute time
    timeV, voltage = np.loadtxt(fileiv, usecols = (1,3), unpack = True, skiprows = 3)
    
    timeL, luminanceL = np.loadtxt(file_luminance, usecols = (0,1), unpack  =  True, skiprows = 4)
    
    iVmin = voltage.argmin()
    tVmin = timeV[iVmin]
    Vmin = voltage[iVmin]
    LtVmin = luminanceL[np.abs(timeL-tVmin).argmin()]
    
    iLmax = luminanceL.argmax()
    tLmax = timeL[iLmax]
    Lmax = luminanceL[iLmax]
    VtLmax = voltage[np.abs(timeV - tLmax).argmin()]
    
    # Define a steady state by time
    tss = tVmin # define steady-state time either by tVmin or by a specific time
    Vss = voltage[np.abs(tss - timeV).argmin()]
    Lss = luminanceL[np.abs(tss - timeL).argmin()]
    
    
    # Initial timesatamp at the starting of the measurement
    t0 = np.loadtxt(fileiv, max_rows = 1, dtype = np.datetime64)
    # Timestamp of the gonio measurement
    t1 = np.array([], dtype = np.datetime64)
    
    # Get the time vector from the gonio files
    for f in files:
        t1 = np.concatenate([t1, [np.loadtxt(f, max_rows = 1, dtype = np.datetime64)]])
    
    # Get the absolute tVmin/tLmax and tss and match it with the closes gonio file
    abs_tVmin = t0 + int(tVmin * 1e6)
    file_Vmin = files[np.abs(t1 - abs_tVmin).argmin()].stem + '.sri'
    
    abs_tLmax = t0 + int(tLmax * 1e6)
    file_Lmax = files[np.abs(t1 - abs_tLmax).argmin()].stem + '.sri'
    
    abs_tss = t0 + int(tss * 1e6)
    file_ss = files[np.abs(t1 - abs_tss).argmin()].stem + '.sri'

    with open(file_log, 'w') as f:
        f.write('# First row, Vmin L@Vmin tVmin and the closest gonio file in V, cd/m2 and s\n')
        f.write('# Second row, V@Lmax Lmax tLmax and the closest gonio file\n')
        f.write('# Third row, V@ss Lss tss and the closest gonio file\n')
        f.write('{:.2f}\t{:.2f}\t{:.2f}\t{}\n'.format(Vmin, LtVmin, tVmin, file_Vmin))
        f.write('{:.2f}\t{:.2f}\t{:.2f}\t{}\n'.format(VtLmax, Lmax, tLmax, file_Lmax))
        f.write('{:.2f}\t{:.2f}\t{:.2f}\t{}\n'.format(Vss, Lss, tss, file_ss))

    # Plots
    _ = plot_spectral_evolution(file_spectra, times = [10,60,120,600,3600, 3600*5, 99999], path = folder / plots_folder / 'spectral-evolution.png')
    
    # ------------------ Voltage luminance plot --------------------------
    tunit = 60
    fig, ax = plt.subplots()
    
    if tunit == 1:
        unit = 's'
    elif tunit == 60:
        unit = 'min'
    elif tunit == 3600:
        unit = 'h'
    
    
    timeV, voltageV = np.loadtxt(fileiv, usecols = (1,3), unpack = True, skiprows = 3)
    lineV, = ax.plot(timeV/tunit, voltageV, 'o-C0', label = 'voltage', markevery = 0.05)
    
    vtimes, vluminances = np.loadtxt(file_luminance, usecols = (0,1), unpack = True)
    
    ax2 = ax.twinx()
    # ax2.set_ylim(0,1750)
    lineL, = ax2.plot(vtimes/tunit, vluminances, 's-C1', label = 'luminance', mfc = 'white', markevery  = len(vluminances) // 25)
    
    ax.set_xlabel(f'Time ({unit})')
    ax.set_ylabel('Voltage (V)')
    ax2.set_ylabel('Luminance (cd m$^{-2}$)')
    file_id = fileiv.stem.split('_')[1]
    ax.set_title(f'Transient for {file_id}')
    
    Vmin, LtVmin, tVmin = np.loadtxt(file_log, skiprows = 3, max_rows=1, usecols = (0,1,2), unpack = True)
    Vss, Lss, tss = np.loadtxt(file_log, skiprows = 5, max_rows=1, usecols = (0,1,2), unpack = True)
    
    
    ax.annotate('V$_{min}$', xy=(tVmin/tunit, Vmin), xytext=(tVmin/tunit, Vmin + 2),
            arrowprops=dict(facecolor='black',arrowstyle  = '->'),
            )
    
    ax.axvline(tss/tunit, ls = '--', color = 'black', lw = 0.5)
    
    lines = [lineV, lineL]
    
    ax.legend(lines, [l.get_label() for l in lines])
    
    fig.savefig(folder / plots_folder / 'ivl_transient.png',  bbox_inches = 'tight')
    
    plt.close()

def create_times_vector(t0, files):
    """Creates a time vector based on the input files adn the t0 as the starting time"""
    vtimes = np.zeros((len(files), ))
    
    for i, file in enumerate(files):
        # Ellapsed time since start of the whole measurement
        t1 = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
        rTimes = np.loadtxt(file, skiprows = 4, usecols = 0)
        eTime = rTimes[rTimes.shape[0] // 2]
        
        # Add the time of the goniometer at the middle of the measurement, in us
        t2 = t1 + int(eTime * 1e6) # The time needs to be add un us

        # Calculate the relative time between the starting of the measurement (t0) and the current measruement
        vtimes[i] = np.float64(t2 - t0) / 1e6
    
    return vtimes

def fit_data_time_series(folder, thickness, error_landscape_file, folder_cal = calibration_dir, ethickness = None):
    
    folder = Path(folder)
    
    # Initial timestamp at the starting of the measurement
    fileiv =  list(folder.glob('*.txt'))[0]
    t0 = np.loadtxt(fileiv, max_rows = 1, dtype = np.datetime64)
    
    # Load the error landscape file
    simEL = load_simdata(folder_cal / error_landscape_file, wl_limit=(450,800))
    
    #Create and additional folder for the error maps
    error_landscape_dir = folder / fits_folder / 'error_landscapes'
    error_landscape_dir.mkdir() if not  error_landscape_dir.exists() else None
         
    # Files
    files = list((folder / processed_folder).glob('*_time-series.sri'))
    
    # Get the times relative to the time t0
    vtimes = create_times_vector(t0, list((folder / raw_folder).glob('*_time-series.dat')))
    
    output = np.zeros((len(files),7))
    
    # Thickness
    print(f'Fittings assuming a thickness of {thickness:.0f} nm')
    
    for i, file in enumerate(files):
        # Calculates the EZP based on the ethickness range
        if ethickness != None:
            error_low, pos_min_low, _, _ = error_landscape(file, thickness - max(ethickness,5), simEL, folder = error_landscape_dir)
            error_high, pos_min_high, _, _ = error_landscape(file, thickness + max(ethickness,5), simEL, folder = error_landscape_dir)
            error_min_high = error_high.min()
            error_min_low = error_low.min()
        else:
           pos_min_high, pos_min_low, error_min_high, error_min_low = [np.nan]*4
           
        error, pos_min, _, _ = error_landscape(file, thickness, simEL, folder = error_landscape_dir)
        # output[i,:] = [vtimes[i], pos_min, error.min(), 0, 0, 0, 0]
        output[i,:] = [vtimes[i], pos_min, error.min(), pos_min_low, error_min_low, pos_min_high, error_min_high]
    
    if ethickness == None:
        ethickness = np.nan
        
    # Saving the data    
    header = str(t0) +'\n'
    header += f'{thickness:.0f} +/- {ethickness:.0f}\n'
    header += ('{:^10s}\t'*3).format('Rel.Time(s)','EZP(bestFit)', 'RMSE')
    np.savetxt(folder / fits_folder / 'time-evolution_emission-zone.dat', output, header = header, fmt = '%10.6g')
    

def plot_fit_vs_experimental(folder, time_to_plot, error_landscape_file, name = 'comparison', ext = '.png', folder_cal = calibration_dir, alternative_EZP = [0.3, 0.7]):
    
    folder = Path(folder)
    
    file_EZP = folder / fits_folder / 'time-evolution_emission-zone.dat'
    with open(file_EZP) as f:
        f.readline()
        thickness = float(f.readline().split(' ')[1])
    
    vtimes, vpos = np.loadtxt(file_EZP, usecols = (0,1), unpack=True)
    
    # Files to plot and EZP for that file
    index2plot = (abs(time_to_plot - vtimes).argmin())
    files = list((folder / processed_folder).glob('*_time-series.sri'))
    pos = vpos[index2plot]
    file = files[index2plot]
    
    # Thickness
    print(f'Plotting data from the file "{file.stem}"')
    print(f'Plots assuming a thickness of {thickness:.0f} nm')
   
    # Load the error landscape file
    simEL = load_simdata(folder_cal / error_landscape_file, wl_limit=(450,800))
    
    plt.ioff()
    
    
    fname = folder / plots_folder / (name +'_fitting_best' + ext)
    
    wl_exp, iNormExpSRI, ipos_sim, wl_sim, NormSimSRI = compare_data(file, thickness, simEL, [pos], fname = fname, ext = ext)
    
    # Plot
    fcomparison = folder /  plots_folder / (name +'_fitting_bad' + ext)
    
    compare_data(file, thickness, simEL, alternative_EZP, fname = fcomparison, ext = ext)

    return simEL['angles'], wl_exp, iNormExpSRI, ipos_sim[0], wl_sim, NormSimSRI[0]
    
    
def plot_EZP_vs_time(folder, ext = '.png'):
    
    folder = Path(folder)
    fileiv =  list(folder.glob('*.txt'))[0]
    file_luminance = folder / 'time-luminance.dat'
    file_EZP = folder / fits_folder / 'time-evolution_emission-zone.dat'
    
    with open(file_EZP) as f:
        f.readline()
        thickness = float(f.readline().split(' ')[1])

# #------------------ Voltage luminance position plot --------------------------
    gridspec_kw = dict(height_ratios = (1,5))
    fig, [sax, ax] = plt.subplots(figsize = (4,3+3/5), nrows = 2, sharex=True, gridspec_kw=gridspec_kw)

    ax2 = ax.twinx()
    ax3 = ax.twinx()
    ax3.spines["right"].set_position(("axes", 1.25))
    ax.set_xlim(-5,200)
    ax.set_ylim(2.5,12)
    # ax.set_xlim(0, 180)
    ax3.set_ylim(0.2,0.8)
    fig.suptitle(f'd = {thickness:.0f} nm')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Voltage (V)')
    ax2.set_ylabel('Luminance (cd m$^{-2}$)')
    ax3.set_ylabel('Rel. emitter position')
    
    tunit = 60

    timeV, voltageV = np.loadtxt(fileiv, usecols = (1,3), unpack = True, skiprows = 3)
    lineV, = ax.plot(timeV/tunit, voltageV, 'o-C0', label = 'voltage', markevery = 0.05)

    timeL, luminanceL = np.loadtxt(file_luminance,usecols = (0,1), unpack = True)

    lineL, = ax2.plot(timeL/tunit, luminanceL, 's-C1', label = 'luminance', mfc = 'white', markevery  = len(luminanceL) // 25)

    timeP, positionP, errorP = np.loadtxt(file_EZP, usecols = (0,1,2), unpack = True)
    lineP, = ax3.plot(timeP/tunit, positionP, 'x-C2', label = 'position')

    sax.plot(timeP/tunit, errorP, 'x-C2')
    sax.set_ylabel('RMSE')
    
    lines = [lineV, lineL, lineP]

    
    ax.legend(lines, [l.get_label() for l in lines], bbox_to_anchor=(1.05, 1.4), loc='upper left')

    fig.savefig( folder /  plots_folder / ('transient+position' + ext), bbox_inches = 'tight')

    
    