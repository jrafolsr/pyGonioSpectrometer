# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:31:51 2020

@author: JOANRR
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.io  import loadmat
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import seaborn as sns
from .data_processing import load_sri_file
from os.path import join as pjoin
import os
from pathlib import Path



def load_simdata(file, wl_limit = None, angle_max = None, pos_lim = None):
    """
    Reads the *.mat file output structure from Mattias' Error Landscape generator and returns a dictionary structure with the fields wl, angle, ipos, dAL and data.
    
    Parameters
    ----------
    file : str or path
        Path to the *.mat file containing the error landscape
    wl_limit: tuple, optional
        2-element tuple with the lower and upper values of the wavelengths to load.
    angle_max : int or float, optional
        Maximum angle to load.
    pos_lim :  tuple, optional
        2-element tuple with the lower and upper values of the simulated emission zone position to load.
    
    Returns
    -------
    output:  dict
        A dictionary with the fields wl, angle, ipos, dAL and data
    
    """
    tdata = loadmat(file)
    
    wl  = tdata['wl_q'][0].astype(np.float)
    angles = tdata['angle_q'][0].astype(np.float)
    dAL = tdata['dAL_sim'][0].astype(np.float)
    ipos = tdata['ipos_sim'][0].astype(np.float)
    # The first of data is the ipos, the second is the thickness
    data =  tdata['qnormSpecRadInt_sim_2D']
    
    if pos_lim is not None:
        imin = gci(pos_lim[0], ipos)
        imax = gci(pos_lim[1], ipos)
        ipos = ipos[imin:imax]
        data = data[imin:imax, :]
    
    if wl_limit is not None:
        wl_filter = (wl >= wl_limit[0]) & (wl <=  wl_limit[1])
    else:
        wl_filter = [True] * len(wl)
    
    if angle_max is not None:
        angle_filter = angles <= angle_max
    else:
        angle_filter = [True] * len(angles)
        
    wl = wl[wl_filter]
    angles = angles[angle_filter]
    
    for i in range(len(ipos)):
        for j in range(len(dAL)):
            data[i,j] = data[i,j][wl_filter, :]
            data[i,j] = data[i,j][:, angle_filter]

    
    output = dict(wl  = wl,\
                  angles = angles,\
                  dAL = dAL,\
                  ipos = ipos,\
                  data =  data)
    

    return output



def load_simdata_python(file, wl_limit = None, angle_max = None, pos_lim = None):
    """
    Reads the *.mat file output structure from SETFOS-Python Error Landscape generator and returns a dictionary structure with the fields wl, angle, ipos, dAL and data.
    
    Parameters
    ----------
    file : str or path
        Path to the *.mat file containing the error landscape
    wl_limit: tuple, optional
        2-element tuple with the lower and upper values of the wavelengths to load.
    angle_max : int or float, optional
        Maximum angle to load.
    pos_lim :  tuple, optional
        2-element tuple with the lower and upper values of the simulated emission zone position to load.
    
    Returns
    -------
    output:  dict
        A dictionary with the fields wl, angle, ipos, dAL and data
    
    """
    tdata = loadmat(file)
    
    wl  = tdata['wl_q'][0].astype(np.float)
    angles = tdata['angle_q'][0].astype(np.float)
    dAL = tdata['dAL_sim'][0].astype(np.float)
    ipos = tdata['ipos_sim'][0].astype(np.float)
    # The first of data is the ipos, the second is the thickness
    data =  tdata['qnormSpecRadInt_sim_2D']
    
    if pos_lim is not None:
        imin = gci(pos_lim[0], ipos)
        imax = gci(pos_lim[1], ipos)
        ipos = ipos[imin:imax]
        data = data[imin:imax, :]
    
    if wl_limit is not None:
        print("OBS! The wavelength filtering is not done properly, it works soso, fix it soon.")
        wl_slice = slice(abs(wl_limit[0] - wl).argmin(), abs(wl_limit[1] - wl).argmin())
    else:
        wl_slice = slice(None)
    
    if angle_max is not None:
        print("OBS! The angle filtering is not done properly, it works soso, fix it soon.")
        angle_slice = slice(abs(angle_max - angles).argmin())
    else:
        angle_slice = slice(None)
        
    angles = angles[angle_slice]
    wl = wl[wl_slice]
    new_data = data[:, :, wl_slice, angle_slice]

    output = dict(wl  = wl,\
                  angles = angles,\
                  dAL = dAL,\
                  ipos = ipos,\
                  data =  new_data)
    

    return output

def interpolate_expdata(sri, wl_i, angles_i, wl_o, angles_o):
    """
    Interpolates the input (experimental sri with wl_i and angles_i to the output wl_o and angles_o and normalizes to max wl of the zero angle
                            
    Parameters
    ----------
    sri : 2D numpy.array
        The array with the spectral radiant intensity (wavelength as rows) for each of the collected angles (as columns). It assumes that there are 3 zero angles, as it is the output of process_data function.
    wl_i : 1D numpy.array
        Array with the input wavelengths, should match the length of the first dimension of the sri.
    angles_i : 1D numpy.array
        Array with the input angles, should match the length of the second dimension of the sri.    
    wl_o : 1D numpy.array
        Array with the output wavelengths to which the data will be interpolated.
    angles_o : 1D np.array
         Array with the output angles to which the data will be interpolated.
         
    Returns
    -------
    NormSRI : numpy.array
        Interpolated and normalized (to zero angles) spectral radiant intensity. It has len(wl_o) rows and len(angles_o) columns. If points outside the range are requested it will extrapole while giving a warning, it is expected to be fine as long as it is not a large extrapolation.
    """
    
    # Average the zero angles, which will be the three central values
    Nzeros = int(len(angles_i - 1)/2)
    zero_slice = slice(Nzeros-1, Nzeros+2)
    zero_angles = angles_i[zero_slice].mean()
    zero_sri = sri[:, zero_slice].mean(axis = -1, keepdims =  True)
    
    neg_slice = slice(0,Nzeros-1) # The negative angles slice object
    pos_slice = slice(Nzeros+2,None) # The negative angles slice object
    
    # Creating the new angles_i vector and sri matrice, once the 3-zeros averaged
    angles_i = (np.hstack([angles_i[neg_slice],zero_angles, angles_i[pos_slice]]))    
    sri = np.hstack([sri[:, neg_slice],zero_sri, sri[:, pos_slice]])
    
    # Warnings if the data has to be extrapolated
    if wl_o[0] < wl_i[0] or wl_o[-1] > wl_i[-1]:
        print('WARNING: The output wavelength range ({0:.0f} - {1:.0f} nm) will be extrapolated from the input range  ({2:.0f} - {3:.0f} nm)'.format(wl_o[0],wl_o[-1],wl_i[0],wl_i[-1]))
    if angles_o[0] < angles_i[0] or angles_o[-1] > angles_i[-1]:
        print('WARNING: The output angle range ({0:.0f} - {1:.0f} °) will be extrapolated from the input range  ({2:.0f} to {3:.0f} °)'.format(angles_o[0],angles_o[-1],angles_i[0],angles_i[-1]))

    # Create the interpolation function for the angles
    f = interp1d(angles_i, sri, kind = 'quadratic', axis = -1, fill_value='extrapolate')
    # Take both hemispheres and do the average
    a_sri = (f(angles_o) + f(-1.*angles_o)) / 2
                          
    # Create the interpolation function for the wavelengths
    g = interp1d(wl_i, a_sri, kind = 'quadratic', axis = 0, fill_value= 0, bounds_error=False) 
    i_sri = g(wl_o)
    
    # Normalize by the angle zero, assumin the matrix is sorted to the absolute angle!!!
    NormSRI = i_sri / i_sri[:,0].max()
    
    return NormSRI

def gci(value, vector):
    """
    Get the closest index of the vector corresponding to the closest value
    
    Parameters
    ----------
    value: float or int,
        The value to search
    vector: numpy.array,
        Vector into which find the value
    
    Returns
    -------
    idx: int
        The index corresponding to the closest queried value.
    
    """
    # Handle exceptions in case the input parameters are out of the range
    if value > vector.max() or value < vector.min():
        raise Exception(f'The input value "{value:.0f}" is out of the range {vector.min()}-{vector.max()}')
    else:
        idx = (np.abs(vector - value)).argmin()
       
    return idx


def error_landscape(file, thickness, simEL, weights = None, plot = False, folder = ''):
    """ 
    Calculates the error landscape for a given thickness with respect the emitter positon within the device. It outputs a *.el file containing two columns: the EZ position and the error at each position.
    
    Parameters
    ----------
    file : stri or path
        Path with the experiments data.
    thickness : floator int
        Thickness of the experimental data.
    simEL : dict
        Dictionary with all the simulation data to compare the experimental data.
    weigths : numpy.array, optional
        You can input a vector to weight the angles in the error calculation, its length must correspond to the length of the simEL['angles']. The default is None.
    plot : bool, optional
        If True, it plots a colormap of the angular spectral radiant intensity for the exp and sim data together with a colormap of the error at the best fit. It also plots the error vs the position of the emitter. The default is False.
    folder : str, optional
        String containing the name of the folder into which save the outputs. The default is ''.
    
    Returns
    -------
    Error_Landscape : 1D numpy.array
        The error for each position of the emitter, corresponding to the vector providel in simEL['ipos']
    pos_min : float
        The position with the minimum error.
    iNormExpSRI : numpy.array
        The normalized and interpolated exp angular SRI.
    NormSimSRI[ipos_min] : numpy.array
        The normalized sim SRI for the minimum error position.
        
    """
    
    wavelengths, angles, sri = load_sri_file(file)
    
    if weights is None:
        weights  = 1.0
    else:
        print('INFO: The error will be weighted.')
    
    wl_sim = simEL['wl']
    angles_sim = simEL['angles']
    dAL_sim = simEL['dAL']
    ipos_sim = simEL['ipos']
    simData = simEL['data']

    # Output names
    foutput = os.path.basename(file)[:-4]
    if folder == '':
        folder = os.path.dirname(file)
    
    iNormExpSRI = interpolate_expdata(sri, wavelengths, angles, wl_sim, angles_sim)
    
    # Find the simulation data according to the input thickness   
    ithickness = gci(thickness, dAL_sim)
    
#     print(ithickness, dAL_SIM[ithickness])
    # The matrix EL_SIM['qnormSpecRadInt_sim_2D'] contains the simulated data
    # for thickness ranging 50 to 600 every 5 nm (111 columns)
    # and the rel. emitter position from x0 to x1 every dx, (n rows)
    # So, we choose the column according to the experimental thickness:
    NormSimSRI = simData[:,ithickness] # Columns are the simulated thicknesses
    thickness_sim = dAL_sim[ithickness] # simulated thickness
    #     print(NormSimSRI[0].shape)
    N_pos = len(ipos_sim) 
    Error_Landscape = np.zeros(ipos_sim.shape)

    for i in range(N_pos):
        # Substract the exp. amd sim. data and do the mean of the abs error in both axis, wavelengths and angles
        Error_Landscape[i] = (np.sqrt((iNormExpSRI - NormSimSRI[i]) ** 2).mean(axis = 0)* weights).mean()
    
    # Take the minimum error, which will give the best fit to the emitter position
    min_error = Error_Landscape.min()
    ipos_min = Error_Landscape.argmin() # The i-th element corresponding to the min
    pos_min = ipos_sim[ipos_min] # The actual position of the emitter with the min error
    
    np.savetxt(pjoin(folder, foutput + '.el'),\
               np.vstack([ipos_sim, Error_Landscape]).T, header = 'Rel.ipos\t Error',\
               fmt = '%4.2f %10.6f')
    
    if plot:
        fig, [ax,ax1,ax2] = plt.subplots(ncols=3, figsize = (12,4), sharex=True, sharey=True)

        ax.pcolormesh(angles_sim,wl_sim,iNormExpSRI, cmap = 'inferno', vmin=0, vmax = 1)
        ax1.pcolormesh(angles_sim,wl_sim,NormSimSRI[ipos_min], cmap = 'inferno', vmin=0, vmax = 1)
        # Error plot
        errplot = ax2.pcolormesh(angles_sim,wl_sim,np.abs(iNormExpSRI - NormSimSRI[ipos_min]), cmap = 'inferno')
        cbar = plt.colorbar(errplot)
        cbar.set_label('Error')
        
        ax.set_ylabel('Wavelength (nm)')
        fig.suptitle('Spectral Radiant Intensity ($W\cdot sr^{-1}\cdot nm^{-1}$)')
        ax.set_title('Experimental')
        ax1.set_title(f'Simulated pos = {pos_min:.2f}')
        ax.set_xlabel('Angle (°)')
        ax1.set_xlabel('Angle (°)')
        ax2.set_xlabel('Angle (°)')
        
        fig.savefig(pjoin(folder, foutput + '_aSRI_heatmap.png'), bbox_inches = 'tight')
        
        fig, ax = plt.subplots()
        ax.plot(ipos_sim, Error_Landscape, 'o-', label = f'd = {thickness_sim:.0f}\nmin_pos = {ipos_sim[ipos_min]:.2f}\nmin_error = {min_error:.2g}')
        ax.set_xlabel('Rel. position of the emitter')
        ax.set_ylabel('$\Delta_{sim}^{meas}$')
        ax.set_title("Error landscape for emitter position")
        ax.legend()
        
        fig.savefig(pjoin(folder, foutput + '_ipos_errorlandscape.png'), bbox_inches = 'tight')
        
        N = len(angles_sim)
        offset = 0.25 * (N-1)
        
        c =  sns.cubehelix_palette(N, start=.5, rot=-.75, reverse = True)
        fig, ax = plt.subplots(figsize = (4,4))

        for i in range(0,N, 1):
            ax.plot(wl_sim, offset - i*0.25 + NormSimSRI[ipos_min][:,i], '--', color = c[i], lw  =1)
            
            # label = f'{angles_sim[i]:}°'

            ax.plot(wl_sim, offset - i*0.25 + iNormExpSRI[:,i], '-', color = c[i], lw = 1)
            
            if i % 2 == 0:
                ax.text(wl_sim.min(),  offset + 0.04 - 0.25*i, f'{angles_sim[i]:.0f}°', fontsize = 'x-small')
                
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('SRI (a.u.)')
        ax.set_title(f'Experimental and simulated emission spectra\n dAL = {thickness}, pos = {pos_min:.2f}')
        # ax.legend()
        
        fig.savefig(pjoin(folder, foutput + '_fitting.png'), bbox_inches = 'tight')
        
    return Error_Landscape, pos_min, iNormExpSRI, NormSimSRI[ipos_min]


def fit_forward_position(file, thickness, simEL, folder = ''):
    """
    Calculates the error landscape for a given thickness with respect the emitter positon within the device assuming there is only the forward emission available.
        
    Parameters
    ----------
    file: str path,
        File path with a *.evolution data file.
    thickness: float
        Thickness of the experimental data.
    simEL: dict,
        Dictionary with all the simulation data to compare the exp data.
    folder : str, optional
        String containing the name of the folder into which save the outputs. The default is ''.
    
    Returns
    -------
    times
    position_w_time
    min_error_w_time
        
    """
    
    data = np.loadtxt(file)
    wl = data[0,2:]
    times = data[2:,0]
    tsri = data[2:,2:] 
    
    wl_sim = simEL['wl']
    dAL_sim = simEL['dAL']
    ipos_sim = simEL['ipos']
    simData = simEL['data']
    
    # Output names
    foutput = os.path.basename(file)[:-10]
    if folder == '':
        folder = os.path.dirname(file)
    
    # Normalize by the maximum (across columns, as it is the wl)
    tsri = tsri / tsri.max(axis = 1, keepdims = True)
    # Create an interpolation function
    f = interp1d(wl, tsri, axis = 1, kind='linear')
    # Interpolate to the simulated wavelengths
    iNormExpSRI= f(wl_sim)

    # Find the simulation data according to the input thickness   
    ithickness = gci(thickness, dAL_sim)
    
# #     print(ithickness, dAL_SIM[ithickness])
    # The matrix EL_SIM['qnormSpecRadInt_sim_2D'] contains the simulated data
    # for thickness ranging 50 to 600 every 5 nm (111 columns)
    # and the rel. emitter position from x0 to x1 every dx, (n rows)
    # So, we choose the column according to the experimental thickness:
    NormSimSRI = simData[:,ithickness] # Columns are the simulated thicknesses
    thickness_sim = dAL_sim[ithickness] # simulated thickness
    #     print(NormSimSRI[0].shape)
    N_times = times.shape[0]
    N_pos = ipos_sim.shape [0]
    
    position_w_time = np.zeros((N_times, ))
    min_error_w_time = np.zeros((N_times, ))
    for i in range(N_times):
        error  = np.zeros((N_pos, ))
        for j in range(N_pos):
        # Substract the exp. amd sim. data and do the mean of the abs error and take the minimum value,column zero from NormSimSRI corresponds to the angle zero
            error[j] = np.sqrt((iNormExpSRI[i,:] - NormSimSRI[j][:,0]) ** 2).mean()
        
        k = error.argmin()
        position_w_time[i] = ipos_sim[k]
        min_error_w_time[i] = error[k]
        # if i % 20 == 0:
        #     fig, ax = plt.subplots()
        #     text = f'Time = {times[i]/60:.2f} min\nd = {thickness:.0f}({thickness_sim:.0f})\n'+ '$\delta_{pos}$ = ' + f'{position_w_time[i]:.2f}'
        #     ax.text(0.95,0.95, text, va = 'top', ha = 'right', transform=ax.transAxes)
        #     ax.plot(wl_sim, NormSimSRI[error.argmin()][:,0],f'--C{i%10}')
        #     ax.plot(wl_sim, iNormExpSRI[i,:], f'-C{i%10}')
        
    header =f'Simulated(experimental) thickness: {thickness:.0f}({thickness_sim:.0f}\n' + 'Time(s)\tRel.ipos\t Error'
    np.savetxt(pjoin(folder, foutput + '+EZ.evolution'),\
                np.vstack([times, position_w_time, min_error_w_time]).T,\
                    header = header,\
                    fmt = '%4.2f\t%10.6f\t%10.6f')
            
    return times, position_w_time, min_error_w_time



def fit_thickness(file, simEL, plot = False, weights = None):
    """
    Finds the thickness corresponding to the minimum error using the error_landscape function.
    
    Parameters
    ----------
    file: str or path
        *.sri file path with the experiments data
    simEL: dict
        Dictionary with all the simulation data to compare the exp data.
    plot: bool,
        If True, it plots the thickness error landscape. The default is False.
    weigths : numpy.array, optional
        You can input a vector to weight the angles in the error calculation, its length must correspond to the length of the simEL['angles']. The default is None.
    
    Returns
    -------
    fitted_thickness : float
        The thicknes giving the minimum error.
    error_pos: 1D numpy.array
        The error vector, corresponding the the minimim error obtained at each thickness.   
    
    """
    dAL_sim = simEL['dAL']

    error_pos = np.zeros(dAL_sim.shape)
    
    for i,d in enumerate(dAL_sim):
        Error_Landscape,_,_,_= error_landscape(file, d, simEL, plot = False, weights = weights)
        error_pos[i] = Error_Landscape.min()
    
    # The index with the minimum error from all the thicknesses tried
    fitted_thickness = dAL_sim[error_pos.argmin()]
    
    if plot:
        fig, ax = plt.subplots()
        ax.plot(dAL_sim, error_pos,'o-', label =  f'Min. Err. Thickness = {fitted_thickness:.0f} nm')
        ax.set_xlabel("Thickness of the AL (nm)")
        ax.set_ylabel('$\Delta_{sim}^{meas}$')
        ax.set_title("Error landscape for AL thickness")
        ax.legend()
        
    return fitted_thickness, error_pos

def min_error_profile(weights, simEL_positions, exp_data, fitting = True):
    """
    Calculates the error with respect the experimental data assuming a linear combination of emmitters at different positions and with different weights.

    Parameters
    ----------
    weights : TYPE
        DESCRIPTION.
    simEL_positions : TYPE
        DESCRIPTION.
    exp_data : TYPE
        DESCRIPTION.
    fitting : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """   
    # Simulated data as a linear combination of emitters in multiple positions
    lc_SimData = np.zeros(exp_data.shape) # linear combination SimData
    
    for i, w in enumerate(weights):
        lc_SimData += w * simEL_positions[i]
    
    lc_SimData /= lc_SimData[:,0].max()
    # Abs error
    # error = ((np.abs(exp_data - lc_SimData)).mean(axis = -1)).mean(axis = -1)
    # Quadratic error
    error = ((np.sqrt((exp_data - lc_SimData) ** 2)).mean(axis = 0)).mean(axis = -1)
#     print (error)
    
    if fitting:
        return error
    else:
        return error, lc_SimData

def compare_data(file, thickness, simEL, positions, fname = None, ext = '.png', legend_text = None):
    """
    Generates as many plots as requested comparing the expermimental data with the simulated SRI for the given EZ positions in the parameters positions.
        
    Parameters
    ----------
    file : str or path
        Path with the experiments data.
    thickness : floator int
        Thickness of the experimental data.
    simEL : dict
        Dictionary with all the simulation data to compare the experimental data.
    positions : list
        List with the position do you want to plot.
    fname : str, optional
        Filename prefix (can be a path too). If None is passed, it takes the input file as prefix. The default is None.
    ext : str, optional
        Figure file format, it accepts '.png' or '.svg'. The default is '.png'.
    legend_text: str, optional
        Text to add in the legend. The default is None.
    Returns
    -------
    iNormExpSRI: numpy.array
        The normalized and interpolated exp angular SRI
    ipos_sim: numpy.array
        Vector with the the emission positions return (the closest the the queried values)
    NormSimSRI: numpy.array
        The normalized sim SRI for the queried positions
        
    """
    file = Path(file)
    
    wavelengths, angles, sri = load_sri_file(file)
    
    wl_sim = simEL['wl']
    angles_sim = simEL['angles']
    dAL_sim = simEL['dAL']
    ipos_sim = simEL['ipos']
    simData = simEL['data']
    
    iNormExpSRI = interpolate_expdata(sri, wavelengths, angles, wavelengths, angles_sim)
    
    # Find the simulation data according to the input thickness   
    ithickness = gci(thickness, dAL_sim)
    
#     print(ithickness, dAL_SIM[ithickness])

    # So, we choose the column according to the experimental thickness:
    NormSimSRI = simData[:,ithickness] # Columns are the simulated thicknesses
    thickness_sim = dAL_sim[ithickness] # simulated thickness
    
    # Convert position sto a list in case only one is passed
    if (type(positions) == int)  or type(positions) == float:
        positions = [positions]
        
    # Find the indices of positions in the ipos_sim vector
    ipositions = [gci(p, ipos_sim) for p in positions]
    
    print('INFO: Plotting data for', ipos_sim[ipositions])

    
    N = len(angles_sim)
    offset = 0.25 * (N-1)
    
    c =  sns.cubehelix_palette(N, start=.5, rot=-.75, reverse = True)
    
    for k in ipositions:
        fig, ax = plt.subplots(figsize = (3,3))
    
        for i in range(0,N, 1):
            kwargs = dict(color = c[i],lw  = 1)
            
            pos = ipos_sim[k]
            ax.plot(wl_sim, offset - i*0.25 + NormSimSRI[k][:,i],'--',  **kwargs)
            
            # label = f'{angles_sim[i]:}°'
    
            ax.plot(wavelengths, offset - i*0.25 + iNormExpSRI[:,i], ls = '-', **kwargs)
            
            if i % 2 == 0:
                ax.text(wl_sim.min(),  offset + 0.04 - 0.25*i, f'{angles_sim[i]:.0f}°', fontsize = 'x-small')
                
        ax.set_xlabel('Wavelength (nm)')
        ax.yaxis.set_ticklabels([])
        ax.set_ylabel('SRI (a.u.)')
        
        if legend_text == None:
            text = 'd$_{AL} = $' + f'{thickness:.0f} ({thickness_sim:.0f}) nm\n'+'EZP = ' +f'{pos:.2f}'
        else:
            text = legend_text
            
        ax.text(0.98,0.98, text, va = 'top', ha = 'right', transform=ax.transAxes, fontsize = 'small')
        
        if fname is None:
            fname = file.stem
            
        fig.savefig(fname.parent / (fname.stem + f'_comparison_delta={pos:.02f}' + ext), bbox_inches = 'tight')
        
    return iNormExpSRI, ipos_sim[ipositions], NormSimSRI[ipositions]

# def arbitrary_profile_fitting(input_file, model_data, thickness, step = 1, w0 = None, method = 'trf'):

#     simEL = model_data

#     wl_sim = simEL['wl']
#     angle_sim = simEL['angles']
    
#     # Preparing the model data, only taking every n-step
#     step = 1
#     initial = step - 1
#     final = -initial if initial != 0 else None
#     slicing = slice(initial, final, step)
#     positions_sim = simEL['ipos'][slicing]

#     # loading, interpolated and normalizing hte experimental data to the simulated ql and angle vectors
#     wl, angle, sri = load_sri_file(input_file)
#     exp_data = interpolate_expdata(sri, wl, angle, wl_sim, angle_sim)
    
#     # Choose the index corresponding to the input thickness
#     ithickness = gci(thickness, simEL['dAL'])
#     thickness_sim = simEL['dAL'][ithickness]
    
#     print(f'The input thickness is {thickness:.2f}. Tghe closest modelled data corresponds to a thickness of {thickness_sim:.0f} nm\n')
#     # Taking the index corresponding to the input thickness and the specified slice for the positions
#     sim_data = simEL['data'][slicing, ithickness]
    
#     if w0 is None:
#         w0 = np.zeros(positions_sim.shape)
#         w0 = w0 + 1/len(w0)
        

#     res_lsq  = least_squares(min_error_profile, w0, bounds=(0, 1),\
#                      args = (sim_data, exp_data),\
#                          method = method ,max_nfev  = 9000,ftol = 1e-10, verbose = 1)

#     w = res_lsq.x

#     # Initial error
#     error0, SimData0 = min_error_profile(w0, sim_data, exp_data, fitting = False)
    
#     # Profile
#     errorf, SimDataf = min_error_profile(w, sim_data, exp_data, fitting = False)

#     # Calculate the error with a single delta fit
#     _, pos_min, _, sim_data_best_delta  = error_landscape(input_file, thickness, simEL)

#     error_s, SimData_s = min_error_profile([1], [sim_data_best_delta], exp_data, fitting = False)

# text =  f'Thickness is {thickness:.0f}\n'
# text += f'Single delta position = {pos_min}\n'
# text += f'Error with a delta = {error_s}\n' 
# text += f'Error with profile = {errorf}\n'
# text += f'Error with homogenous profile {error0}\n'


# print(text)



# fig, [ax, ax1] = plt.subplots(ncols=2, figsize = (10,4), gridspec_kw={'width_ratios': [1.5, 1]})

# for i,a in enumerate(angles_sim):
#     offset = (len(angles_sim) - i - 1)*0.25
#     ax.plot(wl_sim, offset + SimData_s[:,i], ':C' + str(i%10), lw = 1)
#     ax.plot(wl_sim, offset + SimDataf[:,i], '--C' + str(i%10), lw = 1)
#     ax.plot(wl_sim, offset + exp_data[:,i], '-C' + str(i%10), label = f'{a:}°')
#     ax.set_xlabel('Wavelength (nm)')
#     ax.set_ylabel('Norm. SRI (a.u.)')

#     ax.set_title('Fittings')
# ax.legend()

# ax1.set_xlabel('Rel.position of the emitter')
# ax1.set_ylabel('Weight')
# ax1.set_title("Fitted profile")
# ax1.set_xlim(0,1)

# ax1.bar(positions, w, width = 0.02)
# ax1.bar(positions, w0, width = 0.005)
# ax1.axvline(pos_min, c= 'black', ls = '--')
# plt.subplots_adjust(wspace = 0.3)



# text = '$\Delta_{\delta}$ = ' + f'{error_delta:.4g}\n' 
# text += '$\Delta_{profile}$ = ' + f'{errorf:.4g}\n' 
# text += '$\Delta_{h}$ = ' + f'{error0:.4g}\n' 


# ax[2].text(0.95,0.95,text, ha = 'right', va = 'top', transform = ax[2].transAxes, fontsize = 'small')


# name_el_file = os.path.basename(file_sri)[0:-4] + '_' +  method

# fig.suptitle(f'Thickness = {thickness:.0f} nm (' + name_el_file + ')')
# fig.savefig(pjoin(path, name_el_file.png))