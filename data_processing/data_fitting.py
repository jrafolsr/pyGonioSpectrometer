# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:31:51 2020

@author: JOANRR
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.io  import loadmat
from scipy.interpolate import interp1d


def load_simdata(file, wl_limit = None, angle_max = None):
    """
    Reads the *.mat file output structure from Mattias' Error Landscape generator and returns a dictionary structure with the fields wl, angle, ipos, dAL and data.
    
    Parameters:
    ----------
    file: str with the path to the *.mat file containing the error landscape
    wl_limit (opt): tuple with the lower and ipper values of the wavelengths to take
    angle_max (opt): angle maximum to consider
    
    Returns:
    output:  a dictionary structure with the fields wl, angle, ipos, dAL and data
    
    """
    tdata = loadmat(file)
    
    wl  = tdata['wl_q'][0]
    angles = tdata['angle_q'][0]
    dAL = tdata['dAL_sim'][0]
    ipos = tdata['ipos_sim'][0]
    data =  tdata['qnormSpecRadInt_sim_2D']
    
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

def load_simdata_new(file, wl_limit = None, angle_max = None):
    """Reads the *.mat file output structure from Mattias' Error Landscape generator and returns a dictionary structure with the fields wl, angle, ipos, dAL and data.
    
    Parameters:
    ----------
    file: str with the path to the *.mat file containing the error landscape
    wl_limit (opt): tuple with the lower and ipper values of the wavelengths to take
    angle_max (opt): angle maximum to consider
    
    Returns:
    output:  a dictionary structure with the fields wl, angle, ipos, dAL and data
    
    """
    tdata = loadmat(file)
    
    wl  = tdata['wl'][0]
    angles = tdata['angles'][0]
    dAL = tdata['dAl'][0]
    ipos = tdata['ipos'][0]
    data =  tdata['data']
    
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


def interpolate_expdata(sri, wl_i, angles_i, wl_o, angles_o):
    """Interpolates the input (experimental sri with wl_i and angles_i to the output wl_o and angles_o and normalizes to max wl of the zero angle
        Parameters
        ----------
        sri: 2D np.array
            The array with the spectral radiant intensity (wavelength as rows) for each of the collected angles (as columns). It assumes that there are 3 zero angles, as it is the output of process_data function.
        wl_i: 1D np.array
            Array with the input wavelengths, should match the length of the first dimension of the sri.
        angles_i: 1D np.array
            Array with the input angles, should match the length of the second dimension of the sri.    
        wl_o: 1D np.array
            Array with the output wavelengths to which the data will be interpolated.
        angles_i: 1D np.array
             Array with the output angles to which the data will be interpolated.
             
        Returns:
        --------
        NormSRI : interpolated and normalized (to zero angles) spectral radiant intensity. It has len(wl_o) rows and len(angles_o) columns. If points outside the range are requested it will extrapole while giving a warning, it is assumed to be fine as long as it is not a large extrapolation.
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
    """Get the closest index of the vector corresponding to the closest value"""
    # Handle exceptions in case the input parameters are out of the range
    if value > vector.max() or value < vector.min():
        raise Exception(f'The input value "{value:.0f}" is out of the range {vector.min()}-{vector.max()}')
    else:
        idx = (np.abs(vector - value)).argmin()
       
    return idx


def load_sri_file(file):
    temp = np.loadtxt(file)
    wavelengths = temp[2:, 0]
    angles = temp[1,1:]
    sri = temp[2:, 1:]
    return wavelengths, angles, sri

def error_landscape(file, thickness, simEL, weights = None, plot = False):
    """ Calculates the error landscape for a given thickness with respect the emitter positon within the device.
        
        Parameters
        ----------
        file: file path with the experiments data
        thickness: thickness of the experimental data
        simEL: dict structure with all the simulation data to compare the exp data
        weigths: None. You can input a vector to weight the angles in the error calculation, its length must correspond to the length of the simEL['angles']
        plot: if True, it plots a colormap of theangular spectral radiant intensity for the exp and sim data together with a colormap of the error at the best fit. It also plots the error vs the position of the emitter.
        
        Returns
        -------
        Error_Landscape: the error for each position of the emitter, corresponding to the vector providel in simEL['ipos']
        pos_min: the position with the minimum error iNormExpSRI, NormSimSRI[ipos_min]
        iNormExpSRI: the normalized and interpolated exp angular SRI
        NormSimSRI[ipos_min]: the normalized sim SRI for the minimum error position
        
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
    
    iNormExpSRI = interpolate_expdata(sri, wavelengths, angles, wl_sim, angles_sim)
    
    # Find the simulation data according to the input thickness   
    ithickness = gci(thickness, dAL_sim)
    
#     print(ithickness, dAL_SIM[ithickness])
    # The matrix EL_SIM['qnormSpecRadInt_sim_2D'] contains the simulated data
    # for thickness ranging 50 to 600 every 5 nm (111 columns)
    # and the rel. emitter position from x0 to x1 every dx, (n rows)
    # So, we choose the column according to the experimental thickness:
    NormSimSRI = simData[:,ithickness] # Columns are the simulated thicknesses

    #     print(NormSimSRI[0].shape)
    N_pos = len(ipos_sim) 
    Error_Landscape = np.zeros(ipos_sim.shape)

    for i in range(N_pos):
        # Substract the exp. amd sim. data and do the mean of the abs error in both axis, wavelengths and angles
        Error_Landscape[i] = (np.abs(iNormExpSRI - NormSimSRI[i]) * weights).mean(axis = 0).mean()
    
    # Take the minimum error, which will give the best fit to the emitter position
    ipos_min = Error_Landscape.argmin() # The i-th element corresponding to the min
    pos_min = ipos_sim[ipos_min] # The actual position of the emitter with the min error
    
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
        
        fig.savefig(file[:-4] + '_aSRI_heatmap.png', bbox_acnhor = 'tight')
        
        fig, ax = plt.subplots()
        ax.plot(ipos_sim, Error_Landscape, 'o-', label = f'd = {thickness:.0f}\nmin_pos = {ipos_sim[ipos_min]:.2f}')
        ax.set_xlabel('Rel. position of the emitter')
        ax.set_ylabel('$\Delta_{sim}^{meas}$')
        ax.set_title("Error landscape for emitter position")
        ax.legend()
        
        fig.savefig(file[:-4] + '_ipos_errorlandscape.png', bbox_inches = 'tight')
        
        N = len(angles_sim)
        offset = 0.25 * (N-1)
        
        fig, ax = plt.subplots(figsize = (6,4))

        for i in range(0,N, 1):
            ax.plot(wl_sim, offset - i*0.25 + NormSimSRI[ipos_min][:,i], f'--C{i}')
            ax.plot(wl_sim, offset - i*0.25 + iNormExpSRI[:,i], f'-C{i}', label = f'{angles_sim[i]:}°')
            ax.text(450,  offset + 0.04 - 0.25*i, f'{angles_sim[i]:.0f}°')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('SRI (a.u.)')
        ax.set_title(f'Experimental and simulated emission spectra\n dAL = {thickness}, pos = {pos_min:.2f}')
        # ax.legend()
        
        fig.savefig(file[:-4] + '_fitting.png', bbox_inches = 'tight')
        
    return Error_Landscape, pos_min, iNormExpSRI, NormSimSRI[ipos_min]


def fit_thickness(file, simEL, plot = False, weights = None):
    """Finds the thickness corresponding to the minimum error using the error_landscape function.
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
        ax.plot(dAL_sim, error_pos,'o-', label =  f'Min. Err. Thickness = {fitted_thickness} nm')
        ax.set_xlabel("Thickness of the AL (nm)")
        ax.set_ylabel('$\Delta_{sim}^{meas}$')
        ax.set_title("Error landscape for AL thickness")
        ax.legend()
        
    return fitted_thickness

def min_error_profile(weights, simEL_positions, exp_data, fitting = True):
    """Calculates the error with respect the experimental data assuming a linear combination of emmitters at different positions and with different weights."""   
    # Simulated data as a linear combination of emitters in multiple positions
    lc_SimData = np.zeros(exp_data.shape) # linear combination SimData
    
    for i, w in enumerate(weights):
        lc_SimData += w * simEL_positions[i]
    
    # Abs error
    # error = ((np.abs(exp_data - lc_SimData)).mean(axis = -1)).mean(axis = -1)
    # Quadratic error
    error = ((np.sqrt((exp_data - lc_SimData)**2)).mean(axis = -1)).mean(axis = -1)
#     print (error)
    
    if fitting:
        return error
    else:
        return error, lc_SimData
    