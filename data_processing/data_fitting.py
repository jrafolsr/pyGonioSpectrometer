# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:31:51 2020

@author: JOANRR
"""

from os.path import join as pjoin
import numpy as np
from matplotlib import pyplot as plt
from scipy.io  import loadmat
import os 


fcal = pjoin(os.path.dirname(__file__), 'calibration_files')


# This is manual, so far, needs to be a function
# EL_SIM2 = loadmat(pjoin(fcal, 'ErrorLandscape_2D-OLED_Hres_v1'))
# EL_SIM2['qnormSpecRadInt_sim_2D'] = EL_SIM2.pop('ans')
# EL_SIM2['ipos_sim'] = np.arange(0.01, 1.0, 0.01)
# EL_SIM2['dAL_sim'] = np.arange(50, 205, 5)
# EL_SIM2['wl_q'] = np.arange(380,781, 1)
# EL_SIM2['angle_q'] = np.arange(0,90, 10)

# WL_Q = EL_SIM2['wl_q']
# ANGLE_Q = EL_SIM2['angle_q']
# dAL_SIM = EL_SIM2['dAL_sim']
# i_POS = EL_SIM2['ipos_sim']
# SIM_DATA = EL_SIM2['qnormSpecRadInt_sim_2D']

EL_SIM = loadmat(pjoin(fcal, 'Errorlandscape_sim'))
WL_Q = EL_SIM['wl_q'][0]
ANGLE_Q = EL_SIM['angle_q'][0]
dAL_SIM = EL_SIM['dAL_sim'][0]
i_POS = EL_SIM['ipos_sim'][0]
SIM_DATA = EL_SIM['qnormSpecRadInt_sim_2D']



def interpolate_expdata(sri, wl_i, angles_i, wavelengths, angles):   
    # Manually add the 0 intensity at 90°
    angles = (np.hstack([-90, angles, 90]))
    N = len(wavelengths) # Original wavelengths vector length
    
    sri = np.hstack([np.zeros((N,1)), sri, np.zeros((N,1))])
    
    # I assume the data is already sorted
    a_sri = np.zeros((N, len(angles_i)))
    
    for i in range(N):
        # Take both hemispheres and to the average
        a_sri[i,:] = (np.interp(angles_i, angles, sri[i,:])+\
                            np.interp(-1.*angles_i, angles, sri[i,:])) / 2
                          
    # Interpolate to wavelengths from simulation wl_q
    i_sri = np.zeros((len(wl_i), len(angles_i)))  
    
    for i in range(len(angles_i)):
        i_sri[:,i] = np.interp(wl_i, wavelengths, a_sri[:,i])
    
    # Normalize by the angle zero, assumin the matrix is sorted to the absolute angle!!!
    NormSRI = i_sri / i_sri[:,0].max()
    
    return NormSRI

def get_thickness_index(thickness):
    """Gets index of the simulated d_al with the closest value to the input thickness """
    # Handle exceptions in case the input parameters are out of the range
    if thickness > dAL_SIM.max() or thickness < dAL_SIM.min():
        raise Exception(f'The input thickness "{thickness:.0f}" is out of the range {dAL_SIM.min()}-{dAL_SIM.max()}')
    else:
        ithickness = np.argmin(np.abs(dAL_SIM - thickness))
       
    return ithickness

def get_position_index(position):
    """Gets index of the simulated i_pos with the closest value to the input position """
    # Handle exceptions in case the input parameters are out of the range
    if position > i_POS.max() or position < i_POS.min():
        raise Exception(f'The input position "{position}" is out of the range {i_POS.min()}-{i_POS.max()}')
    else:
        ipos =  np.argmin(np.abs(i_POS - position))
    
    return ipos

def load_sri_file(file):
    temp = np.loadtxt(file)
    wavelengths = temp[2:, 0]
    angles = temp[1,1:]
    sri = temp[2:, 1:]
    return wavelengths, angles, sri

def error_landscape(file, thickness, plot = False):
    
    wavelengths, angles, sri = load_sri_file(file)
    iNormExpSRI = interpolate_expdata(sri, WL_Q, ANGLE_Q, wavelengths, angles)
    
    # Find the simulation data according to the input thickness   
    ithickness = get_thickness_index(thickness)
    
#     print(ithickness, dAL_SIM[ithickness])
    # The matrix EL_SIM['qnormSpecRadInt_sim_2D'] contains the simulated data
    # for thickness ranging 50 to 600 every 5 nm (111 columns)
    # and the rel. emitter position from x0 to x1 every dx, (n rows)
    # So, we choose the column according to the experimental thickness:
    NormSimSRI = SIM_DATA[:,ithickness] # Columns are the simulated thicknesses

    #     print(NormSimSRI[0].shape)
    N_pos = len(i_POS)
    
    Error_Landscape = np.zeros(i_POS.shape)

    for i in range(N_pos):
        # Substract the exp. amd sim. data and do the mean of the abs error in both axis, wavelengths and angles
        Error_Landscape[i] = (np.abs(iNormExpSRI - NormSimSRI[i])).mean().mean()
    
    # Take the minimum error, which will give the best fit to the emitter position
    ipos_min = Error_Landscape.argmin() # The i-th element corresponding to the min
    pos_min = i_POS[ipos_min] # The actual position of the emitter with the min error
    
    if plot:
        fig, [ax,ax1,ax2] = plt.subplots(ncols=3, figsize = (12,4), sharex=True, sharey=True)

        ax.pcolormesh(ANGLE_Q,WL_Q,iNormExpSRI, cmap = 'inferno', vmin=0, vmax = 1)
        ax1.pcolormesh(ANGLE_Q,WL_Q,NormSimSRI[ipos_min], cmap = 'inferno', vmin=0, vmax = 1)
        # Error plot
        errplot = ax2.pcolormesh(ANGLE_Q,WL_Q,np.abs(iNormExpSRI - NormSimSRI[ipos_min]), cmap = 'inferno')
        cbar = plt.colorbar(errplot)
        cbar.set_label('Error')
        
        ax.set_ylabel('Wavelength (nm)')
        fig.suptitle('Spectral Radiant Intensity ($W\cdot sr^{-1}\cdot nm^{-1}$)')
        ax.set_title('Experimental')
        ax1.set_title(f'Simulated pos = {pos_min:.2f}')
        ax.set_xlabel('Angle (°)')
        ax1.set_xlabel('Angle (°)')
        ax2.set_xlabel('Angle (°)')
        
        fig, ax = plt.subplots()
        ax.plot(i_POS, Error_Landscape, label = f'd = {thickness:.0f}\nmin_pos = {i_POS[ipos_min]:.2f}')
        ax.set_xlabel('Rel.position of the emitter')
        ax.set_ylabel('$\Delta_{sim}^{meas}$')
        ax.set_title("Error landscape for emitter position")
        ax.legend()
        
        fig, ax = plt.subplots()
        
        for i in range(0,len(ANGLE_Q), 1):
            ax.plot(WL_Q, NormSimSRI[ipos_min][:,i],'--C' + str(i))
            ax.plot(WL_Q, iNormExpSRI[:,i],'-C' + str(i), label = f'{ANGLE_Q[i]:}°')
            ax.set_xlabel('Wavelength (nm)')
            ax.set_ylabel('SRI (a.u.)')
            ax.set_title(f'Experimental and simulated emission spectra\n dAL = {thickness}, pos = {pos_min:.2f}')
        ax.legend()
    
    return i_POS, Error_Landscape, WL_Q, ANGLE_Q, pos_min, iNormExpSRI,NormSimSRI[ipos_min]


def fit_thickness(file, plot = False):
    """Finds the thickness corresponding to the minimum error using the error_landscape function.
        dmin : minimum thickness to start the fit with
        dmax : maximum thickness to end the fit with """
    error_pos = np.zeros(dAL_SIM.shape)
    
    for i,d in enumerate(dAL_SIM):
        pos_emitter, Error_Landscape,_,_,_,_,_= error_landscape(file, d)
        error_pos[i] = Error_Landscape.min()
    
    # The index with the minimum error from all the thicknesses tried
    fitted_thickness = dAL_SIM[error_pos.argmin()]
    
    if plot:
        fig, ax = plt.subplots()
        ax.plot(dAL_SIM, error_pos,'o-', label =  f'Min. Err. Thickness = {fitted_thickness} nm')
        ax.set_xlabel("Thickness of the AL (nm)")
        ax.set_ylabel('$\Delta_{sim}^{meas}$')
        ax.set_title("Error landscape for AL thickness")
        ax.legend()
        
    return fitted_thickness

def min_error_profile(weights, emitter_positions, exp_data, thickness, fitting = True):
    """Calculates the error with respect the experimental data assuming a linear combination of emmitters at different positions and with different weights."""   
    
    if len(weights) != len(emitter_positions):
        raise Exception('Length of the positions and weights do no match')
    
    # Simulated data as a linear combination of emitters in multiple positions
    lc_SimData = np.zeros(exp_data.shape) # linear combination SimData
    
    for w, pos in zip(weights, emitter_positions):
        lc_SimData += w * SIM_DATA[get_position_index(pos), get_thickness_index(thickness)]
    
    error = ((np.abs(exp_data - lc_SimData)).mean(axis = -1)).mean(axis = -1)
#     error = ((np.sqrt((exp_data - SimData)**2)).mean(axis = -1)).mean(axis = -1)
#     print (error)
    
    if fitting:
        return error
    else:
        return error, lc_SimData