#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:20:54 2024

@author: pi
"""
from scipy.interpolate import interp1d
import numpy as np
from pathlib import Path

local_folder = Path(__file__).parent

def check_a_multiple_b(a, b, precision=4, tolerance=1e-4):
    return np.isclose(0.0e-16, np.round(a, precision) % np.round(b, precision), tolerance)

def find_symmetry(angles, intensity):
    """
    Uses the integrated radiant intensity (no need to be calibrated, from the spectra) to calculate the symmetry point.

    Parameters
    ----------
    wavelengths : 1D np.array
        Array with the input wavelengths, [nm] should match the length of the first dimension of the sri. 
    angles : 1D np.array
        Array with the input angles, should match the length of the second dimension of the sri.
    spectra : 2D np.array
        The array with the spectral radiant intensity in aribitray units [a.u.] (wavelength as rows) for each of the collected angles (as columns).

    Returns
    -------
    angle_offset: float
        The offset by which the data should be corrected.
    """
    angles = np.array(angles)
    intensity = np.array(intensity)

    sorting = np.argsort(angles)
    angles = angles[sorting]
    intensity = intensity[sorting]


    angles = np.round(angles, 2)
    
    # Get rid of the annoying three measurments at zero, the interpolation doesn't like it!
    N = angles.shape[0] // 2
    fixed_angles = np.concatenate([angles[0:N-1], [angles[N-1:N+2].mean()], angles[N+2:]])
    fixed_intensity = np.concatenate([intensity[0:N-1], [intensity[N-1:N+2].mean()], intensity[N+2:]])
    fixed_intensity /= fixed_intensity.max()
    
    # Generate a fine angle data using the interpolation
    func = interp1d(fixed_angles, fixed_intensity, kind = 'quadratic')
    
    step = 0.1125
    fine_angles = np.round(np.arange(angles[0], angles[-1]+step/2, step), 2)
    fine_intensity = func(fine_angles)
    
    # Assume the offset is not going to be greater than 80, else, the data would be pretty fucked
    fine_offset = np.arange(angles.min()+5.4, angles.max()-5.4+step/2, step)
    N = len(fine_offset)
    offsets = np.zeros((N, ))
    errors = np.zeros((N, ))
    
    # Find the offset that gives the minimum difference between the 0-90 deg data and the -90 - 0 deg
    
    for i, fa in enumerate(fine_offset):
        offset = fa
        x_new = fine_angles + offset

        # get the positive angles data
        ff_pos = x_new >= 0
#        pos_angles = x_new[ff_pos]
        positive_intensity = fine_intensity[ff_pos] 
        # and the negative
        ff_neg = x_new <= 0
#        neg_angles = np.flip(x_new[ff_neg])
        negative_intensity = np.flip(fine_intensity[ff_neg])

        # Subtract only the common angles (on is going to be larger than the other)
        M_max = min(positive_intensity.shape[0], negative_intensity.shape[0])

        offsets[i] = offset
        errors[i] = np.mean(np.abs(positive_intensity[:M_max] - negative_intensity[:M_max]))

    symmetry_point = np.round(offsets[errors.argmin()], 4)
    
    return symmetry_point


def read_config_file(file):
    config_lines = []
    with open(file, 'r') as f:
        for l in f.readlines():
            config_lines.append(l)
                
    polarization_state  = None
    IRF_basename = None
    abs_calfactor = None
    calibration_folder = local_folder / 'calibration_files'
           
    for l in config_lines:
        if 'Polarizer state' in l:
            polarization_state = l.split(':')[-1].strip()
        elif 'Instrument response function' in l:
            IRF_basename = l.split(':')[-1].strip()
        elif 'Luminance calibration factor' in l:
            abs_calfactor = float(l.split(':')[-1].strip())
    
    if any(var is None for var in [polarization_state, IRF_basename, abs_calfactor]):
        print("AT least one of these variables is not defined: Polarizer state, Instrument response function, Luminance calibration factor")
        
    IRF_file = calibration_folder / (IRF_basename + '_'+ polarization_state + '-pol.txt')   

    
    return polarization_state, IRF_file, abs_calfactor