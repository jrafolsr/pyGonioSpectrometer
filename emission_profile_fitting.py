# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:21:07 2020

@author: JOANRR
"""
from pyGonioSpectrometer.data_processing import load_sri_file,error_landscape, min_error_profile, interpolate_expdata, load_simdata, gci
import numpy as np
from scipy.optimize import least_squares
from matplotlib import pyplot as plt

# Fit profile
filename = r'C:\Users\JOANRR\Documents\p-n-junction-shifting-with-temperature\calibration_files\specGon_testdata_Mattias_newformat.sri'

el_file = r'.\data_processing\calibration_files\Errorlandscape_sim.mat'

simEL = simEL = load_simdata(el_file)

wl_sim = simEL['wl']
angles_sim = simEL['angles']

thickness = 380
step = 1

initial = step - 1
final = -initial if initial != 0 else None
slicing = slice(initial, final, step)
positions = simEL['ipos'][slicing]

wavelengths, angles, sri = load_sri_file(filename)
exp_data = interpolate_expdata(sri,wavelengths, angles, wl_sim, angles_sim)

ithickness = gci(thickness, simEL['dAL'])
sthickness = simEL['dAL'][ithickness]
print(f'Taking simulated data for thickness = {sthickness:.0f} nm')

simData_fixdAL = simEL['data'][slicing, ithickness]




w0 = np.zeros(positions.shape)
w0 = w0 + 1/len(w0)


res_lsq  = least_squares(min_error_profile, w0, bounds=(0, 1),\
                     args = (simData_fixdAL, exp_data),\
                         method = 'dogbox', verbose = 1)#, jac = '3-point')

w = res_lsq.x

# Initial error
error0, SimData0 = min_error_profile(w0, simData_fixdAL, exp_data,fitting = False)
    
# Profile
errorf, SimDataf = min_error_profile(w, simData_fixdAL, exp_data, fitting = False)

# Single delta fit
_, pos_min, _, simData_single  = error_landscape(filename, thickness, simEL)

error_s, SimData_s = min_error_profile([1], [simData_single], exp_data, fitting = False)

text =  f'Thickness is {thickness:.0f}\n'
text += f'Single delta position = {pos_min}\n'
text += f'Error with a delta = {error_s}\n' 
text += f'Error with profile = {errorf}\n'
text += f'Error with homogenous profile {error0}\n'


print(text)



fig, [ax, ax1] = plt.subplots(ncols=2, figsize = (8,4))

for i,a in enumerate(angles_sim):
    offset = (len(angles_sim) - i) *.25 
    ax.plot(wl_sim,offset +  SimData_s[:,i], '-.C' + str(i), lw = 1)
    ax.plot(wl_sim,  offset +  SimDataf[:,i], '--C' + str(i), lw = 1)
    ax.plot(wl_sim,  offset + exp_data[:,i], '-C' + str(i), label = f'{a:}°')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Norm. SRI (a.u.)')
    ax.set_title('Fittings')

ax.legend()

ax1.set_xlabel('Rel.position of the emitter')
ax1.set_ylabel('Weight')
ax1.set_title("Fitted profile")
ax1.set_xlim(0,1)
#     ax1.set_ylim(0,0.075)
# ax.bar([.69], [1], width = 0.02)
ax1.bar(positions, w, width = 0.02)
ax1.bar(positions, w0, width = 0.005)
ax1.axvline(pos_min, c= 'black', ls = '--')
plt.subplots_adjust(wspace = 0.3)

fig.suptitle(f'Thickness = {thickness:.0f} nm')