# -*- coding: utf-8 -*-
"""
GonioLogger -- spectro-goniometer acquisition class for the OPEGLAB setup.

Wraps a Flame/Ocean spectrometer (`SpectraMeasurement`) and the Raspberry-Pi
motor/shutter controller (`RaspberryMotorController`) to:

    * take dark and signal spectra, with automatic integration-time tuning;
    * run a full angular scan over the forward hemisphere (or half of it),
      writing one spectrum per angle to a `*_gonio.dat` file;
    * track a forward "luminance" time-series in a separate `*_luminance.dat`
      file (the angle-0 spectrum at each pass).

Data layout of the output files: row 1 is the wavelength vector, every
subsequent row is `[time, parameter, counts...]`. Counts are background-
subtracted unless noted otherwise.

CHANGED (vs the original tasks.py):
    * Plotting uses a non-GUI Agg backend and saves one PNG snapshot per scan
      instead of an interactive Tk window (the Tk window crashed when run
      alongside the logger thread -- 'Tcl_AsyncDelete').
    * Acquisition methods no longer swallow exceptions and continue with bad
      data; they propagate so the calling script can park the hardware.
    * Dead/commented-out code removed; docstrings and section banners added.

Originally created 2020-06-12 by OPEGLAB; cleanup pass kept behaviour-
equivalent except where flagged with `# CHANGED:`.
"""
# =============================================================================
# Imports and module constants
# =============================================================================
from pyGonioSpectrometer.instrumentation import list_spectrometers, SpectraMeasurement, RaspberryMotorController
from datetime import datetime
from time import sleep, monotonic
import numpy as np
from pathlib import Path
# CHANGED: force a non-GUI backend BEFORE pyplot is imported. This module runs
# concurrently with the logger daemon thread; any Tk-backed figure can have its
# Tcl image objects finalized by the GC off the main thread -> 'Tcl_AsyncDelete'
# hard crash. Agg has no Tk, so that failure class is impossible.
# HW/ENV NOTE: this makes ALL matplotlib in this process non-interactive (no
# popup windows anywhere in this Python session). Fine for the measurement
# script; if you ever import this module in an interactive Spyder kernel that
# already opened figures, run it in a fresh process instead.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# CHANGED: create figures explicitly off pyplot (no global figure manager, no
# Tk) so saving is fully thread-safe regardless of backend/import order.
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from pyGonioSpectrometer.gui_init import WAIT_TIME
from pyGonioSpectrometer import check_a_multiple_b, read_config_file
import seaborn as sns
import traceback

SATURATION_COUNTS = 65535   # Hard saturation limit of the spectrometer (counts)
UPPER_LIM = 58000           # Above this, reduce the integration time (counts)
LOWER_LIM = 10000           # Below this, increase the integration time (counts)

PIXEL_SIZE = 4e-6           # McScience substrate pixel area (m^2)

N_WAVELENGTHS = 2028        # Number of pixels/wavelengths returned by the spectrometer

# =============================================================================
# Output folders
# =============================================================================
DEFAULT_OUTPUT_FOLDER = Path.home() / 'Documents/data/goniologger'
if not DEFAULT_OUTPUT_FOLDER.is_dir():
    DEFAULT_OUTPUT_FOLDER.mkdir()

# CHANGED: scratch folder for the per-scan PNG snapshots (safe to bulk-delete).
TEMP_PLOT_FOLDER = DEFAULT_OUTPUT_FOLDER / 'temp'
if not TEMP_PLOT_FOLDER.is_dir():
    TEMP_PLOT_FOLDER.mkdir()


# =============================================================================
# GonioLogger
# =============================================================================
class GonioLogger():
    def __init__(self, filename, folder = DEFAULT_OUTPUT_FOLDER, angle_step = 10.8, angle_max = 86.4, max_time_per_angle = 5000, integration_time = 100, n_spectra = 1, max_n_spectra= 5, max_intensity_angle = np.round(0.0, 4), suffix_luminance_file = ''):
        """
        Parameters
        ----------
        filename : str
            A string containing the output filename. The default is 'gonio_measurement'.
        folder : str
                Path, relative or absolute, to the directory where to save the data. It should exists, the program does not check if it does.
        angle_step : int or float, optional
            Angular step that the motor will perform in degrees. Preferably a divisor of angle_max and an integer, the program does not check for these conditions to be fullfilled. The default is 5.4.
        angle_max : int or float, optional
            Maximum angle to scan with the goniometer in degrees. The default is 86.4.
        integration_time: int or float
            Sets the integration time in ms.
        n_spectra : int
            Number of spectra that will be averaged.
        max_time_per_angle : int, optional
            Maximum time in miliseconds allowed per angle step (number of spectra x integration time). The default is 5000 ms.
        """
        self.filename = filename
        # CHANGED: running index for the per-scan PNG snapshots in TEMP_PLOT_FOLDER.
        self._plot_counter = 0
        self.calibration_folder = Path(__file__).parent / 'calibration_files'
        self.file_config = Path(__file__).parent / 'local-config.txt'
        if not self.file_config.exists():
            raise ValueError('A file-config must exist!')
        self.metadata = ''
        with open(self.file_config) as fc:
            for l in fc.readlines():
                self.metadata+= '# ' + l
                if not l.endswith('\n'):
                    self.metadata+='\n'
        # Check if angle_step is a multiple of 1.8
        if not check_a_multiple_b(angle_step, 1.8):
            raise ValueError('OBS! The angle_step must be a multiple of 1.8 deg!')
        
        # Check if angle_max is a multiple of angle_step
        if not check_a_multiple_b(angle_max, angle_step):
            raise ValueError('OBS! The angle_max must be a multiple of angle_step deg!')
        
        itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
        timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
        self.luminance_filename = folder / (timestamp  + filename + ('_'+suffix_luminance_file if len(suffix_luminance_file) else '')  + '_luminance.dat')
        with open( self.luminance_filename , 'w') as f:
            f.write(self.metadata)
            f.write(f'# Timestamp at the beginning of the measurement: {itimestamp}\n')
            f.write(f'# The first row is the wavelength vector\n')    
            f.write('# Rel.time\t Integration time\t Intensities\n')
        
        self.folder = folder
        self.angle_step = angle_step
        self.angle_max = angle_max
        self.max_time_per_angle = max_time_per_angle
        self.integration_time = integration_time
        self.n_spectra = n_spectra
        self.max_n_spectra = max_n_spectra
        self.lower_lim=LOWER_LIM
        self.upper_lim=UPPER_LIM
        
        self.background = np.zeros((N_WAVELENGTHS, ))*np.nan
        self.intensities = np.zeros((N_WAVELENGTHS, ))*np.nan
        self.intensities_raw = np.zeros((N_WAVELENGTHS, ))*np.nan

        self.max_intensity_angle = np.round(max_intensity_angle,4) # Define the angle that offers the maximum intensity. We assume 0.0 (forward) as the nitial guess (lambertian)
        self.max_total_counts = 0

        self.current_angular_scan = []
        # Load configuration from file_config
        self.polarization_state, self.IRF_file, self.abs_calfactor = read_config_file(self.file_config)
        if not self.IRF_file.exists():
            raise ValueError(f'IRF file "{self.IRF_file}" does not exist!')


        # Take the first spectromter from the list (it assumes it will be only one)
        if len(list_spectrometers()):
            self.name_spectrometer = list_spectrometers()[0]
        else:
            raise ValueError("No spectrometer in the system, check the connections and try again.")
            
        try: 
            self.flame = SpectraMeasurement(self.name_spectrometer,\
                                           integration_time = self.integration_time,\
                                           n_spectra = self.n_spectra,
                                           max_n_spectra= self.max_n_spectra)
            self.wavelengths = self.flame.get_wavelengths()
            
            # Initializer for data processing
            self.filter_wavelengths = (self.wavelengths >=350) & (self.wavelengths <= 950)
            self.wavelengths_cutted = self.wavelengths[self.filter_wavelengths]
            self.photopic_eye_response = None # I am going to load it later
            
            
            self.write_to_file(np.nan, integration_time, self.wavelengths, self.luminance_filename, fmt_data='%8.2f')
        
        except Exception:
            # CHANGED: fail fast. The original printed and continued, which left
            # self.flame / self.wavelengths undefined and produced confusing
            # AttributeError crashes much later. A spectrometer that won't
            # initialize is fatal -- surface it here.
            traceback.print_exc()
            raise
            
            
        # Make sure the angle step divides the max angle, otherwise raise.
        # CHANGED: removed dead commented-out divmod/quotient variants.
        if not check_a_multiple_b(self.angle_max, self.angle_step):
            raise Exception("The angle step is not a divisor of the max. angle!")
        
    
        # Initialize counters and flags and gonio
        self.gonio = RaspberryMotorController()
        
        
    def take_dark_spectra(self,  debug = False):
        """Acquire the dark/background spectrum with the shutter closed and
        store it in `self.background`. The shutter is reopened afterwards.

        CHANGED: on a read error the shutter is reopened and the exception is
        re-raised (the original swallowed it, leaving a stale/NaN background and
        letting the measurement continue with bad data).
        """
        if debug: print(f'INFO: Taking N = {self.n_spectra:d} dark spectra with integration time of {self.integration_time:.0f}')
        if not self.gonio.shutter_is_closed:
            self.gonio.close_shutter()
        
        try:
            # Take the dark spectra with the shutter closed
            self.background = self.flame.get_averaged_intensities()
        except Exception:
            print('ERROR in take_dark_spectra(): could not read the spectrometer.')
            traceback.print_exc()
            if self.gonio.shutter_is_closed:
                self.gonio.open_shutter()   # leave the setup in a known state
            raise

        if self.gonio.shutter_is_closed:
            self.gonio.open_shutter()
            
            
    def update_integration_time(self, angle = None, factor=0.9):
        """Auto-tune the integration time (and number of averaged spectra) so the
        peak counts sit between LOWER_LIM and UPPER_LIM, then refresh the dark
        spectrum. If `angle` is given, the gonio is moved there for the tuning
        and moved back afterwards (use the angle of maximum intensity to avoid
        saturating elsewhere).

        CHANGED: exceptions propagate instead of being swallowed.
        """
        print(f'\nUpdating the integration time and number of spectra')
        if self.gonio.shutter_is_closed:
            self.gonio.open_shutter()
        
        try:
            if angle is not None:
                self.gonio.move_angle(angle)
                sleep(1.0)
            
            self.integration_time, self.n_spectra = self.flame.adjust_integration_time(max_time = self.max_time_per_angle,\
                                                            lower_limit = self.lower_lim,
                                                            upper_limit = self.upper_lim,
                                                            factor=factor)
            if angle is not None:
                self.gonio.move_angle(-angle)
                sleep(1.0)
            # Take and save a new dark spectra
            self.take_dark_spectra()

        except Exception:
            print('ERROR in update_integration_time().')
            traceback.print_exc()
            raise
    
    def take_spectra(self, debug = False):
        """Acquire one (background-subtracted) signal spectrum into
        `self.intensities`; the raw spectrum is kept in `self.intensities_raw`.
        Takes a dark spectrum first if none has been acquired yet.

        CHANGED: exceptions propagate instead of being swallowed.
        """
        if debug: print(f'INFO: Taking N = {self.n_spectra:d} spectra with integration time of {self.integration_time:.0f}')
        
        if self.gonio.shutter_is_closed:
            self.gonio.open_shutter()
        
        try:
            # Take a dark spectra if not taken
            if np.all(np.isnan(self.background)):
                self.take_dark_spectra()
            
            sleep(0.010)            
            self.intensities_raw = self.flame.get_averaged_intensities()           
            self.intensities = self.intensities_raw- self.background
            
        except Exception:
            print('ERROR in take_spectra().')
            traceback.print_exc()
            raise
    
    def set_nspectra(self, n_spectra):
        """Set the number of spectra to average (on both this object and the
        underlying `flame`) and refresh the dark spectrum."""
        self.flame.n_spectra = self.n_spectra = n_spectra
        self.take_dark_spectra()
    
    def save_spectra(self, etime = np.nan):
        """Acquire a forward spectrum and append it to the luminance time-series
        file with the given elapsed time `etime`."""
        self.take_spectra()
        
        self.write_to_file(etime, self.integration_time, self.intensities, self.luminance_filename)
        
        
    def shutdown(self):
        """Close the shutter and release the gonio and spectrometer resources.
        Call once at the very end of a measurement."""
        if not self.gonio.shutter_is_closed:
            self.gonio.close_shutter()
        
        self.gonio.close()
        self.flame.close()
        
        
    def write_to_file(self, etime, parameter, data, file, cut = True, fmt_data = '%8.0f'):
        """Append one row `[etime, parameter, data...]` to `file`.

        Parameters
        ----------
        etime, parameter : float
            Scalars written as the first two columns (e.g. elapsed time and
            angle, or elapsed time and integration time).
        data : np.ndarray
            The spectrum (counts). If `cut` is True it is restricted to the
            350-950 nm window (`self.filter_wavelengths`) before writing.
        file : path-like
            Output file, opened in append mode.
        fmt_data : str
            numpy `savetxt` format for the data columns.
        """
        if cut:
            data = data[self.filter_wavelengths]
            
        t  = np.hstack((etime, parameter, data))
        t = t.reshape((1, t.shape[0]))
        with open(file, 'a') as f:
           np.savetxt(f, t, fmt = ['% 8.2f', '% 8.2f'] + [fmt_data]*data.shape[0]) 
        
    def check_low_and_high_limits(self):
        """Inspect the last raw spectrum against the count limits.

        Returns
        -------
        int
            ``1``  if any pixel is above UPPER_LIM (too bright -> lower the
                   integration time),
            ``-1`` if the peak is below LOWER_LIM (too dim -> raise it),
            ``0``  if within band (or if no spectrum has been taken yet).
        """
        if np.all(np.isnan(self.intensities)):
            all_good = 0 # No spectrum taken yet -> nothing to act on.
        else:
            if np.any(self.intensities_raw > self.upper_lim):
                print('\n! WARNING: Some values close to the saturation limit. Consider lowering the integration time.')
                all_good =  1
            elif self.intensities_raw.max() < self.lower_lim:
                print(f'\n! WARNING: The max. count is less than {self.lower_lim}. Consider increasing the integration time')
                all_good = -1
            else:
                all_good = 0
                
        return all_good
    
    def take_gonio_measurement(self, suffix = '', header = '', disable_gonio = False, plot = True, parameter_1 = np.nan,half_hemisphere = False):
        """
        Performs a complete measurement for the goniospectrometer setup, by taking spectra at every specified angle, for the whole forward hemisphere.
        
        Parameters
        ----------
        angle_max : int or float
            Maximum angle to scan with the goniometer in deg.
        angle_step : int or float
            Angular step that the motor will perform in deg. Preferably a divisor of angle_max and integer, the program does not check for this conditions to be fullfilled.
        gonio: RaspberryMotorController or ArduinoMotorControllerobject
            Object to controll the motor.
        name_spectrometer : seabreeze.spectrometers.Spectrometer class
            The spectrometer resource as the specified class. You can get it from list_spectrometers().
        integration_time : int or float
            Sets the integration time in ms.
        n_spectra: int
            Number of spectra that will be averaged.
        filename : str, optional
            A string containing the output filename. The default is 'gonio_measurement'.
        folder : str, optional
            Path, relative or absolute, to the directory where to save the data. Should exists, the program does not check if it does. The default is '.'.
        disable_gonio : boolean, optional
            Whether to disable the goniometer motor after the measruement. The default is False.
        plot : boolean, optional
            Whether to plot or not the data. It is better to set it to false for a time series measurement, otherwise one will end with to many open windows. The default is True.
    
        """ 

        # Initalizing some variables
        self.current_angular_scan = []
        current_angle = 0.0
        
        
        try:
            if np.all(np.isnan(self.background)):
                print('INFO: Taking dark spectra....')
                self.take_dark_spectra()
            
            
            # Create the object instance taht will control the spectrometer, assuming it is the
    #        first of the list

            n_angles = int(round(self.angle_max / self.angle_step, 0)) + 1 # Gi e wird results if the factr I use is 100
            
            n_steps = (n_angles - 1) * 2
    
            # Prepraring the plot
            if plot:
                colors = sns.color_palette('rainbow', n_colors= n_angles)
                # CHANGED: build the figure explicitly off pyplot (Agg canvas) so no
                # Tk object is ever created -> no off-thread Tcl_AsyncDelete crash.
                # The snapshot is re-saved after every angle to a single numbered
                # PNG per scan; open it in a viewer that auto-reloads for live-ish
                # feedback. Old files are safe to delete at any time.
                fig = Figure()
                FigureCanvasAgg(fig)
                ax = fig.add_subplot(111)
                ax.set_ylabel('Counts')
                ax.set_xlabel('Wavelength (nm)')
                self._plot_counter += 1
                fig._gonio_savepath = TEMP_PLOT_FOLDER / f'{self.filename}_{self._plot_counter:04d}.png'
                plot_measurement(fig, ax, [], [])
            
            # Timestamps for the header and filename
            itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
            timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
            start_time = monotonic()
            
            path = self.folder / (timestamp + self.filename + ('_' + suffix if len(suffix) else '') + '_gonio.dat')
                        
            with open(path, 'a') as f:
                f.write(self.metadata)
                header = header if header.startswith('#') else '# ' + header
                header = header if header.endswith('\n') else header + '\n'
                f.write(header)
                f.write(f'# Timestamp at the beginning of the measurement: {itimestamp}\n')
                f.write(f'# Integration time in (ms): {self.integration_time:.0f}\n')
                f.write(f'# Number of spectra taken: {self.n_spectra:d}\n')
                f.write(f'# Time(s) Angle(deg) Wavelength(nm) -> (row 1) Counts() -> from row 2\n')
            
            
            
            # Saving the data in the new scheme
            self.write_to_file(monotonic()-start_time, np.nan, self.wavelengths, path, fmt_data = '%8.2f')
            # I am going to save the data with the bkg already substracted
#            self.write_to_file(monotonic()-start_time, np.nan, self.background, path)
          
#            if plot: plot_measurement(fig, ax, self.wavelengths, self.background, 'dark', color = 'black')

            sleep(0.05)
            
            # Open the shutter if it is closed
            if self.gonio.shutter_is_closed:
                self.gonio.open_shutter()
    
            sleep(WAIT_TIME)
            
            # Take spectra at zero
            print(f'\rINFO Step # 0, moved  {0.0: >4.1f}°, position {0.0: >+5.1f}° |' + ' '* (n_steps + 3) + f'| {1/(n_steps + 3)*100:3.0f} % Done...', end =''  )
                  
            temp = self.flame.get_averaged_intensities()
            # CHANGED: removed a dead branch here. The original did
            #   self.max_total_counts = temp.max()
            #   if temp.max() > self.max_total_counts:   # never true
            #       ... update max_intensity_angle ...
            # The condition can never fire right after the assignment, so it was
            # dead code. We keep the baseline assignment; the per-angle loop
            # below is what actually tracks the angle of maximum intensity.
            self.max_total_counts = temp.max()   # baseline = counts at 0 deg
            
            temp = temp - self.background   # subtract the background
            
            self.current_angular_scan.append([current_angle, temp])
            
            # Saving the data in the new scheme
            tt = monotonic()-start_time
            self.write_to_file(tt, 0.0, temp, path)
            self.write_to_file(parameter_1 + tt, self.integration_time, temp, self.luminance_filename)
            
            if plot: plot_measurement(fig, ax, self.wavelengths, temp, '0.0°', color = colors[0])
            
    
            # Move to the last position
            out_angle = self.gonio.move_angle(-1.0 * self.angle_max)

    #        # Initialize the error made in each step
    #        error = (angle_max - out_angle)
            total = 0
            current_angle = -out_angle
            print(f'\rINFO: Step # 1, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |#' + ' '* (n_steps + 2) + f'| {2/(n_steps + 3)*100:3.0f} % Done...', end =''  )
            # Wait longer time, as the angle is larger and take the spectra 
            sleep(WAIT_TIME  + 1.0) # Long enough time to make sure that it waits until the ed of the movement
            
            k = 0
            for k in range(n_steps):        
                # Get the whole spectra (wl and I)
                temp = self.flame.get_averaged_intensities()
 
                if temp.max() > self.max_total_counts:
                    print(f'\nMaximum intensity angle updated from {self.max_intensity_angle:.4f} to {np.round(current_angle, 4):.4f}')
                    self.max_intensity_angle = np.round(current_angle, 4)
                    self.max_total_counts = temp.max()
                
                temp = temp - self.background #  Substract the background
                self.current_angular_scan.append([current_angle, temp])
                # Saving the data in the new scheme
                tt = monotonic()-start_time
                self.write_to_file(tt, current_angle, temp, path)
                
                if round(current_angle, 4) == round(0.0, 4):
                    self.write_to_file(parameter_1 + tt, self.integration_time, temp, self.luminance_filename)
            
                # Check for any values higher than saturation
                if np.any(temp > SATURATION_COUNTS): print('WARNING: Some values are saturating. Consider lowering the integration time.')
                
                if plot and k % 2 == 0:
                    color_offset = -1-k if k < n_angles else k+1-n_angles
#                    print(color_offset, k, n_angles)
                    plot_measurement(fig, ax, self.wavelengths, temp, f'{current_angle:.1f}°', color = colors[color_offset])
                
                # Moving the gonio
                out_angle = self.gonio.move_angle(self.angle_step) 
                
                total += abs(out_angle)
    
                current_angle += out_angle * np.sign(self.angle_step) # The sign is for the case in which I only do the half hemisphere.s
                
                print(f'\rINFO: Step #{k+2:2d}, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |' + '#'* (k + 2) + ' '* (n_steps - k + 1) + f'| {(k+3)/(n_steps + 3)*100:3.0f} % Done...' , end =''  )
                if half_hemisphere & (k == n_steps // 2 -1):
                    print(f'\nReversing the angle step at {current_angle: >+5.1f}°, only taking half the hemisphere')
                    self.angle_step *= -1
                sleep(WAIT_TIME)
            
#            Setting back the angle step, to its original value if the half hemisphere has been set
            if half_hemisphere:
                self.angle_step *= -1
                
            # Take last angle spectra
            temp = self.flame.get_averaged_intensities() - self.background
            self.current_angular_scan.append([current_angle, temp])
            # Saving the data in the new scheme
            self.write_to_file(monotonic()-start_time, current_angle, temp, path)
            
            # Plot the data
            if plot: plot_measurement(fig, ax, self.wavelengths, temp, f'{current_angle:.1f}°', color = colors[-1])

            
            # Go back to zero
            back_angle = -1 * abs(current_angle) if not half_hemisphere else abs(current_angle) # Quick fix for a scanning the same half hemisphere
            # Moving back the exact angle we moved to set everything to the initial position, so correct_drift = False
            out_angle = self.gonio.move_angle(back_angle, correct_drift = False)
            current_angle = current_angle - out_angle if not half_hemisphere else current_angle + out_angle
            
            print(f'\rINFO: Step #{k+3:2d}, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |' + '#'* (n_steps + 3) + f'| {(k+4)/(n_steps + 3)*100:3.0f} % Done...' , end ='\n'  )
            
            # Wait longer time, as the angle is larger and take the spectra
            sleep(WAIT_TIME  + 1.0)
            temp = self.flame.get_averaged_intensities()- self.background
            self.current_angular_scan.append([current_angle, temp])
            if plot: plot_measurement(fig, ax, self.wavelengths, temp, f'{current_angle:.1f}°', color = colors[0])
            
            # Saving the data in the new scheme
            tt = monotonic()-start_time
            self.write_to_file(monotonic()-start_time, current_angle, temp, path)
            self.write_to_file(parameter_1 + tt, self.integration_time, temp, self.luminance_filename)
            

            if disable_gonio: self.gonio.disable_gonio()
              
            
            if plot:
                # CHANGED: figure is not pyplot-managed, so plt.close() is a no-op.
                # Just drop the reference to free it.
                del fig
    
        except KeyboardInterrupt:
            print('INFO: The angle-scan has been cancelled by the user. Going back to 0°.')
            if self.gonio != None:
                # Go back to since the spectrogoniometer movement has been cancelled.
                back_angle = -1 * current_angle          
                out_angle = self.gonio.move_angle(back_angle)
                
        except Exception as e:
            print(e)
            traceback.print_exc()
            print('INFO: Some error has ocurred during the angle-scan. Going back to 0°.')
            if self.gonio != None:
                # Go back to since some error has occurred during the gonio measurement
                back_angle = -1 * current_angle
                out_angle = self.gonio.move_angle(back_angle)

        finally:
                pass
    
    def take_gonio_measurement_halfhemisphere(self, suffix = '', header = '', disable_gonio = False, plot = True, parameter_1 = np.nan,full_hemisphere = True):
        """
        Performs a complete measurement for the goniospectrometer setup, by taking spectra at every specified angle, for the whole forward hemisphere.
        
        Parameters
        ----------
        angle_max : int or float
            Maximum angle to scan with the goniometer in deg.
        angle_step : int or float
            Angular step that the motor will perform in deg. Preferably a divisor of angle_max and integer, the program does not check for this conditions to be fullfilled.
        gonio: RaspberryMotorController or ArduinoMotorControllerobject
            Object to controll the motor.
        name_spectrometer : seabreeze.spectrometers.Spectrometer class
            The spectrometer resource as the specified class. You can get it from list_spectrometers().
        integration_time : int or float
            Sets the integration time in ms.
        n_spectra: int
            Number of spectra that will be averaged.
        filename : str, optional
            A string containing the output filename. The default is 'gonio_measurement'.
        folder : str, optional
            Path, relative or absolute, to the directory where to save the data. Should exists, the program does not check if it does. The default is '.'.
        disable_gonio : boolean, optional
            Whether to disable the goniometer motor after the measruement. The default is False.
        plot : boolean, optional
            Whether to plot or not the data. It is better to set it to false for a time series measurement, otherwise one will end with to many open windows. The default is True.
    
        """ 

        # Initalizing some variables
        self.current_angular_scan = []
        current_angle = 0.0
        
        
        try:
            if np.all(np.isnan(self.background)):
                print('INFO: Taking dark spectra....')
                self.take_dark_spectra()
            
            
            # Create the object instance taht will control the spectrometer, assuming it is the
    #        first of the list

            n_angles = int(round(self.angle_max / self.angle_step, 0)) + 1 # Gi e wird results if the factr I use is 100
            
            n_steps = (n_angles - 1) * 2 if full_hemisphere else (n_angles - 1)
    
            # Prepraring the plot
            if plot:
                colors = sns.color_palette('rainbow', n_colors= n_angles)
                # CHANGED: same Agg-only figure as take_gonio_measurement (no Tk).
                fig = Figure()
                FigureCanvasAgg(fig)
                ax = fig.add_subplot(111)
                ax.set_ylabel('Counts')
                ax.set_xlabel('Wavelength (nm)')
                self._plot_counter += 1
                fig._gonio_savepath = TEMP_PLOT_FOLDER / f'{self.filename}_{self._plot_counter:04d}.png'
                plot_measurement(fig, ax, [], [])
            
            # Timestamps for the header and filename
            itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
            timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
            start_time = monotonic()
            
            path = self.folder / (timestamp + self.filename + ('_' + suffix if len(suffix) else '') + '_gonio.dat')
                        
            with open(path, 'a') as f:
                f.write(self.metadata)
                header = header if header.startswith('#') else '# ' + header
                header = header if header.endswith('\n') else header + '\n'
                f.write(header)
                f.write(f'# Timestamp at the beginning of the measurement: {itimestamp}\n')
                f.write(f'# Integration time in (ms): {self.integration_time:.0f}\n')
                f.write(f'# Number of spectra taken: {self.n_spectra:d}\n')
                f.write(f'# Time(s) Angle(deg) Wavelength(nm) -> (row 1) Counts() -> from row 2\n')
            
            
            
            # Saving the data in the new scheme
            self.write_to_file(monotonic()-start_time, np.nan, self.wavelengths, path, fmt_data = '%8.2f')
            # I am going to save the data with the bkg already substracted
#            self.write_to_file(monotonic()-start_time, np.nan, self.background, path)
          
#            if plot: plot_measurement(fig, ax, self.wavelengths, self.background, 'dark', color = 'black')

            sleep(0.05)
            
            # Open the shutter if it is closed
            if self.gonio.shutter_is_closed:
                self.gonio.open_shutter()
    
            sleep(WAIT_TIME)
            
            # Take spectra at zero
            print(f'\rINFO Step # 0, moved  {0.0: >4.1f}°, position {0.0: >+5.1f}° |' + ' '* (n_steps + 3) + f'| {1/(n_steps + 3)*100:3.0f} % Done...', end =''  )
                  
            temp = self.flame.get_averaged_intensities()
            # CHANGED: removed a dead branch here. The original did
            #   self.max_total_counts = temp.max()
            #   if temp.max() > self.max_total_counts:   # never true
            #       ... update max_intensity_angle ...
            # The condition can never fire right after the assignment, so it was
            # dead code. We keep the baseline assignment; the per-angle loop
            # below is what actually tracks the angle of maximum intensity.
            self.max_total_counts = temp.max()   # baseline = counts at 0 deg
            
            temp = temp - self.background   # subtract the background
            
            self.current_angular_scan.append([current_angle, temp])
            
            # Saving the data in the new scheme
            tt = monotonic()-start_time
            self.write_to_file(tt, 0.0, temp, path)
            self.write_to_file(parameter_1 + tt, self.integration_time, temp, self.luminance_filename)
            
            if plot: plot_measurement(fig, ax, self.wavelengths, temp, '0.0°', color = colors[0])
            
    
            # Move to the last position
            out_angle = self.gonio.move_angle(-1.0 * self.angle_max)

    #        # Initialize the error made in each step
    #        error = (angle_max - out_angle)
            total = 0
            current_angle = -out_angle
            print(f'\rINFO: Step # 1, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |#' + ' '* (n_steps + 2) + f'| {2/(n_steps + 3)*100:3.0f} % Done...', end =''  )
            # Wait longer time, as the angle is larger and take the spectra 
            sleep(WAIT_TIME  + 1.0) # Long enough time to make sure that it waits until the ed of the movement
            
            k = 0
            for k in range(n_steps):        
                # Get the whole spectra (wl and I)
                temp = self.flame.get_averaged_intensities()
 
                if temp.max() > self.max_total_counts:
                    print(f'\nMaximum intensity angle updated from {self.max_intensity_angle:.4f} to {np.round(current_angle, 4):.4f}')
                    self.max_intensity_angle = np.round(current_angle, 4)
                    self.max_total_counts = temp.max()
                
                temp = temp - self.background #  Substract the background
                self.current_angular_scan.append([current_angle, temp])
                # Saving the data in the new scheme
                tt = monotonic()-start_time
                self.write_to_file(tt, current_angle, temp, path)
                
                if round(current_angle, 4) == round(0.0, 4):
                    self.write_to_file(parameter_1 + tt, self.integration_time, temp, self.luminance_filename)
            
                # Check for any values higher than saturation
                if np.any(temp > SATURATION_COUNTS): print('WARNING: Some values are saturating. Consider lowering the integration time.')
                
                if plot and k % 2 == 0:
                    color_offset = -1-k if k < n_angles else k+1-n_angles
#                    print(color_offset, k, n_angles)
                    plot_measurement(fig, ax, self.wavelengths, temp, f'{current_angle:.1f}°', color = colors[color_offset])
                
                # Moving the gonio
                out_angle = self.gonio.move_angle(self.angle_step)
                
                total += abs(out_angle)
    
                current_angle += out_angle
                
                print(f'\rINFO: Step #{k+2:2d}, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |' + '#'* (k + 2) + ' '* (n_steps - k + 1) + f'| {(k+3)/(n_steps + 3)*100:3.0f} % Done...' , end =''  )
                
                sleep(WAIT_TIME)
            
            # Take last angle spectra
            temp = self.flame.get_averaged_intensities() - self.background
            self.current_angular_scan.append([current_angle, temp])
            # Saving the data in the new scheme
            self.write_to_file(monotonic()-start_time, current_angle, temp, path)
            
            # Plot the data
            if plot: plot_measurement(fig, ax, self.wavelengths, temp, f'{current_angle:.1f}°', color = colors[-1])

            
            # Go back to zero
            back_angle = -1 * abs(current_angle)
            # Moving back the exact angle we moved to set everything to the initial position, so correct_drift = False
            out_angle = self.gonio.move_angle(back_angle, correct_drift = False)
            current_angle -= out_angle
            
            print(f'\rINFO: Step #{k+3:2d}, moved {out_angle: >4.1f}°, position {current_angle: >+5.1f}° |' + '#'* (n_steps + 3) + f'| {(k+4)/(n_steps + 3)*100:3.0f} % Done...' , end ='\n'  )
            
            # Wait longer time, as the angle is larger and take the spectra
            sleep(WAIT_TIME  + 1.0)
            temp = self.flame.get_averaged_intensities()- self.background
            self.current_angular_scan.append([current_angle, temp])
            if plot: plot_measurement(fig, ax, self.wavelengths, temp, f'{current_angle:.1f}°', color = colors[0])
            
            # Saving the data in the new scheme
            tt = monotonic()-start_time
            self.write_to_file(monotonic()-start_time, current_angle, temp, path)
            self.write_to_file(parameter_1 + tt, self.integration_time, temp, self.luminance_filename)
            

            if disable_gonio: self.gonio.disable_gonio()
               
    
        except KeyboardInterrupt:
            print('INFO: The angle-scan has been cancelled by the user. Going back to 0°.')
            if self.gonio != None:
                # Go back to since the spectrogoniometer movement has been cancelled.
                back_angle = -1 * current_angle          
                out_angle = self.gonio.move_angle(back_angle)
                
        except Exception as e:
            print(e)
            traceback.print_exc()
            print('INFO: Some error has ocurred during the angle-scan. Going back to 0°.')
            if self.gonio != None:
                # Go back to since some error has occurred during the gonio measurement
                back_angle = -1 * current_angle
                out_angle = self.gonio.move_angle(back_angle)

        finally:
                pass
            
    def take_radiance_and_luminance(self):
        """Convert the current spectrum to absolute radiance and luminance.

        FIXME (flagged, NOT fixed -- needs your input): as written this method
        is broken and would raise. `self.IRF_file` is a pathlib.Path (the IRF
        *filename*), not the loaded instrument-response array, so the line
        `self.intensities * IRF` below multiplies an array by a Path. The IRF
        must be loaded from disk first, e.g.
            IRF = np.loadtxt(self.IRF_file)            # shape/columns?
        but I don't know the file's format/length (does it already match
        `self.wavelengths`, or need interpolation onto it?). This method is not
        called by any current script, so behaviour is unchanged -- tell me the
        IRF file layout and I'll wire it up. `abs_calfactor` units also need a
        bench check against a known reference.
        """
        self.take_spectra()
        # Load calibration files
        IRF = self.IRF_file   # FIXME: this is a Path, not the loaded IRF data (see docstring)
        # EyeResponse = loadmat(path_eye_response)['CIE1988Photopic']
        if self.photopic_eye_response is None:
            EyeResponse = np.loadtxt(self.calibration_folder / 'CIE1988photopic.txt')
            self.photopic_eye_response = (683.002 * np.interp(self.wavelengths_cutted, EyeResponse[:,0], EyeResponse[:,1]))
        
        SpecRadInt = self.intensities * IRF * PIXEL_SIZE / self.abs_calfactor / (self.integration_time/1000)
        # Cut innecessary wavelengths assuming the range 450 - 800 nm to more than enough
        self.SpecRadInt = SpecRadInt[self.filter_wavelengths]
        
        self.SpecLumInt = self.SpecRadInt * self.photopic_eye_response
        
        self.radiance = np.trapz(self.SpecRadInt, self.wavelengths_cutted) / PIXEL_SIZE
        self.luminance = np.trapz(self.SpecLumInt, self.wavelengths_cutted) / PIXEL_SIZE

        return self.radiance , self.luminance
    
def plot_measurement(fig, ax, x, y, label = None, color = None):
    """
    Plots the spectra given the figure and axis handlers together with the x, and y data.
    CHANGED: instead of an interactive Tk window (plt.show/plt.pause, which crashes
    when run alongside a worker thread), this re-saves the figure to the PNG path
    stashed on the figure as `fig._gonio_savepath`. One numbered file per scan,
    refreshed after every angle.
    """
    if color is None:
        ax.plot(x,y, label = label)
    else:
        ax.plot(x,y, label = label, color = color)
    
    if label is not None:
        ax.legend(ncol = 2, fontsize = 'x-small')
    
    # CHANGED: write the snapshot to disk (no GUI). If no path was set, do nothing.
    savepath = getattr(fig, '_gonio_savepath', None)
    if savepath is not None:
        try:
            fig.savefig(savepath, dpi = 100)
        except Exception as e:
            # Never let a plotting/IO hiccup take down a long measurement.
            print(f'WARNING: could not save plot snapshot ({e})')
    
if __name__ == '__main__':
    try:
        # Define folder and filename
        g = GonioLogger('test', angle_step=10.8)
        g.take_gonio_measurement()
        
    except KeyboardInterrupt:
        g.shutdown()
        
    #    g.take_dark_spectra()