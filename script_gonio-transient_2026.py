#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time-resolved gonio + luminance measurement with periodic full angular scans.

Same protocol as time-gonio_measurement_script_snapshots_updated.py, but the
LEC source is now driven by the SweepMyLEC from `pyInstruments.ivl_snapshots`
(the same class used by script_IVL_snapshots_with_T.py), and the steady-state
rate condition is delegated to `SweepMyLEC.check_rate_condition` instead of a
hand-rolled linregress on the voltage trace.

REQUIRES the optional-photodiode / optional-temperature patch to SweepMyLEC
(resource_pd=None and run_logger_thread(..., save_temperature=False)). Without
it the logger thread will try to read a photodiode and a temp.dat file that do
not exist on the gonio rig. See the message accompanying this file.

Deviations from the previous (_updated) script are tagged `# CHANGED:`.

@author: JOANRR (cleanup pass)
"""

# =============================================================================
# Imports
# =============================================================================
# CHANGED: use the shared pyInstruments SweepMyLEC instead of class_sweeper_current.
from pyInstruments.ivl_snapshots import SweepMyLEC
from pyGonioSpectrometer.tasks_updated import GonioLogger
from time import sleep, monotonic
from datetime import datetime
from pathlib import Path
import numpy as np
from timer import Timer
import shutil
# CHANGED: linregress is no longer imported here -- the rate fit now lives in
# SweepMyLEC.check_rate_condition.


# =============================================================================
# Experiment configuration
# =============================================================================

# --- Output folder / file id ------------------------------------------------
folder = Path('/home/pi/Documents/data/joan/2026/20260603_polariton-postICEL-series')
file_id = 'Ag14D2'

# --- Spectrometer / gonio acquisition ---------------------------------------
integration_time = 200      # Initial guess for the integration time (ms)
n_spectra = 2               # Spectra to average for a quick initial luminance
max_n_spectra = 3           # Cap on averaged spectra (bounds acquisition time)
angle_step = 5.4            # Angular step (deg), must be a multiple of 1.8 deg
angle_max = 86.4            # Max angle (deg), must be a multiple of angle_step
max_intensity_angle = 0.0   # Initial guess of the angle of max intensity

max_time_per_fwd_luminance = 3000   # ms, integration budget for the luminance series
max_time_per_angle = 3000           # ms, integration budget during a full scan

min_time_before_update_integration_time = 30    # s, grace before retuning
factor = 0.5                # Integration-time adjustment factor
gonio_transient_time_step = 120     # s, base cadence for full gonio scans

plot_snapshots = False      # True -> one Agg PNG per scan in <goniologger>/temp/

# --- Electrical staircase (constant current) --------------------------------
term = 'REAR'               # Keithley sensing terminal
list_currents = list(np.round(10 ** (np.arange(-1.25, 0.51, 0.25)), 3))  # mA
N_cycles = 3                # Number of staircase cycles

staircase_downscan = True   # Append the reversed staircase (up then down)
max_time_per_current_step = 10.1 * 60   # s, hard cap on the hold per current
min_time_per_current_step = 10 * 60     # s, min hold before accepting SS

# CHANGED: rate_condition is now passed straight to check_rate_condition. In CC
# mode that method fits dV/dt and reports it in mV/min, so the units/meaning are
# identical to before (|dV/dt| <= rate_condition, mV/min).
rate_condition = 1          # mV/min
rate_window_s = 20          # s, regression window for check_rate_condition

# --- I-V sweep (CURRENTLY HARD DISABLED, kept for future use) ----------------
Vstart = 0
Vend = 4.5
step = 0.25
sweep_downscan = True
sweep_nplc = 0.01

# --- Misc timing ------------------------------------------------------------
time_after_sweep = 5        # s, dwell at the same step after the final scan
waiting_calc = lambda current: max_time_per_current_step   # per-step time budget

# --- Instrument resources ---------------------------------------------------
resource_keithley = 'GPIB0::24::INSTR'      # LEC source
# CHANGED: no photodiode on the gonio rig -> None (needs the optional-PD patch).
# If a photodiode SMU IS present, set its GPIB address here instead.
resource_pd = None

# --- Derived configuration (do not touch) -----------------------------------
folder = folder / file_id

initial_current_mA = list_currents[0]
list_currents = [el / 1000 for el in list_currents]   # mA -> A
if staircase_downscan:
    list_currents = list_currents + list_currents[-2::-1]

length_staircase = len(list_currents) / 2 if staircase_downscan else len(list_currents)

sweep_voltage = np.concatenate([
    np.arange(Vstart, 1.75, 0.25),
    np.arange(1.75, 2.75, 0.125),
    np.arange(2.75, Vend + step, step),
])
if sweep_downscan:
    sweep_voltage = np.concatenate([sweep_voltage, sweep_voltage[-2::-1]])


# =============================================================================
# Helpers
# =============================================================================

def forward_radiance(gonio):
    """Rough forward 'radiance' proxy (a.u.) from the last luminance spectrum.
    Returns NaN if no spectrum has been taken yet."""
    if np.all(np.isnan(gonio.intensities)):
        return np.nan
    slc = slice(450, -450)
    return np.trapz(gonio.intensities[slc], gonio.wavelengths[slc]) \
        / gonio.integration_time * 1000


def recent_voltage_mean(m, window_s):
    """Mean logged voltage over the last `window_s` seconds (NaN if too few
    points). Used only for the shorted-device guard now that the dV/dt fit is
    done inside check_rate_condition."""
    ff = m.time_arr >= m.time_arr[-1] - window_s
    if np.count_nonzero(ff & ~np.isnan(m.voltage_arr)) > 0:
        return np.nanmean(m.voltage_arr[ff])
    return np.nan


def take_full_scan(gonio, m, current, tag):
    """Take one full angular scan, handling the integration-time bookkeeping.
    Retunes at the angle of max intensity, switches to max_time_per_angle for
    the scan, then restores max_time_per_fwd_luminance."""
    gonio.update_integration_time(angle=gonio.max_intensity_angle)
    print(f'\nINFO: Moved to angle {gonio.max_intensity_angle:.2f} deg for the '
          f'integration-time adjustment.')

    gonio.max_time_per_angle = max_time_per_angle
    print('\nINFO: Taking a full goniometer scan!\n')
    header = (f'# Current in A:\t {current:.6e}\n'
              f'# Ellapsed time in s:\t {m.main_timer.ellapsed_time():.2f}\n')
    suffix = (f'I={current * 1000:04.2f}mA_'
              f't={m.main_timer.ellapsed_time():05.0f}s_{tag}')
    gonio.take_gonio_measurement(plot=plot_snapshots, suffix=suffix,
                                 header=header,
                                 parameter_1=m.main_timer.ellapsed_time())
    gonio.max_time_per_angle = max_time_per_fwd_luminance


def log_interval(m, gonio_timer, dVdt, condition_flag):
    """Append one row to the intervals log for the current step just finished."""
    with open(logfile, 'a') as f:
        f.write('{:.2f}\t{:.2f}\t{:.4f}\t{:.6e}\t{:.2f}\t{}\n'.format(
            m.main_timer.ellapsed_time(),
            gonio_timer.ellapsed_time(),
            m.voltage_arr[-1],
            m.current_arr[-1] * 1000,
            dVdt,
            condition_flag))


# =============================================================================
# Folder / log setup
# =============================================================================

if not folder.exists():
    folder.mkdir()

logfile = folder / (file_id + '_intervals.log')
with open(logfile, 'w') as f:
    f.write('Total_time[s]\tEllapsed_time[s]\tVoltage[V]\tCurrent[mA]\t'
            'dVdt[mV/min]\tCondition_Status\n')


# =============================================================================
# Bring up instruments
# =============================================================================

# --- Gonio + spectrometer ---------------------------------------------------
goniospectrometer = GonioLogger(
    file_id, folder,
    angle_step=angle_step, angle_max=angle_max,
    max_time_per_angle=max_time_per_fwd_luminance,
    integration_time=integration_time, n_spectra=n_spectra,
    max_n_spectra=max_n_spectra,
    suffix_luminance_file='time-series',
)
goniospectrometer.max_intensity_angle = max_intensity_angle

shutter_is_closed = bool(int(input('Is the shutter closed?\n\t0 --> No\n\t1 -> Yes\n')))
if not shutter_is_closed:
    goniospectrometer.gonio.shutter_is_closed = False
    goniospectrometer.gonio.close_shutter()
goniospectrometer.take_dark_spectra()

# --- Timers (gonio cadence only) --------------------------------------------
gonio_timer = Timer(min_time_step=0.5, max_time_step=300, fix_step=False)
gonio_timer_fullscan = Timer(min_time_step=gonio_transient_time_step,
                             max_time_step=900,
                             intervals=(600, 1800, 3600, 7200), fix_step=False)

# --- Keithley (LEC source) --------------------------------------------------
# CHANGED: new SweepMyLEC constructor takes (resource, resource_pd, ...).
m = SweepMyLEC(resource_keithley, resource_pd, output_folder=folder)
kwargs_current = dict(nplc=1, aver=False, Ncount=1, fw=False, term=term)
m.configure_I(initial_current_mA / 1000, **kwargs_current)

timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss_")
logger_filename = timestamp + file_id + '_voltage-current.txt'

# CHANGED: save_temperature=False -- no PID writes temp.dat on this rig.
m.run_logger_thread(logger_filename, save_temperature=False)

# CHANGED: the old class exposed `keithley_is_on`; this one does not. Give the
# logger a moment to spin up and start filling the data arrays instead. (Until
# >10 points are in the window, check_rate_condition just returns (False, NaN).)
sleep(1.0)

# Hold times reached at steady state during the up-scan, reused on the down-scan.
half_staircase_times = dict()


# =============================================================================
# Main protocol
# =============================================================================

dVdt = np.nan   # persists across iterations for the post-step log/print

try:
    for j in range(N_cycles):
        print(f'\n----------------- Cycle {j:d}/{N_cycles:d} of the staircase -----------------')

        for i, current in enumerate(list_currents):
            staircase_direction = i <= length_staircase
            rate_condition_flag = False

            gonio_timer.initialize()
            gonio_timer_fullscan.initialize()

            # CHANGED: queue-based, thread-safe setpoint change. The logger picks
            # it up, resets m.sub_timer and clears the data arrays; the short
            # sleep lets that happen before we start polling the rate.
            m.request_bias_update(current)
            sleep(0.2)

            sleeping_time = (waiting_calc(current) if staircase_direction
                             else half_staircase_times.get(f'{current * 1000:.4f}',
                                                            waiting_calc(current)))
            print(f'\nGoing to sleep for {sleeping_time} sec')

            # ----- Hold this current step -----------------------------------
            while gonio_timer.ellapsed_time() < sleeping_time:

                # CHANGED: rate condition delegated to the source object. In CC
                # mode `dVdt` comes back in mV/min and `condition` is already
                # |dVdt| <= rate_condition. verbose=False so we keep our own
                # status print (which also shows the gonio counts).
                condition, dVdt = m.check_rate_condition(
                    rate_condition, check_last_n_seconds=rate_window_s, verbose=False)

                # Shorted-device guard (aborts the run; caught below).
                Vmean = recent_voltage_mean(m, rate_window_s)
                if np.isfinite(Vmean) and gonio_timer.ellapsed_time() > 10 and Vmean < 2.4:
                    raise ValueError(f'\nThe device is probably shorted, '
                                     f'stopping the measurement with V = {Vmean:.4f} V')

                radiance = forward_radiance(goniospectrometer)
                print(f'\r{gonio_timer.ellapsed_time(): 6.2f} s  '
                      f'{m.voltage_arr[-1]: 6.2f} V '
                      f'{m.current_arr[-1] * 1000: 8.4f} mA  '
                      f'{dVdt: 10.3f} mV/min, counts = {radiance:.0f} a.u.',
                      end='\r')

                # Steady-state acceptance (up-scan only, after the min hold).
                if (condition
                        and gonio_timer.ellapsed_time() > min_time_per_current_step
                        and staircase_direction):
                    print(f'\nINFO: Steady state condition reached with {dVdt:.2f} '
                          f'mV/min at {gonio_timer.ellapsed_time():.2f}'
                          f'({m.main_timer.ellapsed_time():.2f}) s')
                    half_staircase_times[f'{current * 1000:.4f}'] = gonio_timer.ellapsed_time()
                    rate_condition_flag = True
                    break

                # Periodic full angular scan.
                if gonio_timer_fullscan.istime2measure():
                    take_full_scan(goniospectrometer, m, current, tag='time-series')

                # Otherwise, periodic forward-luminance point + limit check.
                elif gonio_timer.istime2measure():
                    goniospectrometer.save_spectra(m.main_timer.ellapsed_time())
                    limit_hit = goniospectrometer.check_low_and_high_limits()
                    if gonio_timer.ellapsed_time() > min_time_before_update_integration_time:
                        if limit_hit == 1:
                            if goniospectrometer.integration_time > 1:
                                goniospectrometer.update_integration_time(factor=factor)
                        elif limit_hit == -1:
                            if goniospectrometer.integration_time < max_time_per_fwd_luminance:
                                goniospectrometer.update_integration_time(factor=factor)
                        # else: within band, nothing to do.

                sleep(0.005)

            # ----- Step finished --------------------------------------------
            if not rate_condition_flag:
                print(f'\nINFO: Steady state ended by time budget {sleeping_time:.2f} s')

            log_interval(m, gonio_timer, dVdt, rate_condition_flag)

            # Final steady-state full scan for this current.
            take_full_scan(goniospectrometer, m, current, tag='ss_time-series')

            # ----- (Disabled) fast I-V sweep --------------------------------
            # NEEDS WORK FOR LATER USAGE
            # sweeper_filename = file_id + '_run1'
            # m.sweep_and_log(np.append(sweep_voltage, m.voltage_arr[-1]), current,
            #                 logger_filename, logger_aver=False, logger_Ncount=1,
            #                 logger_nplc=1, sweep_fileid=sweeper_filename,
            #                 sweep_ranging='AUTO', sweep_reset=False,
            #                 sweep_delay='AUTO', sweep_nplc=sweep_nplc)

            time0 = monotonic()
            print(f'\nGoing to sleep for {time_after_sweep:.2f} s')
            while (monotonic() - time0) < time_after_sweep:
                print(f'\r{(monotonic() - time0): 6.2f} s  '
                      f'{m.voltage_arr[-1]: 6.2f} V '
                      f'{m.current_arr[-1] * 1000: 8.4f} mA  '
                      f'{dVdt: 10.3f} mV/min', end='\r')
                sleep(0.01)

    print('\nCheckpoint 1: staircase loop complete')
    m.bias_setpoint = list_currents[0]   # CHANGED: was m.current_setpoint
    sleep(0.05)
    print('Checkpoint 2: stopping logger')
    m.stop_logger()

except KeyboardInterrupt:
    print('\nINFO: Terminating program (KeyboardInterrupt)')
except Exception as e:
    print(e)
    print('\nINFO: Terminating program (exception above)')


# =============================================================================
# Cleanup
# =============================================================================

sleep(0.05)
print('\nCheckpoint 3: parking instruments')

# CHANGED: guard the logger stop so a hiccup can't skip the hardware shutdown.
try:
    m.stop_logger()
except Exception as e:
    print(f'(non-fatal) m.stop_logger raised: {e}')

m.outpoff()
if resource_pd is not None:
    m.pd_outpoff()
m.device.close()           # CHANGED: was m.keithley.close(); the new attr is `device`.

goniospectrometer.shutdown()

timestamp = datetime.now().strftime("%Y-%m-%d")
shutil.copy(__file__, folder / (timestamp + '_measurement_script.py'))

print('INFO: Program terminated.')