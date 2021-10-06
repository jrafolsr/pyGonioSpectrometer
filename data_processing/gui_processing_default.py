# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 11:50:58 2021

@author: JOANRR
"""
from pyGonioSpectrometer.data_processing import __file__ as calibration_dir

from pathlib import Path

# Calibration files
calibration_dir = Path(calibration_dir).parent / 'calibration_files'

calibrations_dict ={'arduino_gonio_1': dict(path_IRF = calibration_dir / 'IRF_FlameBlueFiber_wLens2945K.txt',\
                                       abs_calfactor = 7.10e6),\
               'raspberry_gonio2':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference.txt',\
                                       abs_calfactor = 6.08E6),\
               'raspberry_gonio2 + polarizer':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference+LPVISE100a.txt', abs_calfactor = 6.08E6),\
               'raspberry_gonio2_old':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens2945K.txt',\
                                       abs_calfactor = 1.44E6)}
    
default_calibration = None


error_landscapes_dict = {'ErrorLandscape_newPL(ref ITO/SY/Al)' : str(calibration_dir /'ErrorLandscape_newPL.mat'),\
                         'EL_model_ITO/PEDOT(40nm)/SY/Al' :  str(calibration_dir / 'EL_model_ITO-PEDOT-SY-Al.mat'),\
                         'ErrorLandscape_newPL_p-pol' : str(calibration_dir /'ErrorLandscape_newPL_p-pol.mat'),\
                         'ErrorLandscape_newPL_s-pol' : str(calibration_dir /'ErrorLandscape_newPL_s-pol.mat'),\
                         'ErrorLanscape_newPL_281nm_parallel': str(calibration_dir /'ErrorLanscape_newPL_281nm_parallel.mat'),\
                         'ErrorLanscape_newPL_281nm_parallel_p-pol' : str(calibration_dir /'ErrorLanscape_newPL_281nm_parallel_ppol.mat'),\
                         'ErrorLanscape_newPL_281nm_parallel_s-pol' : str(calibration_dir /'ErrorLanscape_newPL_281nm_parallel_spol.mat')   }

default_error_landscape = 'ErrorLandscape_newPL(ref ITO/SY/Al)'