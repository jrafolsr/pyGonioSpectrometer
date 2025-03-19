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
                    'arduino_gonio_1 (newLENS)': dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference_u-pol.txt',\
                                       abs_calfactor = 7.73e6),\
               'arduino_gonio_1 (newLENS p-pol)': dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference_p-pol.txt',\
                                       abs_calfactor = 7.73e6 * 0.4500),\
         
               'arduino_gonio_1 (newLENS s-pol)': dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference_s-pol.txt',\
                                       abs_calfactor = 7.73e6 * 0.3784),\
         
               'raspberry_gonio2':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference.txt',\
                                       abs_calfactor = 6.08E6),\
               'catagoniometer (raspberry_gonio3)' : dict(path_IRF = calibration_dir / 'IRF_USB2000+ShortBlueFiber_wLens-F260SMA-A_PMA12-reference.txt',\
                                                          abs_calfactor = 1),\
               'raspberry_gonio2 + polarizer':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference+LPVISE100a.txt', abs_calfactor = 6.08E6),\
               # 'raspberry_gonio4 + polarizer':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference+LPVISE100a.txt', abs_calfactor = 1),\
               'raspberry_gonio4+Flame1-F110+u-pol':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F110SMA-A+100um-slit_PMA12-reference_u-pol.txt',\
                                       abs_calfactor = 3.275E6),\
               'raspberry_gonio4+Flame1-F110+p-pol':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F110SMA-A+100um-slit_PMA12-reference_p-pol.txt',\
                                       abs_calfactor = 1.550E6),\
               'raspberry_gonio4+Flame1-F110+s-pol':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F110SMA-A+100um-slit_PMA12-reference_s-pol.txt',\
                                       abs_calfactor = 1.221E6),   
               'raspberry_gonio4+Flame1-F260+u-pol':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference_u-pol.txt',\
                                       abs_calfactor = 6.197E6),\
               'raspberry_gonio4+Flame1-F260+p-pol':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference_p-pol.txt',\
                                       abs_calfactor = 2.789E6),\
               'raspberry_gonio4+Flame1-F260+s-pol':dict(path_IRF = calibration_dir / 'IRF_Flame1OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference_s-pol.txt',\
                                       abs_calfactor = 2.345E6)}

default_calibration = None

extra_dir = Path(r'C:\Users\JOANRR\Documents\18_p-n-junction_silent-dipoles\simulations')

error_landscapes_dict = {'ErrorLandscape_newPL(ref ITO/SY/Al)' : str(calibration_dir /'ErrorLandscape_newPL.mat'),\
                         'EL_model_ITO/PEDOT(40nm)/SY/Al' :  str(calibration_dir / 'EL_model_ITO-PEDOT-SY-Al.mat'),\
                         'EL_model_ITO/SY/Ca(20nm)/Al' :  str(calibration_dir / 'EL_model_ITO-SY-Ca20-Al.mat'),\
                         'ErrorLandscape_newPL_p-pol' : str(calibration_dir /'ErrorLandscape_newPL_p-pol.mat'),\
                         'ErrorLandscape_newPL_s-pol' : str(calibration_dir /'ErrorLandscape_newPL_s-pol.mat'),\
                         'ErrorLanscape_newPL_281nm_parallel': str(calibration_dir /'ErrorLanscape_newPL_281nm_parallel.mat'),\
                         'ErrorLanscape_newPL_281nm_parallel_p-pol' : str(calibration_dir /'ErrorLanscape_newPL_281nm_parallel_ppol.mat'),\
                         'ErrorLanscape_newPL_281nm_parallel_s-pol' : str(calibration_dir /'ErrorLanscape_newPL_281nm_parallel_spol.mat'),\
                         'half-rod_ITO-SY-Al_100-400_isotropic_ErrorLandscape': str(extra_dir/'half-rod_ITO-SY-Al_100-400_isotropic_ErrorLandscape.mat'),\
                         'half-rod_ITO-SY-Al_100-400_isotropic_ErrorLandscape_p-pol' : str(extra_dir /'half-rod_ITO-SY-Al_100-400_isotropic_ErrorLandscape_ppol.mat'),\
                         'half-rod_ITO-SY-Al_100-400_isotropic_ErrorLandscape_s-pol' : str(extra_dir /'half-rod_ITO-SY-Al_100-400_isotropic_ErrorLandscape_spol.mat'),\
                         'MR-TADF_ITO-PEDOT-TADF-Al_s-pol_n=1.8' : str(r'C:\Users\JOANRR\Documents\16_Shi-projects\20230322_MR-TADF\MRTADF_150-180nm_isotropic_ErrorLandscape_spol.npz'),\
                         'MR-TADF_ITO-PEDOT-TADF-Al_p-pol_n=1.8' : str(r'C:\Users\JOANRR\Documents\16_Shi-projects\20230322_MR-TADF\MRTADF_150-180nm_isotropic_ErrorLandscape_ppol.npz'),\
                         'MR-TADF_ITO-PEDOT-TADF-Al_s-pol_n=1.6' : str(r'C:\Users\JOANRR\Documents\16_Shi-projects\20230322_MR-TADF\MRTADF_160-175nm_isotropic_n=1.6_ErrorLandscape_spol.npz'),}
    
for f in calibration_dir.glob('*.npz'):
    error_landscapes_dict[f.stem] = str(calibration_dir / f)
    

default_error_landscape = 'ErrorLandscape_newPL(ref ITO/SY/Al)'