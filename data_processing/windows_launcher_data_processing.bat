@echo off
ECHO Activating the conda environment...
call conda.bat activate
ECHO Launching the program...
ECHO Launching the browser...
"C:\Program Files\Mozilla Firefox\firefox" http://127.0.0.1:8055/ -ftimeout=10 
python gui_processing.py -p 8055
call conda.bat deactivate