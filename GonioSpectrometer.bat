@echo off
ECHO Activating the conda environment...
call conda.bat activate
ECHO Launching the program...
ECHO Launching the browser...
"C:\Program Files\Mozilla Firefox\firefox" http://127.0.0.1:8051/ -ftimeout=10 
python GonioSpec_GUI.py
call conda.bat deactivate