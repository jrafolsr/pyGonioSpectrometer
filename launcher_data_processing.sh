#!/bin/bash
echo "Launching program"
python3 /home/pi/Documents/PythonScripts/pyGonioSpectrometer/data_processing/gui_processing.py -p 8055 &
echo "Launching chromium"
sleep 3s
chromium-browser --URL http://127.0.0.1:8055 --remote-debugging-port=9222 &
$SHELL

