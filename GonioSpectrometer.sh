#!/bin/bash
echo "Launching program"
python3 /home/pi/Documents/PythonScripts/pyGonioSpectrometer/gui.py 2>/dev/null&
echo "Launching chromium"
sleep 3s
chromium-browser --URL http://127.0.0.1:8052 --remote-debugging-port=9222 &
$SHELL

