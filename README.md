# pyGoniospectrometer
[![MIT license](http://img.shields.io/badge/license-MIT-yellowgreen.svg)](http://opensource.org/licenses/MIT)

## Description
This project contains the code to control a home-made goniospectrometer instrument. A goniospectrometer can be used to measure the angular spectral emission of, for example, organic light-emitting diodes (OLEDs) or organic light-emitting electrochemical cells (LECs). The data you can obtain can be later used together with optical models to extract rellevant information from your devices (EQE, emission zone position, diple orientation, etc).

This instrument is relatively easy and cheap to built. 

## What will you need?

This is the minimum hardware required:

- A PC, a laptop or a Raspberry (e.g. I used a [Pi400](https://www.raspberrypi.org/products/raspberry-pi-400/) see the [raspberry branch](https://github.com/jrafolsr/pyGonioSpectrometer/tree/raspberry) of the repository for more info).
- An Ocean Optics USB controlled spectrometer (e.g. I used a Flame spectrometer).
- An optical fiber suitable for your purposes.
- A hobbyst stepper motor with 1.8°/step to be used as the goniometer (e.g. I used a NEMA 17 Bipolar 1.8deg 0.4A 12V).
- A second stepper motor to be used as the shutter (*Please note that it will be removed in future versions*).
- A (Big) Easy Driver from Sparkfun.
- An Arduino (just in case you don't choose the Raspberry Pi option!).
- Some wires for the connections.

The following hardware can be useful but it is not required:
- Some space in and optical table.
- Components to mount the optical fiber in the optical table (posts, mounting base, kinematic mount, adapters, etc...)
- A collimating lens (e.g. from Thorlabs).
- Acces to a 3D printer, e.g. to print the device holder and a coupler to the motor shaft.

## The code
The code relies on the great [python-seabreeze](https://python-seabreeze.readthedocs.io/en/latest/) library to control the spectrometer through the ```SpectraMeasurement```class.

The stepper motor is controlled through either the ```ArduinoMotorController``` or the ```RaspberryMotorController``` class. Those objects take care of the goniometer angular movement.

The ```gonio_measurement.py``` or  ```gonio_measurement_time.py``` are scripts to perform a full goniospectrometer measurement.

The script ```GonioSpec_GUI.py``` launches a simple graphical user interface using the [dash](https://dash.plotly.com/) framework to perform a similar job as ```gonio_measurement.py```.

There are currently two versions of the instrument separated (at the moment) via two brances: ```master``` and ```raspberry```:

1. The instrument in the ```master``` branch assumes that the control for the stepper motor is done via an Arduino and that this one is controlled and synchronized with the spectrometer through a computer. The communication is done via the serial port of the Arduino (the USB) and uses the [pyvisa](https://pyvisa.readthedocs.io/en/latest/) library to manage it.

2. **(Recommended)** The instrument in the ```raspberry``` branch is simplified and you just need a Raspberry Pi (I used the Pi400), since the control of the stepper can be done with directly with the Raspberry's in-built GPIO pins and you can skip the PC-Arduino tandem.


## Installation & documentation
At the moment, I am working on making the code and the documentation for the setting-up of the instrument more accessible.
There is a [manual](manual/manual.docx) of how to operate the instrument once it has been built.

