# pyGoniospectrometer
[![MIT license](http://img.shields.io/badge/license-MIT-yellowgreen.svg)](http://opensource.org/licenses/MIT)

## Description
This project contains the code to control a home-made goniospectrometer instrument. A goniospectrometer can be used to measure the angular spectral emission of, for example, organic light-emitting diodes (OLEDs) or organic light-emitting electrochemical cells (LECs).

This instrument is relatively easy and relatively cheap to built. 

## What will you need?
, but you still will need the following hardware:

- A PC, laptop or a Raspberry (see the raspberry branch of the repository for more info).
- An Ocean Optics USB controlled spectrometer (e.g. I used a Flame spectrometer).
- An optical fiber suitable for your purposes.
- A hobbyst stepper motor with 1.8°/step (e.g. I used a NEMA 17 Bipolar 1.8deg 0.4A 12V.
- A (Big) Easy Driver from Sparkfun.
- An Arduino (just in case you don't choose the Raspberry Pi option!).
- Acces to a 3D printer.
- Some wires for the connections.

The following hardware can be useful but it is not required:
- Some space in and optical table.
- Components to mount the optical fiber in the optical table (posts, mounting base, kinematic mount, adapters, etc...)
- A collimating lens (e.g. from Thorlabs).
- Acces to a 3D printer, e.g. to print the device holder.

## The code
The code relies on the fantastic [python-seabreeze](https://python-seabreeze.readthedocs.io/en/latest/) library.
There is also a simple graphical user interface coded using the ```dash``` framework.

There are two current versions separated (at the moment) via two brances: ```master``` and ```raspberry```:
1. The instrument in the ```master``` branch assumes that the control for the stepper motor is done via an Arduino and that this one is controlled and synchronized with the spectrometer through a computer. The communication is done via the serial port of the Arduino (the USB) and uses the [pyvisa](https://pyvisa.readthedocs.io/en/latest/) library to manage it.
2. **(Recommended)** The instrument in the ```raspberry``` branch is simplified and you just need a Raspberry Pi (I used the Pi400), since the control of the stepper can be done with directly with the Raspberry's in-built GPIO pins.


## Installation & documentation
At the moment, I am working on making the code and the documentation for the instrument setting-up of the instrument more accessible. 

