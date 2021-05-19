# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 15:59:20 2020

@author: JOANRR

This module contains a class to control the two stepper motors needed for the goniospectrometer setup: 
    1. The goniometer
    2. The shutter (currently done with a stepper but it eventually be simplified with a servo/actuator)

It is a Python wrapper that communicates with the Arduino that controls the stepper and simplifies the control. From the wrapper, we tell the arduino the number of steps to be made, the resolution and the direction.

Explanation to be improved.

"""
from pyvisa import ResourceManager
from time import sleep, time
import RPi.GPIO as gpio
from math import exp

# Declaration of constants
FACTOR = 1000000 # We'll be working with udeg, so the interger operations work fine
MOTOR_ASTEP = int(1.8 * FACTOR) # Step of the motor in udeg, specified by the manufacturer. E.g. 1.8 deg/step
TIMEOUT = 5000 # Max. number of seconds to wait before perform the next movement of the motor, just as a safety feature
RESOLUTIONS = [1, 2, 4, 8, 16]

rm = ResourceManager()

def list_ports():
    """
    Lists the available resources in the computer. Can be use to determine which corresponds to the Arduino.
    
    Returns:
    --------
    lports: list
        List with the available ports in the computer.
    """
    
    lports = rm.list_resources()
    
    if lports == []:
        print('WARNING: No COM ports found. Check the connections.')
    # else:
        # print('The available COM ports are:')
        # for port in lports:
        #     print(port)
    return lports
        


class ArduinoMotorController():
    """
    Creates an object that controls and simplifies the communication and movement of the motors via the Arduino.
    
    Parameters
    ----------
        port : str
            String containing the name of the Arduino port.
            
    Attributes
    ----------
        motor : obj
            A ResourceManager object.
        port : str
            Str with containing the communication port.
            
    Methods
    -------
        close()
            Clears and close the communication port.
        angle2steps(angle, resolution = 16, slow = True)
            Translates the desired angle to a number of steps according to the specified resolution
        move_angle(angle, direction = None, resolution = 16)
            Tells the gonio motor to move the specified angle.
        enable_gonio()
            Enables the gonio motor.
        disable_gonio()
            Disables the gonio motor.
        move_shutter(angle, direction=None, resolution=16)
            Moves the shutter. The Arduino controls and tracks in which direction should be moved.
        
    """
    
    def __init__(self, port):
        self.port = port
        if port is None:
            raise Exception('ERROR: COM port not provided')
        
        self.motor = rm.open_resource(port)
        
        self.motor.read_termination = '\r\n'
        self.motor.write_termination = '\n' # important for the Arduino to understand the commands.
        # Clear the buffer
        self.motor.clear()

        sleep(0.500)
        
    def close(self):
        """
        Closes the resource
        """
        self.motor.clear() 
        self.motor.close()
        
    def angle2steps(self, angle, resolution = 16, slow = True):
        """
        Translates the angle into steps according to the resolution. 
        
        If slow = False, the fastes step/resolution
        will be calculated and it will deliver a series of commands for the arduino to perform. Nice but
        unnecessary feature. To be simplified and replaced. I originally overkilled the functionalities for the 
        requirements needed.
        
        Parameters
        ----------
            angle : float
                Angle that will be converted into steps.
            resolution : int (only 1, 2, 4, 8, 16), optional
                The resolution of the microsteps. The default is 16 (1/16th of the step).
            slow (deprecated): bool, optional
            
        Returns
        -------
            lsteps : list
                List with the number of steps per command. (To be updated to a int value once slow is removed.)
            lresolution: list
                List with the resolutions per command. (To be updated to a int value once slow is removed.)
            out_angle: float
                Calculated angle that the motor will move given by the snumber_of_steps *  microstep.
        """
    
    #     print(f'Input angle: {angle:.2f}°')
        angle = int (angle * FACTOR) # Angle to udeg and to integer, for a proper modules operation
        
        if slow:
            if resolution not in RESOLUTIONS:
                raise Exception('The resolution value is not accepted. Must be one of this: 1, 2, 4, 8 or 16.')
            A_RESOLUTION = MOTOR_ASTEP / FACTOR / resolution
    #         print(f'Using a fixed step of {A_RESOLUTION:.4f}')
            step = angle // int(A_RESOLUTION * FACTOR)
            
            out_angle = step * A_RESOLUTION     
    #         print(f'Output angle: {out_angle:.2f}°')
            lsteps = [step]
            lresolution = [resolution]
            
        """    # Not using this option so far, nice feature but not needed. Maybe in the future.
        else:
    #         print('Using the fastest configuration to move.')
            #Lists with the data to be output
            lsteps = [angle // MOTOR_ASTEP] # List of steps with first step, as many as 1.8 deg steps as it can fit
            langles = [MOTOR_ASTEP / FACTOR] # List angles
            lresolution = [1] # List with the resolutions
            for i in range(1,5):
                angle = angle %  (MOTOR_ASTEP // lresolution[i-1])
    
                lsteps.append(angle // (MOTOR_ASTEP // 2**i))
    
                langles.append(langles[-1] / 2)
    
                lresolution.append(2**i)
                  
            out_angle = np.sum(np.array(langles)*np.array(lsteps))
    #         print(f'Output angle: {out_angle:.2f}°')"""
        
        return lsteps, lresolution, out_angle
    
    def move_angle(self, angle, direction = None, resolution = 16, slow = True):
        """
        Tells the gobiometer to move a certain angle in the specified direction and resolution.
        
        Parameters
        ----------
        angle : float
            Angle to move. CLoackwise if positive, counterclockwise if negative.
        direction : int (only 1 or 2), optional
            Move clockwise (1) or counterclockwise (2). If not specified the direction is inferred from the sign
            of angle (prefered choice). The default is None.
        resolution: int (only 1, 2, 4, 8, 16), optional
            The resolution of the microsteps. The default is 16 (1/16th of the step). The default is 16.
        slow (deprecated) : bool, optional
        
        Returns
        -------
            out_angle : float
                Calculated angle that the motor will move.
        """
        
        if direction is None:
            if angle > 0: direction = 1
            elif angle < 0: direction = 2
            else: direction = 0
            
        angle = abs(angle) # Makes sure the angle is always positive
        msteps, mresolution, out_angle = self.angle2steps(angle, resolution = resolution, slow = slow)
                  
        for step, res in zip(msteps, mresolution):
            
            string2arduino = f'1,{step:d},{direction:d},{res:d}' #EXPLAIN HERE
#            print(string2arduino)
            self.motor.write(string2arduino)
            
            unfinished = True
            itime = time()
            while unfinished:
                try:
                    unfinished = self.motor.bytes_in_buffer == 0
                    # Temporary solution in case Arduino misses the command
                    if (time() - itime) > 5:
                        print('(WARNING: Failed to move the gonio, trying it again)')
                        self.motor.write(string2arduino)
                        itime = time()    
                    sleep(0.01)
                    
                except KeyboardInterrupt:
                    break
                
            self.motor.clear()
#            print('INFO: '  + self.motor.read())
        
        sleep(0.25)
        
        return out_angle
    
    def disable_gonio(self):
        """
        Disables the gonio motor.
        """
        self.motor.write('1,-1,0,0')
        self.motor.read()
        print('INFO: Gonio motor disabled')
        return None
    
    def enable_gonio(self):
        """
        Enables the gonio motor.
        """
        self.motor.write('1,-1,1,0')
        self.motor.read()
        print('INFO: Gonio motor enabled')
        return None
    
    def move_shutter(self, angle = 270, resolution = 1, delay = None, slow = True):
        """
        Opens or closes the shutter. the Arduino controls and tracks the direction.
        
        Parameters
        ----------
        angle : float, optional
            Angle to move the shutter. The default is 270 deg.
        resolution : int (only 1, 2, 4, 8, 16), optional
            The resolution of the microsteps. The default is 1 (full step).
        delay : float, optional
            Add a delay (in seconds) before sending the command to the Arduino (temporal fix). The default is None.
        slow (deprecated) : bool, optional
        
        """
        
        steps, _, _ = self.angle2steps(angle, resolution = resolution, slow = slow)
        step = steps[0]
        
        if delay is not None:
            sleep(delay)
        
        self.motor.write(f'0,{step:d}')
        
        unfinished = True
        itime = time()
        while unfinished:
            try:
                unfinished = self.motor.bytes_in_buffer == 0
                # Temporary solution to the shutter opening issue.
                if (time() - itime) > 3:
                    print('(WARNING: Failed to move shutter, trying it again)')
                    self.motor.write(f'0,{step:d}')
                    itime = time()
                    
                sleep(0.01)
            except KeyboardInterrupt:
                break
        self.motor.clear()
            
#        print('INFO: '  + self.motor.read())
        sleep(0.5)

            
    def disable_shutter(self):
        """
        Does nothing.
        """
        pass
class RaspberryMotorController():
    """
    Creates an motor controller object assuming a Raspberry system.
    
    Parameters
    ----------

            
    Attributes
    ----------
        pinout : 4-length tuple
        stepPIN, directionPIN, enablePIN, shutterPIN : int
            GPIO pin numbers to control gonio and shutter.
        delay : float
            Delay between stepPIn LOW/HIGH signal. Basically, the speed of the rotation.
            
    Methods
    -------
        close()
            Clears and close the communication port.
        angle2steps(angle, resolution = 16, slow = True)
            Translates the desired angle to a number of steps according to the specified resolution
        move_angle(angle, direction = None, resolution = 16)
            Tells thgpio_functiongpio_functione gonio motor to move the specified angle.
        enable_gonio()
            Enables the gonio motor.
        disable_gonio()
            Disables the gonio motor.
        move_shutter(angle, direction=None, resolution=16)
            Moves the shutter. The Arduino controls and tracks in which direction should be moved.
        
    """
    
    def __init__(self, pinout = (27, 17, 26), delay = 0.5, shutter_angle = 180):
        """
        It assumes the shutter is closed from the beginning.
        Parameters
        ----------
        pinout : 4-element tuple or list, optional
            Tuple containing the pins as int for the motor driver in the order: stepPIN, directionPIN, enablePIN.
            The default is (17, 27, 26).
        delay : float, optional
            Delay between the triggers send to the motor drivers in ms. The default is 0.1 ms.
        shutter_angle : int
            Angle corresponding to the open and closed position of the shutter, it must be an int.The default is 180.
        """
        
        gpio.setmode(gpio.BCM)
        
        [gpio.setup(i, gpio.OUT) for i in pinout]
        
        self.pinout = pinout
        self.stepPIN, self.dirPIN, self.enPIN = pinout
        self.angle_error = 0.00
        
        self.delay = delay
        self.shutter_angle = int(180)
        self.shutter_steps = (self.shutter_angle * FACTOR ) // int(MOTOR_ASTEP // 16) # Using the default resolution of 1/16h
        self.shutter_counter = 0
        # Initialize all the pinouts to LOW
        [gpio.output(i, gpio.LOW) for i in self.pinout]
        
        # Initialize the pwm to control the shutter servo
        self.isclosed = False
        self.shutter_is_closed = True
        self.direction = 0 # Controls the change of direction to correct the error drift
        self.steps_counter = 0
        
        # Gonio rotation speed (hardwarewise)
        self.speed = 144 # deg/s
        
        sleep(1)
        
    def start(self):
        """
        Starts the gonio in case it is closed.
        """
        if self.isclosed:
            self.__init__()

    def close(self):
        """
        Swithces off the motor driver and frees teh GPIOs
        """
        self.isclosed = True
        gpio.cleanup()
        
       
        
    def angle2steps(self, angle, resolution = 16):
        """
        Translates the angle into steps according to the resolution. 
        
        If slow = False, the fastes step/resolution
        will be calculated and it will deliver a series of commands for the arduino to perform. Nice but
        unnecessary feature. To be simplified and replaced. I originally overkilled the functionalities for the 
        requirements needed.
        
        Parameters
        ----------
            angle : float
                Angle that will be converted into steps.
            resolution : int (only 1, 2, 4, 8, 16), optional
                The resolution of the microsteps. The default is 16 (1/16th of the step).
            
        Returns
        -------
            steps : int
                Number of steps
            out_angle: float
                Calculated angle that the motor will move given by the snumber_of_steps *  microstep.
        """
    
    #     print(f'Input angle: {angle:.2f}°')
        angle = int (angle * FACTOR) # Angle to udeg and to integer, for a proper modules operation

        if resolution not in RESOLUTIONS:
            raise Exception('The resolution value is not accepted. Must be one of this: 1, 2, 4, 8 or 16.')
        
        A_RESOLUTION = MOTOR_ASTEP / FACTOR / resolution
#         print(f'Using a fixed step of {A_RESOLUTION:.4f}')
        steps = angle // int(A_RESOLUTION * FACTOR)
        
        out_angle = steps * A_RESOLUTION     
#         print(f'Output angle: {out_angle:.2f}°')
              
        return steps, round(out_angle,4)
    
    def move_angle(self, angle, direction = None, resolution = 16, correct_drift = True):
        """
        Tells the goniometer to move a certain angle in the specified direction and resolution.
        
        Parameters
        ----------
        angle : float
            Angle to move. CLoackwise if positive, counterclockwise if negative.
        direction : int (only 0 or 1), optional
            Move clockwise (0) or counterclockwise (1). If not specified the direction is inferred from the sign
            of angle (prefered choice). The default is None.
        resolution: int (only 1, 2, 4, 8, 16), optional
            The resolution of the microsteps. The default is 16 (1/16th of the step). The default is 16.
        
        Returns
        -------
            out_angle : float
                Calculated angle that the motor will move.
        """
        
        if direction is None:
            if angle >= 0: direction = gpio.LOW
            else: direction = gpio.HIGH
            
        gpio.output(self.dirPIN, direction) # Pull the direction pin LOW/HIGH depending on the rotation 
        
        sleep(0.1)
        
        if correct_drift:
            if self.direction != direction:
                # A change in the direction has occurred, so the error needs to be reversed
                self.angle_error *= -1
                self.direction = direction

        # Makes sure the angle is always positive for the step calculation and adds the drift error
        
        angle = round(abs(round(angle, 1)) + self.angle_error, 4)
        
        self.current_angle = angle
        
        steps, out_angle = self.angle2steps(angle, resolution = resolution)      
        
        if correct_drift: self.angle_error = round(angle - out_angle, 4)
        
        
        self.move_steps(steps)
        

        return out_angle
    
    def move_steps(self, steps):
        
        itime = time()
        # Replacing the fix delay with the accelerator generator.
        delays = self.__accelerator__(steps)
        for delay in delays:
#            for i in range(steps):
            gpio.output(self.stepPIN, gpio.HIGH)
#            sleep(delay)
            gpio.output(self.stepPIN, gpio.LOW)
            sleep(delay)
        
        # Waiting time to ensure the movement has finished, based on the hardware speed of the rotation
        etime = time() - itime
        
        if etime < self.current_angle  / self.speed:
            sleep(self.current_angle  /self.speed - etime)
            
        self.steps_counter += (-1)**(self.direction +2) *steps
        
#        print(f'INFO: Total number of steps perfomed: {self.steps_counter:d}')
            
    def __accelerator__(self, steps, tau = 50, max_delay = 2):
        a = 1/tau
        delays = [(max_delay * (exp(-i*a) + exp((i-steps + 1)*a)) +self.delay) * 0.001 for i in range(steps)]
#        delays = [0.001 for i in range(steps)]
        return  delays
    
    def disable_gonio(self):
        """
        Disables the gonio motor.
        """
        gpio.output(self.enPIN, gpio.HIGH)
        print('INFO: Gonio motor disabled')
        return None
    
    def enable_gonio(self):
        """
        Enables the gonio motor.
        """
        gpio.output(self.enPIN, gpio.LOW)
        print('INFO: Gonio motor enabled')
        return None
    

    def open_shutter(self):
        """
        Opens the shutter
        """
        self.direction = (self.shutter_counter + 1) % 2
        gpio.output(self.dirPIN, self.direction) # Pull the direction pin LOW/HIGH depending on the rotation 
        sleep(0.1)
        
        self.shutter_is_closed = False
        self.current_angle = self.shutter_angle
        
        self.move_steps(self.shutter_steps)
        sleep(2)
        
        print('INFO: Shutter opened.')
        

        
    def close_shutter(self):
        """
        Closes the shutter
        Parameters:
        -----------
        open_angle : float
            Angle at which the shutter is open.
        """
        self.direction = self.shutter_counter % 2
        
        gpio.output(self.dirPIN,self.direction ) # Pull the direction pin LOW/HIGH depending on the rotation

        sleep(0.1)
        self.shutter_is_closed = True
        self.shutter_counter += 1
        self.current_angle = self.shutter_angle
        self.move_steps(self.shutter_steps)
        sleep(2)
        
        print('INFO: Shutter closed.')

    
    def move_shutter (self):
        """
        Checks wheter the shutter is open or closed and moves it based on that.
        """
        if self.shutter_is_closed:
            self.open_shutter()
        else:
            self.close_shutter()
            
            
#class RaspberryMotorController():
#    """
#    Creates an motor controller object assuming a Raspberry system.
#    
#    Parameters
#    ----------
#
#            
#    Attributes
#    ----------
#        pinout : 4-length tuple
#        stepPIN, directionPIN, enablePIN, shutterPIN : int
#            GPIO pin numbers to control gonio and shutter.
#        delay : float
#            Delay between stepPIn LOW/HIGH signal. Basically, the speed of the rotation.
#            
#    Methods
#    -------
#        close()
#            Clears and close the communication port.
#        angle2steps(angle, resolution = 16, slow = True)
#            Translates the desired angle to a number of steps according to the specified resolution
#        move_angle(angle, direction = None, resolution = 16)
#            Tells thgpio_functiongpio_functione gonio motor to move the specified angle.
#        enable_gonio()
#            Enables the gonio motor.
#        disable_gonio()
#            Disables the gonio motor.
#        move_shutter(angle, direction=None, resolution=16)
#            Moves the shutter. The Arduino controls and tracks in which direction should be moved.
#        
#    """
#    
#    def __init__(self, pinout = (4, 17, 27, 12), delay = 0.1, shutter_angles = (0,26)):
#        """
#        Parameters
#        ----------
#        pinout : 4-element tuple or list, optional
#            Tuple containing the pins as int for the motor driver in the order: stepPIN, directionPIN, enablePIN, shutterPIN.
#            The default is (4, 17, 27, 13).
#        delay : float, optional
#            Delay between the triggers send to the motor drivers in ms. The default is 0.1 ms.
#        shutter_angles : 2-length tuple
#            Angles corresponding to the open and closed position of the shutter.
#        """
#        
#        gpio.setmode(gpio.BCM)
#        
#        [gpio.setup(i, gpio.OUT) for i in pinout]
#        
#        self.pinout = pinout
#        self.stepPIN, self.dirPIN, self.enPIN, self.shutterPIN = pinout
#        
#        self.delay = delay / 1000
#        
#        # Initialize all the pinouts to LOW
#        [gpio.output(i, gpio.LOW) for i in self.pinout]
#        
#        # Initialize the pwm to control the shutter servo
#        self.pwm = gpio.PWM(self.shutterPIN, 50)
#        
#        self.pwm.start(0)
#        self.isclosed = False
#        self.open_angle, self.closed_angle = shutter_angles
#        sleep(0.25)
#        self.close_shutter()
#        
#        sleep(1)
#        
#    def start(self):
#        """
#        Starts the gonio in case it is closed.
#        """
#        if self.isclosed:
#            self.__init__()
#
#    def close(self):
#        """
#        Swithces off the motor driver and frees teh GPIOs
#        """
#        self.close_shutter()
#        self.isclosed = True
#        self.disable_gonio()
#
#        sleep(0.25)
#        self.pwm.stop()
#        gpio.cleanup()
#        
#       
#        
#    def angle2steps(self, angle, resolution = 16):
#        """
#        Translates the angle into steps according to the resolution. 
#        
#        If slow = False, the fastes step/resolution
#        will be calculated and it will deliver a series of commands for the arduino to perform. Nice but
#        unnecessary feature. To be simplified and replaced. I originally overkilled the functionalities for the 
#        requirements needed.
#        
#        Parameters
#        ----------
#            angle : float
#                Angle that will be converted into steps.
#            resolution : int (only 1, 2, 4, 8, 16), optional
#                The resolution of the microsteps. The default is 16 (1/16th of the step).
#            
#        Returns
#        -------
#            lsteps : list
#                List wgpio_functionith the number of steps per commthonScripts/pyGonioSpectrometer/instrumentation')
#
#In [14]: and. (To be updated to a int value once slow optional arg is removed.)
#            lresolution: list
#                List with the resolutions per command. (To be updated to a int value once slow optional arg is removed.)
#            out_angle: float
#                Calculated angle that the motor will move given by the snumber_of_steps *  microstep.
#        """
#    
#    #     print(f'Input angle: {angle:.2f}°')
#        angle = int (angle * FACTOR) # Angle to udeg and to integer, for a proper modules operation
#        
#
#        if resolution not in RESOLUTIONS:
#            raise Exception('The resolution value is not accepted. Must be one of this: 1, 2, 4, 8 or 16.')
#        
#        A_RESOLUTION = MOTOR_ASTEP / FACTOR / resolution
##         print(f'Using a fixed step of {A_RESOLUTION:.4f}')
#        step = angle // int(A_RESOLUTION * FACTOR)
#        
#        out_angle = step * A_RESOLUTION     
##         print(f'Output angle: {out_angle:.2f}°')
#        lsteps = [step]
#        lresolution = [resolution]
#            
#        
#        return lsteps, lresolution, out_angle
#    
#    def move_angle(self, angle, direction = None, resolution = 16):
#        """
#        Tells the goniometer to move a certain angle in the specified direction and resolution.
#        
#        Parameters
#        ----------
#        angle : float
#            Angle to move. CLoackwise if positive, counterclockwise if negative.
#        direction : int (only 0 or 1), optional
#            Move clockwise (0) or counterclockwise (1). If not specified the direction is inferred from the sign
#            of angle (prefered choice). The default is None.
#        resolution: int (only 1, 2, 4, 8, 16), optional
#            The resolution of the microsteps. The default is 16 (1/16th of the step). The default is 16.
#        
#        Returns
#        -------
#            out_angle : float
#                Calculated angle that the motor will move.
#        """
#        
#        if direction is None:
#            if angle >= 0: direction = gpio.LOW
#            else: direction = gpio.HIGH
#            
#        gpio.output(self.dirPIN, direction) # Pull the direction pin LOW/HIGH depending on the rotation 
#        sleep(0.1)
#        
#        angle = abs(angle) # Makes sure the angle is always positive for the tep calculation
#        msteps, mresolution, out_angle = self.angle2steps(angle, resolution = resolution)
#         # The list style is a remanent of old code or in case I want to implement a faster movement in teh future.
#         
#        for steps, res in zip(msteps, mresolution):
#            for i in range(steps):
#                gpio.output(self.stepPIN, gpio.HIGH)
#                sleep(self.delay)
#                gpio.output(self.stepPIN, gpio.LOW)
#                sleep(self.delay)
#        
#        return out_angle
#    
#    def disable_gonio(self):
#        """
#        Disables the gonio motor.
#        """
#        gpio.output(self.enPIN, gpio.HIGH)
#        print('INFO: Gonio motor disabled')
#        return None
#    
#    def enable_gonio(self):
#        """
#        Enables the gonio motor.
#        """
#        gpio.output(self.enPIN, gpio.LOW)
#        print('INFO: Gonio motor enabled')
#        return None
#    
#    def move_angle_shutter(self, angle = 45, angle_x_duty = 20, offset = 3.7):
#        """
#        Move the shutter servo up to certain angle.
#        
#        Parameters
#        ----------
#        angle : float, optional
#            Angle to move the shutter. The default is 270 deg.
#        angle_x_duty : float, optional
#            Angle per unit of pwm. The default is 21 deg.
#        offset : float, optional
#            Duty cycle at 0 deg. The default is 3 deg.   
#        """
#        
#        duty = angle / angle_x_duty + offset
##        print(duty)
#        if duty > 10:
#            duty = 10
#            print('INFO: Max. angle reached.')
#        elif duty < 3:
#            duty = 3
#            print('INFO: Min. angle reached.')
#            
#        # Move servo to the desired positon
#        self.pwm.ChangeDutyCycle(duty)
#        
#        sleep(0.5) # Wait sufficient time to reach it
#        # Deactivate servo 
#        self.pwm.ChangeDutyCycle(0)
#        
#    def open_shutter(self):
#        """
#        Opens the shutter
#        """
#        self.move_angle_shutter(self.open_angle)
#        self.shutter_is_closed = False
#        
#        print('INFO: Shutter opened.')
#        
#    def close_shutter(self):
#        """
#        Closes the shutter
#        Parameters:
#        -----------
#        open_angle : float
#            Angle at which the shutter is open.
#        """
#        self.move_angle_shutter(self.closed_angle)
#        self.shutter_is_closed = True
#        print('INFO: Shutter closed.')
#   
#    def move_shutter (self):
#        """
#        Checks wheter the shutter is open or closed and moves it based on that.
#        """
#        if self.shutter_is_closed:
#            self.open_shutter()
#        else:
#            self.close_shutter()