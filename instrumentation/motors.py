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
        self.speed = 110 # deg/s

        sleep(2.0)
        
    def close(self):
        """
        Closes the resource
        """
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
        
        return lsteps, lresolution, round(out_angle,4)
    
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
            itime = time()
            self.motor.write(string2arduino)
            
#            unfinished = True
#            itime = time()
#            while unfinished:
#                try:
#                    unfinished = self.motor.bytes_in_buffer == 0
#                    # Temporary solution in case Arduino misses the command
#                    if (time() - itime) > 5:
#                        print('(WARNING: Failed to move the gonio, trying it again)')
#                        self.motor.write(string2arduino)
#                        itime = time()    
#                    sleep(0.01)
#                    
#                except KeyboardInterrupt:
#                    break
#                
##            self.motor.clear()
#            print('INFO: '  + self.motor.read())
    # Waiting time to ensure the movement has finished, based on the hardware speed of the rotation
        etime = time() - itime
        
        sleeping = 0.15 + angle  / self.speed
        if etime <  sleeping:
            sleep(sleeping - etime)

        
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
#        print(step)
        if delay is not None:
            sleep(delay)
        
        self.motor.write(f'0,{step:d}')
        
#        unfinished = True
#        itime = time()
#        while unfinished:
#            try:
#                unfinished = self.motor.bytes_in_buffer == 0
#                # Temporary solution to the shutter opening issue.
#                if (time() - itime) > 3:
#                    print('(WARNING: Failed to move shutter, trying it again)')
#                    self.motor.write(f'0,{step:d}')
#                    itime = time()
#                    
#                sleep(0.01)
#            except KeyboardInterrupt:
#                break
#        self.motor.clear()
            
#        print('INFO: '  + self.motor.read())
        sleep(0.75)

            
    def disable_shutter(self):
        """
        Does nothing.
        """
        pass