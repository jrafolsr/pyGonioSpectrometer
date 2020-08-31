# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 15:59:20 2020

@author: JOANRR
Contains a class to control the two stepper motors movement needed for the spectro-goniometer setup: the goniometer itself and the shutter.

It is a Python wrapper that communicates with the Arduino that controls the stepper. From the wrapper, we tell the arduino the number of steps to be made, the velocity and the directions.

Explanation to be improved.
"""
import pyvisa
from time import sleep


# Some definitions
FACTOR = 1000000 # We'll be working with udeg, si the interger operations work fine
MOTOR_ASTEP = int(1.8*FACTOR) #udeg
TIMEOUT = 5000 # Max. number of seconds to wait until perform  the next movement of the motor, just as a safety feature
RESOLUTIONS = [1, 2, 4, 8, 16]

rm = pyvisa.ResourceManager()

def list_ports():
    
    lports = rm.list_resources()
    if lports == []:
        print('WARNING: No COM ports found. Check the connections.')
    # else:
        # print('The available COM ports are:')
        # for port in lports:
        #     print(port)
    return lports
        


class ArduinoMotorController():
    """AN object that controls the communication and movement of the motors via the Arduino"""
    def __init__(self, port):
        self.port = port
        if port is None:
            raise Exception('ERROR: COM port not provided')
            
        self.motor = rm.open_resource(port)
        # Clear the buffer
        self.motor.clear()
        self.motor.read_termination = '\r\n'
        self.motor.write_termination = '\n'
        sleep(0.500)

        
    def close(self):
        """Close the resource"""
        self.motor.clear() 
        self.motor.close()
        
    def angle2steps(self, angle, slow = True, resolution = 16):
        """Calculates the angle to move using from the bigger to the smaller steps. Not suitable for small and accurate movements perhaps"""
    
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
    
    def move_angle(self, angle, direction = None, slow = True, resolution = 16):
        """Tells the gonio-motor to move certain angle"""
        if direction is None:
            if angle > 0: direction = 1
            elif angle < 0: direction = 2
            else: direction = 0
            
        angle = abs(angle) # Makes sure the angle is always positive
        msteps, mresolution, out_angle = self.angle2steps(angle, slow, resolution)
                  
        for step, res in zip(msteps, mresolution):
            
            string2arduino = f'1,{step:d},{direction:d},{res:d}' #EXPLAIN HERE
            self.motor.write(string2arduino)
            
            unfinished = True
            while unfinished:
                try:
                    unfinished = self.motor.bytes_in_buffer == 0
                    sleep(0.01)
                except KeyboardInterrupt:
                    break
#            print('INFO: '  + self.motor.read())
        
        return out_angle
    
    def disable_gonio(self):
        self.move_angle(0)
        self.motor.read()
        return None
    
    def move_shutter(self, angle = 270, slow = True, resolution = 1, delay = None):
        steps, _, _ = self.angle2steps(angle, slow, resolution)
        step = steps[0]
        if delay is not None:
            sleep(delay)
        
        self.motor.write(f'0,{step:d}')
        
        unfinished = True
        while unfinished:
            try:
                unfinished = self.motor.bytes_in_buffer == 0
                sleep(0.01)
            except KeyboardInterrupt:
                break
            
#        print('INFO: '  + self.motor.read())
        sleep(0.5)
            
    def disable_shutter(self):
        pass