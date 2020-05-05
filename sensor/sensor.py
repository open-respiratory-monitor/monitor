#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 16:59:47 2020

sensor.py

This defines a sensor class which initializes the pressure sensors and
can poll the sensors.

@author: nlourie
"""

#from datetime import datetime
#import numpy as np
#from scipy import signal
import board
import busio
import adafruit_lps35hw
import time
import numpy as np
import sys
import os

# add the main directory to the PATH
main_path = os.path.dirname(os.getcwd())
sys.path.insert(1, main_path)


class sensor(object):
    
    """ 
    Constituents:
        sensor1 = pressure sensor #1
        sensor2 = pressure sensor #2
        
        p1 = pressure @ sensor 1 in cmH20
        p2 = pressure @ sensor 2 in cmH20
        dp = differential pressure (p2 - p1) in cmH20
    """
    
    def __init__(self,main_path, calfile = '/calibration/Flow_Calibration.txt',verbose = False):
        
        # run in verbose mode?
        self.verbose = verbose
        
        # Initialize the i2c bus
        self.i2c = busio.I2C(board.SCL, board.SDA)
        
        # Using the adafruit_lps35hw class to read in the pressure sensor
            # note the address must be in decimal.
            # allowed addresses are: 
                # 92 (0x5c - if you put jumper from SDO to Gnd)
                # 93 (0x5d - default)
        
        # Set up the sensors        
        self.sensor2 = adafruit_lps35hw.LPS35HW(self.i2c, address = 92)
        self.sensor1 = adafruit_lps35hw.LPS35HW(self.i2c, address = 93)
        
        self.sensor1.data_rate = adafruit_lps35hw.DataRate.RATE_75_HZ
        self.sensor2.data_rate = adafruit_lps35hw.DataRate.RATE_75_HZ
        self.sensor1.low_pass_enabled = True
        self.sensor2.low_pass_enabled = True
    

    
        # Define the unit conversion factor
        self.mbar2cmh20 = 1.01972

        # Load the flow calibration polynomial coefficients
        self.main_path = main_path
        self.calfile = calfile
        
        # flow calibration polynomial
        if self.verbose:
            print(f"trying to load calfile at {calfile}")
        self.flowcal = np.loadtxt(self.main_path + self.calfile,delimiter = '\t',skiprows = 1)

        
        # Zero the sensors
        self.rezero()
        
        # Initialize the class values
        self.read()
    def dp2flow(self,dp_cmh20):
        flow_sign = np.sign(dp_cmh20)
        flow = flow_sign*np.polyval(self.flowcal,np.abs(dp_cmh20))
        return flow        
    
    def rezero(self):
        # Zeroes the sensors
        # Just takes the current readings and sets them to zero pressure for
        # each sensor
        
        # Now read out the pressure difference between the sensors
        print('p1_0 = ',self.sensor1.pressure,' mbar')
        print('p1_0 = ',self.sensor1.pressure*self.mbar2cmh20,' cmH20')
        print('p2_0 = ',self.sensor2.pressure,' mbar')
        print('p2_0 = ',self.sensor2.pressure*self.mbar2cmh20,' cmH20')
        
        print('')
        print('Now zero the pressure:')
        # Not sure why sometimes I have to do this twice??
        self.sensor1.zero_pressure()
        self.sensor1.zero_pressure()
        time.sleep(1)
        self.sensor2.zero_pressure()
        self.sensor2.zero_pressure()
        time.sleep(1)
        print('p1_0 = ',self.sensor1.pressure,' mbar')
        print('p1_0 = ',self.sensor1.pressure*self.mbar2cmh20,' cmH20')
        print('p2_0 = ',self.sensor2.pressure,' mbar')
        print('p2_0 = ',self.sensor2.pressure*self.mbar2cmh20,' cmH20')
        print()

    def read(self):
        # Read the pressure sensors and update the values
        self.p1 = self.sensor1.pressure * self.mbar2cmh20
        self.p2 = self.sensor2.pressure * self.mbar2cmh20
        self.dp = self.p2 - self.p1
        
        # Calculate the flow
        self.flow = self.dp2flow(self.dp)
        
if __name__ == "__main__":
    
    sensor = sensor(main_path = main_path,verbose = True)       
        
        
        
        
        
