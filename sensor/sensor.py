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

import time
import numpy as np
import sys
import os

#load the sensor board modules:
try:
    import board
    import busio
    import adafruit_lps35hw
except Exception as e:
    print("could not load sensor board modules: ",e)

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

    def __init__(self,main_path, calfile = '/calibration/Flow_Calibration.txt',dp_thresh = 0.0,verbose = False):

        # run in verbose mode?
        self.verbose = verbose

        self.dp_thresh = dp_thresh

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

        dp = self.p2 - self.p1
        if np.abs(dp) < self.dp_thresh:
            dp = 0.0
        self.dp = dp

        # Calculate the flow
        self.flow = self.dp2flow(self.dp)


class fakesensor(object):
    # if we're in simulation mode, then instead of using a real sensor, just
    # read in old data like it's real data.
    """
    Constituents:
        sensor1 = pressure sensor #1
        sensor2 = pressure sensor #2

        p1 = pressure @ sensor 1 in cmH20
        p2 = pressure @ sensor 2 in cmH20
        dp = differential pressure (p2 - p1) in cmH20
    """

    def __init__(self,main_path, calfile = '/calibration/Flow_Calibration.txt',datafile = '/calibration/Simulated_Data.txt',dp_thresh = 0.0,verbose = False):
        #datafile = '/calibration/1589499917_sensor_raw.txt'
        self.datafile = main_path + datafile
        self.verbose = verbose
        self.time_arr,self.p1_arr,self.p2_arr,self.dp_arr = np.loadtxt(self.datafile,delimiter = '\t',skiprows = 1,unpack = True)

        self.linenum = 0
        
        self.lastline = len(self.time_arr)-1
        
        # Define the unit conversion factor
        self.mbar2cmh20 = 1.01972

        # Load the flow calibration polynomial coefficients
        self.main_path = main_path
        self.calfile = calfile
        
        # flow calibration polynomial
        if self.verbose:
            print(f"trying to load calfile at {calfile}")
        self.flowcal = np.loadtxt(self.main_path + self.calfile,delimiter = '\t',skiprows = 1)
        
        #print statement
        print(f"Reading Simulated Data from File: {self.datafile}")

    def dp2flow(self,dp_cmh20):
        flow_sign = np.sign(dp_cmh20)
        flow = flow_sign*np.polyval(self.flowcal,np.abs(dp_cmh20))
        return flow
        
    #def rezero(self):
        
    
    def read(self):
        
        # read the fake data from the current line
        self.p1 = self.p1_arr[self.linenum]
        self.p2 = self.p2_arr[self.linenum]
        self.dp = self.p2 - self.p1
        
        # increment the line number
        self.linenum += 1
        
        # start again if you hit the end of the file
        if self.linenum >= self.lastline:
            self.linenum = 0
        
        
if __name__ == "__main__":
    try:
        print("Loading real sensor...")
        sensor = sensor(main_path = main_path,verbose = True)
    except Exception as e:
        print("Could not load real sensor: ",e)

    try:
        print("Loading simulated sensor...")
        sensor = fakesensor(main_path = main_path,verbose = True)
    except Exception as e:
        print("Could not load fake sensor data: ",e)



