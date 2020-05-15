#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 16:43:51 2020

Data Handler

This file defines two loops which are started in the mainwindow

Fast Loop: update data
    - P1
    - P2
    - Flow
    - Volume
    - Pi GPIO State

Slow Loop: update calculations
    - Volume calibration (spline fit)
    - Breath parameters


@author: nlourie
"""

from PyQt5 import uic, QtCore, QtGui, QtWidgets
#from PyQt5.QtWidgets import QMessageBox

import numpy as np


#import time
#import traceback
#import signal
#import board
#import busio
#import adafruit_lps35hw
import os
import sys
from datetime import datetime
from scipy import interpolate
from scipy import signal

try:
    import RPi.GPIO as GPIO
except Exception as e:
    print("Could not load Raspberry Pi GPIO module: ",e)

# add the wsp directory to the PATH
main_path = os.path.dirname(os.getcwd())
sys.path.insert(1, main_path)

# import custom modules
from sensor import sensor
from utils import utils

"""
# Define the fast and slow data objects that will be passed to the GUI thread
"""
class fast_data(object):
    def __init__(self):

        # pressures
        self.p1 = np.array([])
        self.p2 = np.array([])
        self.dp = np.array([])

        # flow
        self.flow = np.array([])

        # volume
        self.vol_raw = np.array([])
        self.vol_drift = np.array([]) # the drift volume which is the spline line through the detrended volume
        self.vol = np.array([])

        # time
        self.t_obj = np.array([]) # datetime object
        self.dt = np.array([])    # dt since first sample in vector
        self.t = np.array([])     # ctime in seconds
        self.fs = []



class slow_data(object):
    def __init__(self):

        ## THINGS THAT HOLD SINGLE VALUES ##
        # respiratory parameters from last breath
        self.pip = None
        self.peep = None
        self.pp = None
        self.vt = None
        self.mve_inf = None
        self.mve_meas = None
        self.rr = None
        self.ie = None
        self.c = None
        self.t_last = datetime.utcnow().timestamp()
        self.dt_last = 0.0

        # pi GPIO state
        self.lowbatt = None
        self.charging = None

    def print_data(self):
        print('Respiratory Parameters:')
        print(f'Time Since Last Detected Breath (s) = {self.dt_last}')
        print(f'PIP = {self.pip}')
        print(f'PEEP = {self.peep}')
        print(f'PP = {self.pp}')
        print(f'VT = {self.vt}')
        print(f'MVE Inferred = {self.mve_inf}')
        print(f'MVE Measured = {self.mve_meas}')
        print(f'RR = {self.rr}')
        print(f'I:E = {self.ie}')
        print(f'C = {self.c}')
        print()
        print('Device Parameters:')
        print(f'Plugged In = {self.charging}')
        print(f'Low Battery = {self.lowbatt}')






"""
# Define the fast and slow loops
"""
class fast_loop(QtCore.QThread):

    """
    This class gets and updates the data that gets plotted in realtime

    Needs to take a few things to get started:

        sensor -- an instance of the sensor class

    """
    # define a new signal that will be used to send updated data back to the main thread
    # this signal returns an object that holds the data to ship out to main
    newdata = QtCore.pyqtSignal(object)

    def __init__(self, main_path, update_time = 1000, time_to_display = 10.0,simulation = False,logdata = False,verbose = False):

        QtCore.QThread.__init__(self)

        # save a continuous log of the pressure data?
        self.logdata = logdata

        # run in simulation mode?
        self.simulation = simulation

        # run in verbose mode?
        self.verbose = verbose

        # get the main path
        self.main_path = main_path
        if self.verbose:
            print(f"fastloop: main path = {self.main_path}")

        # Define the instance of the object that will hold all the data
        self.fastdata = fast_data()
        self.slowdata = slow_data()

        # current time loop is executed
        self.t_obj = datetime.utcnow()
        self.t = self.t_obj.timestamp()

        # time to display is the approx time to show on the screen in seconds
        self.time_to_display = time_to_display #s

        # time between samples
        self.ts = update_time #ms

        # real sample rate
        self.ts_real = []

        # sample frequency - starts out as 1/self.ts but then is updated to the real fs
        self.fastdata.fs = 1.0/self.ts

        # length of vectors
        self.num_samples_to_hold = int(self.time_to_display*1000/self.ts )
        if self.verbose:
            print(f"fastloop: num samples to hold = {self.num_samples_to_hold}")
        # this just holds a number which increments every time the loop runs
        # TODO get rid of this
        self.index = 0



        # Set up the sensor
        if self.simulation:
            self.sensor = sensor.fakesensor(main_path = self.main_path, verbose = self.verbose)
        else:
            self.sensor = sensor.sensor(main_path = self.main_path,verbose = self.verbose)

        # set up file to store sensor data
        if logdata:
            print('creating file to store cal data')
            filename = str(int(datetime.utcnow().timestamp()))
            self.sensor_datafile = open(filename + "_sensor_raw.txt","w")
            self.sensor_datafile.write('time \t p1 \t p2 \t dp')



    def add_new_point(self,arr,new_point,maxlen):
        # adds a new data point to the array,
        # and keeps gets rid of the oldest point

        if len(arr) < maxlen:
            arr = np.append(arr,new_point)
        else:
            arr[:-1] = arr[1:]
            arr[-1] = new_point

        return arr

    def __del__(self):
        self.wait()

    def update(self):
        self.index +=1
        self.t_obj = datetime.utcnow()
        self.t = self.t_obj.timestamp()

        if self.verbose:
            # debugging: print an update of what we're doing
            print("\nfastloop: Index =  %d" % self.index)

        # record the update time
        self.update_time = datetime.utcnow()
        self.fastdata.t_obj = self.add_new_point(self.fastdata.t_obj, self.update_time, self.num_samples_to_hold)
        self.fastdata.t =     self.add_new_point(self.fastdata.t, self.update_time.timestamp(), self.num_samples_to_hold)
        self.fastdata.dt = self.fastdata.t - self.fastdata.t[0]

        # if there's at least two elements in the vector, calculate the real delta between samples
         # if there's at least two elements in the vector, calculate the real average delta between samples
        if len(self.fastdata.dt) >= 2:
            self.ts_real = np.abs(np.mean(self.fastdata.dt[1:] - self.fastdata.dt[:-1]))
            self.fastdata.fs = 1.0/self.ts_real



        # read the sensor pressure and flow data
        # there's probably a clenaer way to do this, but oh well...
        self.sensor.read()
        self.fastdata.p1   = self.add_new_point(self.fastdata.p1,   self.sensor.p1,   self.num_samples_to_hold)
        self.fastdata.p2   = self.add_new_point(self.fastdata.p2,   self.sensor.p2,   self.num_samples_to_hold)
        self.fastdata.dp   = self.add_new_point(self.fastdata.dp,   self.sensor.dp,   self.num_samples_to_hold)



        dp_zero = np.mean(self.fastdata.dp[np.abs(self.fastdata.dp)<0.0])
        if np.isnan(dp_zero):
            dp_zero = 0.0

        flow_zero = self.sensor.dp2flow(dp_zero)

        self.fastdata.dp[np.abs(self.fastdata.dp)<0.0] = 0.0

        self.fastdata.flow = self.sensor.dp2flow(self.fastdata.dp) - flow_zero
        
        # apply a median filter
        self.fastdata.flow = signal.medfilt(self.fastdata.flow,11)
        
        # log the data if we're in logdata mode
        if self.logdata:
            self.log_raw_sensor_data()

        # calculate the raw volume
        self.fastdata.vol_raw = signal.detrend(np.cumsum(self.fastdata.flow)/(self.fastdata.fs*60.0))

        try:
            # correct the detrended volume signal using the slowdata spline fit
            self.apply_vol_corr()

        except Exception as e:
            print("fastloop: error in volume spline correction: ",e)
            print("fastloop: could not apply vol spline correction. using raw volume instead...")
            self.fastdata.vol = self.fastdata.vol_raw
            self.fastdata.vol_drift = 0.0*self.fastdata.vol_raw


        # tell the newdata signal to emit every time we update the data
        self.newdata.emit(self.fastdata)

    def log_raw_sensor_data(self):

        self.sensor_datafile.write('%f \t %f \t %f \t %f\n' %(self.fastdata.t[-1],self.fastdata.p1[-1],self.fastdata.p2[-1],self.fastdata.dp[-1]))


    def apply_vol_corr(self):
        # this uses the current volume minima spline calculation to correct the volume by pinning all the minima to zero
        if len(self.fastdata.vol_raw) >= 10:
            i_min = utils.breath_detect_coarse(-1.0*self.fastdata.vol_raw,self.fastdata.fs,minpeak = 0.05)
        else:
            i_min = []

        if self.verbose:
            print(f"fastloop: found {len(i_min)} volume minima at dt = {self.fastdata.dt[i_min]}")
        if len(i_min) >= 2:
            self.fastdata.vol_corr_spline = interpolate.interp1d(self.fastdata.t[i_min],self.fastdata.vol_raw[i_min],kind = 'linear',fill_value = 'extrapolate')
            self.fastdata.vol_drift = self.fastdata.vol_corr_spline(self.fastdata.t)
        else:
            self.fastdata.vol_drift = np.zeros(len(self.fastdata.vol_raw))
        # apply the correction
        self.fastdata.vol = self.fastdata.vol_raw - self.fastdata.vol_drift


    def run(self):
        if self.verbose:
            print("fast loop: starting fast Loop")
        self.timer = QtCore.QTimer()
        self.timer.setInterval(self.ts)
        self.timer.timeout.connect(self.update)
        self.timer.start()
        self.exec() # YOU NEED THIS TO START UP THE THREAD!

        # NOTE: only QThreads have this exec() function, NOT QRunnables
        # If you don't do the exec(), then it won't start up the event loop
        # QThreads have event loops, not QRunnables

        """
        # NOTE: only QThreads have this exec() function, NOT QRunnables
        #       If you don't do the exec(), then it won't start up the event
        #       loop QThreads have event loops, not QRunnables

        # Source: https://doc.qt.io/qtforpython/overviews/timers.html
        Quote:
          In multithreaded applications, you can use the timer mechanism in any
          thread that has an event loop. To start an event loop from a non-GUI
          thread, use exec()

        """


class slow_loop(QtCore.QThread):

    # define a new signal that will be used to send updated data back to the main thread
    # this signal returns an object that holds the data to ship out to main
    newdata = QtCore.pyqtSignal(object)

    # this signal sends a request to the mainloop to get the current data from the fastloop
    request_fastdata = QtCore.pyqtSignal()

    def __init__(self, main_path, update_time = 5000, verbose = False):
        QtCore.QThread.__init__(self)

        # print stuff for debugging?
        self.verbose = verbose

        # get the main path
        self.main_path = main_path
        if self.verbose:
            print(f"slowloop: main path = {self.main_path}")

        # run in verbose mode?
        self.verbose = verbose

        # this just holds a number which increments every time the loop runs
        # TODO get rid of this
        self.index = 0

        # set up a place to store the data that comes in from the fast loop each time this loop runs
        self.fastdata = fast_data()

        # set up a place to store the slow data that is calculated each time this loop runs
        self.slowdata = slow_data()

        # note the time the loop is executed
        self.t_obj = datetime.utcnow()
        self.t = self.t_obj.timestamp()

        # time between samples
        self.ts = update_time #ms

        # loop sample frequency - starts out as 1/self.ts but then is updated to the real fs
        self.fs = 1.0/self.ts

        # times from the last breath
        self.tsi = self.t # start time of inspiration (absolute time)
        self.dtsi = [] # start time of inspiration (dt since start)
        self.dtei = [] # end time of inspiration (dt since start)
        self.dtee = [] # end time of expiration (dt since start)

        # set up the raspberry pi GPIO THESE ARE BCM -- OTHERWISE THERE'S A CONFLICT IF YOU USE BOARD
        # WHY?  I think that the pressure sensor module must be setting it somewhere.
        # BCM is better anyways, it's the numbers of the inputs not the pins!
        self.gpio_map = {'charging':5, 'lowbatt':4}

        try:
            if True:
                print(f"GPIO mode = {GPIO.getmode()}")
            GPIO.setmode(GPIO.BCM)

            for key in self.gpio_map.keys():
                if True:
                    print(f"slowloop: setting GPIO key: {key}")
                GPIO.setup(self.gpio_map[key],GPIO.IN)
                #time.sleep(1)
        except Exception as e:
            print('slowloop: unable to set up Pi GPIO: ',e)

    def __del__(self):
        self.wait()

    def update(self):
        """
        ### This is the slow loop ####

        Here's what we want to do:

            1. Calculate the correction to the volume signal
                # calibrations
                self.slowdata.vol_corr_spline = np.array([])
                self.slowdata.vol_drift_params = np.array([])
            2. fit the mimima and maxima of the volume signal using peak finder
            3. calculate the breath parameters of the last breath:
        """
        self.index +=1
        if self.verbose:
            print("\nslowloop: %d" % self.index)

        # emit the request data signal to get the current fastloop data vectors
        self.request_fastdata.emit()

        # note the time the loop is executed
        self.t_obj = datetime.utcnow()
        self.t = self.t_obj.timestamp()

        # try to fit a spline
        try:
            # find the volume minima
            self.find_vol_min()

            # calculate breath parameters
            #self.calculate_breath_params()
        except Exception as e:
            print("slowloop: vol min calc error: ",e)

        # try to calculate the breath parameters
        try:
            self.calculate_breath_params()

        except Exception as e:
            print("slowloop: error calculating breath parameters: ",e)

        
        # Do these whether or not a breath is detected!
        # how long ago was the last breath started?
        self.slowdata.dt_last = np.round((self.t - self.slowdata.t_last),2)
        
        # now measure the realtime value (it fluctuates but stays near the real value, and is valuable if no breaths are delivered)
            # we don't display a full minute so need to scale answer
        scale = 60.0/self.fastdata.dt[-1]
        flow_in = self.fastdata.flow[self.fastdata.flow > 0]
        flow_out = self.fastdata.flow[self.fastdata.flow < 0]
        mve_meas_in = np.round(np.abs(np.trapz(flow_in)/(self.fastdata.fs*60.0)*scale),1)
        mve_meas_out = np.round(np.abs(np.trapz(flow_out)/(self.fastdata.fs*60.0)*scale),1)
        #average the flow in and out to get a sensible result regardless of flow sensor drift
        self.slowdata.mve_meas = np.round(np.mean([mve_meas_in, mve_meas_out]),1)


        self.check_power()


        # tell main that there's new slow data: emit the newdata signal
        self.newdata.emit(self.slowdata)



    def check_power(self):
        # Recall: self.gpio_map = {'charging':7, 'lowbatt':29}
        # try to check the power status
        try:
            self.slowdata.lowbatt  = GPIO.input(self.gpio_map['lowbatt'])
            self.slowdata.charging = GPIO.input(self.gpio_map['charging'])

        except Exception as e:
            print("slowloop: could not check power state: ",e)

    def find_vol_min(self):
        """
        ## find the min of the volume signal using peak finder ##

        """
        # step 1: find index of min and max
        self.i_min_vol = utils.breath_detect_coarse(-1*self.fastdata.vol, fs = self.fastdata.fs,minpeak = 0.0)

        if self.verbose:
            print(f"slowloop: found {len(self.slowdata.index_of_min)} peaks at dt = {self.fastdata.dt[self.i_min_vol]}")






    def calculate_breath_params(self):
        """

        ## THINGS THAT HOLD SINGLE VALUES ##
        # respiratory parameters from last breath
        self.pip = []
        self.peep = []
        self.pp = []
        self.vt = []
        self.mve = []
        self.rr = []
        self.ie = []
        self.c = []
        """
        if len(self.i_min_vol) >= 2:
            # define the last breath
            self.tsi = self.fastdata.t[self.i_min_vol[-2]]
            self.tee = self.fastdata.t[self.i_min_vol[-1]]
            self.dtsi = self.fastdata.dt[self.i_min_vol[-2]]
            self.dtee = self.fastdata.dt[self.i_min_vol[-1]]
            
            # update the time of the last breath
            self.slowdata.t_last = self.tee
            
            # index range of the last breath
            index_range = np.arange(self.i_min_vol[-2],self.i_min_vol[-1]+1)


            # get tidal volume (in mL) and the end inspiration time (both defined at vol peak over last breath)
            self.slowdata.vt = np.max(self.fastdata.vol[index_range])*1000
            self.i_max_vol_last = index_range[np.argmax(self.fastdata.vol[index_range])]
            self.dtei = self.fastdata.dt[self.i_max_vol_last]

            # get pip: peak pressure over last breath
            self.slowdata.pip = np.max(self.fastdata.p1[index_range])

            # get respiratory rate (to one decimal place)
            self.slowdata.rr = np.round(60.0/(self.dtee - self.dtsi),1)

            # get i:e ratio (to one decimal place)
            dt_exp = self.dtee - self.dtei
            dt_insp = self.dtei - self.dtsi
            self.slowdata.ie = np.round(np.abs(dt_exp / dt_insp),1)

            # get minute volume
            # first infer it from the last breath: RR * VT
            self.slowdata.mve_inf = np.round((self.slowdata.rr * self.slowdata.vt/1000.0),1)
            

            # get peep: average pressure over the 50 ms about the end of expiration
            dt_peep = 0.05
            self.slowdata.peep = np.mean(self.fastdata.p1[(self.fastdata.dt>=self.dtee-(dt_peep/2)) & (self.fastdata.dt<=self.dtee+(dt_peep/2))])

            # get pp: plateau pressure -- defined as mean pressure over 50 ms before the end of inspiration
            dt_pip = 0.05
            self.slowdata.pp = np.mean(self.fastdata.p1[(self.fastdata.dt>=self.dtei-dt_pip) & (self.fastdata.dt<=self.dtei)])

            # get static lung compliance. Cstat = VT/(PP - PEEP), reference: https://www.mdcalc.com/static-lung-compliance-cstat-calculation#evidence
            self.slowdata.c = np.round(self.slowdata.vt/(self.slowdata.pp - self.slowdata.peep),1)

            # round the pip and peep and vt
            self.slowdata.pip = np.round(self.slowdata.pip,1)
            self.slowdata.peep = np.round(self.slowdata.peep,1)
            self.slowdata.pp = np.round(self.slowdata.pp,1)
            self.slowdata.vt = np.round(self.slowdata.vt,1)

            

        else:
            print("slowloop: no breath detected!")
            
        

    def update_fast_data(self,fastdata):
        # this is a slot connected to mainwindow.newrequest
        # this takes the fastdata from the main window and updates the internal value of fastdata
        self.fastdata = fastdata

        if self.verbose:
            print(f"slowloop: received new fastdata from main (updated at {self.fastdata.t[-1]}, fs = {self.fastdata.fs} Hz")




    def run(self):
        if self.verbose:
            print("slowloop: starting slowloop")
        self.timer = QtCore.QTimer()
        self.timer.setInterval(self.ts)
        self.timer.timeout.connect(self.update)
        self.timer.start()
        self.exec() # YOU NEED THIS TO START UP THE THREAD!
        # NOTE: only QThreads have this exec() function, NOT QRunnables
        # If you don't do the exec(), then it won't start up the event loop
        # QThreads have event loops, not QRunnables
        """
        # NOTE: only QThreads have this exec() function, NOT QRunnables
        #       If you don't do the exec(), then it won't start up the event
        #       loop QThreads have event loops, not QRunnables

        # Source: https://doc.qt.io/qtforpython/overviews/timers.html
        Quote:
          In multithreaded applications, you can use the timer mechanism in any
          thread that has an event loop. To start an event loop from a non-GUI
          thread, use exec()

        """