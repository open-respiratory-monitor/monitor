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
        self.flow_raw = np.array([])
        self.flow = np.array([])

        # volume
        self.vol_raw = np.array([])
        self.vol_trend = np.array([])
        self.vol_detrend = np.array([])
        self.v_drift = np.array([]) # the drift volume which is the spline line through the detrended volume
        self.vol = np.array([])

        ## THINGS THAT HOLD ARRAYS ##
        #indices of volume min and max
        self.index_of_min = np.array([])
        self.index_of_max = np.array([])

        # times of the volume min and max
        self.vmin_times = np.array([])
        self.vmax_times = np.array([])

        # volume values at the min and max
        self.vmin_detrend = np.array([]) # used for fitting the spline
        self.vmin = np.array([]) # min after applying the spline. should all be zero! just a diagnostic
        self.vmax  = np.array([]) # peaks after applying the spline. this is used for calculating breath params.

        # calibrations
        self.vol_corr_spline = None
        self.vol_drift_params = None

        # time
        self.t_obj = np.array([]) # datetime object
        self.dt = np.array([])   # dt since first sample in vector
        self.t = np.array([])   # ctime in seconds

        # pi GPIO state
        self.lowbatt = False
        self.charging = True

class slow_data(object):
    def __init__(self):

        ## THINGS THAT HOLD ARRAYS ##
        #indices of volume min and max
        self.index_of_min = np.array([])
        self.index_of_max = np.array([])

        # times of the volume min and max
        self.vmin_times = np.array([])
        self.vmax_times = np.array([])

        # volume values at the min and max
        self.vmin_detrend = np.array([]) # used for fitting the spline
        self.vmin = np.array([]) # min after applying the spline. should all be zero! just a diagnostic
        self.vmax  = np.array([]) # peaks after applying the spline. this is used for calculating breath params.

        # calibrations
        self.vol_corr_spline = None
        self.vol_drift_params = None

        ## THINGS THAT HOLD SINGLE VALUES ##
        # times from the last breath
        self.tsi = [] # start time of inspiration
        self.tei = [] # end time of inspiration
        self.tse = [] # start time of expiration (same as tei)
        self.tee = [] # end time of expiration

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

    def __init__(self, main_path, update_time = 1000, time_to_display = 20.0,verbose = False):

        QtCore.QThread.__init__(self)

        # run in verbose mode?
        self.verbose = verbose

        # get the main path
        self.main_path = main_path
        if self.verbose:
            print(f"fastloop: main path = {self.main_path}")

        # time to display is the approx time to show on the screen in seconds
        self.time_to_display = time_to_display #s

        # time between samples
        self.ts = update_time #ms

        # real sample rate
        self.ts_real = []

        # sample frequency - starts out as 1/self.ts but then is updated to the real fs
        self.fs = 1.0/self.ts

        # length of vectors
        self.num_samples_to_hold = int(self.time_to_display*1000/self.ts )
        if self.verbose:
            print(f"fastloop: num samples to hold = {self.num_samples_to_hold}")
        # this just holds a number which increments every time the loop runs
        # TODO get rid of this
        self.index = 0


        # Define the instance of the object that will hold all the data
        self.fastdata = fast_data()
        self.slowdata = slow_data()

        # Set up the sensor
        self.sensor = sensor.sensor(main_path = self.main_path,verbose = self.verbose)

        """
        # Correction equations:
        # Line to fit the flow drift - will hold polynomial fit parameters
        self.vol_drift_params = []

        # spline curve to fit the volume minima
        self.vol_corr_spline = []
        self.tcal_min = [] # minimum time over which the spline curve applies
        self.tcal_max = [] # maximum time overwhich the spline curve applies
        """

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

        if self.verbose:
            # debugging: print an update of what we're doing
            print("\nfastloop: Index =  %d" % self.index)

        # record the update time
        self.update_time = datetime.utcnow()
        self.fastdata.t_obj = self.add_new_point(self.fastdata.t_obj, self.update_time, self.num_samples_to_hold)
        self.fastdata.t =     self.add_new_point(self.fastdata.t, self.update_time.timestamp(), self.num_samples_to_hold)
        self.fastdata.dt = self.fastdata.t - self.fastdata.t[0]

        # if there's at least two elements in the vector, calculate the real average delta between samples
        if len(self.fastdata.dt) >= 2:
            self.ts_real = np.abs(np.mean(self.fastdata.dt[1:] - self.fastdata.dt[:-1]))
            self.fs = 1.0/self.ts_real



        # read the sensor pressure and flow data
        # there's probably a cleaner way to do this, but oh well...
        self.sensor.read()
        self.fastdata.p1   = self.add_new_point(self.fastdata.p1,   self.sensor.p1,   self.num_samples_to_hold)
        self.fastdata.p2   = self.add_new_point(self.fastdata.p2,   self.sensor.p2,   self.num_samples_to_hold)
        self.fastdata.dp   = self.add_new_point(self.fastdata.dp,   self.sensor.dp,   self.num_samples_to_hold)
        self.fastdata.flow = self.add_new_point(self.fastdata.flow, self.sensor.flow, self.num_samples_to_hold)


        # filter the flowdata
        #self.fastdata.flow = signal.savgol_filter(self.fastdata.flow_raw,75,2)

        #detrend the flow
        self.fastdata.flow = signal.detrend(self.fastdata.flow,type = 'constant')


        # calculate the raw volume
        # volume is in liters per minute! so need to convert fs from (1/s) to (1/m)
            # fs (1/min) = fs (1/s) * 60 (s/min)
        self.fastdata.vol_raw = np.cumsum(self.fastdata.flow)/(self.fs*60.0)

        # apply the volume spline correction
        if self.verbose:
            if len(self.fastdata.t) >= self.num_samples_to_hold:
                print("fastloop: vectors are full")
                print(f"fastloop: fs = {self.fs}")


        try:
            self.find_vol_peaks()


        except Exception as e:
            print("fastloop: could not run peakfinder")
            print(e)
        
        if len(self.fastdata.index_of_min) >= 2:
            self.calculate_vol_drift_spline()
            self.apply_vol_corr()
        else:
            self.fastdata.vol = self.fastdata.vol_raw
            self.fastdata.v_drift = 0.0*self.fastdata.vol_raw

            self.vol = np.copy(self.vol_raw)
            self.v_drift = 0.0*self.vol_raw
            
        # tell the newdata signal to emit every time we update the data
        self.newdata.emit(self.fastdata)

    def find_vol_peaks(self):
        """
        ## find the min and max of the volume signal using peak finder ##
        # times of the volume min and max
        self.slowdata.vmin_times = np.array([])
        self.slowdata.vmax_times = np.array([])

        #TODO delete this
        # Old code from monitor_v7.py
        negative_mean_subtracted_volume = [-1*(v-np.mean(self.vol)) for v in self.vol]
        i_valleys = breath_detect_coarse(negative_mean_subtracted_volume,fs = self.fs,plotflag = False)
        self.i_valleys = i_valleys

        """
        # step 1: find index of min and max
        self.fastdata.index_of_min = utils.breath_detect_coarse(-1*self.fastdata.vol_raw, fs = self.fs)


        if self.verbose:
            print(f"fastdata: found {len(self.fastdata.index_of_min)} peaks!")
        # step 2:
        self.fastdata.vmin_times = self.fastdata.t[self.fastdata.index_of_min]
        self.fastdata.vmin = self.fastdata.vol_raw[self.fastdata.index_of_min]


    def calculate_vol_drift_spline(self):
        # correct the volume signal by ensuring that the lung volume is zero after every breath

        # fit a spline to the detrended volume minima
        self.fastdata.vol_corr_spline = interpolate.interp1d(self.fastdata.vmin_times,self.fastdata.vmin,kind = 'linear',fill_value = 'extrapolate')

    def apply_vol_corr(self):
        # this uses the current volume minima spline calculation to correct the volume by pinning all the minima to zero



        if (self.fastdata.vol_corr_spline is None) or (len(self.fastdata.index_of_min) == 0):
            if self.verbose:
                print("fastloop: no spline fit to apply to volume data")
            self.fastdata.vol = np.copy(self.fastdata.vol_raw)
            self.fastdata.v_drift = 0.0 * (self.fastdata.vol_raw)

            pass

        else:
            # calculate the drift volume using the spline. because we made the spline interpolate it will work outside the correction
            self.fastdata.v_drift = self.fastdata.vol_corr_spline(self.fastdata.t)

            # calculate the corrected volume
            self.fastdata.vol = self.fastdata.vol_raw - self.fastdata.v_drift


            if self.verbose:
                print("fastloop: applied spline volume correction")
                print(f"fastloop: max V = {np.max(self.fastdata.vol)}")

    def run(self):
        if self.verbose:
            print("fastloop: starting fast Loop")
        self.timer = QtCore.QTimer()
        self.timer.setInterval(self.ts)
        self.timer.timeout.connect(self.update)
        self.timer.start()
        self.exec() # YOU NEED THIS TO START UP THE THREAD!



###################################################################################
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


        # time between samples
        self.ts = update_time #ms

        # loop sample frequency - starts out as 1/self.ts but then is updated to the real fs
        self.fs = 1.0/self.ts




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
        if self.verbose:
            print(f"slowloop: slowloop.fastdata.dp = {self.fastdata.dp}")

        # detrend the volume signal
        self.detrend_volume()

        # try to fit a spline
        try:
            # find the volume peaks
            self.find_vol_peaks()

            # calculate the sline through the volume minima
            self.calculate_vol_drift_spline()

            # apply the volume correction
            self.apply_vol_corr()

            # calculate breath parameters
            #self.calculate_breath_params()
        except:
            if self.verbose:
                print("slowloop: no breathdetected")


        # tell main that there's new slow data: emit the newdata signal
        self.newdata.emit(self.slowdata)



    def apply_vol_corr(self):
        # this uses the current volume minima spline calculation to correct the volume by pinning all the minima to zero

        if self.verbose:
            print("slowloop: correcting volume")

        # calculate the drift volume using the spline. because we made the spline interpolate it will work outside the correction
        self.fastdata.v_drift = self.slowdata.vol_corr_spline(self.fastdata.t)

        # calculate the corrected volume
        self.fastdata.vol = self.fastdata.vol_raw - self.fastdata.v_drift


    def calculate_vol_drift_spline(self):
        # correct the volume signal by ensuring that the lung volume is zero after every breath

        # fit a spline to the detrended volume minima
        self.slowdata.vol_corr_spline = interpolate.interp1d(self.slowdata.vmin_times,self.slowdata.vmin_detrend,kind = 'linear',fill_value = 'extrapolate')

    def detrend_volume(self):
        # fit a line through the volume to get rid of any slow drifts. this helps the plots stay nice if there's no breaths

        # step 1: calculate the linear drift and spline fits to minimum
        self.slowdata.vol_drift_params = np.polyfit(self.fastdata.t,self.fastdata.vol_raw,1)

        # step 2: detrend the raw volume
        self.fastdata.vol_detrend = self.fastdata.vol_raw - np.polyval(self.slowdata.vol_drift_params,self.fastdata.t)

    def find_vol_peaks(self):
        """
        ## find the min and max of the volume signal using peak finder ##
        # times of the volume min and max
        self.slowdata.vmin_times = np.array([])
        self.slowdata.vmax_times = np.array([])

        #TODO delete this
        # Old code from monitor_v7.py
        negative_mean_subtracted_volume = [-1*(v-np.mean(self.vol)) for v in self.vol]
        i_valleys = breath_detect_coarse(negative_mean_subtracted_volume,fs = self.fs,plotflag = False)
        self.i_valleys = i_valleys

        """
        fastdata_samplerate = 1 / (self.fastdata.dt[-1] - self.fastdata.dt[-2])
        if self.verbose:
            print(f"slowloop: fastdata_samplerate = {fastdata_samplerate}")


        # step 1: find index of min and max
        self.slowdata.index_of_min = utils.breath_detect_coarse(-1*self.fastdata.vol_detrend, fs = fastdata_samplerate)
        self.slowdata.index_of_max = utils.breath_detect_coarse(self.fastdata.vol_detrend,    fs = fastdata_samplerate)
        if self.verbose:
            print(f"slowloop: found {len(self.slowdata.index_of_min)} peaks!")
        # step 2:
        self.slowdata.vmin_times = self.fastdata.t[self.slowdata.index_of_min]
        self.slowdata.vmax_times = self.fastdata.t[self.slowdata.index_of_max]
        self.slowdata.vmin_detrend = self.fastdata.vol_detrend[self.slowdata.index_of_min]


    def calculate_breath_params(self):
        """

        ## THINGS THAT HOLD SINGLE VALUES ##
        # times from the last breath
        self.tsi = [] # start time of inspiration
        self.tei = [] # end time of inspiration
        self.tse = [] # start time of expiration (same as tei)
        self.tee = [] # end time of expiration

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


    def update_fast_data(self,fastdata):
        # this is a slot connected to mainwindow.newrequest
        # this takes the fastdata from the main window and updates the internal value of fastdata
        if self.verbose:
            print(f"slowloop: received new fastdata from main")
        #self.fastdata = fastdata



    def run(self):
        print("starting slowloop")
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