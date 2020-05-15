#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 18:58:25 2020

utils.py

This is just a collection of useful functions used elsewhere in the code


@author: nlourie
"""
from scipy import signal
import numpy as np

def breath_detect_coarse(flow,fs,minpeak = 0.05,plotflag = False):
    """
    %% This function detects peaks of flow signal
    
    % Inputs:
    %           flow:       flow signal
    %           fs:         sampling frequency
    %           plotflag:   set to 1 to plot
    
    % Output:
    %           peak (location, amplitude)
    
    % Written in matlab by: Chinh Nguyen, PhD
    % Email: c.nguyen@neura.edu.au
    % Source: 
    
    % Updated on: 12 Nov 2015.
    % Ver: 1.0
    
    # Converted to python from Matlab by: Nate Lourie, PhD
    # Email: nlourie@mit.edu
    # Updated on: April, 2020
    
    """
    # detect peaks of flow signal
    minpeakwidth = fs*0.3
    peakdistance = fs*1.5
    #print('peakdistance = ',peakdistance)
    minPeak = minpeak # flow threshold = 0.05 (L/s)
    minpeakprominence = 0.05
    
    peak_index, _  = signal.find_peaks(flow, 
                                    height = minPeak,
                                    distance = peakdistance,
                                    prominence = minpeakprominence,
                                    width = minpeakwidth)
    
    #print('found peaks at index = ',peak_index)
    return peak_index


def zerophase_lowpass(x,lf,fs,force_length = False):
    # filter flow using S-G filter
    
    # filter flow
    N = 2 # the higher the sharper the peak is
    
    if force_length:
        l_lfilter = lf
    else:
    
    
    # length = 1s (cut off frequency = 0.965 Hz)
    # Refer. Schafer, 2011, What is S-G filter.
    # fc_n = (N+1)/(3.2M - 4.6) (normalized unit)
    # fc(Hz) = fc_n*fs/2;
    # N: order: M length;
    # fc=lf(Hz) -> fc_n = 2*lf/fs;
    # M = [(N+1)*fs/(2*lf)+4.6]/3.2;
    
    # l_lfilter = round(fs*lf); % the longer, the smoother
    # l_lfilter = round(fs*lf); % the longer, the smoother
    
        l_lfilter = np.int(np.round(((N+1)*fs/(2*lf)+4.6)/3.2))
        #print('l_lfilter = ',l_lfilter)
        if np.mod(l_lfilter,2) == 0:
            l_lfilter = l_lfilter + 1
            #print('made odd: l_lfilter = ',l_lfilter)
        
    # filter flow signal
    x_filt = signal.savgol_filter(x,polyorder = N,window_length = l_lfilter)
    return x_filt

def butter_lowpass(cutoff, fs, order=5):
    #from here: https://stackoverflow.com/a/25192640
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a