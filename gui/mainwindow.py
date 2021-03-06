#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:36:10 2020


mainwindow.py

This is where the mainwindow of the GUI is defined


@author: nlourie
"""

from PyQt5 import uic, QtCore, QtGui, QtWidgets
import pyqtgraph as pg
import os
import sys
import numpy as np

# add the wsp directory to the PATH
main_path = os.path.dirname(os.getcwd())
sys.path.insert(1, main_path)


# import custom modules
#from sensor import sensor
from data_handler import data_handler
from utils import utils

class MainWindow(QtWidgets.QMainWindow):

    """
    The class taking care for the main window.
    It is top-level with respect to other panes, menus, plots and
    monitors.
    """

    # this is a signal that sends fastdata from the mainloop to the slowloop
    request_from_slowloop = QtCore.pyqtSignal(object)
    request_to_update_cal = QtCore.pyqtSignal(object)

    def __init__(self,  config, main_path, mode = 'normal', verbose = False,simulation = False,logdata = False,*args, **kwargs):

        """
        Initializes the main window for the MVM GUI. See below for subfunction setup description.
        """

        super(MainWindow, self).__init__(*args, **kwargs)
        #uic.loadUi('mainwindow.ui', self)  # Load the .ui file

        # set mode
        if mode.lower() == 'debug':
            fast_update_time = 100
            slow_update_time = 1000
            mode_verbose = True

        else:
            fast_update_time = 10
            slow_update_time = 1000
            mode_verbose = False

        # configuration
        self.config = config

        # run in verbose mode? do this if requested specifically or if running in debug mode
        self.verbose = (verbose) or (mode_verbose)

        # define the top level main path
        self.main_path = main_path
        if self.verbose:
            print(f"main: main path = {self.main_path}")

        # define the slow and fast data classes that hold the data generated by the loops
        self.fastdata = data_handler.fast_data()
        self.slowdata = data_handler.slow_data()

        # Start up the fast loop (data acquisition)
        self.fast_loop = data_handler.fast_loop(main_path = self.main_path, update_time = fast_update_time, simulation = simulation, logdata = logdata,verbose = self.verbose)
        self.fast_loop.start()
        self.fast_loop.newdata.connect(self.update_fast_data)

        # start up the slow loop (calculations)
        self.slow_loop = data_handler.slow_loop(main_path = self.main_path, update_time = slow_update_time, verbose = self.verbose)
        self.slow_loop.start()
        self.slow_loop.newdata.connect(self.update_slow_data)

        # if the slowloop requests new data, send it the current fastdata
        self.slow_loop.request_fastdata.connect(self.slowloop_request)
        self.request_from_slowloop.connect(self.slow_loop.update_fast_data)


        ### GUI stuff ###
        self.setWindowTitle("Open Respiratory Monitor")

        self.graph1 = pg.PlotWidget()
        self.graph2 = pg.PlotWidget()
        self.graph3 = pg.PlotWidget()

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.graph1)
        layout.addWidget(self.graph2)
        layout.addWidget(self.graph3)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)

        # make the window with a graph widget
        self.setCentralWidget(widget)

        # set the plot properties
        self.graph1.setBackground('k')
        self.graph1.showGrid(x = True, y = True)
        self.graph2.setBackground('k')
        self.graph2.showGrid(x = True, y = True)
        self.graph3.setBackground('k')
        self.graph3.showGrid(x = True, y = True)

        # Set the label properties with valid CSS commands -- https://groups.google.com/forum/#!topic/pyqtgraph/jS1Ju8R6PXk
        labelStyle = {'color': '#FFF', 'font-size': '12pt'}
        self.graph1.setLabel('left','P','cmH20',**labelStyle)
        self.graph2.setLabel('left','Flow','L/m',**labelStyle)
        self.graph3.setLabel('left','V','mL',**labelStyle)
        self.graph3.setLabel('bottom', 'Time', 's', **labelStyle)

        # change the plot range
        self.graph1.setYRange(-40,40,padding = 0.1)
        self.graph2.setYRange(-100,100,padding = 0.1)
        self.graph3.setYRange(0,1500,padding = 0.1)

        # make a QPen object to hold the marker properties
        yellow = pg.mkPen(color = 'y',width = 2)
        pink = pg.mkPen(color = 'm', width = 2)
        green = pg.mkPen(color = 'g', width = 2)
        
        # define the curves to plot
        self.data_line1 = self.graph1.plot(self.fastdata.dt,    self.fastdata.p1,       pen = yellow)
        #self.data_line1b = self.graph1.plot(self.fastdata.dt,   self.fastdata.p2, pen = bluepen)
        self.data_line2 = self.graph2.plot(self.fastdata.dt,    self.fastdata.flow,     pen = pink)
        self.data_line3 = self.graph3.plot(self.fastdata.dt,    self.fastdata.vol*1000,      pen = green)
        # update the graphs at regular intervals (so it runs in a separate thread!!)
        # Stuff with the timer
        self.t_update = 10 #update time of timer in ms
        self.timer = QtCore.QTimer()
        self.timer.setInterval(self.t_update)
        self.timer.timeout.connect(self.update_plots)
        self.timer.start()



    ### gui-related functions
    def update_plots(self):

        # update the plots with the new data

        self.data_line1.setData(self.fastdata.dt,   self.fastdata.p1)
        #self.data_line1b.setData(self.fastdata.dt,   self.fastdata.p2)
        self.data_line2.setData(self.fastdata.dt,   self.fastdata.flow)
        self.data_line3.setData(self.fastdata.dt,   self.fastdata.vol*1000) #update the data


    ### slots to handle data transfer between threads ###
    def update_fast_data(self,data):
        if self.verbose:
            print("main: received new data from fastloop!")

        self.fastdata = data

    def update_slow_data(self,data):
        if self.verbose:
            print("main: received new data from slowloop!")
        self.slowdata = data

        #os.system('cls' if os.name == 'nt' else 'clear')
        data.print_data()

    def slowloop_request(self):
        if self.verbose:
            print(f"main: received request for data from slowloop")
        self.request_from_slowloop.emit(self.fastdata)
