#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:35:15 2020


monitor.py

This is the top level monitor script for the open source respiratory monitor


@author: nlourie
"""

import sys
import os
import os.path
from PyQt5 import QtCore, QtWidgets
import yaml



from gui import mainwindow


# add the main directory to the PATH
main_path = os.getcwd()
sys.path.insert(1, main_path)

def main():
    """
    Main function.
    """
    # Load configuration
    
    settings_file = main_path + '/config/default_settings.yaml'
    print('settings file = ',settings_file)
    with open(settings_file) as fsettings:
        config = yaml.load(fsettings, Loader=yaml.FullLoader)
    #print('Config:', yaml.dump(config), sep='\n')
   
    config = []

    app = QtWidgets.QApplication(sys.argv)
    
    if 'verbose' in sys.argv:
        verbose = True
        print('running in verbose mode')
    else:
        verbose = False
        
    if 'sim' in sys.argv:
        simulation = True
        print('running in simulation mode')
    else:
        simulation = False
        
    if 'debug' in sys.argv:
        mode = 'debug'
        print('running in slow debug mode')
    else:
        mode = 'normal'
        
    if 'log' in sys.argv:
        logdata = True
        print('saving continuous log of pressure data')
    else:
        logdata = False
        
    if verbose:    
        print(f"toplevel: main path = {main_path}")
    
    if 'help' in sys.argv:
        print_help()
        
    try:    
        window = mainwindow.MainWindow(config = config, main_path = main_path, mode = mode, simulation = simulation, logdata = logdata, verbose = verbose)
    except Exception as e:
        print("Problem running app: ",e)
        print_help()
    
    
    window.show()
    app.exec_()

def print_help():
    print('Help File for monitor.py:')



if __name__ == "__main__":
    main()
    
    
