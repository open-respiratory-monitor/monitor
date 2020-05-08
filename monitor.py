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


from gui import mainwindow


# add the main directory to the PATH
main_path = os.getcwd()
sys.path.insert(1, main_path)

def main():
    """
    Main function.
    """
    # Load configuration
    config = []


    app = QtWidgets.QApplication(sys.argv)
    print(f"toplevel: main path = {main_path}")
    window = mainwindow.MainWindow(config = config, main_path = main_path, mode = 'normal', verbose = False)

    window.show()
    app.exec_()



if __name__ == "__main__":
    main()