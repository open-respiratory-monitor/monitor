#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:12:25 2020

@author: nlourie
"""
from PyQt5 import uic, QtCore, QtGui, QtWidgets
import time
import sounddevice
import soundfile

class audio_alarm(QtCore.QThread):
    
    def __init__(self,alarm_sound_high_file = 'hamilton_high.wav', alarm_sound_low_file = 'hamilton_low.wav'):
        
        QtCore.QThread.__init__(self)
        self.alarm_sound_high = soundfile.read(alarm_sound_high_file)
        self.alarm_sound_low  = soundfile.read(alarm_sound_low_file)
        self.set_priority('high')
        
    
    def __del__(self):
        self.wait()
    
    def set_priority(self,priority):
        if priority.lower() == 'low':
            self.priority = 'low'
            self.alarm_sound = self.alarm_sound_low
        else:
            self.priority = 'high'
            self.alarm_sound = self.alarm_sound_high
            
        self.sound = self.alarm_sound[0]
        self.sound_fs = self.alarm_sound[1]
        self.sound_length = len(self.sound) * (1/self.sound_fs)*1000 # ms
    
    def sound_once(self):
        sounddevice.play(self.sound, self.sound_fs)
        
    def sound_continuous(self,priority = 'high'):
        
        self.set_priority(priority)
        self.start()
    
    def sound_high_continuous(self):
        self.sound_continuous(priority = 'high')
        
    def sound_low_continuous(self):
        self.sound_continuous(priority = 'low')
    
    def print_pushed(self):
        print('pushed')
    
    def silence(self):
        if self.isRunning():
            self.terminate()
        else:
            pass
    
    def run(self):
        # this starts the Qthread with a timer that just plays the alarm each time it loops
        self.timer = QtCore.QTimer()
        self.timer.setInterval(self.sound_length)
        self.timer.timeout.connect(self.sound_once)
        self.timer.start()
        self.exec() # YOU NEED THIS TO START UP THE THREAD!

if __name__ == '__main__':
    
    class MainWindow(QtWidgets.QMainWindow):
    
        # Define the 
        
        
        
        def __init__(self, *args, **kwargs):
            super(MainWindow, self).__init__(*args, **kwargs)
            
            self.alarm = audio_alarm()
            
            
            layout = QtWidgets.QVBoxLayout()
    
            sound_high_button = QtWidgets.QPushButton("Sound High Alarm")
            sound_high_button.pressed.connect(self.alarm.sound_continuous)
            
            sound_low_button = QtWidgets.QPushButton("Sound Low Alarm")
            sound_low_button.pressed.connect(self.alarm.sound_low_continuous)
            
            silence_button = QtWidgets.QPushButton("Silence")
            silence_button.pressed.connect(self.alarm.silence)
            
            layout.addWidget(sound_low_button)
            layout.addWidget(sound_high_button)
            layout.addWidget(silence_button)
            
            w = QtWidgets.QWidget()
            w.setLayout(layout)
    
            self.setCentralWidget(w)
    
            self.show()
        
        def closeEvent(self,event):
            QtWidgets.QApplication.quit()  
            
    # Standard way to start up the event loop in GUI mode
    app = QtWidgets.QApplication([])
    #app.setQuitOnLastWindowClosed(False) #<-- otherwise it will quit once all windows are closed

    window = MainWindow()
    app.exec_()     
     
     