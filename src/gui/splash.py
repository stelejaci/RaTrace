from PyQt5 import QtGui, QtCore, QtWidgets
import time
from utils.configuration_class import config

ANIMATION_TIME = config.getint('view', 'splash_screen_transition')

app = QtWidgets.QApplication([])
splash_pix = QtGui.QPixmap('../assets/RaTrace_medium.png')
splash = QtWidgets.QSplashScreen(splash_pix)
# add fade to splashscreen
opaqueness = 0.0
step = 0.01
splash.setWindowOpacity(opaqueness)
splash.show()
while opaqueness < ANIMATION_TIME:
    splash.setWindowOpacity(opaqueness)
    time.sleep(step)  # Gradually appears
    opaqueness += step
time.sleep(1)  # hold image on screen for a while
splash.close()  # close the splash screen

# widget = YourWidget()
# widget.show() # This is where you'd run the normal application
# from PyQt5.QtCore import QTimer
# QTimer.singleShot(2000, splash.close)
# app.exec_()