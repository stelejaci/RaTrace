import time
from PyQt5 import QtCore, QtWidgets
import sys,os
from utils.configuration_class import config
from gui import simulation_gui_class
from raytracer import simulation_class
import numpy as np

angles = np.linspace(-30,0,31,endpoint=True)


def run_gui_iteration(app, iter):
    simulation = simulation_class.SimulationClass()
    simulation_GUI = simulation_gui_class.SimulationGuiClass(simulation)
    simulation_GUI.show()

    if config.getboolean('scenes', 'load_scene_at_startup'):
        scene_file = os.path.join(config.get('scenes', 'scenes_folder'), config.get('scenes', 'scene_file'))
        simulation.load_scene(scene_file, param=angles[iter])
        # print(simulation)

    if config.getboolean('simulation', 'start_simulation_at_startup'):
        simulation.run()

    simulation_GUI.canvas.axis_lims = [-0.1, 0.1, -0.1, 0.1]
    simulation_GUI.canvas.graph.axis(simulation_GUI.canvas.axis_lims)

    simulation_GUI.repaint()
    QtCore.QTimer.singleShot(2000, lambda: take_screenshot(app, simulation_GUI, iter))

    if config.getboolean('simulation', 'exit_after_run'):
        # print(f'Closing GUI')
        QtCore.QTimer.singleShot(3000, simulation_GUI.close)


def take_screenshot(app, window, iter):
    pixmap = window.grab()
    pixmap.save(f'G:/screenshot_{iter:04d}.png')
    print(f"Screenshot {iter} saved")



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    for iter in range(len(angles)):
        run_gui_iteration(app, iter)
        # Process events until the current iteration is complete
        while app.activeWindow():
            app.processEvents()
            QtCore.QThread.msleep(100)

    sys.exit(app.exec_())