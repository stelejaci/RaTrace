import time
from PyQt5 import QtCore, QtWidgets
import sys,os
from utils.configuration_class import config
from gui import simulation_gui_class
from raytracer import simulation_class
import numpy as np


# f_lens = [18, 16]
# p_BD   = [30, 40]
f_lens = [18, 16, 14, 12, 10]
p_BD   = [30, 40, 50, 60, 70]
axis_lims = [-20, 100, -20, 20]
# axis_X_lims = [-20, 100]

def run_gui_iteration(app, iter):
    simulation = simulation_class.SimulationClass()
    simulation_GUI = simulation_gui_class.SimulationGuiClass(simulation)
    # simulation_GUI.canvas.axis_lims = axis_lims
    simulation_GUI.show()

    scene_file = '../scenes/scene_10_SCENE_WITH_PARAMETERS_IN_A_LOOP.py'
    # scene_file = os.path.join(config.get('scenes', 'scenes_folder'), config.get('scenes', 'scene_file'))
    simulation.load_scene(scene_file, param=[f_lens[iter], p_BD[iter]])

    # Set the number of rays per source to raytrace. This overrides the number in the config.ini file
    simulation.set_nr_of_rays_per_source(10)

    simulation.run()
    # simulation_GUI.canvas.axis_lims = axis_lims
    # simulation_GUI.update_graphics()
    # simulation_GUI.canvas.graph.set_xlim(axis_X_lims[0], axis_X_lims[1])
    # simulation_GUI.canvas.draw()
    # simulation_GUI.set_axis_and_redraw(axis_lims)
    # simulation_GUI.canvas.axis_lims = axis_lims
    # simulation_GUI.canvas.graph.axis(axis_lims)
    # simulation_GUI.canvas.update_items()

    # Taking a screenshot and closing figure
    QtCore.QTimer.singleShot(2000, lambda: take_screenshot(app, simulation_GUI, iter))
    QtCore.QTimer.singleShot(3000, simulation_GUI.close)
    pass

def take_screenshot(app, window, iter):
    pixmap = window.grab()
    pixmap.save(f'../screenshot_{iter:04d}.png')
    print(f"Screenshot {iter} saved")



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    for iter in range(len(f_lens)):
        run_gui_iteration(app, iter)

        # Process events until the current iteration is complete
        while app.activeWindow():
            app.processEvents()
            QtCore.QThread.msleep(100)

    sys.exit(app.exec_())