from PyQt5 import QtCore, QtWidgets, QtGui
import sys
import os

from utils.configuration_class import config
from gui import simulation_gui_class, splash
from raytracer import simulation_class




app = QtWidgets.QApplication(sys.argv)

simulation = simulation_class.SimulationClass()

simulation_GUI = simulation_gui_class.SimulationGuiClass(simulation)
simulation_GUI.show()

if config.getboolean('scenes', 'load_scene_at_startup'):
    scene_file = os.path.join( config.get('scenes', 'scenes_folder'), config.get('scenes', 'scene_file') )
    simulation.load_scene(scene_file)
    # print(simulation)

# simulation_GUI.canvas.axis_lims = [53.2, 53.275, 114.1, 114.2]
# simulation_GUI.canvas.graph.axis(simulation_GUI.canvas.axis_lims)

sys.exit(app.exec_())

