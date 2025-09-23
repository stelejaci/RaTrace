import sys
from PyQt5 import QtWidgets
from datetime import datetime
import time

from utils.configuration_class import config
from raytracer import simulation_class
from utils import varia




# Start a QApplication in order to be able to use signals
app = QtWidgets.QApplication(sys.argv)

for y in [0,-1,-2,-3]:
    # Instantiate a simulation object
    simulation = simulation_class.SimulationClass()

    # Load a scene file that takes an argument
    simulation.load_scene('../scenes/scene_01_Hello_world.py', param=y)

    # Set the number of rays per source to raytrace. This overrides the number in the config.ini file
    # simulation.set_nr_of_rays_per_source(1000)

    # Finally, raytrace the entire thing
    simulation.run()

    # Write the list of rays into a text file
    varia.output_source_data_in_text_file("G:/test.txt", simulation.sources, add_timestamp=True)

    # If it simulates too quickly (imagine the horror), output files are written within the same second, and overwrites them
    time.sleep(1)