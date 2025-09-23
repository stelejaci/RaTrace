import sys
from PyQt5 import QtWidgets
from datetime import datetime


from utils.configuration_class import config
from raytracer import simulation_class
from utils import varia




# Start a QApplication in order to be able to use signals
app = QtWidgets.QApplication(sys.argv)

# Instantiate a simulation object
simulation = simulation_class.SimulationClass()

# Load a static scene file
simulation.load_scene('../scene_01_Hello_world.py')

# Set the number of rays per source to raytrace. This overrides the number in the config.ini file
# simulation.set_nr_of_rays_per_source(10)

# Finally, raytrace the entire thing
simulation.run()

# Write the list of rays into a text file
current_time = datetime.now()
timestamp_string = current_time.strftime("%Y%m%d_%H%M%S")
export_filename_full = f"G:/test_{timestamp_string}.txt"
varia.output_source_data_in_text_file(export_filename_full, simulation.sources)
