import sys
from PyQt5 import QtWidgets
from raytracer import simulation_class
from utils import varia


# Start a QApplication in order to be able to use signals
app = QtWidgets.QApplication(sys.argv)

# Instantiate a simulation object
simulation = simulation_class.SimulationClass()

# Load a static scene file
scene_name = 'scene_07_THREE_LENS_OBJECTIVE_multiple_parallel_beams_and_lenses'
simulation.load_scene(f"../scenes/{scene_name}.py")

# Set the number of rays per source to raytrace. This overrides the number in the config.ini file
simulation.set_nr_of_rays_per_source(20)

# Finally, raytrace the entire thing
simulation.run()

# Write the list of rays into a text file
export_filename = f"../scenes/{scene_name}.txt"
varia.output_source_data_in_text_file(export_filename, simulation.sources, add_timestamp=True)
