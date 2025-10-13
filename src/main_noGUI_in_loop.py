import sys
from PyQt5 import QtWidgets
import time
from raytracer import simulation_class
from utils import varia

# Start a QApplication in order to be able to use signals
app = QtWidgets.QApplication(sys.argv)

f_lens = [18, 16, 14, 12, 10]
p_BD   = [30, 40, 50, 60, 70]

for iter in range(len(f_lens)):
    # Instantiate a simulation object
    simulation = simulation_class.SimulationClass()

    # Load a scene file that takes an argument
    simulation.load_scene('../scenes/scene_10_SCENE_WITH_PARAMETERS_IN_A_LOOP.py', param=[f_lens[iter], p_BD[iter]])

    # Set the number of rays per source to raytrace. This overrides the number in the config.ini file
    simulation.set_nr_of_rays_per_source(1000)

    # Finally, raytrace the entire thing
    simulation.run()

    # Write the list of rays into a text file
    varia.output_source_data_in_text_file("../test.txt", simulation.sources, add_timestamp=True)

    # If it simulates too quickly (imagine the horror), output files are written within the same second, and overwrites them
    time.sleep(1)