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
    scene_name = 'scene_20_SCENE_WITH_PARAMETERS_IN_A_LOOP'
    simulation.load_scene(f"../scenes/{scene_name}.py", param=[f_lens[iter], p_BD[iter]])

    # Set the number of rays per source to raytrace. This overrides the number in the config.ini file
    simulation.set_nr_of_rays_per_source(100)

    # Finally, raytrace the entire thing
    simulation.run()

    # Write the list of rays into a text file
    export_filename = f"../scenes/{scene_name}_f{f_lens[iter]}mm_pBD{p_BD[iter]}mm.txt"
    varia.output_source_data_in_text_file(export_filename, simulation.sources, add_timestamp=False)

    # If it simulates too quickly (imagine the horror), output files are written within the same second, and overwrites them
    time.sleep(1)