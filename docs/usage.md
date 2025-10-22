# Usage

Just download the files, and run one of the "main_XXX.py" files

There are 4 ways to run RaTrace:
* The primary way to use RaTrace is with the GUI. This way most of the settings can be changed, and new scenes can be loaded:
``` python main_GUI.py```
* With a GUI but automated and the possibility for taking screenshots. The program is closed when finished.
``` python main_GUI_in_loop.py```
* Without GUI, the program finishes when finished.
``` python main_noGUI.py```
* Without GUI but automated, the program finishes when finished.
``` python main_noGUI_in_loop.py```

When using the primary way-of-use via the GUI, the scene that is loaded is described in the config.ini file. 
Also, most of the options that are available in the GUI can also be set in the same config.ini file.

``` 
[scenes]
scene_file = scene.py
scenes_folder = D:/RaTrace/scenes
load_scene_at_startup = 1
start_simulation_after_loading_scene = 1
reset_axis_after_loading_scene = 0

[simulation]
nr_of_rays = 1000
use_phase_information = 1
generate_reflected_rays = 0
items_are_ordered = 0
export_folder = ../
add_timestamp_to_export_file = 1
exit_after_run = 1

[view]
show_axis_and_grid = 1
black_background = 0
nr_of_rays_to_plot = 100
show_elements_properties = 0
intensity_coded_colors = 0
show_pixels = 0
intensity_scaler = 0.1
show_noncolliding_rays = 1
splash_screen_transition = 0
```

The scenes itself are written in Python and are dynamically loaded whenever a new scene is loaded. See the next chapter for examples.