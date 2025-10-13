<p align="center">
  <img src="assets/RaTrace_medium.png", alt="RaTrace", width=400, height=400, style="display: block; margin: 0 auto" />
</p>

<p align="center">
  A 2D raytracer with a easy-to-use graphical user interface, written in Python.
</p>

---

[![Buy Me A Coffee](https://img.shields.io/badge/Buy%20Me%20A%20Coffee-support%20my%20work-FFDD00?style=flat&labelColor=101010&logo=buy-me-a-coffee&logoColor=white)](https://coff.ee/stelejaci)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<p><i>!! Disclaimer !! Mind you that this is very much a work-in-progress. This started as a personal hobby project and I never had the intent to release this publicly.
The code is in many places sub-optimal and buggy, but it works for my intent of use. Ok, most of the time.</i>:blush:</p>

---

## Table of contents
* [Overview](#Overview)
* [Usage](#Usage)
* [Examples](#Examples)
* [Graphical user interface](#Graphical user interface)
* [Syntax](#Syntax)

---

## Overview

<p align="center"> <img src="assets/screenshot_01.png", alt="scene_01_Hello_world", width=800, height=400, style="display: block; margin: 0 auto" /> </p>

<b>Implemented features</b>
* GUI for 2D raytracing
* Scene creation via Python scripts
* Simulation of static scenes, with or without UI
* Automated scripts for looped simulations with different scenes
* Exact raytracing for analytically described elements (spherical, parabolic, flat surfaces)
* Accurate raytracing for segments-based, more "complex" elements
* "Fast" raytrace mode for ordered elements or "slow" mode for full raytracing
* Wavelength dispersion
* Tracking of ray phase information
* Export ray information to a text file
* Color coding rays: wavelength, rainbow, fixed, intensity-scaling
* Support for:
  * Light sources: point source, diffusing plane source, parallel plane source, laser source, virtual rays, double coherent point source
  * Glass elements: spherical lens, ideal lens, glass slab, prism, filter
  * Mirrors: flat, parabolic, semi-transparant
  * Surfaces: black absorber, diffuse scattering plane
  * Targets: display surface, imager

<b>To be implemented features</b>
* Glass elements: plano-convex lens, biprism, microlens array, asphere
* Light source: B/W image source
* Better error handling when scene is faulty
* Diffuse scattering sphere

<b>Known bugs</b>
* First screenshot in looped gui does not set the axis correctly
* Warning concerning colors

<p align="center">
  <img src="assets/scene_01_Hello_world_02.png", alt="scene_01_Hello_world", width=400, height=150, style="display: block; margin: 0 auto" />
</p>

---

## Usage

There are 4 ways to run RaTrace:
* The primary way to use RaTrace is with the GUI. This way most of the settings can be changed, and scenes could loaded:
> python.exe main_GUI.py
* With a GUI but automated and the possibility for taking screenshots. The program is closed when finished.
> python.exe main_GUI_in_loop.py
* Without GUI, the program finishes when finished.
> python.exe main_noGUI.py
* Without GUI but automated, the program finishes when finished.
> python.exe main_noGUI_in_loop.py

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

The scenes itself are Python code and are dynamically loaded whenever a new scene is loaded. See the next chapter for examples.

### Note:
The standard units for calculations in RaTrace are <b>mm</b> and <b>radians</b>. You can import and use other units from the utils.vara module in the following way:
```
from utils.varia import mm, nm, deg
wavelength_mm = 660*nm
angle_radians = 30*deg
```

---

## Examples

I invite u to browse and try out the different example scenes that are available in the scenes folder.

### Simple on-axis lens
Below you find the content of a very simple scene file. 
First, all the modules used in the scene are imported. Next the ```load_scene``` module must always be present and contains all necessary elements, in this case a point light, a ideal lens and a beam dump:
```
import numpy as np
from utils.varia import mm,nm, deg
from utils.optics import N_glass
from light import point_source_class
from elements import ideal_thin_lens_class, black_plate_class

def load_scene():
    light = point_source_class.PointSourceClass(p0=np.array([-15, 0]), n0=np.array([1,0]), wavelength=660*nm, fan_angle=30*deg)
    lens = ideal_thin_lens_class.IdealThinLensClass(p0=np.array([0,0]),  n0=np.array([-1,0]), f=10*mm, diameter=10*mm, N=N_glass)
    beam_dump = black_plate_class.BlackPlateClass(p0=np.array([30, 0]), n0=np.array([-1,0]), length=10*mm)
    info = 'Ideal lens projecting a point source'
    return [light, lens, beam_dump, info]
```

When the simulation is run in the GUI, this results in the following (clipped) screenshot:   

![](scenes/scene_01_SIMPLE_ON_AXIS_LENS_fan_beam_ideal_lens_beamdump.png)

### Chromatic aberration
In this example I ony show what's essential in the load_scene module and skip all boilerplate code before and after this section:
```
beam_red  = plane_source_class.PlaneSourceClass(p0=np.array([-20,0]), n0=np.array([1,0]),  diameter=20*mm, wavelength=660*nm, intensity_distribution='equidistant', plot_color='wavelength')
beam_green = plane_source_class.PlaneSourceClass(p0=np.array([-20,0.2]), n0=np.array([1,0]),  diameter=20*mm, wavelength=520*nm, intensity_distribution='equidistant', plot_color='wavelength')
beam_blue = plane_source_class.PlaneSourceClass(p0=np.array([-20,0.4]), n0=np.array([1,0]),  diameter=20*mm, wavelength=450*nm, intensity_distribution='equidistant', plot_color='wavelength')
lens = spherical_lens_class.SphericalLensClass (p0=np.array([0,  0]), n0=np.array([-1,0]), thickness=5*mm, f=40*mm, diameter=25 * mm, N=[1.7, 0.05])
beam_dump = black_plate_class.BlackPlateClass(  p0=np.array([40, 0]), n0=np.array([-1,0]), length=20 * mm, thickness=1*mm)
```
![](scenes/scene_04_CHROMATIC_ABERRATION_spherical_lens_with_dispersion.png)

### Newtonian telescope
The list of elements in the load_scene module becomes a bit too long to show here, so I refer to the file itself in case you are interested.

![](scenes/scene_08_NEWTONIAN_TELESCOPE_parabolic_mirror_lenses_display_1.png)

---

## Graphical user interface

---

## Syntax

### Light sources

### Glass elements

### Surfaces

### Mirrors

### Targets


