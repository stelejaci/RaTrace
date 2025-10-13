# RaTrace

---

[![Buy Me A Coffee](https://img.shields.io/badge/Buy%20Me%20A%20Coffee-support%20my%20work-FFDD00?style=flat&labelColor=101010&logo=buy-me-a-coffee&logoColor=white)](https://coff.ee/stelejaci)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<p align="center">
  RaTrace is a 2D raytracer with a easy-to-use graphical user interface, written in Python.
</p>

<p align="center">
  <img src="assets/RaTrace_medium.png", alt="RaTrace", width=50, height=50, style="display: block; margin: 0 auto" />
</p>

<p><i>Mind you that this is very much a work-in-progress. This started as a personal hobby project, never to be released to the outside world. The code is in many places sub-optimal, but it works for my intent of use.</i></p>


---

## Table of contents
* [Overview](#Overview)
* [Examples](#Examples)
* [Syntax](#Syntax)

---

## Overview

<p align="center"> <img src="assets/screenshot_01.png", alt="scene_01_Hello_world", width=400, height=200, style="display: block; margin: 0 auto" /> </p>

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
* Support for:
  * Light sources: point, plane, laser
  * Glass elements: spherical lens, ideal lens, slab, prism (biprism, microlens array, asphere, filter)
  * Mirrors: flat, parabolic, semi-transparant
  * Surfaces: black absorber, diffuse scattering (diffuse sphere)
  * Targets: display surface, imager
* ...

<p align="center">
  <img src="assets/scene_01_Hello_world_02.png", alt="scene_01_Hello_world", width=400, height=150, style="display: block; margin: 0 auto" />
</p>

---

## Examples

---
## Syntax

