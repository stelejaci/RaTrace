import numpy as np
from utils.varia import mm, nm
from utils.optics import N_glass
from light import plane_source_class
from elements import spherical_lens_class
from display import display_class

def load_scene():
    beam = plane_source_class.PlaneSourceClass(p0=np.array([-10, 0]), n0=np.array([1,0]), diameter=10*mm, wavelength=450*nm, intensity_distribution='equidistant', plot_color=(0,0,1,0.5))
    lens = spherical_lens_class.SphericalLensClass(p0=np.array([0, 0]), n0=np.array([-1,0]), f=15*mm, diameter=12*mm, thickness=3*mm, N=N_glass)
    display = display_class.DisplayClass(p0=np.array([17, 0]), n0=np.array([-1,0]), length=10 * mm)
    info = ('Spherical aberration, parallel beam via spherical lens on display')
    return [beam, lens, display, info]
