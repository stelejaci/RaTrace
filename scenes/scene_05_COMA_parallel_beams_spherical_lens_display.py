import numpy as np
from utils.varia import mm, nm, µm
from utils.optics import N_glass
from light import plane_source_class
from elements import spherical_lens_class
from display import imager_class

def load_scene():
    beam1 = plane_source_class.PlaneSourceClass(p0=np.array([-10, 0]), n0=np.array([1,.00]), diameter=10*mm, wavelength=660*nm, intensity_distribution='random', plot_color=(1,0,0,1))
    beam2 = plane_source_class.PlaneSourceClass(p0=np.array([-10,-4]), n0=np.array([1,0.3]), diameter=10*mm, wavelength=660*nm, intensity_distribution='random', plot_color=(0,1,0,1))
    lens = spherical_lens_class.SphericalLensClass(p0=np.array([0, 0]), n0=np.array([-1,0]), f=20*mm, diameter=18*mm, thickness=5*mm, N=N_glass)
    display = imager_class.ImagerClass(p0=np.array([23, 4]), n0=np.array([-1,0]), length=10 * mm, pixel_size=5*µm)
    info = ('Coma in off-axis beams, with imager')
    return [beam1, beam2, lens, display, info]
