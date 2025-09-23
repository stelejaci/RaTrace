import numpy as np
import os, sys
sys.path.append(os.path.abspath('..'))

from utils.varia import nm
from utils.optics import N_glass
from light.point_source_class import PointSourceClass
from elements.spherical_lens_class import SphericalLensClass
from elements.black_plate_class import BlackPlateClass


def load_scene():
    light = PointSourceClass(p0=np.array([0, 0]), n0=np.array([1,0]), wavelength=660*nm, fan_angle=np.radians(20), intensity=1, intensity_distribution='equiangular', plot_color='rainbow')
    lens = SphericalLensClass(p0=np.array([100,0]),  n0=np.array([-1,0]), f=50, thickness=15, diameter=50, N=N_glass, blur_angle=0, nr_of_secondary_rays=1, plot_resolution=1)
    beam_dump = BlackPlateClass(p0=np.array([200, 0]), n0=np.array([-1,0]), length=50, thickness=1, plot_color=(0.5,0.5,0.5,1))
    info = 'Simple scene'
    return [light, lens, beam_dump, info]
