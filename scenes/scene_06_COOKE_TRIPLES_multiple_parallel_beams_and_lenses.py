import numpy as np
from utils.varia import mm, nm
from light import plane_source_class
from elements import spherical_lens_class, black_plate_class

def load_scene():
    wavelength = 660*nm
    t1, s1, t2, s2, t3 = 25.56*mm, 7.89*mm, 4.51*mm, 14.14*mm, 19.92*mm
    beam1  = plane_source_class.PlaneSourceClass(p0=np.array([-100,   0]), n0=np.array([1,0.0]), wavelength=wavelength, diameter=30*mm, intensity_distribution='equidistant', plot_color=(1,0,0,1))
    beam2  = plane_source_class.PlaneSourceClass(p0=np.array([-100, -10]), n0=np.array([1,0.1]), wavelength=wavelength, diameter=30*mm, intensity_distribution='equidistant', plot_color=(0,1,0,1))
    beam3  = plane_source_class.PlaneSourceClass(p0=np.array([-100, -20]), n0=np.array([1,0.2]), wavelength=wavelength, diameter=30*mm, intensity_distribution='equidistant', plot_color=(0,0,1,1))
    lens1 = spherical_lens_class.SphericalLensClass(p0=np.array([0,          0]), n0=np.array([-1,0]), R0=  62.63*mm, R1=-348.5*mm, thickness=t1, diameter=60, N=1.745, plot_resolution=1)
    lens2 = spherical_lens_class.SphericalLensClass(p0=np.array([t1+s1,      0]), n0=np.array([-1,0]), R0= -66.47*mm, R1= 39.62*mm, thickness=t2, diameter=50, N=1.647, plot_resolution=1)
    lens3 = spherical_lens_class.SphericalLensClass(p0=np.array([t1+s1+t2+s2,0]), n0=np.array([-1,0]), R0= 119.50*mm, R1=-52.41*mm, thickness=t3, diameter=50, N=1.697, plot_resolution=1)
    beam_dump = black_plate_class.BlackPlateClass(p0=np.array([149.3, 0]), n0=np.array([-1,0]), length=50*mm, thickness=2*mm)
    info = ('Patent 2503751 (1948): Photographic objective of the Cooke triplet type')
    return [beam1, beam2, beam3, lens1, lens2, lens3, beam_dump, info]
