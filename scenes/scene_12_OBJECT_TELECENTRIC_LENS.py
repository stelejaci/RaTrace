import numpy as np
from utils.varia import mm, nm, deg
from utils.optics import N_glass
from light import point_source_class
from elements import black_plate_class, spherical_lens_class

def load_scene():
    lens_diameter, lens_thickness, lens_f = 80*mm, 20*mm, 100*mm
    aperture_diameter_inner, aperture_diameter_outer, x_ap = 5*mm, 30*mm, lens_f
    source_x_far = -200*mm
    source_x_close = -100*mm

    source_far_top      = point_source_class.PointSourceClass(p0=np.array([source_x_far,    0]), n0=np.array([1, 0]), wavelength=660 * nm, fan_angle=5 * deg, intensity_distribution='equiangular', plot_color=(1, 0, 0, 1))
    source_far_bottom   = point_source_class.PointSourceClass(p0=np.array([source_x_far,  -25]), n0=np.array([1, 0]), wavelength=660 * nm, fan_angle=5 * deg, intensity_distribution='equiangular', plot_color=(0, 1, 0, 1))
    source_close_top    = point_source_class.PointSourceClass(p0=np.array([source_x_close,  0]), n0=np.array([1, 0]), wavelength=660 * nm, fan_angle=5 * deg, intensity_distribution='equiangular', plot_color=(0, 0, 1, 1))
    source_close_bottom = point_source_class.PointSourceClass(p0=np.array([source_x_close,-25]), n0=np.array([1, 0]), wavelength=660 * nm, fan_angle=5 * deg, intensity_distribution='equiangular', plot_color=(1, 0, 1, 1))
    lens = spherical_lens_class.SphericalLensClass(p0=np.array([-lens_thickness/2,0]), n0=np.array([-1, 0]), f=lens_f, diameter=lens_diameter, thickness=lens_thickness, N=N_glass, plot_resolution=1*mm)
    aperture_top    = black_plate_class.BlackPlateClass(p0=np.array([x_ap,  aperture_diameter_inner/2+(aperture_diameter_outer-aperture_diameter_inner)/4]), n0=np.array([-1, 0]), length=(aperture_diameter_outer-aperture_diameter_inner)/2, thickness=1 * mm)
    aperture_bottom = black_plate_class.BlackPlateClass(p0=np.array([x_ap, -aperture_diameter_inner/2-(aperture_diameter_outer-aperture_diameter_inner)/4]), n0=np.array([-1, 0]), length=(aperture_diameter_outer-aperture_diameter_inner)/2, thickness=1 * mm)
    beamdump = black_plate_class.BlackPlateClass(p0=np.array([2.0*lens_f,0]), n0=np.array([-1, 0]), length=100 * mm, thickness=5 * mm)
    info = 'Object-telecentric lens'
    return [source_far_top, source_far_bottom, source_close_top, source_close_bottom, lens, aperture_top, aperture_bottom, beamdump, info]
