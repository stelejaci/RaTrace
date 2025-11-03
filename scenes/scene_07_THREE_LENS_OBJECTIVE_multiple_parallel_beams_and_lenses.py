import numpy as np
from utils.varia import mm, nm
from light import plane_source_class
from elements import spherical_lens_class, black_plate_class, plano_spherical_lens_class, aperture_class

def load_scene():
    d1, l, d2, b1, b2, d3, s, d4 = 4.3*mm, 2.3*mm, 2.5*mm, 4.2*mm, 1.3*mm, 0.6*mm, 0.0001*mm, 4.7*mm
    beam1  = plane_source_class.PlaneSourceClass(p0=np.array([-50,   0]), n0=np.array([1,0.00]), wavelength=660*nm, diameter=15*mm, intensity_distribution='equidistant', plot_color=(1,0,0,1))
    beam2  = plane_source_class.PlaneSourceClass(p0=np.array([-50,  -7]), n0=np.array([1,0.15]), wavelength=660*nm, diameter=15*mm, intensity_distribution='equidistant', plot_color=(0,1,0,1))
    beam3  = plane_source_class.PlaneSourceClass(p0=np.array([-50, -14]), n0=np.array([1,0.31]), wavelength=660*nm, diameter=15*mm, intensity_distribution='equidistant', plot_color=(0,0,1,1))
    lens1 = plano_spherical_lens_class.PlanoSphericalLensClass(p0=np.array([-b1-d2-l-d1,0]), n0=np.array([-1,0]),   R=23.4*mm, thickness=d1, diameter=20, N=1.58315, plot_resolution=1*mm)
    lens2 = spherical_lens_class.SphericalLensClass(p0=np.array([-b1-d2,0]),      n0=np.array([-1,0]),   R0=  -81.3*mm, R1= 22.4*mm,   thickness=d2, diameter=20, N=1.58215, plot_resolution=1*mm)
    lens3 = spherical_lens_class.SphericalLensClass(p0=np.array([b2,0]),          n0=np.array([-1,0]),   R0= -322.0*mm, R1= 27.6*mm,   thickness=d3, diameter=20, N=1.58215, plot_resolution=1*mm)
    lens4 = spherical_lens_class.SphericalLensClass(p0=np.array([b2+d3+s,0]),     n0=np.array([-1,0]),   R0=   27.6*mm, R1=-48.3*mm,   thickness=d4, diameter=20, N=1.67110, plot_resolution=1*mm)
    aperture = aperture_class.ApertureClass(p0=np.array([0,0]),  n0=np.array([-1,0]), diameter_inner=14*mm, diameter_outer=26*mm)
    beam_dump   = black_plate_class.BlackPlateClass(  p0=np.array([93.6, 16]), n0=np.array([-1,0]), length= 45*mm, thickness=2.0*mm)
    info = ('Patent 1849681 (1931): Photographic three-lens objective')
    return [beam1, beam2, beam3, lens1, lens2, aperture, lens3, lens4, beam_dump, info]
