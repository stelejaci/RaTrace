import numpy as np
from utils.varia import mm, nm, deg
from utils.optics import N_glass
from light import plane_source_class
from elements import spherical_lens_class, black_plate_class, flat_mirror_class, parabolic_mirror_class
from display import display_class

def load_scene():

    L_tel, px_exit, D_lens, T_wall, f_eye = 1000, 100, 20, 3, 17

    beam1  = plane_source_class.PlaneSourceClass(p0=np.array([-100, 5]), n0=np.array([1,0.000]), wavelength=660*nm, diameter=180*mm, intensity_distribution='equidistant', plot_color=(1,0,0,1))
    beam2  = plane_source_class.PlaneSourceClass(p0=np.array([-100, 0]), n0=np.array([1,0.002]), wavelength=660*nm, diameter=180*mm, intensity_distribution='equidistant', plot_color=(0,1,0,1))
    beam3  = plane_source_class.PlaneSourceClass(p0=np.array([-100,-5]), n0=np.array([1,0.004]), wavelength=660*nm, diameter=180*mm, intensity_distribution='equidistant', plot_color=(0,0,1,1))

    tubewall1 = black_plate_class.BlackPlateClass(p0=np.array([500,  100]), n0=np.array([0,-1]), length=L_tel*mm, thickness=T_wall, plot_color='black')
    tubewall2 = black_plate_class.BlackPlateClass(p0=np.array([(px_exit+D_lens/2+L_tel)/2, -100]), n0=np.array([0,1]), length=L_tel-px_exit-D_lens/2, thickness=T_wall, plot_color='black')
    tubewall3 = black_plate_class.BlackPlateClass(p0=np.array([(px_exit-D_lens/2)/2, -100]), n0=np.array([0,1]), length=px_exit-D_lens/2, thickness=T_wall, plot_color='black')

    parabolicmirror = parabolic_mirror_class.ParabolicMirrorClass(p0=np.array([L_tel,0]), n0=np.array([-1,0]), f=1000*mm, diameter=200+2*T_wall, thickness=10*mm, plot_resolution=10*mm, plot_color=(0,0,0.5,1))
    flatmirror = flat_mirror_class.FlatMirrorClass(p0=np.array([px_exit, 3]), n0=np.array([1,-1]), length=35 * mm, thickness=5 * mm, plot_color=(0,0,0.5,1))

    ocularwall1 = black_plate_class.BlackPlateClass(p0=np.array([(px_exit-D_lens/2), -130]), n0=np.array([1,0]),   length=55*mm, thickness=T_wall, plot_color='black')
    ocularwall2 = black_plate_class.BlackPlateClass(p0=np.array([(px_exit+D_lens/2), -130]), n0=np.array([-11,0]), length=55*mm, thickness=T_wall, plot_color='black')
    lens_ocular = spherical_lens_class.SphericalLensClass(p0=np.array([px_exit,-150]),  n0=np.array([0,1]), f=50*mm, thickness=5, diameter=D_lens, N=N_glass, blur_angle=0*deg, nr_of_secondary_rays=1)

    lens_eye = spherical_lens_class.SphericalLensClass(p0=np.array([px_exit,-180]),  n0=np.array([0,1]), f=f_eye, thickness=1, diameter=5*mm, N=N_glass, blur_angle=0*deg, nr_of_secondary_rays=1)
    retina = display_class.DisplayClass(p0=np.array([100.2, -180-17]), n0=np.array([0,1]), length=2*mm)

    info = 'Newtonian telescope with ocular and eye'

    return [beam1, beam2, beam3, tubewall1, tubewall2, tubewall3, parabolicmirror, flatmirror, lens_ocular, lens_eye, ocularwall1, ocularwall2, retina, info]

