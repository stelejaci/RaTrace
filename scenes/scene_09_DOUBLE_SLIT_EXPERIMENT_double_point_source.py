import numpy as np
from utils.varia import mm,nm, µm, deg
from light import double_point_source_class
from display import imager_class

def load_scene():
    double_beam = double_point_source_class.DoublePointSourceClass(p0=np.array([0,0]), n0=np.array([1,0]), spacing=0.1*mm, wavelength=450*nm, fan_angle=1*deg, intensity_distribution='random', plot_color='rainbow')
    imager = imager_class.ImagerClass(p0=np.array([100, 0])*mm, n0=np.array([-1,0]), length=2*mm, pixel_size=0.001*mm+0*double_beam.wavelength)
    info = 'Double slit experiment'
    return [double_beam, imager, info]
