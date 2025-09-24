import numpy as np
from utils.varia import mm, nm
from light import plane_source_class
from elements import glass_element_class, black_plate_class

def load_scene():
    fan_beam  = plane_source_class.PlaneSourceClass(p0=np.array([-40, 0]), n0=np.array([1,0]), wavelength=450*nm, diameter=10*mm, intensity_distribution='equidistant')
    glass_slab = glass_element_class.GlassParallelPlate(p0=np.array([0,-2]), n0=np.array([-1,-1]), thickness=10*mm, length=30*mm, N=1.5)
    beam_dump = black_plate_class.BlackPlateClass(p0=np.array([40, 0]), n0=np.array([-1,0]), length=20 * mm, thickness=2 * mm)
    info = ('SNELL''s LAW, parallel beam through glass slab on beam dump')
    return [fan_beam, glass_slab, beam_dump, info]
