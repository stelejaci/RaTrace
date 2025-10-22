import numpy as np
from utils.varia import mm, deg
from light import point_source_class, plane_source_class
from elements import semi_transparent_mirror_class, black_plate_class


def load_scene():
    light = plane_source_class.PlaneSourceClass(p0=np.array([-30,0]), n0=np.array([1,0]), diameter=5*mm)
    # light = point_source_class.PointSourceClass(p0=np.array([-30,0]), n0=np.array([1,0]), fan_angle=10*deg)
    mirror = semi_transparent_mirror_class.SemiTransparentMirror(p0=np.array([0,0]), n0=-np.array([1,1]), length=10*mm, transmission=0.5)
    beam_dump_1 = black_plate_class.BlackPlateClass(p0=np.array([10,0]), n0=np.array([-1,0]), length=10 * mm, thickness=0.5 * mm)
    beam_dump_2 = black_plate_class.BlackPlateClass(p0=np.array([0,-10]), n0=np.array([0,1]), length=10 * mm, thickness=0.5 * mm)
    info = 'Semitransparent mirror'
    return [light, mirror, beam_dump_1, beam_dump_2, info]