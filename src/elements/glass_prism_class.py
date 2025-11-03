import numpy as np
from utils.varia import mm, deg
from utils import varia
from utils.optics import N_glass
from utils import geometry
from elements import glass_element_class
import matplotlib.pyplot as plt


class GlassPrism(glass_element_class.GlassElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([1,1]), angle_apex=90*deg, length=10*mm, N=N_glass, generate_reflections=False, is_active=True, is_visible=True):
        self.angle_apex = angle_apex
        self.length = length
        self.N = N

        n0 = geometry.normalize(n0)
        r = geometry.orientation_from_normal(n0)
        pts = np.empty((0, 2))
        pts = np.append(pts, [p0 + r *length/2],                                  axis=0)
        pts = np.append(pts, [p0 - n0*length/2 * np.tan((180*deg-angle_apex)/2)], axis=0)
        pts = np.append(pts, [p0 - r *length/2],                                  axis=0)

        super().__init__(p0=p0, n0=n0, pts=pts, N=N, generate_reflections=generate_reflections, is_active=is_active, is_visible=is_visible)
        self.name = 'Prism'

    def __str__(self):
        txt = 'Prism'
        return txt

    def plot(self, graph):
        if self.is_visible:
            poly_face = plt.Polygon(self.pts, closed=True, facecolor='cyan', edgecolor='none', alpha=varia.alpha_from_N(self.N), zorder=0)
            poly_edge = plt.Polygon(self.pts, closed=True, fill=False, edgecolor='blue', linewidth=1, alpha=0.2, zorder=0)
            graph.add_patch(poly_face)
            graph.add_patch(poly_edge)
            super().plot(graph)
            # if config.getboolean('view', 'show_elements_properties'):
            #     canvas_class.plot_normals(graph, self)