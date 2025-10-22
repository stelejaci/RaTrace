import numpy as np
import matplotlib.pyplot as plt
from utils.varia import mm, X, Y
from light import light_class
from elements import element_class
from utils import geometry


class FlatMirrorClass(element_class.ElementClass):
    def __init__(self, p0=np.array([0, 0]), n0=np.array([-1, 0]), length=10 * mm, thickness=1 * mm, plot_color=(0,0,0.5,1), is_active=True, is_visible=True):
        self.length     = length
        self.thickness  = thickness
        self.plot_color = plot_color

        pts = geometry.construct_plate(p0=p0, n=n0, thickness=thickness, length=length)
        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Flat mirror'


    def propagate_ray(self, ray):
        new_rays = list()
        if (self.n_coll is not None)    and    (geometry.point_is_on_PP_line(ray.p1, self.pts[0], self.pts[1])):
            R = 2 * np.dot(self.n_coll, -ray.r) * self.n_coll + ray.r
            ray_new = light_class.RayClass(p0=ray.p1, r=R, intensity=ray.intensity, wavelength=ray.wavelength, N=ray.N, ray_parent=ray, source_element=self, plot_color=ray.plot_color, is_active=True, is_visible=True)
            new_rays.append(ray_new)
        return new_rays


    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='darkgrey')  #, linewidth=5)
        graph.add_patch(poly)
        graph.plot([self.pts[0][X], self.pts[1][X]], [self.pts[0][Y], self.pts[1][Y]], color='cyan', linewidth=2, linestyle='solid', alpha=1, zorder=5)
        super().plot(graph)