import numpy as np
from utils.varia import mm
from utils import geometry
from elements import element_class
import matplotlib.pyplot as plt


class BlackPlateClass(element_class.ElementClass):
    def __init__(self, p0=np.array([0, 0]), n0=np.array([-1, 0]), length=10 * mm, thickness=1 * mm, plot_color=(0.5,0.5,0.5,1), is_active=True, is_visible=True):
        self.length     = length
        self.thickness  = thickness
        self.plot_color = plot_color

        pts = geometry.construct_plate(p0=p0, n=n0, thickness=thickness, length=length)
        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Black element'

    def propagate_ray(self, ray):
        new_rays = list()
        return new_rays

    def __str__(self):
        txt = f'{self.name} --> Element ID={self.ID}, p0={self.p0}, n0={self.n0}, length={self.length}, thickness={self.thickness}'
        return txt

    def plot(self, graph):
        if self.is_visible:
            poly = plt.Polygon(self.pts, closed=True, facecolor=self.plot_color, linewidth=5)
            graph.add_patch(poly)
        super().plot(graph)

