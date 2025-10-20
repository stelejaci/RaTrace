import numpy as np
from utils.varia import mm
from utils import geometry
from elements import diffuse_element_class
import matplotlib.pyplot as plt




class DiffusePlateClass(diffuse_element_class.DiffuseElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([-1,0]), length=10 * mm, thickness=1 * mm, Kd=0.0, Ks=0.0, alpha=1, nr_of_scattered_rays=1, n_light=None, is_active=True, is_visible=True):
        self.length = length
        self.thickness = thickness
        if n_light is None:
            n_light = n0

        # Derived parameters
        pts = geometry.construct_plate(p0=p0, n=n0, thickness=thickness, length=length)
        super().__init__(pts=pts, p0=p0, n0=n0, Kd=Kd, Ks=Ks, alpha=alpha, nr_of_scattered_rays=nr_of_scattered_rays, n_light=n_light, is_active=is_active, is_visible=is_visible)
        self.name = 'Diffuse plate'

    def __str__(self):
        txt = f'{self.name} --> Element ID={self.ID}, p0={self.p0}, n0={self.n0}, length={self.length}, thickness={self.thickness}, Kd={self.Kd}, Ks={self.Ks}, alpha={self.alpha}, number of scattered rays={self.nr_of_scattered_rays}, n_light={self.n_light}'
        return txt

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='grey')
        graph.add_patch(poly)
        super().plot(graph)
