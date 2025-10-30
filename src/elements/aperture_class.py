import numpy as np
from utils.varia import mm, X, Y
from utils import geometry, optics
from elements import element_class
import matplotlib.pyplot as plt

EPSILON = 1e-6

class ApertureClass(element_class.ElementClass):
    def __init__(self, p0=np.array([0, 0]), n0=np.array([-1, 0]), diameter_inner=10*mm, diameter_outer=20*mm, plot_color=(0.5,0.5,0.5,1), is_active=True, is_visible=True):
        self.diameter_inner = np.min([diameter_inner, diameter_outer])
        self.diameter_outer = np.max([diameter_inner, diameter_outer])
        self.plot_color = plot_color

        pts = optics.construct_aperture(p0=p0, n0=n0, Di=diameter_inner, Do=diameter_outer)
        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Aperture'

    def check_collision(self, ray):
        self.p_coll, t1_coll, t0_coll = None, None, None

        for i_face in [0, 2]:
            [p, t0, t1] = geometry.intersection_of_PP_line_with_PR_line(p00=self.pts[i_face], p01=self.pts[i_face + 1], p10=ray.p0, r1=ray.r)

            if (0 <= t0 <= 1) and t1 > EPSILON:  # t0<0: p lies before line segment | t0>1: p lies behind the line segment | t1<0: p lies in the backward direction of the ray
                if (t1_coll is None) or (t1 < t1_coll):  # There is not yet a closest point, or a new closest point is found
                    t1_coll = t1
                    t0_coll = t0
                    self.p_coll = p
                    self.n_coll = self.n[i_face]
                    self.i_coll = i_face

        return [self.p_coll, t1_coll, t0_coll]

    def propagate_ray(self, ray):
        new_rays = list()
        return new_rays

    def __str__(self):
        txt = f'{self.name} --> Element ID={self.ID}, p0={self.p0}, n0={self.n0}'
        return txt

    def plot(self, graph):
        if self.is_visible:
            graph.plot([self.pts[0][X], self.pts[1][X]], [self.pts[0][Y], self.pts[1][Y]], color=self.plot_color, linewidth=1, linestyle='solid', alpha=1, zorder=5)
            graph.plot([self.pts[2][X], self.pts[3][X]], [self.pts[2][Y], self.pts[3][Y]], color=self.plot_color, linewidth=1, linestyle='solid', alpha=1, zorder=5)
        super().plot(graph)

