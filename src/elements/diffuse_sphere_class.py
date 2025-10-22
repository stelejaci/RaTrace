import numpy as np
from utils.varia import X, Y
from utils import geometry
import matplotlib.pyplot as plt
from utils.configuration_class import config
from gui import canvas_class
from elements import diffuse_element_class


class DiffuseSphereClass(diffuse_element_class):
    def __init__(self, p0, D, Kd, Ks, alpha, angular_spacing):
        self.p0 = p0
        self.D = D
        self.R = self.D/2
        self.generate_mesh(angular_spacing)

        super().__init__(Kd=Kd, Ks=Ks, alpha=alpha)
        self.name = 'Sphere'

    def generate_mesh(self, angular_spacing):
        self.pts = np.array([[0,0],[0,0]])

    def check_collision(self, ray):
        self.p_coll, t_coll, t0 = None, None, None
        # Method derived from: https://www.embibe.com/exams/intersection-between-circle-and-line/, but adjusted for an arbitrary circle in the plane, not through the origin
        # The frame is rotated 90°, such that X is horizontal, because of rounding errors, other wise c and D become too large, losing accuracy further on
        x0, y0, m, R = ray.p0[X], ray.p0[Y], ray.r[X]/ray.r[Y], self.R
        cst = x0-m*y0

        A = 1+m**2
        B = 2*m*(cst-self.p0[X]) - 2*self.p0[Y]
        C = self.p0[Y]**2+(cst-self.p0[X])**2-R**2
        D = B**2-4*A*C

        if D>0:
            y1 = (-B + np.sqrt(D)) / (2 * A)
            y2 = (-B - np.sqrt(D)) / (2 * A)
            x1 = np.sqrt(R**2 - (y1 - self.p0[Y])**2) + self.p0[X]  # x derived from the circle equation, is more acurate than from the line equation
            x2 = np.sqrt(R**2 - (y2 - self.p0[Y])**2) + self.p0[X]
            p11 = np.array([x1,y1])
            p12 = np.array([x2,y2])
            d1 = np.linalg.norm(ray.p0-p11)
            d2 = np.linalg.norm(ray.p0-p12)
            self.p_coll = p11  if  d1<d2  else  p12
            self.n_coll = geometry.normalize(self.p_coll-self.p0)
            t_coll = np.min([d1, d2])

        return [self.p_coll, t_coll, t0]

    def __str__(self):
        txt = 'Sphere'
        return txt

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor=canvas_functions.PLOT_COLOR_BLACK)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            canvas_class.plot_normals(graph, self)




