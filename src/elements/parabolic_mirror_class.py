import numpy as np
import matplotlib.pyplot as plt
from utils.varia import mm, X, Y
from light import light_class
from elements import element_class
from utils import geometry, optics
from utils.configuration_class import config


class ParabolicMirrorClass(element_class.ElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([1,0]), f=100*mm, diameter=100*mm, thickness=10 * mm, plot_color=(0,0,0.5,1), plot_resolution=0.1, is_active=True, is_visible=True):
        self.f          = f
        self.diameter   = diameter
        self.thickness  = thickness
        self.plot_resolution = plot_resolution    # Inter-point distance for plotting the lens
        self.plot_color = plot_color

        [pts, self.p1, self.f0, self.p_corners] = optics.construct_parabolic_mirror(p0=p0, f=self.f, D=self.diameter, T=thickness, resolution=plot_resolution)
        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Parabolic mirror'

        # Rotate the mirror to align with the normal vector n0
        p_rot = np.copy(self.p0)
        super().align_with_n0(p_rot=p_rot, n_ref=np.array([-1,0]))    # Rotate all general element-specific properties, including n0
        angle = geometry.angle_between_vectors(self.n0, np.array([-1,0]))  # The "normal" orientation of n0 is to the left ([-1,0]), so subtract 180° from that angle
        self.p1        = geometry.rotate_with_P_angle(pts=self.p1,        p_rot=p_rot, angle=angle)
        self.f0        = geometry.rotate_with_P_angle(pts=self.f0,        p_rot=p_rot, angle=angle)
        self.p_corners = geometry.rotate_with_P_angle(pts=self.p_corners, p_rot=p_rot, angle=angle)

    # collision checker coded explicitly because it is done analytically instead of based on the points as is normally done
    def check_collision(self, ray):
        self.p_coll, t1_coll, t0 = None, None, None

        # Check collision against the straight sides of the mirror
        for i_segment in range(3):
            [self.p_coll, t0, t1, n_coll] = geometry.intersection_of_PR_line_with_PPN_line(P=ray.p0, R=ray.r, P0=self.pts[-4+i_segment], P1=self.pts[-3+i_segment], N=self.n[-3+i_segment])

            if self.p_coll is not None: # Stop looking for intersections if the ray intersects with one of the segments
                self.n_coll = None      # Prevents the creation of a new ray ni the ray propagator
                break

        # No intersection with the outer segments, continue looking for an intersection with the parabola itself
        if self.p_coll is None:
            [self.p_coll, t0, t1, self.n_coll] = geometry.line_parabola_intersections(f=self.f, p0=self.p0, n0=self.n0, D=self.diameter, p_line=ray.p0, r_line=ray.r)

        return [self.p_coll, t0, t1]

    def propagate_ray(self, ray):
        ray_new = []
        if self.n_coll is not None:
            R = 2 * np.dot(self.n_coll, -ray.r) * self.n_coll + ray.r
            ray_new = light_class.RayClass(p0=ray.p1, r=R, intensity=ray.intensity, wavelength=ray.wavelength, N=ray.N, ray_parent=ray, source_element=self, plot_color=ray.plot_color, is_active = True, is_visible=True)
        return ray_new

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor=self.plot_color, linewidth=5)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            graph.scatter(self.p1[0], self.p1[1], color='black', s=25)
            graph.scatter(self.f0[0], self.f0[1], color='green', s=25)
            graph.scatter(self.p_corners[:,0], self.p_corners[:,1], color='black', s=25)
            graph.scatter(self.f0[X], self.f0[Y], color='red', s=25)
            graph.text(self.p1[X], self.p1[Y], f'p1', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
            graph.text(self.f0[X], self.f0[Y], f'f={self.f:0.2f}mm', color='red', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
        super().plot(graph)
