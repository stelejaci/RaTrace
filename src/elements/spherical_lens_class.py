import numpy as np
from utils import varia
from utils.varia import mm, X, Y
from utils import optics
from utils.optics import N_glass
from utils import geometry
from elements import glass_element_class
import matplotlib.pyplot as plt
from utils.configuration_class import config


class SphericalLensClass(glass_element_class.GlassElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([1,0]), R0=None, R1=None, f=None, thickness=2*mm, diameter=10*mm, N=N_glass, plot_resolution=1, generate_reflections=False):
        # Position p0 is the position at the center of the first surface
        # Front face radius R0 is positive for convex faces (bulging out, "fat")
        # Back  face radius R1 is negative for convex faces (bulging out, "fat")
        self.diameter  = diameter     # Diameter of the lens
        self.thickness = thickness  # Width of the lens along its optical axis
        self.plot_resolution = plot_resolution    # Inter-point distance for plotting the lens

        # If f is given, calculate R0 and R1 (also if they are given), else calculate f from R0 and R1
        # H0 and H1 are the first and second principal points, i.e. the distance to resp. p0 and p1
        [self.f, self.R0, self.R1, self.H0, self.H1] = optics.derive_lens_properties(N=N, f=f, R0=R0, R1=R1, T=self.thickness)

        # Construct the outline of the lens in 0° orientation, rotate later
        [pts, self.p1, self.C0, self.C1, self.p_corners] = optics.construct_lens(p0=p0, R0=self.R0, R1=self.R1, T=self.thickness, D=self.diameter, resolution=self.plot_resolution)

        # Initialise the GlassElementsClass instance
        super().__init__(p0=p0, n0=n0,  N=N, pts=pts, is_active=True, is_visible=True, generate_reflections=generate_reflections)
        self.name = 'Spherical lens'

        # Rotate the lens to align with the normal vector n0
        p_rot = np.copy(self.p0)
        super().align_with_n0(p_rot=p_rot, n_ref=np.array([-1,0]))    # Rotate all general element-specific properties, including n0
        angle = geometry.angle_between_vectors(self.n0, np.array([-1,0]))  # The "normal" orientation of n0 is to the left ([-1,0]), so subtract 180° from that angle
        self.p1        = geometry.rotate_with_P_angle(pts=self.p1,        p_rot=p_rot, angle=angle)
        self.C0        = geometry.rotate_with_P_angle(pts=self.C0,        p_rot=p_rot, angle=angle)
        self.C1        = geometry.rotate_with_P_angle(pts=self.C1,        p_rot=p_rot, angle=angle)
        self.p_corners = geometry.rotate_with_P_angle(pts=self.p_corners, p_rot=p_rot, angle=angle)


    # collision checker coded explicitly because it is done analytically instead of based on the points as is normally done
    def check_collision(self, ray):
        self.p_coll, t1_coll, t0 = None, None, None

        # TODO: Add check for collision with bottom straight pieces
        [pt0, t0] = geometry.line_arc_intersections(C_arc=self.C0, p_ends=self.p_corners[0:2], p_line=ray.p0, r_line=ray.r)
        [pt1, t1] = geometry.line_arc_intersections(C_arc=self.C1, p_ends=self.p_corners[2: ], p_line=ray.p0, r_line=ray.r)

        if t0 is not None and t1 is None:   case = 0    # Intersection with only arc 1
        elif t0 is None and t1 is not None: case = 1    # Intersection with only arc 2
        elif t0 is None and t1 is None:     case = 2    # No intersections
        elif t0 < t1:                       case = 0    # 2 intersections, first is closer
        else:                               case = 1    # 2 intersections, second is closer

        if case == 0:    # Intersection with arc 1 is the one to select
            self.n_coll = geometry.normalize(np.sign(self.R0) * (pt0 - self.C0))  # Normal at the front of the lens
            self.p_coll = pt0
            return [self.p_coll, t0, t0]
        elif case == 1:   # Intersection with arc 2 is the one to select
            self.n_coll  = -geometry.normalize(np.sign(self.R1) * (pt1 - self.C1))  # Normal at the back of the lens
            self.p_coll  = pt1
            return [self.p_coll, t1, t1]
        else:       # No intersections
            return [None, None, None]


    def plot(self, graph):
        if self.is_visible:
            poly_face = plt.Polygon(self.pts, closed=True, facecolor='cyan', edgecolor='none', alpha=varia.alpha_from_N(self.N), zorder=0)
            poly_edge = plt.Polygon(self.pts, closed=True, fill=False, edgecolor='blue', linewidth=1, alpha=0.2, zorder=0)
            graph.add_patch(poly_face)
            graph.add_patch(poly_edge)
            if config.getboolean('view', 'show_elements_properties'):
                graph.scatter(self.p1[0], self.p1[1], color='black', s=25)
                graph.scatter([self.C0[0],self.C1[0]], [self.C0[1],self.C1[1]], color='black', s=20)
                graph.scatter(self.p_corners[:,0], self.p_corners[:,1], color='black', s=20)
                r = geometry.orientation_from_normal(self.n0)
                p_H0, p_H1 = self.p0 - self.n0*self.H0, self.p1 - self.n0*self.H1
                p_H0t, p_H0b = p_H0 + r*self.diameter/2,  p_H0 - r*self.diameter/2
                p_H1t, p_H1b = p_H1 + r*self.diameter/2,  p_H1 - r*self.diameter/2
                p_f0, p_f1 = p_H0 + self.n0*self.f, p_H1 - self.n0*self.f
                graph.plot([p_H0t[X], p_H0b[X]], [p_H0t[Y], p_H0b[Y]], color='black', linewidth=1, linestyle='dashed', alpha=1, zorder=5)
                graph.plot([p_H1t[X], p_H1b[X]], [p_H1t[Y], p_H1b[Y]], color='black', linewidth=1, linestyle='dashed', alpha=1, zorder=5)
                graph.plot([p_f0[X], (3*p_H0t[X]+p_H0b[X])/4], [p_f0[Y], (3*p_H0t[Y]+p_H0b[Y])/4], color='black', linewidth=1, linestyle='--', alpha=0.2, zorder=5)
                graph.plot([p_f1[X], (3*p_H1t[X]+p_H1b[X])/4], [p_f1[Y], (3*p_H1t[Y]+p_H1b[Y])/4], color='black', linewidth=1, linestyle='--', alpha=0.2, zorder=5)
                graph.plot([self.C0[X], self.p_corners[1,X]], [self.C0[Y], self.p_corners[1,Y]], color='black', linewidth=1, linestyle=':', alpha=0.2, zorder=5)
                graph.plot([self.C1[X], self.p_corners[3,X]], [self.C1[Y], self.p_corners[3,Y]], color='black', linewidth=1, linestyle=':', alpha=0.2, zorder=5)
                graph.scatter([p_f0[X], p_f1[X]], [p_f0[Y], p_f1[Y]], color='black', s=20)
                graph.text(self.p1[X], self.p1[Y], f'p1', color='black', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'), ha='center', va='bottom', fontsize=8)
                graph.text(p_f0[X], p_f0[Y], f'f={self.f:0.2f}mm', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
                graph.text(p_f1[X], p_f1[Y], f'f={self.f:0.2f}mm', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
                graph.text(p_H0t[X], p_H0t[Y], f'H0', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
                graph.text(p_H1t[X], p_H1t[Y], f'H1', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
                graph.text(self.C0[X], self.C0[Y], f'C0', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
                graph.text(self.C1[X], self.C1[Y], f'C1', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
            super().plot(graph)

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, n0={self.n0}, N={self.N}, R0={self.R0}, R1={self.R1}, f={self.f}, thickness={self.thickness}, diameter={self.diameter}, number of points={self.nr_of_pts}'
        return s
