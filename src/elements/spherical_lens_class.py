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
    def __init__(self, p0=np.array([0,0]), n0=np.array([1,0]), R0=None, R1=None, f=None, thickness=2*mm, diameter=10*mm, N=N_glass, blur_angle=0, nr_of_secondary_rays=1, plot_resolution=0.1):
        # Position p0 is the position at the center of the first surface
        # Front face radius R0 is positive for convex faces (bulging out, "fat")
        # Back  face radius R1 is negative for convex faces (bulging out, "fat")
        self.diameter  = diameter     # Diameter of the lens
        self.thickness = thickness  # Width of the lens along its optical axis
        self.plot_resolution = plot_resolution    # Inter-point distance for plotting the lens

        # If f is given, calculate R0 and R1 (also if they are given), else calculate f from R0 and R1
        [self.f, self.R0, self.R1] = optics.derive_lens_radii_and_f(N=N, f=f, R0=R0, R1=R1, T=self.thickness)

        # Construct the outline of the lens in 0° orientation, rotate later
        [pts, self.p1, self.C0, self.C1, self.p_corners] = optics.construct_lens(p0=p0, R0=self.R0, R1=self.R1, T=self.thickness, D=self.diameter, resolution=self.plot_resolution)

        # Initialise the GlassElementsClass instance
        super().__init__(p0=p0, n0=n0,  N=N, pts=pts, blur_angle=blur_angle, nr_of_secondary_rays=nr_of_secondary_rays, is_active=True, is_visible=True)
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
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=varia.alpha_from_N(self.N), zorder=0)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            graph.scatter(self.p1[0], self.p1[1], color='blue', s=25)
            graph.scatter([self.C0[0],self.C1[0]], [self.C0[1],self.C1[1]], color='cyan', s=25)
            graph.scatter(self.p_corners[:,0], self.p_corners[:,1], color='blue', s=25)
            p_f0 = (self.p0+self.p1)/2 + self.n0*self.f
            p_f1 = (self.p0+self.p1)/2 - self.n0*self.f
            graph.scatter([p_f0[X], p_f1[X]], [p_f0[Y], p_f1[Y]], color='cyan', s=25)
            graph.text(self.p1[X], self.p1[Y], f'p1', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
            graph.text(p_f0[X], p_f0[Y], f'f={self.f:0.2f}mm', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
            graph.text(p_f1[X], p_f1[Y], f'f={self.f:0.2f}mm', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
        super().plot(graph)

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, n0={self.n0}, N={self.N}, R0={self.R0}, R1={self.R1}, f={self.f}, thickness={self.thickness}, diameter={self.diameter}, blur angle={self.blur_angle}, number of secondary rays={self.nr_of_secondary_rays}, number of points={self.nr_of_pts}'
        return s
