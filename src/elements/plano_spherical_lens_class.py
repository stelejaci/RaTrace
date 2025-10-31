import numpy as np
from utils import varia
from utils.varia import mm, X, Y
from utils import optics
from utils.optics import N_glass
from utils import geometry
from elements import glass_element_class
import matplotlib.pyplot as plt
from utils.configuration_class import config


class PlanoSphericalLensClass(glass_element_class.GlassElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([1,0]), R=None, f=None, thickness=5*mm, diameter=10*mm, N=N_glass, blur_angle=0, nr_of_secondary_rays=1, plot_resolution=0.1):
        # Position p0 is the position at the center of the first surface
        # Front face radius R is positive for convex faces (bulging out, "fat")
        self.diameter  = diameter     # Diameter of the lens
        self.thickness = thickness  # Width of the lens along its optical axis
        self.plot_resolution = plot_resolution    # Inter-point distance for plotting the lens

        # If f is given, calculate R (also if they are given), else calculate f from R
        [self.f, self.R] = optics.derive_planosphericallens_radius_and_f(N=N, f=f, R=R)

        # Construct the outline of the lens in 0° orientation, rotate later
        [pts, self.p1, self.C, self.p_corners] = optics.construct_planosphericallens(p0=p0, R=self.R, T=self.thickness, D=self.diameter, resolution=self.plot_resolution)

        # Initialise the GlassElementsClass instance
        super().__init__(p0=p0, n0=n0,  N=N, pts=pts, blur_angle=blur_angle, nr_of_secondary_rays=nr_of_secondary_rays, is_active=True, is_visible=True)
        self.name = 'Plano-convex lens'

        # Rotate the lens to align with the normal vector n0
        p_rot = np.copy(self.p0)
        super().align_with_n0(p_rot=p_rot, n_ref=np.array([-1,0]))    # Rotate all general element-specific properties, including n0
        angle = geometry.angle_between_vectors(self.n0, np.array([-1,0]))  # The "normal" orientation of n0 is to the left ([-1,0]), so subtract 180° from that angle
        self.p1        = geometry.rotate_with_P_angle(pts=self.p1,        p_rot=p_rot, angle=angle)
        self.C         = geometry.rotate_with_P_angle(pts=self.C,         p_rot=p_rot, angle=angle)
        self.p_corners = geometry.rotate_with_P_angle(pts=self.p_corners, p_rot=p_rot, angle=angle)


    # collision checker coded explicitly because it is done analytically instead of based on the points and segments as is normally done
    def check_collision(self, ray):
        self.p_coll, t1_coll, t0 = None, None, None

        # TODO: Add check for collision with bottom straight pieces
        [pt_arc, t0_arc] = geometry.line_arc_intersections(C_arc=self.C, p_ends=self.p_corners[0:2], p_line=ray.p0, r_line=ray.r)
        [pt_line, t0_line, t1_line] = geometry.intersection_of_PR_line_with_PP_line(p00=ray.p0, r0=ray.r, p10=self.p_corners[2], p11=self.p_corners[3])

        if t0_arc is not None and t0_line is None:   case = 0    # Intersection with front surface
        elif t0_arc is None and t0_line is not None: case = 1    # Intersection with flat back surface
        elif t0_arc is None and t0_line is None:     case = 2    # No intersections
        elif t0_arc < t0_line:                       case = 0    # 2 intersections, first one is closer
        else:                                        case = 1    # 2 intersections, second one is closer

        if case == 0:    # Intersection with arc 1 is the one to select
            self.n_coll = geometry.normalize(np.sign(self.R) * (pt_arc - self.C))  # Normal at the front of the lens
            self.p_coll = pt_arc
            return [self.p_coll, t0_arc, t0_arc]
        elif case == 1:   # Intersection with the flat back surface is the one to select
            self.n_coll  = -self.n0  # Normal at the flat back of the lens
            self.p_coll  = pt_line
            return [self.p_coll, t0_line, t1_line]
        else:       # No intersections
            return [None, None, None]


    def plot(self, graph):
        if self.is_visible:
            poly_face = plt.Polygon(self.pts, closed=True, facecolor='cyan', edgecolor='none', alpha=varia.alpha_from_N(self.N), zorder=0)
            poly_edge = plt.Polygon(self.pts, closed=True, fill=False, edgecolor='blue', linewidth=1, alpha=0.2, zorder=0)
            graph.add_patch(poly_face)
            graph.add_patch(poly_edge)
            if config.getboolean('view', 'show_elements_properties'):
                graph.scatter(self.p1[0], self.p1[1], color='blue', s=25)
                graph.scatter(self.C[0], self.C[1], color='cyan', s=25)
                graph.scatter(self.p_corners[:,0], self.p_corners[:,1], color='blue', s=25)
                p_f0 = (self.p0+self.p1)/2 + self.n0*self.f
                p_f1 = (self.p0+self.p1)/2 - self.n0*self.f
                graph.scatter([p_f0[X], p_f1[X]], [p_f0[Y], p_f1[Y]], color='cyan', s=25)
                graph.text(self.p1[X], self.p1[Y], f'p1', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
                graph.text(p_f0[X], p_f0[Y], f'f={self.f:0.2f}mm', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
                graph.text(p_f1[X], p_f1[Y], f'f={self.f:0.2f}mm', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
        super().plot(graph)

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, n0={self.n0}, N={self.N}, R={self.R}, f={self.f}, thickness={self.thickness}, diameter={self.diameter}, blur angle={self.blur_angle}, number of secondary rays={self.nr_of_secondary_rays}, number of points={self.nr_of_pts}'
        return s
