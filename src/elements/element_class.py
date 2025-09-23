import numpy as np
import os, sys
sys.path.append(os.path.abspath('..'))

from utils import geometry
from utils.varia import mm, µm, nm, deg, X, Y, EPSILON
from utils.configuration_class import config

EPSILON = 1e-6




class ElementClass:
    nr_of_elements = 0

    def __init__(self, p0, n0, pts, is_active=True, is_visible=True):
        self.p0 = p0    # Some principal center point on the first surface
        self.n0 = geometry.normalize(n0)    # Principal normal orientation, connected to p0
        self.pts = pts  # The points describing te surface, could be for raytracing purposes, could be for plotting, or both. Points are described in CLOCKWISE order.
        self.is_active = is_active
        self.is_visible = is_visible

        self.ID = ElementClass.nr_of_elements
        ElementClass.nr_of_elements += 1

        self.i_coll = None  # The index of the closest face that a ray collided with
        self.n_coll = None  # The normal of that face a ray collided with
        self.p_coll = None  # The point on the surface that a ray collided with

        self.construct_element()
        self.calculate_normals_and_orthogonal_vectors()
        self.name = "Element"

    def construct_element(self):
        self.nr_of_pts = np.shape(self.pts)[0]
        self.nr_of_faces = self.nr_of_pts
        self.p_COG = np.mean(self.pts, axis=0)
        self.pts_BB = geometry.define_bounding_box(self.pts)
        self.pts = np.vstack([self.pts, self.pts[0, :]])  # Close the shape, even in the case of a 1D shape like a beam dump, diaphragm, ...

    def calculate_normals_and_orthogonal_vectors(self):
        self.r      = np.empty([0,2])    # The orientations of the faces
        self.n      = np.empty([0,2])    # The normals of the faces
        self.pts_n  = np.empty([0,2])    # The attachment points of the normals

        # Calculate the normals of the surfaces that connect each of the two points
        for i_pt in range(self.nr_of_faces):
            r = geometry.direction_between_two_points(self.pts[i_pt], self.pts[i_pt + 1])   # Orientation of the segment
            n = geometry.normal_from_orientation(r)                                          # Normal of the segment
            pt_n = (self.pts[i_pt] + self.pts[i_pt+1])/2                                    # The middle of the segment
            self.r      = np.append(self.r, [r], axis=0)
            self.n      = np.append(self.n, [n], axis=0)
            self.pts_n  = np.append(self.pts_n, [pt_n], axis=0)

    # The standard collision checker, useful for point-and-segment defines items. Not for analytically described items such as lenses etc.
    def check_collision(self, ray):
        self.p_coll, t1_coll, t0_coll = None, None, None

        for i_face in range(self.nr_of_faces):
            [p, t0, t1] = geometry.intersection_of_PP_line_with_PR_line(p00=self.pts[i_face], p01=self.pts[i_face+1], p10=ray.p0, r1=ray.r)

            if (0 <= t0 <= 1) and t1 > EPSILON:  # t0<0: p lies before line segment | t0>1: p lies behind the line segment | t1<0: p lies in the backward direction of the ray
                if (t1_coll is None) or (t1<t1_coll):   # There is not yet a closest point, or a new closest point is found
                    t1_coll = t1
                    t0_coll = t0
                    self.p_coll = p
                    self.n_coll = self.n[i_face]
                    self.i_coll = i_face

        return [self.p_coll, t1_coll, t0_coll]

    def align_with_n0(self, p_rot, n_ref):
        angle = geometry.angle_between_vectors(self.n0, n_ref)  # The "normal" orientation of n0 is to the left ([-1,0])
        self.n0 = geometry.rotate_with_P_angle(pts=self.n0, p_rot=p_rot, angle=-angle, is_orientation=True)  # Also n0 is rotated by the next rotate_by_angle command, so first anti-rotate this direction
        self.rotate_by_angle(angle=angle, p_rot=p_rot)

    def rotate_by_angle(self, angle, p_rot=np.array([0,0])):
        self.p0     = geometry.rotate_with_P_angle(pts=self.p0,     p_rot=p_rot, angle=angle)
        self.n0     = geometry.rotate_with_P_angle(pts=self.n0,     p_rot=p_rot, angle=angle, is_orientation=True)
        self.pts    = geometry.rotate_with_P_angle(pts=self.pts,    p_rot=p_rot, angle=angle)
        self.p_COG  = geometry.rotate_with_P_angle(pts=self.p_COG,  p_rot=p_rot, angle=angle)
        self.pts_n  = geometry.rotate_with_P_angle(pts=self.pts_n,  p_rot=p_rot, angle=angle)
        self.n      = geometry.rotate_with_P_angle(pts=self.n,      p_rot=p_rot, angle=angle, is_orientation=True)
        self.r      = geometry.rotate_with_P_angle(pts=self.r,      p_rot=p_rot, angle=angle)

    def move_by(self, dp):
        self.p0     = geometry.move_by(self.p0,     dp)
        self.pts    = geometry.move_by(self.pts,    dp)
        self.COG    = geometry.move_by(self.COG,    dp)
        self.pts_n  = geometry.move_by(self.pts_n,  dp)

    def move_to(self, p):
        dp = p - self.p0
        self.move_by(dp)

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, nr_of_points={self.nr_of_pts}'
        return s

    def plot(self, graph):
        print(f' --> Plotting element {self.ID + 1}/{ElementClass.nr_of_elements}: {self.name}')
        if config.getboolean('view', 'show_elements_properties'):
            graph.scatter(self.pts[:, X], self.pts[:, Y], color='black', s=2)
            graph.scatter(self.pts_n[:, X], self.pts_n[:, Y], color='orange', s=2)
            graph.scatter(self.p0[0], self.p0[1], color='black', s=25)
            graph.quiver(self.p0[X], self.p0[Y], 1 * self.n0[X], 1 * self.n0[Y], color='orange', width=0.003)
            graph.text(self.p0[X], self.p0[Y], f'p0', color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
            for i_normal in range(self.nr_of_faces):
                graph.quiver(self.pts_n[i_normal, X], self.pts_n[i_normal, Y], self.n[i_normal, X], self.n[i_normal, Y], color='orange', width=0.001, scale=50)




if __name__ == '__main__':
    p0 = np.array([0, 0])
    pts = np.array([[0,-1],[0,1],[1,1],[1,-1]])+10
    element = ElementClass(p0=p0, pts=pts)
    print(element)

    from light.light_class import RayClass
    ray = RayClass(p0=np.array([-1,0]))
    print(ray)

    [p_coll, t1, t0] = element.check_collision(ray)
    print(f'p_coll={p_coll}, t1={t1}, t0={t0}')