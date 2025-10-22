import numpy as np
from utils import varia
from utils.varia import mm, deg, X, Y
from utils import optics
from utils.optics import N_glass
from utils import geometry
from light import light_class
from elements import glass_element_class
from utils.configuration_class import config

RED, WHITE = '\033[31m', '\033[0m'


class IdealThinLensClass(glass_element_class.GlassElementClass):
    def __init__(self, p0=np.array([10,0]), n0=np.array([-1,0]), f=100*mm, diameter=10*mm, N=N_glass, blur_angle=0, nr_of_secondary_rays=1):
        self.diameter  = diameter     # Diameter of the lens
        self.f = f

        self.r = geometry.orientation_from_normal(n0)
        pts = geometry.points_from_position_direction_length(p0=p0, r=self.r, L=self.diameter, symmetric=True)

        # Initialise the GlassElementsClass instance
        super().__init__(p0=p0, n0=n0,  N=N, pts=pts, blur_angle=blur_angle, nr_of_secondary_rays=nr_of_secondary_rays, is_active=True, is_visible=True)
        self.name = 'Ideal thin lens'

    def check_collision(self, ray):
        self.p_coll, self.n_coll, t1, t0 = None, None, None, None

        for i_segment in range(2):  # Loop over the 2 "sides" of the lens, because one points away, one towards the ray
            [p_coll, t0_tmp, t1_tmp, n_coll] = geometry.intersection_of_PR_line_with_PPN_line(ray.p0, ray.r, self.pts[i_segment], self.pts[i_segment+1], self.n[i_segment])
            if (p_coll is not None)   and   (not np.isclose(t0_tmp,0)):
                self.p_coll, self.n_coll, t0, t1 = p_coll, n_coll, t0_tmp, t1_tmp
                break

        return [self.p_coll, t0, t1]

    def propagate_ray(self, ray):
        new_rays = list()

        # Refract a ray on the ideal lens
        ro = optics.refract_ray_on_ideal_lens(p0=self.p0, n0=self.n0, f=self.f, p=ray.p0, r=ray.r, p_coll=self.p_coll)

        # Create the outgoing ray
        if np.isclose(ray.length, 0):
            print(f'{RED}The ray has zero length, possibly colliding within the ideal lens itself{WHITE}')
        else:
            new_ray = light_class.RayClass(p0=ray.p1, r=ro, intensity=ray.intensity, wavelength=ray.wavelength, ray_parent=ray, N=ray.N, source_element=self, is_active=ray.is_active, is_visible=ray.is_visible, is_virtual=ray.is_virtual, plot_color=ray.plot_color)
            new_rays.append(new_ray)

        return new_rays

    def plot(self, graph):
        graph.plot([self.pts[0][X], self.pts[1][X]], [self.pts[0][Y], self.pts[1][Y]], color='cyan', linewidth=3, linestyle='solid', alpha=1, zorder=5)
        varia.plot_arrow_end_at_P(graph, P=self.pts[0], r= self.r[0], s=np.sign(self.f)*self.diameter/50, angle=30*deg, col='cyan')
        varia.plot_arrow_end_at_P(graph, P=self.pts[1], r=-self.r[0], s=np.sign(self.f)*self.diameter/50, angle=30*deg, col='cyan')
        if config.getboolean('view', 'show_elements_properties'):
            p_f0 = self.p0 + self.n0 * self.f
            p_f1 = self.p0 - self.n0 * self.f
            graph.plot([p_f0[X], p_f1[X]], [p_f0[Y], p_f1[Y]], color='cyan', linewidth=1, linestyle='solid', alpha=1, zorder=5)
            graph.scatter([p_f0[X], p_f1[X]], [p_f0[Y], p_f1[Y]], color='cyan', s=25)
            graph.text(p_f0[X], p_f0[Y], f'f={self.f:0.2f}mm', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
            graph.text(p_f1[X], p_f1[Y], f'f={self.f:0.2f}mm', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=8)
            # p_txt = self.p0 + 1.05 * self.r[0] * self.diameter/2
            # graph.text(p_txt[X], p_txt[Y], f'{self.name} {self.ID}', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=10, rotation=90)
        super().plot(graph)

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, n0={self.n0}, N={self.N}' # , R0={self.R0}, R1={self.R1}, f={self.f}, thickness={self.thickness}, diameter={self.diameter}, blur angle={self.blur_angle}, number of secondary rays={self.nr_of_secondary_rays}, number of points={self.nr_of_pts}'
        return s
