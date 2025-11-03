import numpy as np
from utils import geometry, optics
from utils.varia import mm, X, Y
from elements import element_class
from light import light_class
from utils.configuration_class import config


class SemiTransparentMirror(element_class.ElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([-1,0]), length=10*mm, transmission=0.5, is_active=True, is_visible=True):
        self.length = length
        self.transmission = transmission

        self.r = geometry.orientation_from_normal(n0)
        pts = geometry.points_from_position_direction_length(p0=p0, r=self.r, L=self.length, symmetric=True)

        # Initialise the ElementClass instance
        super().__init__(p0=p0, n0=n0,  pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Semi transparent mirror'

    def check_collision(self, ray):
        self.p_coll, self.n_coll, t1, t0 = None, None, None, None

        for i_segment in range(2):  # Loop over the 2 "sides" of the lens, because one points away, one towards the ray
            [p_coll, t0_tmp, t1_tmp, n_coll] = geometry.intersection_of_PR_line_with_PPN_line(ray.p0, ray.r, self.pts[i_segment], self.pts[i_segment+1], self.n[i_segment])
            if (p_coll is not None)   and   (not np.isclose(t0_tmp,0)):
                self.p_coll, self.n_coll, t0, t1 = p_coll, n_coll, t0_tmp, t1_tmp
                break

        return [self.p_coll, t1, t0]

    def propagate_ray(self, ray):
        new_rays = list()

        # print(f'Propagating on SSM:')
        # print(ray)
        if (self.n_coll is not None)    and    (geometry.point_is_on_PP_line(ray.p1, self.pts[0], self.pts[1])):
            r_refl = optics.calculate_reflected_orientation(ray, self.n_coll)
            ray_reflected = light_class.RayClass(p0=ray.p1, r=r_refl, intensity=(1-self.transmission)*ray.intensity, wavelength=ray.wavelength, N=ray.N, ray_parent=ray, source_element=self, plot_color=ray.plot_color, is_active = True, is_visible=True)
            new_rays.append(ray_reflected)
            ray_transmitted = light_class.RayClass(p0=ray.p1, r=ray.r, intensity=self.transmission*ray.intensity, wavelength=ray.wavelength, N=ray.N, ray_parent=ray, source_element=self, plot_color=ray.plot_color, is_active = True, is_visible=True)
            new_rays.append(ray_transmitted)
        return new_rays

    def __str__(self):
        s = f'' # f'Glass parallel plate --> ID= {self.ID}, p0={self.p0}, n0={self.n0}, thickness={self.thickness}, length={self.length}, N={self.N}, nr_of_points={self.nr_of_pts}'
        return s

    def plot(self, graph):
        graph.plot([self.pts[0][X], self.pts[1][X]], [self.pts[0][Y], self.pts[1][Y]], color='blue', linewidth=3, linestyle='solid', alpha=1, zorder=5)
        if config.getboolean('view', 'show_elements_properties'):
            p_txt = self.p0 + 1.05 * self.r[0] * self.length / 2
            graph.text(p_txt[X], p_txt[Y], f'{self.name} {self.ID}', color='blue', horizontalalignment='center', verticalalignment='bottom', fontsize=10, rotation=90)
        super().plot(graph)
