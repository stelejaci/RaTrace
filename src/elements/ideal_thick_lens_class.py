import numpy as np
import os, sys
sys.path.append(os.path.abspath('..'))

from utils import varia
from utils.varia import mm, µm, nm, deg, X, Y
from utils import optics
from utils.optics import N_air, N_glass
from utils import geometry
from light import light_class
from elements import element_class
from elements import glass_element_class
import matplotlib.pyplot as plt
from utils.configuration_class import config
from gui import canvas_class




# ToDo rewrite to new structure, inherit from glass_element_class too
class ThickLensClass(element_class.ElementClass):
    def __init__(self, p0, f, D, N, blur_angle, nr_of_secondary_rays, p1=None, n=None, length=None):
        self.p0 = p0    # Mid-point on the front lens surface
        self.p1 = p1    # Mid-point on the back lens surface
        self.length = length
        self.n = n     # Towards mid-FOV
        self.f = f
        self.D = D
        self.N = N
        self.blur_angle = blur_angle
        self.nr_of_secondary_rays = nr_of_secondary_rays

        # The secondary principal point can be given explicitly, then the length and principal normal can be calculated from that.
        # OR the length of the lens is given, as well as the principal normal, then p1 is calculated
        if self.p1 is not None:
            self.length = geometry.distance_between_2_points(self.p0, self.p1)
            self.n = -geometry.distance_between_2_points(self.p0, self.p1)
        elif self.length is not None  and  self.n is not None:
            self.n = geometry.normalize(n)     # Towards mid-FOV
            self.p1 = self.p0 - self.length * self.n
        else:
            print('Lens definition not well defined, check parameter p1, OR parameters length and normal')

        # We go from an abstract description of the lens (p0,p1,n,r) to an explicit description with 4 points and normals, and redefine some parameters
        self.r   = geometry.perpendicular_direction_to_vector(-self.n)
        self.pts = geometry.create_points_from_position_direction_length(p0=self.p0, r=self.r, L=self.D, sort_left_to_right=False, symmetric=True)
        self.pts = np.append(self.pts, [self.pts[1]-self.n*self.length], axis=0)
        self.pts = np.append(self.pts, [self.pts[0]-self.n*self.length], axis=0)

        # Initialise the ElementsClass instance
        super().__init__(p0=self.p0, pts=self.pts)
        self.name = 'Thick lens'

        # For generating randomly scattered angles, one needs the Gaussian-based probability cumulative sum for sampling
        # self.generate_ray_displacement_matrices(blur_angle=blur_angle, nr_of_secondary_rays=self.nr_of_secondary_rays)
        (self.scattering_angles, self.scattering_pcs) = varia.generate_gaussian_pcs(FWHM=self.blur_angle, verbose=False)

    def propagate_ray(self, ray):
        new_rays = list()

        # The ray hits the sides of the lens
        if self.i_coll in [1,3]:
            return new_rays

        # Calculate the crossing of the ray with the lens' optical axis
        # So, calculate the intersection between the pt-rico-line of the ray and the pt-rico line of the lens' optical axis
        [p_on_optical_axis, t1, t2] = geometry.intersection_of_PR_line_with_PR_line(self.p0, self.n_coll, ray.p0, ray.r)

        if np.dot(ray.r, self.n_coll) < 0:
            image_distance = optics.calculate_image_distance(t1, self.f)
            if   t1>0 and image_distance>0:
                direction = geometry.direction_between_two_points(ray.p1, self.p0-image_distance*self.n_coll)
            elif t1>0 and image_distance<0:
                direction = geometry.direction_between_two_points(self.p0-image_distance*self.n_coll, ray.p1)
            elif t1<0 and image_distance>0:
                direction = geometry.direction_between_two_points(ray.p1, self.p0-image_distance*self.n_coll)
            elif t1<0 and image_distance<0:
                direction = geometry.direction_between_two_points(self.p0-image_distance*self.n_coll, ray.p1)    # This is not reached. Is the "-" correct?
            elif np.isnan(t1):  # Ray on and parallel to optical axis --> propagate undisturbed
                direction = ray.r
            elif t1==0:         # Ray goes through the lens center --> propagate undisturbed
                direction = ray.r
            else:
                # print('How did this happen, are optical elements overlapping maybe???')
                direction = np.array([1,0])
        elif np.dot(ray.r, self.n_coll) > 0:
            image_distance = optics.calculate_image_distance(-t1, self.f)
            if   t1<0 and image_distance>0:
                direction = varia.direction_between_two_points(ray.p1, self.p0+image_distance*self.n_coll)
            elif t1<0 and image_distance<0:
                try:
                    direction = geometry.direction_between_two_points(self.p0+image_distance*self.n_coll, ray.p1)
                except:
                    pass
            else:
                # print('How did this happen, are optical elements overlapping maybe???')
                direction = np.array([1,0])

        # Create one new (passive, non-active) ray that propagates between the principle planes
        if self.length>0:
            propagating_ray = light_class.RayClass(p0=ray.p1, r=-self.n_coll, intensity=ray.intensity, wavelength=ray.wavelength, ray_parent=ray, N=self.N, source_element=self, is_active=False, plot_color=ray.plot_color)
            propagating_ray.p1 = ray.p1 - self.length * self.n_coll
            propagating_ray.length = self.length
            new_rays.append(propagating_ray)

        # Generate scattered rays, following a random Gaussian distribution
        scattered_rays = generate_scattered_rays(p0=propagating_ray.p1, direction=direction, intensity=ray.intensity, wavelength=ray.wavelength,
                                                      ray_parent=propagating_ray, source_element=self, plot_color=ray.plot_color,
                                                      scattering_angles=self.scattering_angles, probability_cumsum=self.scattering_pcs,
                                                      randomize=True, blur_angle=self.blur_angle, nr_of_secondary_rays=self.nr_of_secondary_rays)
        new_rays = new_rays + scattered_rays

        return new_rays

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, n0={self.n0}, N={self.N}, R0={self.R0}, R1={self.R1}, f={self.f}, thickness={self.thickness}, diameter={self.diameter}, blur angle={self.blur_angle}, number of secondary rays={self.nr_of_secondary_rays}, number of points={self.nr_of_pts}'
        return s

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=misc.alpha_from_N(self.N), zorder=5)
        graph.add_patch(poly)
        graph.plot(self.pts[[0, 1], 0], self.pts[[0, 1], 1], c='c', zorder=6)
        graph.plot(self.pts[[2, 3], 0], self.pts[[2, 3], 1], c='c', zorder=6)
        graph.plot(self.pts[[1, 2], 0], self.pts[[1, 2], 1], c=canvas_functions.PLOT_COLOR_BLACK, linewidth=3, zorder=6)
        graph.plot(self.pts[[3, 0], 0], self.pts[[3, 0], 1], c=canvas_functions.PLOT_COLOR_BLACK, linewidth=3, zorder=6)
        if config.getboolean('view', 'show_elements_properties'):
            canvas_class.plot_normals(graph, self)





# def generate_scattered_rays(p0, direction, intensity, wavelength, ray_parent, source_element, plot_color, scattering_angles, probability_cumsum, randomize=True, blur_angle=0*deg, nr_of_secondary_rays=1):
#     rays = list()
#
#     # If no blur angle is given, propagate with only one ray in the nominal, refracted direction
#     if blur_angle == 0:
#         scattered_ray = light_functions.RayClass(p0=p0, r=direction, intensity=intensity, wavelength=wavelength, ray_parent=ray_parent, N=N_air, source_element=source_element, is_active=True, plot_color=plot_color)
#         rays.append(scattered_ray)
#     # If the lens has some blurring, generate a bunch of scattering directions and forthcoming scattered rays
#     else:
#         # Generate the scattering angles
#         angles_sampling = misc.generate_samples_from_pcs(x=scattering_angles, probability_cumsum=probability_cumsum, nr_of_samples=nr_of_secondary_rays, randomize=randomize, verbose=False)
#
#         # Generate the scattered rays itself
#         for angle in angles_sampling:
#             scattered_ray_direction = misc.rotate_direction_over_angle(direction, angle)
#             # M = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
#             # scattered_ray_direction = np.matmul(M, direction)
#             # scattered_ray_direction = misc.normalize(scattered_ray_direction)
#
#             scattered_ray_intensity = intensity / nr_of_secondary_rays
#             scattered_ray = light_functions.RayClass(p0=p0, r=scattered_ray_direction, intensity=scattered_ray_intensity, wavelength=wavelength, ray_parent=ray_parent, N=N_air, source_element=source_element, is_active=True, plot_color=plot_color)
#
#             rays.append(scattered_ray)
#
#     return rays


if __name__ == "__main__":
    pass