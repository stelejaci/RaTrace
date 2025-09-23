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
from utils.configuration_class import config
from gui import canvas_class
import matplotlib.pyplot as plt



class GlassElementClass(element_class.ElementClass):
    def __init__(self, p0, n0, pts, N=N_glass, blur_angle=0, nr_of_secondary_rays=1, is_active=True, is_visible=True):
        self.N = N
        self.blur_angle = blur_angle
        self.nr_of_secondary_rays = nr_of_secondary_rays

        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Glass element'

        # For generating randomly scattered angles, one needs the Gaussian-based probability cumulative sum for sampling
        # self.generate_ray_displacement_matrices(blur_angle=blur_angle, nr_of_secondary_rays=self.nr_of_secondary_rays)
        (scattering_angles, scattering_pcs) = varia.generate_gaussian_pcs(FWHM=self.blur_angle, verbose=False)
        self.scattering_angles_sampling = varia.generate_samples_from_pcs(x=scattering_angles, probability_cumsum=scattering_pcs, nr_of_samples=nr_of_secondary_rays, randomize=False, verbose=False)

    def propagate_ray(self, ray):
        # Select the proper refractive indices
        Ni = ray.N  # The incoming N is always that of the ray
        if np.dot(ray.r, self.n_coll)<0:
            No = self.N         # The element normal points towards the ray, meaning the ray enters the element, thus No=N
        else:
            No = N_air          # The element normal points away from the ray, meaning the ray leaves the element, thus No=N_air

        # Calculate dispersion
        Ni = optics.calculate_refraction_index(ray.N, ray.wavelength)
        No = optics.calculate_refraction_index(No, ray.wavelength)

        ray_new = optics.refract_ray(ray, self.n_coll, Ni=Ni, No=No)
        ray_new.source_element = self

        # Virtual rays should not scatter in/on glass
        if ray.is_virtual:
            return ray_new

        # Refract the ray at the collision surface, but only at the exiting side of a glass element
        if No>Ni:   # Entering glass
            return ray_new
        else:       # Leaving glass
            rays_scattered = self.generate_scattered_rays(ray_new)
            return rays_scattered

    def generate_scattered_rays(self, ray):
        rays_scattered = []
        for angle in self.scattering_angles_sampling:
            scattered_ray_direction = geometry.rotate_direction_over_angle(ray.r, angle)
            # M = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
            # scattered_ray_direction = np.matmul(M, direction)
            # scattered_ray_direction = misc.normalize(scattered_ray_direction)
            scattered_ray_intensity = ray.intensity / self.nr_of_secondary_rays
            scattered_ray = light_class.RayClass(p0=ray.p0, r=scattered_ray_direction, intensity=scattered_ray_intensity, wavelength=ray.wavelength, ray_parent=ray, N=ray.N, source_element=ray.source_element, is_active=True, is_visible=True, plot_color=ray.plot_color)
            rays_scattered.append(scattered_ray)
        return rays_scattered

    def __str__(self):
        s = f'{self.name} --> Element ID= {self.ID}, p0={self.p0}, n0={self.n0}, N={self.N}, blur angle={self.blur_angle}, number of secondary rays={self.nr_of_secondary_rays}, number of points={self.nr_of_pts}'
        return s

    def plot(self, graph):
        if config.getboolean('view', 'show_elements_properties'):
            p_txt = np.array([(np.min(self.pts[:,X]) + np.max(self.pts[:,X]))/2, np.max(self.pts[:][Y])])
            # p_txt = (self.p0  + self.p1)/2
            graph.text(p_txt[X], p_txt[Y], f'{self.name} {self.ID}', color='cyan', horizontalalignment='center', verticalalignment='bottom', fontsize=10, rotation=0)
        super().plot(graph)

class FresnellPrismClass(GlassElementClass):
    def __init__(self, p0, pitch, angle, thickness, length, N, is_active=True):
        self.nr_of_lenslets = int(np.ceil(length/pitch))
        self.pitch = pitch
        self.thickness = thickness
        self.angle = angle
        self.length = length
        self.N = N
        
        pts = np.empty((0,2))
        # pts = np.append(pts, [p0[X]+length/2], axis=0)
        pts = np.append(pts, [np.array([p0[X]+length/2, p0[Y]])], axis=0)
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y]])], axis=0)
        for i in range(self.nr_of_lenslets):
            p1 = np.array([self.pitch*(i+0.1),  self.thickness])
            p2 = np.array([self.pitch*(i+1.0),  self.thickness-self.pitch*(1-0.1)*np.tan(self.angle)])
            pts = np.append(pts, [pts[1] + p1], axis=0)
            pts = np.append(pts, [pts[1] + p2], axis=0)
        # pts = np.append(pts, [np.array([p2[X],0])], axis=0)
        
        super().__init__(p0, pts, N=N, is_active=is_active)
        self.name = 'fresnell prism'

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=geometry.alpha_from_N(self.N), zorder=5)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            canvas_class.plot_normals(graph, self)


class GlassParallelPlate(GlassElementClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([-1,0]), thickness=1*mm, length=10*mm, N=N_glass, is_active=True, is_visible=True):
        self.thickness = thickness
        self.length = length

        pts = geometry.construct_plate(p0=p0, n=n0, thickness=thickness, length=length)
        super().__init__(p0=p0, n0=n0, pts=pts, N=N, is_active=is_active)
        self.name = 'glass parallel plate'

    def __str__(self):
        s = f'Glass parallel plate --> ID= {self.ID}, p0={self.p0}, n0={self.n0}, thickness={self.thickness}, length={self.length}, N={self.N}, nr_of_points={self.nr_of_pts}'
        return s

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=varia.alpha_from_N(self.N), zorder=5)
        graph.add_patch(poly)
        super().plot(graph)


class GlassBiprism(GlassElementClass):
    def __init__(self, p0, thickness, angle_apex, length, N, is_active=True):
        self.thickness = thickness
        self.angle_apes = angle_apex
        self.length = length
        self.N = N

        pts = np.empty((0, 2))
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y]])],                                                    axis=0)
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y] + self.thickness])],                                   axis=0)
        pts = np.append(pts, [np.array([p0[X],          p0[Y] + self.thickness + length/(2*np.tan(angle_apex/2))])], axis=0)
        pts = np.append(pts, [np.array([p0[X]+length/2, p0[Y] + self.thickness])],                                   axis=0)
        pts = np.append(pts, [np.array([p0[X]+length/2, p0[Y]])],                                                    axis=0)

        super().__init__(p0, pts, N=N, is_active=is_active)
        self.name = 'biprism'

    def __str__(self):
        txt = 'Biprism'
        return txt

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=geometry.alpha_from_N(self.N), zorder=5)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            canvas_class.plot_normals(graph, self)


class GlassPrism(GlassElementClass):
    def __init__(self, p0, angle, length, N, is_active=True):
        self.angle = angle
        self.length = length
        self.N = N

        pts = np.empty((0, 2))
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y]])],                        axis=0)
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y] + length*np.sin(angle)])], axis=0)
        pts = np.append(pts, [np.array([p0[X]+length/2, p0[Y]])],                        axis=0)

        super().__init__(p0, pts, N=N, is_active=is_active)
        self.name = 'prism'

    def __str__(self):
        txt = 'Prism'
        return txt

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=geometry.alpha_from_N(self.N), zorder=5)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            canvas_class.plot_normals(graph, self)


class MicroLensArray(GlassElementClass):
    def __init__(self, p0, thickness, pitch, peaktopeak, shape, length, N, is_active=True):
        self.thickness = thickness
        self.pitch = pitch
        self.peaktopeak = peaktopeak
        self.length = length
        self.N = N
        self.is_active = is_active
        self.shape = shape
        self.nr_of_samples_per_pitch = 1000
        self.no_lens_zone = 0*mm

        pts = np.empty((0, 2))
        pts = np.append(pts, [np.array([p0[X]+length/2, p0[Y] + self.thickness])], axis=0)
        pts = np.append(pts, [np.array([p0[X]+length/2, p0[Y]])],                  axis=0)
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y]])],                  axis=0)
        pts = np.append(pts, [np.array([p0[X]-length/2, p0[Y] + self.thickness])], axis=0)

        if self.shape == 'sinusoidal':
            nr_of_samples = round( self.length / self.pitch * self.nr_of_samples_per_pitch )
            dx = np.linspace(self.no_lens_zone, self.length-self.no_lens_zone, num=nr_of_samples, endpoint=True)
            theta = dx / self.pitch * 2*np.pi
            dy = self.peaktopeak/2 * np.sin(theta)
            for ind in range(len(dx)):
                pts = np.append(pts, [pts[3] + np.array([dx[ind],dy[ind]])],  axis=0)
        elif self.shape == 'circular':
            # See Matlab code: D:\Documents\SIM - 2D scanner simulation\aa Matlab proof-of-concepts & testing\ah Microlens array\run_construct_MLA_surface_amplitude_given.m
            nr_of_lenslets = int(np.floor(self.length/self.pitch))
            A = self.peaktopeak # Height of the lenslets
            W = self.pitch  # Width of each lenslet
            yM = 1/(2*A)*(W/2)**2 - A/2  # Radius of a lenslet minus the amplitude A
            R = yM + A
            dx = W/self.nr_of_samples_per_pitch
            x1 = np.linspace(dx, W, self.nr_of_samples_per_pitch)  # X coordinates over the lenslet, with the origin at the start of the lenslet. Skip the first position to avoid double points.
            x2 = x1 - W/2  # X coordinates over the lenslet, with the origin at the center of the lenslet

            for i_lenslet in range(nr_of_lenslets):
                x_start_lenslet = i_lenslet*W - self.length/2 # Starting position of each lenslet w.r.t. the MLA's left side

                for i_pt in range(len(x2)):
                    y = np.sqrt(R**2-x2[i_pt]**2)
                    dy = y-yM
                    pts = np.append(pts, [np.array([p0[X] + x_start_lenslet + x1[i_pt], p0[Y] + self.thickness + dy])], axis=0)
        else:
            print('Microlens shape not yet defined')

        super().__init__(p0, pts, N=N, is_active=is_active)
        self.name = 'microlens array'

    def __str__(self):
        txt = 'Microlens array'
        return txt

    def plot(self, graph):
        poly = plt.Polygon(self.pts, closed=True, facecolor='cyan', alpha=geometry.alpha_from_N(self.N), zorder=5)
        graph.add_patch(poly)
        if config.getboolean('view', 'show_elements_properties'):
            canvas_class.plot_normals(graph, self)
