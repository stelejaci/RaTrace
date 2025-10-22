import numpy as np
from utils.varia import deg, X, Y
from utils.optics import N_air
from utils import geometry
from elements import element_class
from light import light_class
from utils.configuration_class import config


class DiffuseElementClass(element_class.ElementClass):
    def __init__(self, p0, n0, pts, Kd=0.0, Ks=0.0, alpha=1, nr_of_scattered_rays=1, n_light=None, is_active=True, is_visible=True):
        self.Kd = Kd        # Intensity of the diffuse reflection
        self.Ks = Ks        # Intensity of the specular reflection
        self.alpha = alpha  # Sharpness (spotsize) of the specular reflection
        self.nr_of_scattered_rays = nr_of_scattered_rays
        self.n_light = n_light  # The orientation roughly towards the light. Only for plotting the BSDF functionality.

        # Initialise the ElementsClass instance
        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Diffuse element'

        if n_light is None:
            self.n_principle = self.n0  # If no orientation is given, take the normal for simplicity's sake

    def propagate_ray(self, ray):
        # There are two methods for handling ray scattering at a surface
        # First method: intensity weighing
        #   - Generate rays with randomly distributed angles around the surface's normal (-90°->90°)
        #   - Scale the incoming intensity with the surface's reflectivity at the collision point.
        #   - Drawback: Lots of rays are generated and raytraced that go into the void, not efficient
        # Second method: importance sampling
        #   - Generate random rays with an angular distribution around the normal that follow the surface's reflectivity profile
        #   - Keep the intensity equal for all scattered rays.
        #   - Use the cumsum etc for getting the correct distribution of ray directions.
        #   - Benefit: more efficient method
        new_rays = list()

        # Virtual rays should not propagate further from diffuse surfaces
        if ray.is_virtual: return new_rays

        # The first method with intensity weighing
        angle_n = geometry.angle_from_vector(self.n_coll)
        angles_ro = angle_n - np.pi / 2 * (np.random.rand(self.nr_of_scattered_rays) * 2 - 1)
        for angle_ro in angles_ro:
            ro = geometry.vector_from_angle(angle_ro)
            reflectivity = self.calculate_reflectivity(ray.r, ro, self.n_coll, verbose=False)
            intensity = reflectivity * ray.intensity / self.nr_of_scattered_rays
            scattered_ray = light_class.RayClass(p0=ray.p1, r=ro, intensity=intensity, wavelength=ray.wavelength, ray_parent=ray, N=N_air, source_element=self, plot_color=ray.plot_color, is_active=True, is_visible=True)
            new_rays.append(scattered_ray)
        return new_rays

    def calculate_reflectivity(self, ri, ro, n, verbose=False):
        # Info on Lambertian surfaces, and diffuse reflection:
        # https://www.oceanopticsbook.info/view/surfaces/lambertian-brdfs --> Most comprehensive explanation to understand the misconception of "equal intensity"
        # --> The intensity of individually scattered rays (which go in ALL directions) goes with the cosine law (for a Lambertian surface)
        # --> Otherwise stated: the probability of a ray in a CERTAIN direction (in a range of all directions) goes according to the cosine law
        # --> However, when looking at the with a radiometer with a certain viewing space angle, total intensity integrates over all intensities coming from the entire visible surface patch (covered by the space angle)
        # --> And since the area covered goes with the inverse cosine of the viewing angle, the measured (!) intensity is constant from all directions
        # So, we let the rays scatter in ALL directions, but let their intensity scale with the cosine law of the angle between the normal and viewing direction V (not L !!)
        # ( Note: a laser only illuminates a very narrow area of the surface, so integration would also be over that very limited area )
        # https://www.researchgate.net/publication/273019368_Performance_analysis_of_a_car-to-car_visible_light_ommunication_system
        # https://www.researchgate.net/publication/265544226_THE_NON-FLUORESCENT_FLAT_PLATE_SOLAR_CONCENTRATOR

        if not (self.Kd) and (not self.Ks):
            total_reflection = 1.0
        else:
            L = -ri     # Direction to the laser source
            V =  ro     # V is the viewing direction, so it should point away from the surface.

            diffuse_reflection = np.absolute(self.Kd * np.dot(V, n))    # Is it dot(V,n) or dot(L,n) ??? --> dot(V,n), want dat is de intensiteit (probability) van de ray in de gescatterde richting
            diffuse_reflection = np.max([diffuse_reflection, 0])

            R = 2 * np.dot(n, L) * n - L                                # The reflected vector of ri over n
            specular_reflection = self.Ks * np.power(np.dot(V, R), self.alpha)

            total_reflection = specular_reflection + diffuse_reflection
            total_reflection = total_reflection / (self.Kd + self.Ks)  # Normalize to 1

            if verbose:  print(f'R={total_reflection}  /  D={diffuse_reflection}  /  S={specular_reflection}')
        return total_reflection

    def __str__(self):
        txt = f'{self.name} --> Element ID={self.ID}, p0={self.p0}, n0={self.n0}, Kd={self.Kd}, Ks={self.Ks}, alpha={self.alpha}, number of scattered rays={self.nr_of_scattered_rays}, n_light={self.n_light}'
        return txt

    def plot(self, graph):
        super().plot(graph)
        if config.getboolean('view', 'show_elements_properties'):
            self.plot_BSDF_at_P(P=self.p0, R_max=1, col='black', graph=graph)

    def plot_BSDF_at_P(self, P, R_max, col, graph):
        # One can only plot the BSDF with n0 as the (inverse) incoming light vector, since light can come from anywhere
        ri = self.n_light
        angle_sampling = np.linspace(-1, 1, 1000, endpoint=True) * 90*deg
        angle_main = geometry.angle_from_vector(self.n0)
        angles = angle_main + angle_sampling
        ro = geometry.vector_from_angle(angles)
        R = np.zeros(len(angles))
        for i_r in range(len(ro)):
            r = ro[i_r]
            R[i_r] = self.calculate_reflectivity(-ri, r, self.n0, verbose=False) # ri gets reversed in the calculation, so we must first un-reverse it
        graph.plot(P[X]+R_max*R*ro[:,X],  P[Y]+R_max*R*ro[:, Y],  c=col, linewidth=2)
