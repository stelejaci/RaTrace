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

        # Generate reflection_tables
        self.reflection_table_angular_spacing_deg = 0.01    # This is in degrees
        self.generate_diffuse_and_specular_cumsum_tables()

        # Initialise the ElementsClass instance
        super().__init__(p0=p0, n0=n0, pts=pts, is_active=is_active, is_visible=is_visible)
        self.name = 'Diffuse element'

        if n_light is None:
            self.n_principle = self.n0  # If no orientation is given, take the normal for simplicity's sake

    def generate_diffuse_and_specular_cumsum_tables(self):
        # The entire range of angles between -90° and +90°.
        self.angles_deg = np.arange(start=-90, stop=90+self.reflection_table_angular_spacing_deg, step=self.reflection_table_angular_spacing_deg)
        # The arange command suffers from floating point errors and results in numbers such as 44.9999999996 instead of 45.0. This gives the correct precision to the numbers
        self.angles_deg = np.round(self.angles_deg / self.reflection_table_angular_spacing_deg) * self.reflection_table_angular_spacing_deg

        self.diffuse_reflection  = self.Kd * np.absolute(np.cos(self.angles_deg*deg))           # The angle here is between r_out and n
        self.specular_reflection = self.Ks * np.power(np.cos(self.angles_deg*deg), self.alpha)  # The angle here is between r_out and reflected vector (over n) of r_in
        self.diffuse_reflection_cumsum  = np.cumsum(self.diffuse_reflection)
        self.specular_reflection_cumsum = np.cumsum(self.specular_reflection)

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

        # The second method with importance sampling
        for i_ray in range(self.nr_of_scattered_rays):
            ro = self.generate_random_deflection_direction()
            intensity = ray.intensity / self.nr_of_scattered_rays
            scattered_ray = light_class.RayClass(p0=ray.p1, r=ro, intensity=intensity, wavelength=ray.wavelength, ray_parent=ray, N=N_air, source_element=self, plot_color=ray.plot_color, is_active=True, is_visible=True)
            new_rays.append(scattered_ray)
        return new_rays

    def generate_random_deflection_direction(self, ri, n):
        angle_n  = geometry.angle_from_vector(self.n_coll)  # The angle of the normal
        angle_ri = geometry.angle_from_vector(-ri)          # The angle of the incoming ray, but reversed
        angle_rr = 2*angle_n - angle_ri                 # The angle of the perfectly reflected ray

        self.generate_total_reflection_cumsum_table_for_angle(angle_rr=angle_rr)

        cumsum_max = np.max(self.total_reflection_cumsum)
        cumsum_random = np.random.rand() * cumsum_max
        angle_random_deg = np.interp(cumsum_random, self.total_reflection_cumsum, self.angles_deg)
        angle_ro = angle_n - angle_random_deg*deg
        ro = geometry.vector_from_angle(angle_ro)

        return ro

    def generate_total_reflection_cumsum_table_for_angle(self, angle_rr):
        # All angles in this method are in degrees
        # Rounding the angle to one of the entries in the angles table-entries
        angle_rr = angle_rr/deg
        angle_rr = np.round(angle_rr/self.reflection_table_angular_spacing_deg) * self.reflection_table_angular_spacing_deg
        nr_of_entries = len(self.diffuse_reflection_cumsum)

        # print(angle_rr)
        i_angle_rr = np.where(self.angles_deg==angle_rr)[0][0]
        i_angle_0  = np.where(self.angles_deg==0*deg   )[0][0]

        self.total_reflection_cumsum = np.copy(self.diffuse_reflection_cumsum)
        self.total_reflection = np.copy(self.diffuse_reflection)

        if i_angle_rr > 0:
            self.specular_reflection_cumsum_shifted = np.concatenate( (np.zeros(i_angle_rr-i_angle_0), self.specular_reflection_cumsum[0:nr_of_entries-(i_angle_rr-i_angle_0)]), axis=0 )
            self.total_reflection_cumsum += self.specular_reflection_cumsum_shifted

        if False:
            fig, ax = plt.subplots(3,2)
            ax[0,0].plot(self.angles_deg, self.diffuse_reflection)
            ax[1,0].plot(self.angles_deg, self.specular_reflection)
            ax[2,0].plot(self.angles_deg, self.total_reflection)
            ax[0,1].plot(self.angles_deg, self.diffuse_reflection_cumsum)
            ax[1,1].plot(self.angles_deg, self.specular_reflection_cumsum)
            ax[2,1].plot(self.angles_deg, self.total_reflection_cumsum)
            ax[0,0].grid()
            ax[1,0].grid()
            ax[2,0].grid()
            ax[0,1].grid()
            ax[1,1].grid()
            ax[2,1].grid()
            plt.show()

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
