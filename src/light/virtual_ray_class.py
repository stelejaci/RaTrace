import numpy as np
from utils import varia
from utils.varia import nm, deg, X, Y
from utils.optics import N_air
from light  import light_class
from utils import geometry
from utils.configuration_class import config

WAVELENGTH_DEFAULT = 660*nm


class VirtualRayClass(light_class.LightSourceClass):
    def __init__(self, origin=None, destination=None, target_p0=None, target_n0=None, wavelength=WAVELENGTH_DEFAULT, intensity=1, plot_color=None):
        self.origin = origin
        self.destination = destination
        self.target_p0 = target_p0
        self.target_n0 = target_n0

        super().__init__(wavelength=wavelength, p0=np.array([0,0]), n0=np.array([1,0]), intensity=intensity*10, plot_color=plot_color, is_virtual=True)
        self.name = 'virtual ray'

        self.target_p_coll = None

    def generate_ray(self, origin_element):
        super().reset() # Explicitly reset, otherwise rays keep adding when simulating again

        from display import imager_class
        if isinstance(origin_element, imager_class.ImagerClass)  and  ('centroid' in self.origin)  and  (self.destination is not None):
            self.p0 = origin_element.COG_WCF
            self.n0 = geometry.normalize(self.destination - self.p0)
            ray = light_class.RayClass(p0=self.p0, r=self.n0, intensity=self.intensity, wavelength=self.wavelength, N=N_air, ray_parent=None, source_element=self, plot_color = self.plot_color, is_virtual=True)
            self.nr_of_original_rays += 1
        self.add_ray(ray)

    def process(self):
        if not self.rays: return  # Don't bother plotting if there are no rays yet, otherwise the source origin will be plotted on the wrong spot

        # Calculate intersection of the last ray with the target plane
        last_ray = self.rays[-1]
        self.target_r0 = geometry.orientation_from_normal(self.target_n0)
        [self.target_p_coll, t1, t2]  = geometry.intersection_of_PR_line_with_PR_line(pt1=last_ray.p0, r1=last_ray.r, pt2=self.target_p0, r2=self.target_r0)

    def __str__(self):
        txt = f'{self.name} {self.ID+1}/{light_class.LightSourceClass.nr_of_sources} --> ID={self.ID}, wavelength={self.wavelength/nm}nm, intensity={self.intensity}, plot color={self.plot_color}'
        return txt

    def plot(self, graph):
        if not self.rays: return  # Don't bother plotting if there are no rays yet, otherwise the source origin will be plotted on the wrong spot

        # Plotting the origin of the source
        print(f' --> Plotting origin of source {self.ID+1}/{light_class.LightSourceClass.nr_of_sources}')
        col = varia.load_colormap(color=self.plot_color, N_rays=1)
        col = col[0]
        graph.scatter(self.p0[X], self.p0[Y], c=col, s=10)

        # Plot the origin and principal radiating direction
        graph.quiver(self.p0[X], self.p0[Y], 1*self.n0[X], 1*self.n0[Y], color=col, width=0.002, scale=20)
        graph.scatter(self.target_p_coll[X], self.target_p_coll[Y], c=col, s=20, zorder=7)

        if config.getboolean('view', 'show_elements_properties'):
            graph.text(self.p0[X]-0.01, self.p0[Y], f'p=[{self.p0[X]:.3f}, {self.p0[Y]:.3f}]', color='black', horizontalalignment='right', verticalalignment='center', fontsize=8)
            if self.target_p_coll is not None:
                graph.text(self.target_p_coll[X]-0.01, self.target_p_coll[Y], f'p=[{self.target_p_coll[X]:.3f}, {self.target_p_coll[Y]:.3f}]', color='black', horizontalalignment='right', verticalalignment='center', fontsize=8)

        # Plot the rays itself
        super().plot(graph)

