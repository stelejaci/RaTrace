import numpy as np
from light import light_class
from utils import varia
from utils.varia import mm, µm, nm, deg, X, Y
from utils.optics import N_air, N_glass
from utils import geometry

WAVELENGTH_DEFAULT = 660*nm




class DiffusePlaneSourceClass(light_class.LightSourceClass):
    def __init__(self, p0=[0,0], n0=[1,0], diameter=10*mm, wavelength=WAVELENGTH_DEFAULT, fan_angle=30*deg, intensity=1, plot_color=None):
        self.diameter = diameter
        self.fan_angle = fan_angle

        super().__init__(wavelength=wavelength, p0=p0, n0=n0, intensity=intensity, plot_color=plot_color)
        self.name = 'Diffuse plane source'

    def generate_rays(self, N_rays):
        # First, reset all ray information
        super().reset()

        # Set the number of rays originating from the source
        self.nr_of_original_rays = N_rays

        # Setting the plotting colors
        cols = varia.load_colormap(color=self.plot_color, N_rays=N_rays, wavelength=self.wavelength)

        # Generate the rays on random places with random angles
        t_range = np.random.rand(N_rays) * 2 - 1    # t = [-1,+1]
        t_range = np.sort(t_range)
        angle_range = (np.random.rand(N_rays)-0.5)*self.fan_angle   # +/- half angle
        main_angle = geometry.angle_from_vector(self.n0)

        for ind in range(len(t_range)):
            p0 = self.p0 + self.r * self.diameter/2 * t_range[ind]
            intensity = self.intensity/N_rays
            col = cols[round(ind / max([(self.nr_of_original_rays - 1),1]) * (len(cols) - 1))]
            angle = main_angle + angle_range[ind]
            r = geometry.vector_from_angle(angle)
            ray = light_class.RayClass(p0=p0, r=r, intensity=intensity, wavelength=self.wavelength, N=N_air, ray_parent=None, source_element=self, plot_color=col)
            self.add_ray(ray)

    def __str__(self):
        n_str  = varia.format_1D_array(arr=self.n0, fmt='7.4f')
        p0_str = varia.format_1D_array(arr=self.p0, fmt='7.4f')
        txt = f'Diffuse plane source nr {self.ID+1}/{light_class.LightSourceClass.nr_of_sources} --> ID={self.ID}, p0={p0_str}, n0={n_str}, nr of rays={self.nr_of_rays}, wavelength={self.wavelength/nm}nm, intensity={self.intensity}, diameter={self.diameter/mm}mm,  plot color={self.plot_color}'
        return txt

    def plot(self, graph):
        # Plotting the origin of the source, e.g. a plane, a point, a laser source, ...
        print(f' --> Plotting origin of source {self.ID+1}/{light_class.LightSourceClass.nr_of_sources}')
        col = varia.colormap_wavelength(N=1, wavelength=self.wavelength)
        col = col[0]
        graph.plot(self.p0[X] + self.r[X] * self.diameter / 2 * np.array([-1, 1]),  self.p0[Y] + self.r[Y] * self.diameter / 2 * np.array([-1, 1]), color=col,  linewidth=3, linestyle='solid')
        # Plot the principal radiating direction
        graph.quiver(self.p0[X], self.p0[Y], 10*self.n0[X], 10*self.n0[Y], color=col, width=0.2, units='xy')
        # Plot the rays
        super().plot(graph)






if __name__ == '__main__':
    print('\nTest 1')
    rays = []
    rays.append(light_class.RayClass())
    rays.append(light_class.RayClass())
    print(rays[0])
    print(rays[1])

    print('\nTest 2')
    PlaneSource = PlaneSourceClass(p0=np.array([0,0]), n=np.array([1,0]), wavelength=450*nm, intensity=10, diameter=10*mm, intensity_distribution='equidistant', plot_color='rainbow')
    PlaneSource.generate_rays(N_rays=3)
    print(PlaneSource)
    for ray in PlaneSource.rays:
        print(ray)
