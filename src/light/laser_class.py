import numpy as np
from utils import varia
from utils.varia import mm, µm, nm, deg, X, Y
from utils import optics, geometry
from utils.optics import  N_air, N_glass
from light import light_class
from utils.configuration_class import config



class LaserClass(light_class.LightSourceClass):
    def __init__(self, p0=np.array([0,0]), n0=np.array([0,-1]), beamwaist=100*µm, beamwaist_distance=100*mm, sampling_width_at_BW=300*µm, wavelength=450*nm, M2=1, intensity=1, plot_color='wavelength'):
        self.beamwaist = beamwaist
        self.beamwaist_distance = beamwaist_distance
        self.sampling_width_at_BW = sampling_width_at_BW
        self.M2 = M2

        super().__init__(wavelength=wavelength, p0=p0, n0=n0, intensity=intensity, plot_color=plot_color)
        self.name = 'laser source'

        # Derived optical parameters
        self.w0 = self.beamwaist / 2    # Half beamwaist
        self.ZR = optics.GLB_calculate_ZR(w0=self.w0, M2=self.M2, wavelength=self.wavelength)   # Rayleigh distance

    def generate_rays(self, N_rays):
        # Equi-angular sampling of rays is not the correct way to go. The better and more correct approach is to sample equi-distant points along the horizontal axis,
        # and then give the rays the appropriate intensity according to the Gaussian distribution. The max horizontal sampling range (distance) is two times the
        # beam width at near edge, because at near edge the line is the closest and the widest, so appears the widest too. This is effectively an implementation
        # of a "random sampling with weighed intensities". Another approach is "importance sampling with equal intensities", this might converge faster. However,
        # the "importance" is different at different points in space (where the rays collide with objects), and it is thus difficult to know where and how to sample
        # at initialisation.
        # h_sampling = 20*mm  # The absolute height at which the beam is sampled, in positive direction, so this is near edge of the FOV
        x_sampling = (2*np.random.rand(N_rays)-1) * self.sampling_width_at_BW/2  # [-SBW -> +SBW], times 2 for beam width, times 2 for twice the beam width sampling range
        x_sampling[0] = 0   # Force the first ray of the beam to be the principal ray. It doesn't matter if this is the first, because all the others are random anyway
        self.nr_of_original_rays = N_rays
        for i_ray in range(N_rays):
            r = np.array([x_sampling[i_ray], -self.beamwaist_distance])  # The end point of the ray, and we suppose that the rays are pointing down initially
            r = geometry.normalize(r)  # Orientation of the ray

            angle = geometry.angle_between_vectors(self.n0, np.array([0,-1]))  # Until now, the rays are considered pointing down
            r = geometry.rotate_with_P_angle(pts=r, angle=angle, is_orientation=True)  # Also n0 is rotated by the next rotate_by_angle command, so first anti-rotate this direction

            # We give the ray originating from the laser None intensity, as the intensity at a certain point is not constant, but should be calculated when the collision point is known
            ray = light_class.RayClass(p0=self.p0, r=r, intensity=None, wavelength=self.wavelength, N=N_air, ray_parent=None, source_element=self, plot_color=self.plot_color)
            self.add_ray(ray)

    def calculate_intensity_at_P(self, P):
        return optics.GLB_calculate_intensity_at_P(w0=self.w0, ZR=self.ZR, intensity_total=self.intensity, p0=self.p0, n0=self.n0, beamwaist_distance=self.beamwaist_distance, P=P)

    def __str__(self):
        txt = f'Laser --> Element ID={self.ID}: p0={self.p0}, n0={self.n0}, w0={self.w0}, M2={self.M2}, ZR={self.ZR:.2f}mm, beamwaist={self.beamwaist/µm}µm, beamwaist distance={self.beamwaist_distance}mm, wavelength={self.wavelength/nm}nm'
        return txt

    def plot(self, graph):
        col = varia.colormap_wavelength(N=1, wavelength=self.wavelength)
        graph.scatter(self.p0[X], self.p0[Y], c=col[0], s=20)
        super().plot(graph)
        pt_BW = self.p0 + self.n0 * self.beamwaist_distance
        graph.plot([self.p0[X], pt_BW[X]], [self.p0[Y], pt_BW[Y]], color=col[0], linewidth=1, linestyle='dashed')
        if config.getboolean('view', 'show_elements_properties'):
            graph.quiver(self.p0[X], self.p0[Y], 100 * self.n0[X], 100 * self.n0[Y], color=col[0], width=0.5, units='xy')
            I_max = self.calculate_intensity_at_P(pt_BW)
            self.plot_gaussian_profile_at_P(P=pt_BW, I_max=I_max, sampling_width=self.sampling_width_at_BW, col=col, graph=graph)
            if self.rays  and  self.rays[0].p1 is not None:
                z = geometry.distance_between_2_points(self.p0, self.rays[0].p1)
                SW = optics.GLB_calculate_width_at_z(self.w0, ZR=self.ZR, z=z)
                self.plot_gaussian_profile_at_P(P=self.rays[0].p1, I_max=I_max, sampling_width=SW, col=col, graph=graph)

    def plot_gaussian_profile_at_P(self, P, I_max, sampling_width, col, graph):
        r = geometry.orientation_from_normal(self.n0)
        x_sampling = np.linspace(-1, 1, 1000, endpoint=True) * sampling_width/2
        pts = np.empty([len(x_sampling), 2])
        for i_pt in range(len(x_sampling)):
            pt = P + x_sampling[i_pt] * r
            I = 1 * self.calculate_intensity_at_P(pt) / I_max
            pts[i_pt,:] = pt - I * self.n0
        graph.plot(P[X] - np.array([-1,1])*self.r[X]*sampling_width/2, P[Y] - np.array([-1,1])*self.r[Y]*sampling_width/2, color='black', linewidth=1, linestyle='solid')
        graph.plot(P[X] - np.array([0,1])*self.n0[X], P[Y] - np.array([0,1])*self.n0[Y], color='black', linewidth=1, linestyle='solid')
        graph.plot(pts[:, X], pts[:, Y], color=col[0], linewidth=2, linestyle='solid')
        graph.scatter(P[X], P[Y], c=col[0], s=20)

