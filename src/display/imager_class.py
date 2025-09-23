import numpy as np
import math
import os, sys
sys.path.append(os.path.abspath('..'))

from utils import varia
from utils.varia import mm, µm, nm, deg, X, Y
from utils import geometry
from utils.configuration_class import config
from display import display_class

IAS_min   = 0.15
IAS_width = 25
IAS_slope = 12




# Todo: getting rid of the pixels (???)

class ImagerClass(display_class.DisplayClass):
    def __init__(self, p0=np.array([0, 0]), n0=np.array([-1, 0]), length=10 * mm, pixel_size=5*µm, is_active=True, is_visible=True):
        self.pixel_size = pixel_size

        # Instantiate the display class
        super().__init__(p0=p0, n0=n0, length=length, is_active=is_active, is_visible=is_visible)
        self.name = 'imager'

        # Derived parameters & reset
        self.number_of_pixels = int(np.floor(self.length / self.pixel_size))
        self.reset()

    def reset(self):
        self.pixels    = self.generate_pixels()
        self.intensity = np.zeros((self.number_of_pixels,))
        self.phase     = np.zeros((self.number_of_pixels,))
        self.image     = np.zeros((10,self.number_of_pixels,))
        self.COG_IMF   = np.empty(1)
        self.COG_WCF   = np.empty(2)
        # self.pixels_x = np.linspace(0, self.length, self.number_of_pixels)
        self.pixels_x = np.array([pixel.x for pixel in self.pixels])
        super().reset()

    def generate_pixels(self):
        pixels = []

        x_range = np.linspace(self.pts[0][X], self.pts[1][X], self.number_of_pixels, endpoint=True)
        y_range = np.linspace(self.pts[0][Y], self.pts[1][Y], self.number_of_pixels, endpoint=True)

        for x, y in zip(x_range, y_range):
            pixel = PixelClass(pt=[x,y], r=self.r[0], n=self.n0, size=self.pixel_size, ID=len(pixels))
            pixels.append(pixel)

        return pixels

    # def set_ROI(self, ROI_position, ROI_width):
    #     nr_of_active_pixels = 0
    #     for pixel in self.pixels:
    #         if np.abs( pixel.pt_center[X] - ROI_position[X] ) > np.abs(ROI_width/2*self.r[X]):
    #             pixel.is_active = False
    #         else:
    #             nr_of_active_pixels += 1
    #     print('Enabling imager ROI: ' + str(nr_of_active_pixels) + '/' + str(self.number_of_pixels) + ' pixels active')
    #     return

    # TODO: image_multiplicator implementation, is now 1
    def process_image(self):
        # Checking which rays fall into which pixel, then adding intensity to that pixel
        E_tot = 0
        for ray in self.cast_rays:
            ray.imager_pixel_ID = int(np.round(ray.p1_element_rel*(self.number_of_pixels - 1)))  # The rays belongs to, or is cast into a pixel with this ID
            incoming_intensity = ray.intensity * self.calculate_imager_angular_sensitivity(-ray.r)  # Minus direction because the ray comes in the opposite direction of the imager normal
            phase = ray.phase_end  if  config.getboolean('simulation', 'use_phase_information')  else 0     # Take into account (or not) the ray's phase
            E_ray = np.sqrt(incoming_intensity) * np.exp(1j * np.array(phase))
            self.pixels[ray.imager_pixel_ID].E += E_ray    # Add all complex intensities to each pixel

        # Transfer pixel intensity values into a 1D-array
        for i_px in range(len(self.pixels)):
            self.pixels[i_px].intensity = np.abs(  self.pixels[i_px].E) ** 2
            self.pixels[i_px].phase     = np.angle(self.pixels[i_px].E)

            self.intensity[i_px] = self.pixels[i_px].intensity
            self.phase[i_px]     = self.pixels[i_px].phase

        # For better visual interpretation, make a 2D image
        self.image = np.tile(self.intensity, (10, 1))

        self.process_centroid()

    def process_centroid(self):
        # Calculate peak and summed intensities
        self.peak_intensity = np.max(self.intensity)
        self.summed_intensity = np.sum(self.intensity)

        # Masking the image above a certain intensity threshold
        self.peak_cutoff_threshold = 0.00 * self.peak_intensity  # The intensity threshold for background removal
        self.intensity_masked = np.copy(self.intensity)
        self.intensity_masked[self.intensity_masked < self.peak_cutoff_threshold] = 0.0  # Background-substracted image

        # Centre-of-gravity (or COG) of the pulse, in imager frame coordinates
        x = np.linspace(0, self.number_of_pixels - 1, self.number_of_pixels)
        self.summed_intensity_masked = np.sum(self.intensity_masked)
        self.COG_IMF = np.dot(self.intensity_masked, x) / self.summed_intensity_masked

        fraction = self.COG_IMF / (self.number_of_pixels - 1)  # Fraction of the total imager length the COG is positioned
        pts = varia.sort_x_left_to_right(self.pts)  # Sort the imager points in +X order because it might be inverted, corrupting the following calculation
        self.COG_WCF[X] = self.pts[0,X] + fraction * (self.pts[1,X] - self.pts[0,X])  # Position of the COG in world frame coordinates
        self.COG_WCF[Y] = self.pts[0,Y] + fraction * (self.pts[1,Y] - self.pts[0,Y])

        # Calculating the FWHM
        [self.FWHM, self.FWHM_pts] = varia.calculate_FWHM(self.intensity)

    def calculate_imager_angular_sensitivity(self, r, verbose=False):
        # Empirical formula that approximates the experimental results, see screenshot from Matlab simulation
        # Also, see report:  " REP - LS - SG18003 VITA5000 imager efficiency vs incident angle.docx "
        # TODO: Pre-compute the erf-function

        r = geometry.normalize(r)
        angle = geometry.angle_between_vectors(self.n0, r) * 180 / math.pi
        multiplicator = (1 / 2) * ((1 + IAS_min) + (1 - IAS_min) * math.erf((-angle + IAS_width) / IAS_slope))
        return 1 + 0 * multiplicator

    # def plot_imager_efficiency_map(self, p, scale=1, col='black'):
    #     print("Plotting imager efficiency map")
    #
    #     angle_normal = misc.angle_from_vector(self.n)
    #     self.efficiency_map_angle = np.arange(angle_normal + 90 * deg, angle_normal - 90 * deg, -1 * deg)
    #     self.efficiency_map_points = np.zeros((len(self.efficiency_map_angle), 2))
    #     self.efficiency_map = np.zeros((len(self.efficiency_map_angle, )))
    #
    #     for i_angle in range(len(self.efficiency_map_angle)):
    #         angle = self.efficiency_map_angle[i_angle]
    #         V = np.array([np.cos(angle), np.sin(angle)])
    #         self.efficiency_map[i_angle] = self.calculate_imager_angular_sensitivity(V, verbose=False)
    #         self.efficiency_map_points[i_angle, X] = p[X] + V[X] * scale * self.efficiency_map[i_angle]
    #         self.efficiency_map_points[i_angle, Y] = p[Y] + V[Y] * scale * self.efficiency_map[i_angle]
    #     plt.plot(self.efficiency_map_points[:, X], self.efficiency_map_points[:, Y], color=col, linewidth=2)
    #     return

    def __str__(self):
        txt = f'Imager --> Element ID={self.ID}, p0={self.p0}, n0={self.n0}, length={self.length}, pixel size={self.pixel_size}, number of pixels={self.number_of_pixels}'
        return txt

    def plot(self, graph):
        super().plot(graph)
        if config.getboolean('view', 'show_pixels'):
            for pixel in self.pixels:
                pixel.plot(graph)
            graph.scatter(self.pixels[ 0].p0[X],    self.pixels[0].p0[Y],   s=20, facecolor='g', marker='o')
            graph.scatter(self.p0[X],               self.p0[Y],             s=20, facecolor='g', marker='o')
            graph.scatter(self.pixels[-1].p0[X],    self.pixels[0-1].p0[Y], s=20, facecolor='g', marker='o')
            graph.text(self.pixels[0].p0[X]  - 0.2*self.n0[X], self.pixels[0].p0[Y]  - 0.2*self.n0[Y],  f'pixel 0: p=[{self.pixels[0].p0[X]:.3f}, {self.pixels[0].p0[Y]:.3f}]', color='g', horizontalalignment='left', verticalalignment='center', fontsize=8)
            graph.text(self.p0[X]            - 0.2*self.n0[X], self.p0[Y]            - 0.2*self.n0[Y],  f'pixel {self.number_of_pixels/2-1}: p=[{self.p0[X]:.3f}, {self.p0[Y]:.3f}]', color='g', horizontalalignment='left', verticalalignment='center', fontsize=8)
            graph.text(self.pixels[-1].p0[X] - 0.2*self.n0[X], self.pixels[-1].p0[Y] - 0.2*self.n0[Y],  f'pixel {self.number_of_pixels-1}: p=[{self.pixels[-1].p0[X]:.3f}, {self.pixels[-1].p0[Y]:.3f}]', color='g', horizontalalignment='left', verticalalignment='center', fontsize=8)


class PixelClass:
    def __init__(self, pt=[0, 0], r=[1, 0], n=[0, 1], size=5, ID=0):
        self.p0 = pt    # Coordinate at the center of the pixel
        self.n = n
        self.r = r
        self.size = size
        self.ID = ID
        self.intensity = 0
        self.E = 0
        self.ray_ID_list = list()
        self.x = ID*size    # U-coordinate, i.e. in imager frame

        self.pts = geometry.points_from_position_direction_length(p0=self.p0, r=self.r, L=self.size, sort_left_to_right=True, symmetric=True)

        self.is_active = True

    def plot(self,graph):
        graph.scatter(self.p0[X], self.p0[Y], s=5, facecolor='g', marker='o')
        graph.scatter(self.pts[1,X], self.pts[1,Y], s=20, facecolor='g', marker='+', linewidth=1)
        if self.ID == 0:
            graph.scatter(self.pts[0,X], self.pts[0,Y], s=20, facecolor='g', marker='+', linewidth=1)
