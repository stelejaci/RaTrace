import numpy as np
import matplotlib.pyplot as plt
import math
from utils.varia import nm, deg, X, Y
from utils.optics import N_air


class FilterClass:
    def __init__(self, p0, r, D, CWL_0=450*nm, T_max_0=0.9, FWHM_0=30.5*nm, slope_0=0.45):
        self.p0 = p0
        self.D = D
        self.r = r/np.linalg.norm(r)
        self.is_active = True
        self.name = 'filter'

        self.T_max_0 = T_max_0
        self.CWL_0   = CWL_0
        self.FWHM_0  = FWHM_0
        self.slope_0 = slope_0
        
        # Derived parameters
        self.n = misc.calculate_perpendicular_direction_to_vector(self.r)
        self.pts = misc.create_points_from_position_direction_length(self.p0, self.r, self.D, symmetric=True, sort_left_to_right=True)

    def check_collision(self, ray):
        if not self.is_active:
            return [None, None, None]
    
        # Calculate intersection between the 2-pt lens line and pt-and-rico ray line
        [p, t0, t1] = misc.calculate_intersection_of_PP_line_with_PR_line(p00=self.pts[0], p01=self.pts[1], p10=ray.p0, r1=ray.r)
    
        # If the ray does not hit the black body in forward direction, return none
        if t0 > 1 or t0 < 0 or t1 < 0:  # t0<0: p lies before line segment | t0>1: p lies behind the line segment | t1<0: p lies in the backward direction of the ray
            return [None, None]
    
        return [p, t1, t0]

    def propagate_ray(self, ray):
        # print('Tracing rays towards filter')
        
        filter_transmission = self.calculate_filter_transmission(ray.r, ray.wavelength)
        I = filter_transmission * ray.I
        ray_new = light_functions.RayClass(p0=ray.p1, r=ray.r, I=I, wavelength=ray.wavelength, ID_parent=ray.ID, N=N_air, col=ray.col)
        # rays.append(ray_new)
        
        return ray_new
    
    def calculate_filter_transmission(self, r, wavelength):
        r = misc.normalize(r)
        AOI = misc.angle_between_vectors(self.n, r) * 180 / np.pi

        # In the empirical formulae, derived in Matlab, we used dimension-less values for CWL, FWHM and wavelength, e.g. 450 instead of 450nm for the CWL
        # THESE EMPRIRICAL FORMULAE ARE ONLY VALID FOR THE 450nm FILTER
        CWL   = self.CWL_0/nm   - 1 * np.power(AOI/9,  2)
        FWHM  = self.FWHM_0/nm  - 1 * np.power(AOI/23, 2)
        T_max = self.T_max_0 - 0 * np.power(AOI/1,  2)
        slope = self.slope_0 - 1 * np.power(AOI/65, 2)
        
        T = T_max/4 * (1 + math.erf((FWHM/2 + (wavelength/nm - CWL)) * slope)) * (1 + math.erf((FWHM/2 - (wavelength/nm - CWL)) * slope))
  
        return T
        
    def plot(self, col='purple'):
        print("Plotting filter")

        if self.is_active:
            plt.plot([self.pts[0][X], self.pts[1][X]], [self.pts[0][Y], self.pts[1][Y]], color=col, linewidth=4)

    def plot_transmission_map(self, p, scale=1, col='black', wavelength=450*nm, verbose=False, ):
        angle_normal = misc.angle_from_vector(self.n)
        transmission_map_angle = np.arange(angle_normal + 90*deg, angle_normal - 90*deg, -1*deg)
        transmission_map_points = np.zeros((len(transmission_map_angle), 2))
        transmission_map = np.zeros((len(transmission_map_angle, )))

        for i_angle in range(len(transmission_map_angle)):
            angle = transmission_map_angle[i_angle]
            V = np.array([np.cos(angle), np.sin(angle)])
            transmission_map[i_angle] = -self.calculate_filter_transmission(V, wavelength=wavelength)
            transmission_map_points[i_angle, X] = p[X] + V[X] * scale * transmission_map[i_angle]
            transmission_map_points[i_angle, Y] = p[Y] + V[Y] * scale * transmission_map[i_angle]
        plt.plot(transmission_map_points[:, X], transmission_map_points[:, Y], color=col, linewidth=2)
        
        return
