import numpy as np
from utils import varia
from utils.varia import mm, X, Y
from utils import geometry
from utils.configuration_class import config
from elements import element_class

IAS_min   = 0.15
IAS_width = 25
IAS_slope = 12



class DisplayClass(element_class.ElementClass):
    nr_of_displays = 0

    def __init__(self, p0=np.array([0,0]), n0=np.array([-1,0]), length=10*mm, is_active=True, is_visible=True):
        # Basic parameters
        self.length = length
        self.name = 'Display'

        # Derived parameters
        r = geometry.orientation_from_normal(n0)
        pts = geometry.points_from_position_direction_length(p0, r, self.length, symmetric=True, sort_left_to_right=False)

        DisplayClass.reset(self)
        super().__init__(p0=p0, n0=n0,  pts=pts, is_active=is_active, is_visible=is_visible)
        DisplayClass.nr_of_displays += 1

    def reset(self):
        self.cast_rays = list()

    def check_collision(self, ray):
        if not self.is_active:
            return [None, None, None]
    
        # Calculate intersection between the 2-pt lens line and pt-and-rico ray line
        [p, t0, t1] = geometry.intersection_of_PP_line_with_PR_line(p00=self.pts[0], p01=self.pts[1], p10=ray.p0, r1=ray.r)
    
        # If the ray does not hit the imager in forward direction, return none
        if t0 > 1 or t0 < 0 or t1 < 0:  # t0<0: p lies before line segment | t0>1: p lies behind the line segment | t1<0: p lies in the backward direction of the ray
            return [None, None, None]
        
        return [p, t1, t0]

    def propagate_ray(self, ray):
        # Add the ray to the list of rays cast onto the imager, because the raytrace_ray function adds ray information only after the check_collision function
        self.cast_rays.append(ray)
        return

    def process_cast_rays(self):
        print('start processing')
        self.cast_rays_x     = np.empty(len(self.cast_rays))
        self.cast_rays_col   = np.empty((len(self.cast_rays),4))
        self.cast_rays_ID    = np.empty(len(self.cast_rays))
        self.cast_rays_phase = np.empty(len(self.cast_rays))

        for i_cast_ray in range(len(self.cast_rays)):
            self.cast_rays_ID[i_cast_ray]       = self.cast_rays[i_cast_ray].ID
            self.cast_rays_x[i_cast_ray]        = self.cast_rays[i_cast_ray].p1_element_rel * self.length
            self.cast_rays_phase[i_cast_ray]    = self.cast_rays[i_cast_ray].phase_end

            # Store the ID and cast position of the ray onto the display
            self.cast_rays[i_cast_ray].display_ID = self.ID
            self.cast_rays[i_cast_ray].display_x  = self.cast_rays_x[i_cast_ray]

            # In case of one single ray simulated, this results in a list or so --> STILL TO DEBUG !!!
            col = self.cast_rays[i_cast_ray].plot_color
            if isinstance(col, list):
                col = col[0]
            col = varia.load_colormap(color=col ,N_rays=1, wavelength=self.cast_rays[i_cast_ray].wavelength)
            self.cast_rays_col[i_cast_ray] = col[0]

        print('end of processing')

    def __str__(self):
        txt = f'Display --> Element ID={self.ID}, p0={self.p0}, n0={self.n0}, length={self.length}'
        return txt

    def plot(self, graph):
        if self.is_visible:
            graph.plot([self.pts[0,X], self.pts[1,X]], [self.pts[0,Y], self.pts[1,Y]], color='green', linewidth=2)
            if config.getboolean('view', 'show_elements_properties'):
                graph.scatter(self.p0[X], self.p0[Y], s=10, facecolor='g')
                p_txt = self.p0 - self.n0*0.5
                graph.text(p_txt[X], p_txt[Y], f'{self.name} {self.ID}', color='green', horizontalalignment='left', verticalalignment='bottom', fontsize=10, rotation=0)
            super().plot(graph)
