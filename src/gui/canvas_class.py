import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
from utils.varia import X, Y
from utils.configuration_class import config


PLOT_COLOR_BLACK = 'grey'
MAX_NR_OF_SCATTERPOINTS_PLOTTED = 100000
rng = np.random.default_rng()


class CanvasClass(FigureCanvasQTAgg):

    def __init__(self, simulation=None, width=5, height=4, dpi=100):
        self.simulation = simulation
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.graph = fig.add_subplot(111)
        super().__init__(fig)
        fig.subplots_adjust(left=0.055, right=0.99, bottom=0.03, top=0.99)
        self.axis_lims = []

    def set_nr_of_rays_to_plot(self, value):
        self.nr_of_rays_to_plot = value

    def clear(self):
        if config.getboolean('scenes', 'reset_axis_after_loading_scene'):
            self.axis_lims = []
        # self.graph.cla()

    def update_entire_scene(self):
        # Store the current axis limits before plotting, to set again later on
        if self.axis_lims:
            self.axis_lims = self.graph.axis()

        # Clear the figure and plot all the items
        self.update_items()

        # Set the correct axis limits
        if not self.axis_lims:
            self.axis_lims = self.graph.axis()
        self.graph.axis(self.axis_lims)

        # Set the background to black if required
        col = 'black'   if    config.getboolean('view','black_background')    else    'white'
        self.graph.set_facecolor(col)

        # Some graph aesthetics
        SHOW_AXIS_AND_GRID = config.getboolean('view', 'show_axis_and_grid')
        self.graph.grid(visible=SHOW_AXIS_AND_GRID, alpha=0.2)
        self.graph.tick_params(axis='both', which='both', bottom=SHOW_AXIS_AND_GRID, top=SHOW_AXIS_AND_GRID, left=SHOW_AXIS_AND_GRID, right=SHOW_AXIS_AND_GRID)
        self.graph.tick_params(axis='both', which='both', labelbottom=SHOW_AXIS_AND_GRID, labelleft=SHOW_AXIS_AND_GRID)
        self.graph.spines['top'].set_visible(SHOW_AXIS_AND_GRID)
        self.graph.spines['right'].set_visible(SHOW_AXIS_AND_GRID)
        self.graph.spines['bottom'].set_visible(SHOW_AXIS_AND_GRID)
        self.graph.spines['left'].set_visible(SHOW_AXIS_AND_GRID)
        # self.graph.text((self.axis_lims[0]+self.axis_lims[1])/2, self.axis_lims[3]-0.1*(self.axis_lims[3]-self.axis_lims[2]), f'{self.simulation.info}', fontsize=14, horizontalalignment='center', verticalalignment='bottom')
        self.draw()

    def update_items(self):
        self.graph.cla()
        for source in self.simulation.sources:
            source.plot(self.graph)
        for element in self.simulation.elements:
            element.plot(self.graph)
        for display in self.simulation.displays:
            display.plot(self.graph)
        self.graph.axis('equal')
        self.draw()

    # # @timeit
    # def plot_full_ray(self, source, i_ray):
    #     self.plot_ray_segment(source, i_ray)
    #     for ID_child in source.rays[i_ray].ID_children:
    #         self.plot_full_ray(source, ID_child)
    #
    #
    # def plot_flat_mirror(self, mirror):
    #     poly = plt.Polygon(mirror.pts, closed=True, facecolor='grey', alpha=1.0, zorder=5)
    #     self.graph.add_patch(poly)
    #     self.graph.plot(mirror.pts[[0,1],0], mirror.pts[[0,1],1], c=PLOT_COLOR_BLACK, zorder=6)
    #     if config.getboolean('view', 'show_elements_properties'):
    #         self.plot_normals(self.graph, mirror)
    #
    # def plot_black_body(self, surface):
    #     poly = plt.Polygon(surface.pts, closed=True, facecolor=PLOT_COLOR_BLACK)
    #     self.graph.add_patch(poly)
    #     if config.getboolean('view', 'show_elements_properties'):
    #         self.plot_normals(self.graph, surface)
    #
    # def plot_displays(self):
    #     for display in self.simulation.displays:
    #         display.plot(self.graph)
    #         # if isinstance(display, display_functions.ImagerClass):
    #         #     self.plot_imager(display)
    #         # else:
    #         #     self.plot_display(display)
    #
    # def plot_plate_BSDF(self, surface):
    #         p0 = ( surface.pts[0] + surface.pts[1] )/2
    #         n  = surface.n[0]
    #         ri = self.simulation.sources[0].r
    #         p1_n  = p0 + 10*n
    #         p1_ri = p0 - 10*ri  # ri is from the source, towards the plate
    #
    #         # angles = np.linspace(-np.pi/2, np.pi/2, 11)
    #         angles = np.linspace(-np.pi/2, np.pi/2, 361)
    #         R  = np.zeros(len(angles))
    #         ro = np.zeros([len(angles),2])
    #         angle_n = geometry.angle_from_vector(n)
    #         for i_angle in range(len(angles)):
    #             angle_o = angle_n - angles[i_angle]
    #             ro[i_angle,:] = geometry.vector_from_angle(angle_o)
    #             R[i_angle]  = surface.calculate_reflectivity(-ri, ro[i_angle], n)
    #         self.graph.plot(p0[X]+10*R*ro[:,X], p0[Y]+10*R*ro[:,Y], c=PLOT_COLOR_BLACK, linewidth=5)
    #
    # def plot_sphere(self, surface):
    #     circle = plt.Circle(surface.p0, radius=surface.R, facecolor='0.5', edgecolor='0.5')
    #     self.graph.add_patch(circle)



class CanvasDisplayClass(FigureCanvasQTAgg):
    def __init__(self, simulation=None, width=5, height=4, dpi=100):
        self.simulation = simulation
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.graph = fig.add_subplot(111)
        super().__init__(fig)
        fig.subplots_adjust(left=0.15, right=0.99, bottom=0.10, top=0.99)
        self.greyscale_mode = False
        self.zoom_on_centroid = True

    def update_graphs(self, graph_type_ind):
        for display in self.simulation.displays:
            self.update_graph(display=display, graph_type_ind=graph_type_ind)

    def update_graph(self, display, graph_type_ind):
        self.graph.cla()
        self.xlims, self.ylims = np.array([0, 1]),  np.array([0, 1])

        if not display.cast_rays:
            return

        # self.graph_types = ["1D scattered", "2D scattered", "2D greyscale", "Centroid", "Phase plot"]
        if graph_type_ind == 0:
            if len(display.cast_rays) <= MAX_NR_OF_SCATTERPOINTS_PLOTTED:
                size_pts = 100/len(display.cast_rays)
                self.graph.scatter(display.cast_rays_x, 0*display.cast_rays_x, s=size_pts, c=display.cast_rays_col)
                self.xlims = np.array([0,display.length])
                self.ylims = np.array([-1,1])
                # self.graph.set_ylim(-1,1)
                self.graph.tick_params(axis='x', which='both', top=False, bottom=config.getboolean('view', 'show_axis_and_grid'))
                self.graph.grid(config.getboolean('view', 'show_axis_and_grid'), which='both', axis='x')
                self.graph.set_yticks([])
        elif graph_type_ind == 1:
            if len(display.cast_rays) <= MAX_NR_OF_SCATTERPOINTS_PLOTTED:
                size_pts = 100/len(display.cast_rays)
                random_y_array = 0.5+np.random.random(*display.cast_rays_x.shape)
                self.graph.scatter(display.cast_rays_x, random_y_array, s=size_pts, c=display.cast_rays_col)
                self.xlims = np.array([0,display.length])
                self.ylims = np.array([0,2])
                self.graph.grid(config.getboolean('view', 'show_axis_and_grid'), which='both', axis='x')
                self.graph.set_yticks([])
        elif graph_type_ind == 2:
            if self.greyscale_mode:
                self.graph.imshow(display.image, cmap='gray', aspect='auto')
            else:
                self.graph.imshow(display.image, cmap='jet', aspect='auto')
            self.graph.tick_params(axis='x', which='both', top=False, bottom=config.getboolean('view', 'show_axis_and_grid'))
            self.ylims = np.array([0,1.1*display.peak_intensity])
            if self.zoom_on_centroid and display.FWHM is not None:
                self.xlims = display.COG_IMF + 3 * display.FWHM * np.array([-1,1])
            else:
                self.xlims = np.array([0,display.number_of_pixels-1])
        elif graph_type_ind == 3:
            self.graph.plot(display.intensity, color='blue', marker='.', linestyle='-', linewidth=1, markersize=4)
            if display.FWHM is not None:
                self.graph.plot((display.FWHM_pts[0][X], display.FWHM_pts[1][X]), (display.FWHM_pts[0][Y], display.FWHM_pts[1][Y]), color='red', marker='.', linestyle='--', linewidth=1, markersize=4)
                self.graph.plot((display.FWHM_pts[0][X]+display.FWHM_pts[1][X])/2*np.array([1,1]), display.peak_intensity*np.array([0, 1]), color='red', marker='.', linestyle='--', linewidth=1, markersize=4)
                self.graph.text(self.xlims[1] - 0.1*display.FWHM, self.ylims[1]-0.05*display.peak_intensity, f'COG = {display.COG_IMF:.2f}', color=(0,1,0), ha='right')
                self.graph.text(self.xlims[1] - 0.1*display.FWHM, self.ylims[1]-0.10*display.peak_intensity, f'FWHM = {display.FWHM:.1f}', color=(1,0,0), ha='right')
            self.graph.plot(display.COG_IMF * np.array([1, 1]), display.peak_intensity * np.array([0, 1]), color=(0, 1, 0), marker='.', linestyle='--', linewidth=1, markersize=4)
            if self.zoom_on_centroid  and  display.FWHM is not None:
                self.xlims = display.COG_IMF + 3 * display.FWHM * np.array([-1, 1])
            else:
                self.xlims = np.array([0, display.number_of_pixels - 1])
            self.ylims = np.array([0, 1.1 * display.peak_intensity])
            self.graph.grid(config.getboolean('view', 'show_axis_and_grid'), which='both')
        elif graph_type_ind == 4:
            if len(display.cast_rays) <= np.inf + MAX_NR_OF_SCATTERPOINTS_PLOTTED:
                size_pts = 100 / np.sqrt(len(display.cast_rays))
                self.graph.scatter(display.cast_rays_x, display.cast_rays_phase, s=size_pts, c=display.cast_rays_col)
                self.graph.scatter(display.pixels_x + display.pixel_size/2, display.phase, s=10+0*size_pts, facecolor=(1,1,1,0), edgecolor=(0,0,0,1), marker='o')
                self.xlims = np.array([0, display.length])
                self.ylims = np.array([-np.pi, np.pi])
                self.graph.tick_params(axis='x', which='both', top=False, bottom=config.getboolean('view', 'show_axis_and_grid'))
                self.graph.grid(config.getboolean('view', 'show_axis_and_grid'), which='both', axis='x')
                # self.graph.set_yticks([])

        self.graph.set_xlim(self.xlims)
        self.graph.set_ylim(self.ylims)

        self.draw()
        print('Display graph updated')
