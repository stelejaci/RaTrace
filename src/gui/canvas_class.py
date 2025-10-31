import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np
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
            x_mm  = np.array([pixel.x for pixel in display.pixels])
            x_ind = np.arange(display.number_of_pixels)
            self.graph.plot(x_mm, display.intensity, color='blue', marker='.', linestyle='-', linewidth=1, markersize=4)
            self.ylims = np.array([0, 1.1 * display.peak_intensity])
            self.xlims = np.array([0, display.length])
            self.graph.grid(config.getboolean('view', 'show_axis_and_grid'), which='both')
        elif graph_type_ind == 3:
            if self.greyscale_mode:
                self.graph.imshow(display.image, cmap='gray', aspect='auto')
            else:
                self.graph.imshow(display.image, cmap='jet', aspect='auto')
            self.graph.tick_params(axis='x', which='both', top=False, bottom=config.getboolean('view', 'show_axis_and_grid'))
            self.ylims = np.array([0, 1.1 * display.peak_intensity])
            if self.zoom_on_centroid and display.FWHM is not None:
                self.xlims = display.COG_IMF + 3 * display.FWHM * np.array([-1, 1])
            else:
                self.xlims = np.array([0, display.number_of_pixels - 1])
        elif graph_type_ind == 4:
            if len(display.cast_rays) <= np.inf + MAX_NR_OF_SCATTERPOINTS_PLOTTED:
                size_pts = 100 / np.sqrt(len(display.cast_rays))
                self.graph.scatter(display.cast_rays_x, display.cast_rays_phase, s=size_pts, c=display.cast_rays_col)
                self.graph.scatter(display.pixels_x + display.pixel_size/2, display.phase, s=10+0*size_pts, facecolor=(1,1,1,0), edgecolor=(0,0,0,1), marker='o')
                self.xlims = np.array([0, display.length])
                self.ylims = np.array([-np.pi, np.pi])
                self.graph.tick_params(axis='x', which='both', top=False, bottom=config.getboolean('view', 'show_axis_and_grid'))
                self.graph.grid(config.getboolean('view', 'show_axis_and_grid'), which='both', axis='x')

        self.graph.set_xlim(self.xlims)
        self.graph.set_ylim(self.ylims)

        self.draw()
        print('Display graph updated')
