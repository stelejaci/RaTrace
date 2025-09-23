import sys
import matplotlib
import numpy as np
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import pyqtSignal
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import os, sys
sys.path.append(os.path.abspath('..'))
from datetime import datetime

from gui import canvas_class
from utils.configuration_class import config
from display import display_class, imager_class
from utils import varia


NR_OF_RAYS_TO_PLOT      = config.getint('view', 'nr_of_rays_to_plot')
INTENSITY_SCALER        = config.getfloat('view', 'intensity_scaler')
NR_OF_RAYS_LIST         = [1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,5000,10000,20000,30000,50000,100000,200000,300000,500000,1000000,2000000,3000000,5000000,10000000]
NR_OF_RAYS_TO_PLOT_LIST = [1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,5000,10000,20000,30000,50000,100000]
INTENSITY_SCALER_LIST   = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
RED, WHITE              = '\033[31m', '\033[0m'



class SimulationGuiClass(QtWidgets.QWidget):
    def __init__(self, simulation=None):
        super().__init__()
        self.window_title = 'RaTrace'
        self.setWindowTitle(self.window_title)
        self.setWindowIcon(QtGui.QIcon('../assets/RaTrace_thumbnail.png'))
        self.setGeometry(0, 50, 1800, 1000)

        self.simulation = simulation
        self.nr_of_rays_to_plot = NR_OF_RAYS_TO_PLOT

        self.setup_UI()

    def setup_UI(self):
        self.canvas = canvas_class.CanvasClass(simulation=self.simulation, width=5, height=4, dpi=100)
        self.canvas_toolbar = NavigationToolbar2QT(self.canvas, self)

        self.tab_widget = QtWidgets.QTabWidget()
        self.tab_widget.setEnabled(True)
        self.tab_setup      = SetupTabWidget(      simulation=self.simulation )
        self.tab_simulation = SimulationTabWidget( simulation=self.simulation )
        self.tab_view       = ViewTabWidget(       simulation=self.simulation )
        self.tab_display    = DisplayTabWidget(    simulation=self.simulation )
        self.tab_widget.addTab(self.tab_setup, "Setup")
        self.tab_widget.addTab(self.tab_simulation, "Simulation")
        self.tab_widget.addTab(self.tab_view, "View")
        self.tab_widget.addTab(self.tab_display, "Display")

        # En/disable certain controls if a display or imager is present/absent
        self.display_present = True if self.simulation.displays else False
        self.tab_widget.setTabEnabled(3, self.display_present)
        self.imager_present = any(isinstance(display, display_class.ImagerClass) for display in self.simulation.displays)
        self.tab_view.show_pixels_checkbox.setEnabled(self.imager_present)

        layout_left = QtWidgets.QVBoxLayout()
        layout_left.addWidget(self.canvas,  98)
        layout_left.addWidget(self.canvas_toolbar)

        layout_main = QtWidgets.QHBoxLayout()
        layout_main.addLayout(layout_left,60)
        layout_main.addWidget(self.tab_widget,20)
        self.setLayout(layout_main)

        self.simulation.scene_file_loaded_signal.connect(self.clear_and_update_graphics)
        self.simulation.simulation_run_done_signal.connect(self.update_graphics)

        self.tab_view.redraw_scene_signal.connect(self.update_graphics)
        # self.tab_setup.config_file_selected_signal.connect(self.update_window_title)
        # self.canvas.update_source_progress_signal.connect(self.tab_view.update_progress_bar)

        # Connect each of the sources' signals to the progress bar
        for source in self.simulation.sources:
            source.update_source_progress_signal.connect(self.tab_view.update_progress_bar)
        # self.tab_setup.config_file_selected_signal.connect(self.simulation.load_configuration)
        # self.canvas.update_components()

    @QtCore.pyqtSlot()
    def clear_and_update_graphics(self):
        self.setWindowTitle(f'{self.window_title}  [{self.simulation.scene_file}]  -  {self.simulation.info}')
        filename = os.path.split(self.simulation.scene_file)[1]
        self.tab_setup.scene_file_text.setPlainText(filename)
        self.display_present = True if self.simulation.displays else False
        self.tab_widget.setTabEnabled(3, self.display_present)
        self.imager_present = any(isinstance(display, imager_class.ImagerClass) for display in self.simulation.displays)
        self.tab_view.show_pixels_checkbox.setEnabled(self.imager_present)
        self.canvas.clear()
        self.tab_display.canvas_display.graph.cla()
        self.update_graphics()

    @QtCore.pyqtSlot()
    def update_graphics(self):
        print('Updating canvas and display graph ...')
        self.canvas.set_nr_of_rays_to_plot(self.tab_view.nr_of_rays_to_plot)
        self.canvas.update_entire_scene()
        self.tab_display.update_graphs()

    # @QtCore.pyqtSlot(bool)
    # def update_color_code_intensities(self, value):
    #     self.canvas.set_color_code_ray_intensities(value)

    # @QtCore.pyqtSlot(bool)
    # def update_show_elements_normals(self, value):
    #     self.canvas.set_show_elements_normals(value)

    # @QtCore.pyqtSlot(str)
    # def update_window_title(self, filename_full):
    #     self.setWindowTitle(f'{self.window_title}  [{filename_full}]')


class SetupTabWidget(QtWidgets.QWidget):
    # config_file_selected_signal = pyqtSignal(str)

    def __init__(self, simulation):
        super().__init__()
        self.simulation = simulation

        self.scene_file_label = QtWidgets.QLabel("Scene file:")
        self.scene_file_browse_button = QtWidgets.QPushButton("Browse ...")
        self.scene_file_browse_button.clicked.connect(self.browse_scene_file)
        self.scene_file_reload_button = QtWidgets.QPushButton("Reload")
        self.scene_file_reload_button.clicked.connect(self.reload_scene_file)
        layout_scene_file   = QtWidgets.QHBoxLayout()
        layout_scene_file.addWidget(self.scene_file_label)
        layout_scene_file.addWidget(self.scene_file_browse_button)
        layout_scene_file.addWidget(self.scene_file_reload_button)

        self.scene_file_text = QtWidgets.QPlainTextEdit("Selected and loaded scene file")
        self.scene_file_text.setMaximumHeight(self.scene_file_label.sizeHint().height()*2)
        self.scene_file_text.setReadOnly(True)

        self.start_simulation_after_loading_scene_checkbox = QtWidgets.QCheckBox("Start simulation after loading scene")
        self.start_simulation_after_loading_scene_checkbox.setChecked(config.getboolean('scenes', 'start_simulation_after_loading_scene'))
        self.start_simulation_after_loading_scene_checkbox.stateChanged.connect(self.update_start_simulation_after_loading_scene)

        self.reset_axis_after_loading_scene_checkbox = QtWidgets.QCheckBox("Reset axis after loading scene")
        self.reset_axis_after_loading_scene_checkbox.setChecked(config.getboolean('scenes', 'reset_axis_after_loading_scene'))
        self.reset_axis_after_loading_scene_checkbox.stateChanged.connect(self.update_reset_axis_after_loading_scene)

        self.model_params_label = QtWidgets.QLabel("Model parameters")
        self.model_params_tree = QtWidgets.QTreeWidget()
        self.model_params_tree.setColumnCount(2)
        self.model_params_tree.setHeaderLabels(["Element parameters", "Value"])
        tree_items = []
        tree_item_1 = QtWidgets.QTreeWidgetItem(["Lens"])
        tree_item_1_1 = QtWidgets.QTreeWidgetItem(["Focal distance", "1.8"])
        tree_item_1_2 = QtWidgets.QTreeWidgetItem(["Position", "[64.5, 0.0, 178.4]"])
        tree_item_1.addChild(tree_item_1_1)
        tree_item_1.addChild(tree_item_1_2)
        tree_items.append(tree_item_1)
        tree_item_2 = QtWidgets.QTreeWidgetItem(["Imager"])
        tree_item_2_1 = QtWidgets.QTreeWidgetItem(["Pixel size", "6.0"])
        tree_item_2_2 = QtWidgets.QTreeWidgetItem(["Position", "[74.5, 0.0, 188.4]"])
        tree_item_2.addChild(tree_item_2_1)
        tree_item_2.addChild(tree_item_2_2)
        tree_items.append(tree_item_2)
        self.model_params_tree.insertTopLevelItems(0, tree_items)
        layout_controls = QtWidgets.QVBoxLayout()
        layout_controls.addLayout(layout_scene_file)
        layout_controls.addWidget(self.scene_file_text)
        layout_controls.addWidget(self.start_simulation_after_loading_scene_checkbox)
        layout_controls.addWidget(self.reset_axis_after_loading_scene_checkbox)
        layout_controls.addStretch(1)
        layout_controls.addWidget(self.model_params_label)
        layout_controls.addWidget(self.model_params_tree,15)
        self.setLayout(layout_controls)

    def browse_scene_file(self):
        self.simulation.scene_file, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', config.get('scenes', 'scenes_folder'), 'Scene files (*.py)')
        print(f'Scene file selected: {self.simulation.scene_file}')
        self.simulation.load_scene(self.simulation.scene_file)
        config.set('scenes', 'scenes_folder', os.path.dirname(self.simulation.scene_file))
        config.set('scenes', 'scene_file',    os.path.basename(self.simulation.scene_file))

    def reload_scene_file(self):
        print(f'Scene file reloaded: {self.simulation.scene_file}')
        self.simulation.load_scene(self.simulation.scene_file)

    def update_start_simulation_after_loading_scene(self, value):
        config.set('scenes', 'start_simulation_after_loading_scene', str(int(value == 2)))

    def update_reset_axis_after_loading_scene(self, value):
        config.set('scenes', 'reset_axis_after_loading_scene', str(int(value == 2)))


class SimulationTabWidget(QtWidgets.QWidget):
    def __init__(self, simulation=None):
        super().__init__()
        self.simulation = simulation

        self.nrofrays_label = QtWidgets.QLabel(f'Number of rays per light source: {self.simulation.nr_of_rays_per_source}')
        self.nrofrays_slider = QtWidgets.QSlider()
        self.nrofrays_slider.setRange(0,len(NR_OF_RAYS_LIST)-1)
        self.nrofrays_slider.setSingleStep(1)
        self.nrofrays_slider.setPageStep(3)
        slider_pos = np.where(np.array(NR_OF_RAYS_LIST) == self.simulation.nr_of_rays_per_source)[0][0]
        self.nrofrays_slider.setValue(slider_pos)
        self.nrofrays_slider.setOrientation(QtCore.Qt.Horizontal)
        self.nrofrays_slider.setTickPosition(QtWidgets.QSlider.TickPosition.TicksAbove)
        self.nrofrays_slider.valueChanged.connect(self.update_nr_of_rays_per_source)
        self.sampling_space_label = QtWidgets.QLabel("Sampling space")
        self.sampling_3D_radioButton = QtWidgets.QRadioButton("3D")
        self.run_simulation_button = QtWidgets.QPushButton("Start simulation")
        self.run_simulation_button.clicked.connect(self.simulation.run)
        self.use_phase_information_checkbox = QtWidgets.QCheckBox("Use phase information")
        self.use_phase_information_checkbox.setChecked(config.getboolean('simulation', 'use_phase_information'))
        self.use_phase_information_checkbox.stateChanged.connect(self.update_use_phase_information)
        self.items_are_ordered_checkbox = QtWidgets.QCheckBox("Items are ordered (faster, but items must be correctly ordered)")
        self.items_are_ordered_checkbox.setChecked(config.getboolean('simulation', 'items_are_ordered'))
        self.items_are_ordered_checkbox.stateChanged.connect(self.update_items_are_ordered)
        # self.items_are_ordered_checkbox.setEnabled(False)
        self.generate_reflected_rays_checkbox = QtWidgets.QCheckBox("Generate reflected rays")
        self.generate_reflected_rays_checkbox.setChecked(config.getboolean('simulation', 'generate_reflected_rays'))
        self.generate_reflected_rays_checkbox.stateChanged.connect(self.update_generate_reflected_rays)
        self.generate_reflected_rays_checkbox.setEnabled(False)
        self.export_ray_data_button = QtWidgets.QPushButton("Export ray data ...")
        self.export_ray_data_button.clicked.connect(self.export_ray_data_button_clicked)
        self.add_timestamp_to_export_file_checkbox = QtWidgets.QCheckBox("Add timestamp to file name")
        self.add_timestamp_to_export_file_checkbox.stateChanged.connect(self.update_add_timestamp_to_export_file)
        self.add_timestamp_to_export_file_checkbox.setChecked(config.getboolean('simulation', 'add_timestamp_to_export_file'))

        font = QtGui.QFont()
        font.setPointSize(16)
        self.run_simulation_button.setFont(font)
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setProperty("value", 0)
        self.simulation.raytracer.raytracing_progress_signal.connect(self.update_progress_bar)
        font.setPointSize(12)
        self.export_ray_data_button.setFont(font)

        layout_simulation = QtWidgets.QVBoxLayout()
        layout_simulation.addWidget(self.nrofrays_label)
        layout_simulation.addWidget(self.nrofrays_slider)
        layout_simulation.addWidget(self.use_phase_information_checkbox)
        layout_simulation.addWidget(self.items_are_ordered_checkbox)
        layout_simulation.addWidget(self.generate_reflected_rays_checkbox)
        layout_simulation.addStretch(1)
        layout_simulation.addWidget(self.run_simulation_button)
        layout_simulation.addWidget(self.progressbar)
        layout_export = QtWidgets.QHBoxLayout()
        layout_export.addWidget(self.export_ray_data_button)
        layout_export.addWidget(self.add_timestamp_to_export_file_checkbox)
        layout_simulation.addLayout(layout_export)
        self.setLayout(layout_simulation)

    @QtCore.pyqtSlot(float)
    def update_progress_bar(self, progress_pct):
        self.progressbar.setProperty("value", 100*progress_pct)
        # print(f'Updating progress bar {100*progress_pct:.0f}%')

    def update_nr_of_rays_per_source(self, value):
        self.simulation.nr_of_rays_per_source = NR_OF_RAYS_LIST[value]
        self.nrofrays_label.setText(f'Nr of rays: {self.simulation.nr_of_rays_per_source}')
        config.set('simulation', 'nr_of_rays', str(self.simulation.nr_of_rays_per_source))

    def update_use_phase_information(self, value):
        config.set('simulation', 'use_phase_information', str(int(value==2)))

    def update_items_are_ordered(self, value):
        config.set('simulation', 'items_are_ordered', str(int(value==2)))

    def update_generate_reflected_rays(self, value):
        config.set('simulation', 'generate_reflected_rays', str(int(value==2)))

    def update_add_timestamp_to_export_file(self, value):
        config.set('simulation', 'add_timestamp_to_export_file', str(int(value==2)))

    def export_ray_data_button_clicked(self):
        # Only save ray data when there is available, i.e. when a simulation is run
        if not self.simulation.sources[0].rays:
            print(f'No ray data available yet to write to file')
            return

        # Request the export file name by the user
        folder = config.get('simulation', 'export_folder')
        filename = os.path.join(folder, 'ray_data.txt')
        filename = QtWidgets.QFileDialog.getSaveFileName(None, f'Select file to export ray data to', filename, 'Text files (*.txt)')[0]
        print(f'File selected: {filename}')

        # Update the config file
        config.set('simulation', 'export_folder', os.path.dirname(filename))

        # Finally, write the ray data to the export file
        try:
            add_timestamp_to_export_file = config.getboolean('simulation', 'add_timestamp_to_export_file')
            varia.output_source_data_in_text_file(filename, self.simulation.sources, add_timestamp_to_export_file)
            print(f'Ray data exported to: {filename}')
        except:
            print(f'{RED}Ray data could NOt be exported to: {filename}{WHITE}')



class ViewTabWidget(QtWidgets.QWidget):
    redraw_scene_signal = pyqtSignal()

    def __init__(self, simulation=None):
        super().__init__()
        self.simulation = simulation
        self.nr_of_rays_to_plot = NR_OF_RAYS_TO_PLOT
        self.intensity_scaler = INTENSITY_SCALER

        self.nrofrays_label = QtWidgets.QLabel(f'Nr of plotted rays: {NR_OF_RAYS_TO_PLOT}')
        self.nrofrays_slider = QtWidgets.QSlider()
        self.nrofrays_slider.setRange(0,len(NR_OF_RAYS_TO_PLOT_LIST)-1)
        self.nrofrays_slider.setSingleStep(1)
        self.nrofrays_slider.setPageStep(3)
        self.nrofrays_slider.setOrientation(QtCore.Qt.Horizontal)
        slider_pos = np.where(np.array(NR_OF_RAYS_TO_PLOT_LIST)==NR_OF_RAYS_TO_PLOT)[0][0]
        self.nrofrays_slider.setValue(slider_pos)
        self.nrofrays_slider.setTickPosition(QtWidgets.QSlider.TickPosition.TicksAbove)
        self.nrofrays_slider.valueChanged.connect(self.update_nr_of_rays_per_source_to_plot)

        self.intensity_scaler_label = QtWidgets.QLabel(f'Intensity scaler: {INTENSITY_SCALER}')
        self.intensity_scaler_slider = QtWidgets.QSlider()
        self.intensity_scaler_slider.setRange(0,len(INTENSITY_SCALER_LIST)-1)
        self.intensity_scaler_slider.setSingleStep(1)
        self.intensity_scaler_slider.setPageStep(3)
        self.intensity_scaler_slider.setOrientation(QtCore.Qt.Horizontal)
        slider_pos = np.where(np.array(INTENSITY_SCALER_LIST)==INTENSITY_SCALER)[0][0]
        self.intensity_scaler_slider.setValue(slider_pos)
        self.intensity_scaler_slider.setTickPosition(QtWidgets.QSlider.TickPosition.TicksAbove)
        self.intensity_scaler_slider.valueChanged.connect(self.update_intensity_scaler)

        self.viewport_label = QtWidgets.QLabel("Viewport")
        self.viewport_dropdown = QtWidgets.QComboBox()
        self.viewport_dropdown.addItems(["Overview", "Surface", "Laser beam", "Lens", "Imager"])

        self.auto_redraw_checkbox = QtWidgets.QCheckBox("Auto redraw when simulation is done")
        self.auto_redraw_checkbox.setChecked(True)
        self.show_axis_and_grid_checkbox = QtWidgets.QCheckBox("Show axis and grid")
        self.show_axis_and_grid_checkbox.setChecked(config.getboolean('view', 'show_axis_and_grid'))
        self.show_axis_and_grid_checkbox.stateChanged.connect(self.update_show_axis_and_grid)
        self.intensity_coded_colors_checkbox = QtWidgets.QCheckBox("Intensity-coded ray colors")
        self.intensity_coded_colors_checkbox.setChecked(config.getboolean('view', 'intensity_coded_colors'))
        self.intensity_coded_colors_checkbox.stateChanged.connect(self.update_intensity_coded_colors)
        self.show_elements_properties_checkbox = QtWidgets.QCheckBox("Show element properties")
        self.show_elements_properties_checkbox.setChecked(config.getboolean('view', 'show_elements_properties') )
        self.show_elements_properties_checkbox.stateChanged.connect(self.update_show_elements_properties_checkbox)
        self.background_black_checkbox = QtWidgets.QCheckBox("Black background")
        self.background_black_checkbox.setChecked(config.getboolean('view', 'black_background'))
        self.background_black_checkbox.stateChanged.connect(self.update_background_black_checkbox)
        self.show_pixels_checkbox = QtWidgets.QCheckBox("Show pixels (slow!)")
        self.show_pixels_checkbox.setChecked(config.getboolean('view', 'show_pixels'))
        self.show_pixels_checkbox.stateChanged.connect(self.update_show_pixels_checkbox)
        self.show_non_colliding_rays_checkbox = QtWidgets.QCheckBox("Show non-colliding rays")
        self.show_non_colliding_rays_checkbox.setChecked(config.getboolean('view', 'show_noncolliding_rays'))
        self.show_non_colliding_rays_checkbox.stateChanged.connect(self.update_show_noncolliding_rays)

        self.redraw_scene_button = QtWidgets.QPushButton("Redraw scene")
        self.redraw_scene_button.clicked.connect(self.redraw_scene)
        font = QtGui.QFont()
        font.setPointSize(16)
        self.redraw_scene_button.setFont(font)
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setProperty("value", 0)
        layout_view = QtWidgets.QVBoxLayout()
        layout_view.addWidget(self.nrofrays_label)
        layout_view.addWidget(self.nrofrays_slider)
        layout_view.addWidget(self.intensity_scaler_label)
        layout_view.addWidget(self.intensity_scaler_slider)
        layout_view.addStretch()
        layout_view.addWidget(self.viewport_label)
        layout_view.addWidget(self.viewport_dropdown)
        layout_view.addStretch()
        layout_view.addWidget(self.auto_redraw_checkbox)
        layout_view.addWidget(self.show_axis_and_grid_checkbox)
        layout_view.addWidget(self.intensity_coded_colors_checkbox)
        layout_view.addWidget(self.show_elements_properties_checkbox)
        layout_view.addWidget(self.background_black_checkbox)
        layout_view.addWidget(self.show_pixels_checkbox)
        layout_view.addWidget(self.show_non_colliding_rays_checkbox)
        layout_view.addStretch(1)
        layout_view.addWidget(self.redraw_scene_button)
        layout_view.addWidget(self.progressbar)
        self.setLayout(layout_view)

    def update_nr_of_rays_per_source_to_plot(self, value):
        self.nr_of_rays_to_plot = NR_OF_RAYS_TO_PLOT_LIST[value]
        self.nrofrays_label.setText(f'Nr of plotted rays: {self.nr_of_rays_to_plot}')
        config.set('view', 'nr_of_rays_to_plot', str(self.nr_of_rays_to_plot))

    def update_intensity_scaler(self, value):
        self.intensity_scaler = INTENSITY_SCALER_LIST[value]
        self.intensity_scaler_label.setText(f'Intensity scaler: {self.intensity_scaler}')
        config.set('view', 'intensity_scaler', str(self.intensity_scaler))

    def update_intensity_coded_colors(self, value):
        config.set('view', 'intensity_coded_colors', str(int(value==2)))
        self.redraw_scene()

    def update_show_axis_and_grid(self, value):
        config.set('view', 'show_axis_and_grid', str(int(value==2)))
        self.redraw_scene()

    def update_show_elements_properties_checkbox(self, value):
        config.set('view', 'show_elements_properties', str(int(value==2)))     # a checkbox has either state 0 or 2
        self.redraw_scene()

    def update_show_pixels_checkbox(self, value):
        config.set('view', 'show_pixels', str(int(value==2)))
        self.redraw_scene()

    def update_show_noncolliding_rays(self, value):
        config.set('view', 'show_noncolliding_rays', str(int(value==2)))
        self.redraw_scene()

    def update_background_black_checkbox(self, value):
        config.set('view', 'black_background', str(int(value==2)))
        if config.getboolean('view', 'black_background'):
            canvas_class.PLOT_COLOR_BLACK = 'darkgrey'
        else:
            canvas_class.PLOT_COLOR_BLACK = 'black'
        self.redraw_scene()

    def redraw_scene(self):
        print(f'Redrawing scene with {self.nr_of_rays_to_plot} plotted rays per source')
        # self.redraw_scene_signal.emit(self.nr_of_rays_to_plot)
        self.redraw_scene_signal.emit()

    @QtCore.pyqtSlot(float)
    def update_progress_bar(self, progress_pct):
        self.progressbar.setProperty("value", 100*progress_pct)




class DisplayTabWidget(QtWidgets.QWidget):
    def __init__(self, simulation=None):
        super().__init__()
        self.simulation = simulation
        self.canvas_display = canvas_class.CanvasDisplayClass(simulation=self.simulation, width=5, height=4, dpi=100)
        # self.canvas_display.figure.tight_layout()
        self.canvas_display_toolbar = NavigationToolbar2QT(self.canvas_display, self)

        self.graph_types = ["Scatterplot 1D", "Scatterplot 2D", "Image 2D", "Image centroid", "Phase plot"]
        self.graph_type_ind = 1
        self.graphtype_dropdown = QtWidgets.QComboBox()
        self.graphtype_dropdown.addItems(self.graph_types)
        self.graphtype_dropdown.setCurrentIndex(self.graph_type_ind)
        self.graphtype_dropdown.currentIndexChanged.connect(self.graphtype_dropdown_changed)
        self.greyscale_checkbox = QtWidgets.QCheckBox("Greyscale image")
        self.greyscale_checkbox.stateChanged.connect(self.greyscale_checkbox_changed)
        self.greyscale_checkbox.setChecked(self.canvas_display.greyscale_mode)
        self.zoom_on_centroid_checkbox = QtWidgets.QCheckBox("Zoom on centroid")
        self.zoom_on_centroid_checkbox.setChecked(self.canvas_display.zoom_on_centroid)
        self.zoom_on_centroid_checkbox.stateChanged.connect(self.zoom_on_centroid_checkbox_changed)
        self.export_imager_data_button = QtWidgets.QPushButton("Export imager data ...")
        self.export_imager_data_button.clicked.connect(self.export_imager_data_button_clicked)
        self.add_timestamp_to_export_file_checkbox = QtWidgets.QCheckBox("Add timestamp to file name")
        self.add_timestamp_to_export_file_checkbox.stateChanged.connect(self.update_add_timestamp_to_export_file)
        self.add_timestamp_to_export_file_checkbox.setChecked(config.getboolean('simulation', 'add_timestamp_to_export_file'))

        font = QtGui.QFont()
        font.setPointSize(12)
        self.export_imager_data_button.setFont(font)

        layout_main = QtWidgets.QVBoxLayout()
        layout_main.addWidget(self.canvas_display, 50)
        layout_main.addWidget(self.canvas_display_toolbar, 5)
        layout_sub1 = QtWidgets.QHBoxLayout()
        layout_sub1.addWidget(self.graphtype_dropdown)
        layout_sub1.addWidget(self.greyscale_checkbox)
        layout_sub1.addWidget(self.zoom_on_centroid_checkbox)
        layout_main.addLayout(layout_sub1)
        layout_main.addStretch(50)
        layout_export = QtWidgets.QHBoxLayout()
        layout_export.addWidget(self.export_imager_data_button)
        layout_export.addWidget(self.add_timestamp_to_export_file_checkbox)
        layout_main.addLayout(layout_export)
        self.setLayout(layout_main)

    def update_graphs(self):
        print(f'Updating graphs')
        self.canvas_display.update_graphs(graph_type_ind = self.graph_type_ind)

    def graphtype_dropdown_changed(self, text):
        self.graph_type_ind = self.graphtype_dropdown.currentIndex()  # Get selected index
        print("Graph type index changed to: " + str(self.graph_type_ind) + " (" + self.graph_types[self.graph_type_ind] + ")")
        self.update_graphs()

    def greyscale_checkbox_changed(self, value):
        self.canvas_display.greyscale_mode = value
        print("Image mode changed: " + str(value))
        self.update_graphs()

    def zoom_on_centroid_checkbox_changed(self, value):
        self.canvas_display.zoom_on_centroid = value
        print("Zoom on centroid: " + str(value))
        self.update_graphs()

    def update_add_timestamp_to_export_file(self, value):
        config.set('simulation', 'add_timestamp_to_export_file', str(int(value==2)))

    def export_imager_data_button_clicked(self):
        # Only save ray data when there is available, i.e. when a simulation is run
        if not self.simulation.sources[0].rays:     # Only check the first source for now
            print(f'No cast ray data available yet to write to file')
            return
        elif not isinstance(self.simulation.displays[0], imager_class.ImagerClass):   # Only check the first display for now
            print(f'No imager present in the simulation')
            return

        # Request the export file name by the user
        folder = config.get('simulation', 'export_folder')
        filename = os.path.join(folder, 'imager_data.txt')
        filename = QtWidgets.QFileDialog.getSaveFileName(None, f'Select file to export imager data to', filename, 'Text files (*.txt)')[0]
        print(f'File selected: {filename}')

        # Update the config file
        folder = os.path.dirname(filename)
        config.set('simulation', 'export_folder', folder)

        # Add timestamp data to the filename, is needed
        if config.getboolean('simulation', 'add_timestamp_to_export_file'):
            current_time = datetime.now()
            timestamp_string = current_time.strftime("%Y%m%d_%H%M%S")
            filename_base = os.path.splitext(filename)[0]
            filename = filename_base + '_' + timestamp_string + '.txt'
            print(f'Adding timestamp information to filename: {filename}')

        # Finally, write the ray data to the export file
        try:
            varia.output_imager_data_in_text_file(filename, self.simulation.displays)
            print(f'Imager data exported to: {filename}')
        except:
            print(f'{RED}Imager data could NOt be exported to: {filename}{WHITE}')


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    simulation_GUI = SimulationGuiClass()
    simulation_GUI.show()
    sys.exit(app.exec_())
