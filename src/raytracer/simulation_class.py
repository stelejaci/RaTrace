from PyQt5.QtCore import pyqtSignal, QObject, pyqtSlot
import importlib.util
import time
import pathlib
import inspect
import os
from raytracer  import raytrace_class
from elements   import element_class
from light      import light_class
from display    import display_class, imager_class
from utils.configuration_class import config

INITIAL_NR_OF_RAYS = config.getint('simulation', 'nr_of_rays')
RED, WHITE = '\033[31m', '\033[0m'




class SimulationClass(QObject):
    simulation_run_done_signal = pyqtSignal()
    scene_file_loaded_signal = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.raytracer = raytrace_class.RaytracerClass()
        self.nr_of_rays_per_source = INITIAL_NR_OF_RAYS
        self.reset_simulation()

    def reset_simulation(self):
        self.reset_items()
        element_class.ElementClass.nr_of_elements  = 0
        light_class.LightSourceClass.nr_of_sources = 0
        display_class.DisplayClass.nr_of_displays  = 0

    def reset_items(self):
        self.sources  = list()
        self.elements = list()
        self.displays = list()
        self.info     = None

    @pyqtSlot(str)
    def load_scene(self, scene_file, param=None):
        path = pathlib.Path(scene_file)
        if not path.suffix:
            print(f'No scene file given: {scene_file}')
            return

        # First, reset the simulation
        self.reset_simulation()

        # Loading the scene file as a python module
        self.scene_file = os.path.normpath(scene_file)
        module_name = self.scene_file.split('/')[-1].split('.')[0]
        spec = importlib.util.spec_from_file_location(module_name, self.scene_file)
        module = importlib.util.module_from_spec(spec)

        try:
            spec.loader.exec_module(module)
        except:
            print(f'Scene file does not exist: {scene_file}')
            self.scene_loaded_successfully = False
            return

        # Loading the scene itself, by running the "load_scene" method in the file
        try:
            func = module.load_scene
            sig = inspect.signature(func)
            if len(sig.parameters) == 0:    # If load_scene takes NO parameter, as for single static simulations, run without arguments
                items = module.load_scene()     # Run the "load_scene" method in each scene file
            else:       # If load_scene DOES take a parameter, as for multiple dynamic simulations in a loop, run with the param argument
                items = module.load_scene(param=param)
            self.scene_loaded_successfully = True
        except (RuntimeError, TypeError, NameError, OSError, AttributeError) as err:
            print(f'{RED}Error in the scene file: {self.scene_file}{WHITE}')
            print(f'{RED}Error message: {err}{WHITE}')
            self.scene_loaded_successfully = False
            return

        if self.classify_items(items):  # Give every item in the different item lists its proper classification and ID, and check if the items are set correctly
            print(f'Scene succesfully loaded: {self.scene_file}')
            print(self)
            self.scene_file_loaded_signal.emit()
        else:
            print(f'{RED}Scene NOT loaded, one of the items is not valid{WHITE}')

        if config.getboolean('scenes', 'start_simulation_after_loading_scene'):
            self.run()

        return

    # ToDo: items have numbers following the order of creation. Renumber the items following their order of listing
    def classify_items(self, items):
        self.reset_items()
        for item in items:
            if isinstance(item, light_class.LightSourceClass):
                item.ID = len(self.sources)
                self.sources.append(item)
            elif isinstance(item, display_class.DisplayClass):
                item.ID = len(self.displays)
                self.displays.append(item)
            elif isinstance(item, element_class.ElementClass):
                item.ID = len(self.elements)
                self.elements.append(item)
            elif isinstance(item, str):
                self.info = item
            else:
                print(f'{RED}Item {item} is not recognised as a valid type{WHITE}')
                return False
        return True

    def run(self):
        if not self.scene_loaded_successfully:
            return

        # Explicitly reset (sources and) displays, otherwise the cast rays list keeps adding
        self.initialise_sources()
        self.initialise_displays()

        # Raytrace all beams over all elements and displays, treat a display as an element.
        for source in self.sources:
            if source.is_virtual: continue
            start_time = time.time()
            self.raytracer.raytrace_source(source=source, elements=self.elements + self.displays)   # Treat displays as a kind of beam dump too
            print(f'Processing time for ray tracing: {time.time()-start_time:.1f} s')

        # Process display information
        for display in self.displays:
            display.process_cast_rays()
            if isinstance(display, imager_class.ImagerClass):
                display.process_image()

        # Process virtual rays
        for source in self.sources:
            if not source.is_virtual: continue
            if ('imager 0' in source.origin)  and  self.displays:
                source.generate_ray(self.displays[0])
            # self.raytracer.raytrace_source(source=source, elements=self.elements + self.displays)   # Treat displays as a kind of beam dump too
            self.raytracer.raytrace_source(source=source, elements=self.elements)   # Treat displays as a kind of beam dump too
            source.process()

        self.simulation_run_done_signal.emit()

    def set_nr_of_rays_per_source(self, nr_of_rays_per_source):
        self.nr_of_rays_per_source = nr_of_rays_per_source

    def initialise_sources(self):
        for source in self.sources:
            if source.is_virtual: continue
            source.reset()      # This was grayed out, why? If not done this, rays will add up
            source.generate_rays(N_rays=self.nr_of_rays_per_source)

    def initialise_displays(self):
        for display in self.displays:
            display.reset()

    def __str__(self):
        txt = ''
        txt += 'Simulation\n'
        txt += u' \u21b3 Sources\n'
        for source in self.sources:
            txt += u'    \u21b3 ' + source.__str__() + '\n'
        txt += u' \u21b3 Elements\n'
        for element in self.elements:
            txt += u'    \u21b3 ' + element.__str__() + '\n'
        txt += u' \u21b3 Displays\n'
        for display in self.displays:
            txt += u'    \u21b3 ' + display.__str__() + '\n'

        return txt


# if __name__ == "__main__":
#     tmp = light_functions.SourceClass()
#     print('done')