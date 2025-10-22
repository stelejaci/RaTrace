from PyQt5 import QtCore, QtWidgets
import sys,os
from gui import simulation_gui_class
from raytracer import simulation_class


f_lens = [18, 16, 14, 12, 10]
p_BD   = [30, 40, 50, 60, 70]


def run_gui_iteration(app, iter):
    simulation = simulation_class.SimulationClass()
    simulation_GUI = simulation_gui_class.SimulationGuiClass(simulation)
    simulation_GUI.show()

    # Load a scene file that takes an argument
    scene_name = 'scene_20_SCENE_WITH_PARAMETERS_IN_A_LOOP'
    simulation.load_scene(f"../scenes/{scene_name}.py", param=[f_lens[iter], p_BD[iter]])

    # Set the number of rays per source to raytrace. This overrides the number in the config.ini file
    simulation.set_nr_of_rays_per_source(10)

    # Finally, raytrace the entire thing
    simulation.run()

    # Taking a screenshot and closing figure
    screenshot_filename = f'../scenes/{scene_name}_iter{iter:04d}.png'
    QtCore.QTimer.singleShot(2000, lambda: take_screenshot(app, simulation_GUI, screenshot_filename))
    QtCore.QTimer.singleShot(3000, simulation_GUI.close)


def take_screenshot(app, window, filename):
    pixmap = window.grab()
    pixmap.save(filename)
    print(f"Screenshot saved: {filename}")



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    for iter in range(len(f_lens)):
        run_gui_iteration(app, iter)

        # Process events until the current iteration is complete
        while app.activeWindow():
            app.processEvents()
            QtCore.QThread.msleep(100)

    sys.exit(app.exec_())
