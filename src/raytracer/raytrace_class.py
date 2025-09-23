import time

import numpy as np
from PyQt5.QtCore import pyqtSignal, QObject
import os, sys
sys.path.append(os.path.abspath('..'))

from utils import varia
from utils.varia import mm, µm, nm, deg, X, Y
from light import light_class
from elements import element_class
from utils.configuration_class import config


EPSILON = 1e-6*mm


class RaytracerClass(QObject):
    raytracing_progress_signal = pyqtSignal(float)

    def __init__(self):
        super().__init__()

    def raytrace_source(self, source, elements):
        print(f'Raytracing source {source.ID+1}/{light_class.LightSourceClass.nr_of_sources} ({source.nr_of_rays} rays)')

        nr_of_rays = source.nr_of_rays

        for i_ray in range(nr_of_rays):
            ray = source.rays[i_ray]
            self.raytrace_ray_and_children(ray, source, elements)
            # Update the progress bar every 1%, but use an epsilon approach, otherwise it skips steps because of floating point errors slipping through
            if np.abs(np.mod((i_ray+1)/nr_of_rays, 0.01) - 0.01) < EPSILON:
                varia.print_progress(i_ray, nr_of_rays)
                self.raytracing_progress_signal.emit((i_ray+1)/nr_of_rays)

    def raytrace_ray_and_children(self, ray, source, elements):
            self.raytrace_ray(ray, source, elements)
            for ID_child in ray.ID_children:
                self.raytrace_ray_and_children(source.rays[ID_child], source, elements)

    def raytrace_ray(self, ray, source, elements):
        if not ray.is_active:
            return
        if config.getboolean('simulation', 'items_are_ordered'):
            self.raytrace_ray_with_sorted_elements(ray,source,elements)
        else:
            self.raytrace_ray_full(ray, source, elements)

    def raytrace_ray_with_sorted_elements(self, ray, source, elements):
        # Keep track of the elements over which to loop, omit elements that rays are already passed through, if needed
        i_element_start = 0  # Start from the beginning, loping over all elements
        if isinstance(ray.source_element, element_class.ElementClass)  and  config.getboolean('simulation', 'items_are_ordered'):
            i_element_start = ray.source_element.ID   # If the checkbox is set, take the last collission element as the starting element

        # Loop over the elements, look for a collision
        for element in elements[i_element_start:]:
            [p1, length, p1_element_rel] = element.check_collision(ray)  # Retrieve the collision point of the ray on the element and length of the terminated ray

            # The ray is not hitting the element
            if p1 is None:
                continue

            # The ray hits an element, in that case stop looping over the elements.
            # This speeds-up the algorithm because it does not look for other collisions.
            # Since the elements should be ordered by the user, this is valid is most situations
            if length > EPSILON:  # length>eps: avoid collision with the current line segment
                ray.p1 = p1
                ray.element_hit = element
                ray.p1_element_rel = p1_element_rel
                break

        # IF a colliding point and element is found, generate (a) propagating ray(s)
        if ray.element_hit:
            new_rays = ray.element_hit.propagate_ray(ray)
            ray.is_active = False

            # Propagating the ray could result in a bunch of new rays, a single ray, or no new rays at all
            if isinstance(new_rays, list):
                for new_ray in new_rays:
                    new_ray.ID = source.nr_of_rays
                    source.add_ray(new_ray)
                    ray.ID_children.append(new_ray.ID)
            elif new_rays is not None:
                new_rays.ID = source.nr_of_rays
                source.add_ray(new_rays)
                ray.ID_children.append(new_rays.ID)
            else:
                # Not needed, but mentioned explicitly, just to clarify that in this case no new rays are returned
                # E.g. when rays are casted onto a display
                pass

    def raytrace_ray_full(self, ray, source, elements):
        # Loop over ALL the elements, look for a collision
        for element in elements:
            [p1, length, p1_element_rel] = element.check_collision(ray)  # Retrieve the collision point of the ray on the element and length of the terminated ray

            # The ray is not hitting the element
            if p1 is None:
                continue

            if length > EPSILON:    # length>eps: avoid collision with the current line segment
                if (ray.length is None) or (length < ray.length):
                    ray.p1 = p1
                    ray.element_hit = element
                    ray.p1_element_rel = p1_element_rel

        # IF a colliding point and element is found, generate (a) propagating ray(s)
        if ray.element_hit:
            new_rays = ray.element_hit.propagate_ray(ray)
            ray.is_active = False

            # Propagating the ray could result in a bunch of new rays, a single ray, or no new rays at all
            if isinstance(new_rays, list):
                for new_ray in new_rays:
                    new_ray.ID = source.nr_of_rays
                    source.add_ray(new_ray)
                    ray.ID_children.append(new_ray.ID)
            elif new_rays is not None:
                new_rays.ID = source.nr_of_rays
                source.add_ray(new_rays)
                ray.ID_children.append(new_rays.ID)
            else:
                # Not needed, but mentioned explicitly, just to clarify that in this case no new rays are returned
                # E.g. when rays are casted onto a display
                pass
        # print(ray)
        # time.sleep(1)
