from surface import Surface
import matplotlib.pyplot as plt
from math import *


class CylinderSurface(Surface):

    def __init__(self, r, h):
        self.r = r
        self.h = h

    def is_in(self, photon):
        point_interaction = photon.get_last_position()
        is_in_surface = True
        is_in_surface *= (point_interaction[2] < self.h) and (point_interaction[2] > 0)
        is_in_surface *= point_interaction[0] ** 2 + point_interaction[1] ** 2 < self.r ** 2

        return is_in_surface

    def plot_surface(self, ax):

        dt = 0.01
        t = 0
        x = []
        y = []
        z1 = []
        z2 = []
        r = self.r
        h = self.h
        while t < 2*pi:
            x.append(r * cos(t))
            y.append(r * sin(t))
            z1.append(0)
            z2.append(h)
            t += dt

        ax.plot(x, y, z1)
        ax.plot(x, y, z2)
