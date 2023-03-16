
from photon import Photon
from source import Source
import numpy as np
import matplotlib.pyplot as plt


class HalfDiskSource(Source):

    def __init__(self, source_radius):
        self.r = [0, source_radius]
        self.phi = [-np.pi / 2, np.pi / 2]

    def born_photon(self, energy):
        point_born = [0, 0, 0]
        r_born = np.sqrt(np.random.uniform(0, 1)) * self.r[1]
        phi_born = np.random.uniform(self.phi[0], self.phi[1])
        point_born[0] = r_born * np.cos(phi_born)
        point_born[1] = r_born * np.sin(phi_born)
        point_born[2] = 0
        return Photon(point_born, energy)

    def plot_source(self, ax, d3='3d'):
        phi = np.linspace(-np.pi / 2, np.pi / 2, num=50)
        X = self.r[1] * np.cos(phi)
        Y = self.r[1] * np.sin(phi)
        Z = np.zeros(50)
        x = [0, 0]
        y = [-self.r[1], self.r[1]]
        z = [0, 0]
        if d3 == '3d':
            ax.plot(X, Y, Z, 'b')
            ax.plot(x, y, z, 'b')
        else:
            ax.plot(X, Y, 'b')
            ax.plot(x, y, 'b')
