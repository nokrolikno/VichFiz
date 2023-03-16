

from photon import Photon
from numpy import random
from source import Source
import numpy as np
import collections
from multiprocessing import Pool
import time


def timer(f):
    def tmp(*args, **kwargs):
        t = time.time()
        res = f(*args, **kwargs)
        print("function run time: %f" % (time.time()-t))
        return res

    return tmp


class Modeling:

    def __init__(self, surface, source, n, start_energy, min_energy):
        self.photones = []
        self.dell_photones = []
        self.surface = surface
        self.source = source
        self.n = n
        self.start_energy = start_energy
        self.min_energy = min_energy

    def get_surface(self):
        return self.surface

    def get_source(self):
        return self.source

    @timer
    def set_photones(self):

        for _ in range(self.n):
            photon = self.source.born_photon(self.start_energy)
            if self.surface.is_in(photon):
                self.photones.append(photon)

    def set_point_interaction(self):
        photones = []
        for i in range(len(self.photones)):
            photon = self.photones[i]

            if photon.next_interaction(self.min_energy):

                if self.surface.is_in(photon):
                    photones.append(photon)
                else:
                    photon.delete_last_position()
                    self.dell_photones.append(photon)

            else:
                photon.delete_last_position()
                if len(photon.get_points_interaction()) > 1:
                    self.dell_photones.append(photon)

        self.photones = photones

    def get_delete_photones(self):
        return self.dell_photones

    def get_photones(self):
        return self.photones

    def start_of_modeling(self):

        print("start modeling")

        print('photon making')

        self.set_photones()

        print('calculation of the following interactions')

        while len(self.photones) != 0:
            self.set_point_interaction()

        print('finish modeling')
