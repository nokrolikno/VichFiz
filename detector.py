from modeling import Modeling
from photon import Photon
import numpy as np
from math import pi, e, log, exp, sqrt


class Detector:

    def __init__(self, z, x=0, y=0):
        self.position = [x, y, z]
        self.rate = []
        self.energy_ranges = {}

    def get_position(self):
        return self.position

    def get_hist_rate(self):
        x_list = []
        y_list = []
        width_list = []
        for range_energy, rate in self.energy_ranges.items():
            width = range_energy[1] - range_energy[0]
            x_list.append(range_energy[0])
            y_list.append(rate)
            width_list.append(width)

        return [x_list, y_list, width_list]

    def set_energy_ranges(self, n_bins_hist, max_energy, min_energy):

        energy_list = np.linspace(min_energy, max_energy, n_bins_hist)
        self.energy_ranges = {(energy_list[i], energy_list[i + 1]): 0 for i in range(len(energy_list) - 1)}

    def get_key_energy(self, energy):
        keys_energy = self.energy_ranges.keys()
        for key in keys_energy:
            if key[0] <= energy <= key[1]:
                return key

    def flow_rate(self, modeling, n_bins_hist, energy_max, energy_min):
        print('calculation of contributions to the detector')
        photones = modeling.get_delete_photones()
        self.set_energy_ranges(n_bins_hist+1, energy_max, energy_min)
        position = self.position
        for photon in photones:
            weight = photon.get_weight()
            mu = photon.get_mu()
            a = photon.get_energy_photon()
            points = photon.get_points_interaction()
            sigma_total = photon.get_sigma_total()
            for i in range(1, len(photon.get_points_interaction())):
                a1 = a[i - 1]
                a2 = a1 * a1
                mu1 = mu[i]
                mu2 = mu1 * mu1

                r_pd2 = (points[i][0] - position[0])**2 +\
                        (points[i][1] - position[1])**2 +\
                        (points[i][2] - position[2])**2

                sigma_relationship = 1 / 4 / pi * (1 + mu2 + a2 * (1 - mu2) / (1 + a1 * (1 - mu1))) /\
                                     (1 + a1*(1 - mu1))**2 /\
                                     ((1 + a1) / a2 * (2 * (1 + a1) / (1 + 2*a1) - log(1 + 2 * a1, e) / a1) +
                                      log(1 + 2 * a1, e) / (2 * a1) - (1 + 3 * a1) / (1 + 2 * a1)**2)

                nu = weight[i] * exp(-sigma_total[i] * sqrt(r_pd2)) / r_pd2 * sigma_relationship
                self.energy_ranges[self.get_key_energy(a1)] += nu / len(photones)
