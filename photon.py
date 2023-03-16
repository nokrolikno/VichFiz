
from math import *
from numpy import *
import numpy as np
import matplotlib.pyplot as plt


def load_data(name):
    data = []
    with open(name) as f:
        for line in f:
            data.append([float(x) for x in line.split()])
    return data


def plot_energy(x, y, _xlim):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, '--.')
    ax.grid(True)
    plt.xlim(_xlim)
    plt.show()


class Photon:

    data = np.array(load_data('energy.txt'))
    sigma_total_interp = data[:, 2]
    sigma_compton_interp = data[:, 1]
    energy_photon_interp = data[:, 0]

    mc2 = 0.511

    def __init__(self, point_born, energy_photon):
        self.points_interaction = []
        self.energy_photon = []
        self.weight = [1.0]
        self.sigma_compton = []
        self.sigma_total = []
        self.mu = []
        self.points_interaction.append(point_born)
        self.energy_photon.append(energy_photon)
        sigma_compton, sigma_total = self.interpolate_linear(energy_photon)
        self.sigma_compton.append(sigma_compton)
        self.sigma_total.append(sigma_total)

        cos_tetha = 2 * random.random() - 1
        self.mu.append(cos_tetha)
        cos_ksi, sin_ksi = self.get_ksi()
        self.omega1 = cos_tetha * cos_ksi
        self.omega2 = cos_tetha * sin_ksi
        self.omega3 = sqrt(1 - cos_tetha * cos_tetha)

        length = - math.log(random.random(), math.e) / sigma_total

        point_interaction = [0, 0, 0]
        point_interaction[0] = length * self.omega1 + self.points_interaction[-1][0]
        point_interaction[1] = length * self.omega2 + self.points_interaction[-1][1]
        point_interaction[2] = length * self.omega3 + self.points_interaction[-1][2]
        self.points_interaction.append(point_interaction)

    @staticmethod
    def get_ksi():
        while True:
            a = 1 - 2 * random.random()
            b = 1 - 2 * random.random()
            d = a * a + b * b
            if d <= 1:
                break

        return [a / sqrt(d), b / sqrt(d)]

    def interpolate_linear(self, energy):

        left = 0
        right = len(self.energy_photon_interp) - 1
        energy_photon_left = self.energy_photon_interp[left]
        energy_photon_right = self.energy_photon_interp[right]

        if energy <= energy_photon_left:
            return [self.sigma_compton_interp[0], self.sigma_total_interp[0]]

        if energy >= energy_photon_right:
            return [self.sigma_compton_interp[-1], self.sigma_total_interp[-1]]

        while right - left > 1:
            i = (right - left) // 2 + left
            if energy < self.energy_photon_interp[i]:
                right = i
                energy_photon_right = self.energy_photon_interp[right]
            else:
                left = i
                energy_photon_left = self.energy_photon_interp[left]

        sigma_compton_0 = self.sigma_compton_interp[left]
        sigma_compton_1 = self.sigma_compton_interp[right]

        sigma_total_0 = self.sigma_total_interp[left]
        sigma_total_1 = self.sigma_total_interp[right]

        sigma_compton = sigma_compton_0 + (sigma_compton_1 - sigma_compton_0) *\
                        (energy - energy_photon_left) / (energy_photon_right - energy_photon_left)

        sigma_total = sigma_total_0 + (sigma_total_1 - sigma_total_0) * \
                        (energy - energy_photon_left) / (energy_photon_right - energy_photon_left)

        return [sigma_compton, sigma_total]

    def next_interaction(self, min_energy):
        sigma_total = self.sigma_total[-1]
        sigma_compton = self.sigma_compton[-1]
        weight = self.weight[-1] * sigma_compton / sigma_total

        if sigma_compton / sigma_total <= random.random() or weight < e-10 or min_energy > self.get_last_energy():
            return False

        self.weight.append(weight)

        a_old = self.get_last_energy() / self.mc2
        while True:
            g1 = random.random()
            g2 = random.random()
            a = a_old * (1 + 2 * a_old * g1) / (1 + 2 * a_old)
            if g2 * (1 + 2 * a_old + 1 / (1 + 2 * a_old)) < self.p(a, a_old):
                break

        mu = 1 - 1 / a + 1 / a_old
        self.mu.append(mu)
        cos_ksi, sin_ksi = self.get_ksi()

        energy_photon = a_old * self.mc2 / (1 + a_old * (1 - mu))
        sigma_compton, sigma_total = self.interpolate_linear(energy_photon)
        self.energy_photon.append(energy_photon)
        self.sigma_total.append(sigma_total)
        self.sigma_compton.append(sigma_compton)

        omega1_ = self.omega1
        omega2_ = self.omega2
        omega3_ = self.omega3

        omega3 = omega3_ * mu + sqrt((1 - mu * mu) * (1 - omega3_ * omega3_)) * cos_ksi
        omega2 = (omega2_ * (mu - omega3_ * omega3) +
                  omega1_ * sin_ksi * sqrt((1 - mu * mu) * (1 - omega3_ * omega3_))) / (1 - omega3_ * omega3_)
        omega1 = (omega1_ * (mu - omega3_ * omega3) -
                  omega2_ * sin_ksi * sqrt((1 - mu * mu) * (1 - omega3_ * omega3_))) / (1 - omega3_ * omega3_)

        self.omega1 = omega1
        self.omega2 = omega2
        self.omega3 = omega3

        length = - math.log(random.random(), math.e) / sigma_total
        point_interaction = [0, 0, 0]
        point_interaction[0] = length * omega1 + self.points_interaction[-1][0]
        point_interaction[1] = length * omega2 + self.points_interaction[-1][1]
        point_interaction[2] = length * omega3 + self.points_interaction[-1][2]
        self.points_interaction.append(point_interaction)

        return True

    def p(self, x, a_old):
        return x / a_old + a_old / x + (1 / a_old - 1 / x) * (2 + 1 / a_old - 1 / x)

    def get_trajectory(self):
        x = []
        y = []
        z = []
        for point in self.points_interaction:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        trajectory = [x, y, z]
        return trajectory

    def get_last_position(self):
        return self.points_interaction[-1]

    def get_last_energy(self):
        return self.energy_photon[-1]

    def delete_last_position(self):
        self.points_interaction.pop()

    def get_points_interaction(self):
        return self.points_interaction

    def get_sigma_total(self):
        return self.sigma_total

    def get_mu(self):
        return self.mu

    def get_energy_photon(self):
        return self.energy_photon

    def get_weight(self):
        return self.weight
