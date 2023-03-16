from photon import Photon
from source import Source
from surface import Surface
from cylinder_surface import CylinderSurface
from half_disk_source import HalfDiskSource
from modeling import Modeling
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
from detector import Detector
import time


figure = 2

def timer(f):
    def tmp(*args, **kwargs):
        t = time.time()
        res = f(*args, **kwargs)
        print("function run time: %f" % (time.time()-t))
        return res

    return tmp


def plot_born_photon(modeling):
    plt.rcParams['legend.fontsize'] = 10
    fig, ax = plt.subplots()
    modeling.get_source().plot_source(ax, '2d')
    for photon in modeling.get_delete_photones():
        trajectory = photon.get_trajectory()
        ax.scatter(trajectory[0][0], trajectory[1][0])
    plt.savefig('born_photones_{:.1e}_{}.pdf'.format(n, source_radius), bbox_inches='tight')


@timer
def plot_trajectory(modeling, detector=None):
    print('plot trajectory')
    plt.rcParams['legend.fontsize'] = 10
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    if detector is not None:
        ax.scatter(*detector.get_position(), color='r', s=100)

    modeling.get_source().plot_source(ax)
    modeling.get_surface().plot_surface(ax)

    for photon in modeling.get_delete_photones():
        trajectory = photon.get_trajectory()
        ax.plot(trajectory[0], trajectory[1], trajectory[2], '--.')

    plt.savefig('trajectory_{:.1e}_{}.png'.format(n, surface_height), bbox_inches='tight')


@timer
def plot_hist(detector):
    print('plot hist')
    global figure
    plt.figure(figure)
    figure += 1
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    x_list, y_list, width_list = detector.get_hist_rate()
    ax.bar(x_list, y_list, width_list, color='b', align='edge', edgecolor='r', alpha=0.7)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
    ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
    plt.xlabel('Энергия, МэВ')
    plt.ylabel('Плотность потока')
    plt.title('Зависимость вклада от энергии фотона')
    position = detector.get_position()
    plt.suptitle('Координаты детектора ({}, {}, {}), количество частиц {:.1e}'
                 .format(round(position[0], 2), round(position[1], 2), round(position[2], 2), n))
    plt.savefig('plot_hist_{}_{}_{}_{}.png'.format(figure-2, round(position[0], 2), round(position[1], 2), round(position[2], 2), bbox_inches='tight'))


n = 500
start_energy = 3.5
min_energy = 0.1
surface_radius = 25
surface_height = 10
source_radius = 20
n_bins_hist = 20


def main():

    start_time = time.process_time()
    surface = CylinderSurface(surface_radius, surface_height)
    source = HalfDiskSource(source_radius)

    detectors = [
        Detector(surface_height),
        Detector(0, x=surface_radius, y=0),
        Detector(surface_height / 3, x=surface_radius, y=0),
        Detector(2 * surface_height / 3, x=surface_radius, y=0),
        Detector(surface_height, x=surface_radius),
        Detector(surface_height / 3),
        Detector(2 * surface_height / 3),
        Detector(surface_height, x=surface_radius / 2),
        Detector(surface_height, y=surface_radius / 2),
        Detector(surface_height, x=-surface_radius / 2),
        Detector(surface_height, y=-surface_radius / 2),
    ]

    modeling = Modeling(surface, source, n // len(detectors), start_energy, min_energy)
    modeling_many = Modeling(surface, source, n, start_energy, min_energy)
    modeling.start_of_modeling()
    modeling_many.start_of_modeling()
    for detector in detectors:
        detector.flow_rate(modeling_many, n_bins_hist, start_energy, min_energy)

    print("FULL MODELING TIME: {:g} s".format(time.process_time() - start_time))

    for detector in detectors:
        plot_trajectory(modeling, detector)

    plot_born_photon(modeling)
    for detector in detectors:
        plot_hist(detector)
    plt.show()


if __name__ == '__main__':
    main()
