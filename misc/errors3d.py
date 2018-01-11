import sys

sys.path.insert(0, '../build/Release')

import common3d
import eikonal as eik
import itertools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import speedfuncs3d

from matplotlib.colors import LogNorm
from matplotlib.widgets import Slider, RadioButtons

matplotlib.interactive(True)

eps = np.finfo(float).eps

line_styles = ['o-', 'o--', 'o-.', 'o:']
line_cycler = itertools.cycle(line_styles)

def rms(x): return np.sqrt(x.dot(x)/x.size)

class PlotState():
    def __init__(self, marcher, Ms):
        self._marcher = marcher
        self._max_slice = int(Ms[-1] - 1)
        self._init_slice = int(self._max_slice/2)
        self._current_slice = self._init_slice
        self._axis = 'z'

    @property
    def marcher(self):
        return self._marcher

    @marcher.setter
    def marcher(self, value):
        if type(value) == str:
            self._marcher = common3d.get_marcher_by_name(value)
        else:
            self._marcher = value

    @property
    def min_slice(self):
        return 0

    @property
    def max_slice(self):
        return self._max_slice

    @property
    def init_slice(self):
        return self._init_slice

    @property
    def current_slice(self):
        return self._current_slice

    @current_slice.setter
    def current_slice(self, value):
        assert(value >= 0)
        assert(value <= self._max_slice)
        self._current_slice = int(value)

    @property
    def axis(self):
        return self._axis

    @axis.setter
    def axis(self, value):
        assert(value in {'x', 'y', 'z'})
        self._axis = value

if __name__ == '__main__':
    marchers = common3d.mid1_marchers
    s = speedfuncs3d.s0
    f = speedfuncs3d.f0

    minpow = 1
    maxpow = 5
    Ms = np.power(2, np.arange(minpow, maxpow + 1)) + 1

    diff = {marcher: np.zeros((Ms[-1], Ms[-1], Ms[-1])) for marcher in marchers}

    E = dict()
    E['inf'] = {marcher: np.zeros(Ms.shape) for marcher in marchers}
    E['rms'] = {marcher: np.zeros(Ms.shape) for marcher in marchers}

    for marcher, (i, M) in itertools.product(marchers, enumerate(Ms)):
        u = common3d.get_exact_soln(f, M)
        U = common3d.compute_soln(marcher, s, M)
        e_inf = common3d.relerr(u, U, np.inf)
        e_rms = rms((u - U).flatten())
        E['inf'][marcher][i] = e_inf
        E['rms'][marcher][i] = e_rms
        if i == len(Ms) - 1:
            diff[marcher] = np.maximum(np.abs(u - U), eps)
        print("%s (M = %d, e_inf = %g, e_rms = %g)" %
              (common3d.get_marcher_name(marcher), M, e_inf, e_rms))

    plt.figure(1)

    plt.subplot(121)
    for marcher in marchers:
        plt.loglog(Ms, E['inf'][marcher], next(line_cycler),
                   label=common3d.get_marcher_name(marcher))
    plt.title('l-inf error')

    plt.subplot(122)
    for marcher in marchers:
        plt.loglog(Ms, E['rms'][marcher], next(line_cycler),
                   label=common3d.get_marcher_name(marcher))
    plt.title('rms error')
    legend = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.show(1)

    fig, ax = plt.subplots()

    plt.subplots_adjust(left=0.25, bottom = 0.1, right=0.95, top=0.95)

    plot_state = PlotState(marchers[0], Ms)

    cmin = min(diff[m][diff[m] > eps].min() for m in marchers)
    cmax = max(diff[m].max() for m in marchers)
    vmin = 10**np.floor(np.log10(np.abs(cmin)))
    vmax = 10**np.ceil(np.log10(np.abs(cmax)))

    def get_slice(plot_state):
        if plot_state.axis == 'x':
            return diff[plot_state.marcher][:, plot_state.current_slice, :]
        elif plot_state.axis == 'y':
            return diff[plot_state.marcher][plot_state.current_slice, :, :]
        else:
            return diff[plot_state.marcher][:, :, plot_state.current_slice]

    img = plt.imshow(get_slice(plot_state), norm=LogNorm(vmin=vmin, vmax=vmax))
    plt.clim(cmin, cmax)
    plt.xticks([])
    plt.yticks([])
    plt.colorbar()

    ax_axis_radio_buttons = plt.axes([0.025, 0.525, 0.2, 0.4])
    axis_radio_buttons = RadioButtons(
        ax_axis_radio_buttons, ('x', 'y', 'z'), active=2)
    def axis_radio_buttons_update(label):
        plot_state.axis = label
        img.set_data(get_slice(plot_state))
        fig.canvas.draw_idle()
    axis_radio_buttons.on_clicked(axis_radio_buttons_update)

    ax_marcher_radio_buttons = plt.axes([0.025, 0.125, 0.2, 0.4])
    marcher_radio_buttons = RadioButtons(
        ax_marcher_radio_buttons,
        tuple(common3d.get_marcher_name(marcher) for marcher in marchers),
        active=0)
    def marcher_radio_buttons_update(marcher_name):
        plot_state.marcher = marcher_name
        img.set_data(get_slice(plot_state))
        fig.canvas.draw_idle()
    marcher_radio_buttons.on_clicked(marcher_radio_buttons_update)

    ax_slider = plt.axes([0.1, 0.025, 0.8, 0.05])
    slice_slider = Slider(
        ax_slider,
        'Slice',
        plot_state.min_slice,
        plot_state.max_slice,
        plot_state.init_slice,
        valfmt='%0.0f')
    def slider_update(val):
        plot_state.current_slice = int(slice_slider.val)
        img.set_data(get_slice(plot_state))
        fig.canvas.draw_idle()
    slice_slider.on_changed(slider_update)

    plt.show()
