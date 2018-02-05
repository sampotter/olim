import sys

sys.path.insert(0, '../build/Release')

import argparse
import eikonal as eik
import glob
import h5py
import numpy as np
import os.path

import wx
import wx.lib.agw.aui as aui

import matplotlib as mpl
import matplotlib.figure
import matplotlib.backends.backend_wxagg as mpl_backend

from common3d import compute_soln, get_exact_soln, get_marcher_name, marchers, \
    time_marcher
from itertools import product
from speedfuncs3d import get_speed_func_name, get_soln_func, speed_funcs
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ColumnSorterMixin, \
    ListCtrlAutoWidthMixin

def get_dataset_names(f):
    mnames = list(f.keys())
    snames = list(f[mnames[0]].keys())
    return [mstr + '/' + sstr for mstr, sstr in product(mnames, snames)]

class CheckListCtrl(wx.ListCtrl, CheckListCtrlMixin, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT |
                             wx.SUNKEN_BORDER)
        CheckListCtrlMixin.__init__(self)
        ListCtrlAutoWidthMixin.__init__(self)

class LogLogPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        global fig, ax, canvas, legend
        fig = mpl.figure.Figure()
        ax = fig.add_subplot(111)
        legend = ax.legend()
        legend.set_visible(False)
        canvas = mpl_backend.FigureCanvas(self, wx.ID_ANY, fig)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(canvas, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

class DatasetSelectPanel(wx.Panel, ColumnSorterMixin):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self._list = CheckListCtrl(self)
        self._list.InsertColumn(0, 'OLIM')
        self._list.InsertColumn(1, 'Type')
        self._list.InsertColumn(2, 'Function')

        self.itemDataMap = dict()
        self._dataset_names = list(get_dataset_names(hdf5_file))
        self._plots = dict()
        for dataset_name in self._dataset_names:
            mname, sname = dataset_name.split('/')
            if mname == 'basic_3d':
                olim, type_ = 'basic_3d', 'finite_diff'
            else:
                olim, type_ = mname.split('_')
            index = self._list.InsertItem(sys.maxsize, olim)
            self._list.SetItem(index, 1, type_)
            self._list.SetItem(index, 2, sname)
            self._list.SetItemData(index, index)
            self.itemDataMap[index] = (olim, type_, sname)
        self._list.OnCheckItem = self.OnCheckItem
        
        self._list.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self._list.SetColumnWidth(1, wx.LIST_AUTOSIZE)
        self._list.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(self._list, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

        ColumnSorterMixin.__init__(self, 3)

    def GetListCtrl(self):
        return self._list
    
    def OnCheckItem(self, index, flag):
        dataset_name = self._dataset_names[index]
        if not flag:
            self._plots[dataset_name].remove()
            del self._plots[dataset_name]
            if len(self._plots) == 0:
                assert legend.get_visible()
                legend.set_visible(False)
        else:
            t = hdf5_file[dataset_name + '/t']
            rms = hdf5_file[dataset_name + '/rms']
            line, = ax.loglog(t, rms, 'o-', label=dataset_name)
            self._plots[dataset_name] = line
            if len(self._plots) == 1:
                assert not legend.get_visible()
                legend.set_visible(True)
        ax.legend()
        canvas.draw()

class PlotFrame(wx.Frame):
    def __init__(self, parent, id=-1, title="test",
                 pos=wx.DefaultPosition, size=(800, 600),
                 style=wx.DEFAULT_FRAME_STYLE):
        wx.Frame.__init__(self, parent, id, title, pos, size, style)
        self._mgr = aui.AuiManager()
        self._mgr.SetManagedWindow(self)
        self._panels = dict()
        self._panels['loglog'] = LogLogPanel(self)
        self._panels['dataset_select'] = DatasetSelectPanel(self)
        self._mgr.AddPane(
            self._panels['loglog'],
            aui.AuiPaneInfo().CenterPane())
        self._mgr.AddPane(
            self._panels['dataset_select'],
            aui.AuiPaneInfo().Left().Caption("Dataset Selection"))
        self._mgr.Update()
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, event):
        self._mgr.UnInit()
        event.Skip()

class PlotApp(wx.App):
    def __init__(self, size=(800, 600), redirect=False,
                 filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self, redirect, filename, useBestVisual, clearSigInt)
        self._frame = PlotFrame(None, size=size)
        self.SetTopWindow(self._frame)
        self._frame.Show()

    def OnInit(self):
        return True
        
def rms(x): return np.sqrt(x.dot(x)/x.size)

def get_ns(args):
    minpow = args.minpow
    maxpow = args.maxpow
    steps = args.step
    ns = np.logspace(minpow, maxpow, steps*(maxpow - minpow) + 1, base=2)
    return (2*np.round(ns/2)).astype(int) + 1

def get_dataset_name(Marcher, s):
    mname = get_marcher_name(Marcher)
    sname = get_speed_func_name(s)
    return '%s/%s' % (mname.replace(' ', '_'), sname)

def get_rms(Marcher, s, ns):
    return np.array([
        rms((get_exact_soln(get_soln_func(s), n) -
             compute_soln(Marcher, s, n)).flatten())
        for n in ns])

def get_times(Marcher, s, ns):
    return np.array([time_marcher(Marcher, s, n) for n in ns])

INITIAL_SIZE = (1000, 800)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--minpow', type=int, default=3)
    parser.add_argument('--maxpow', type=int, default=7)
    parser.add_argument('--step', type=int, default=2)
    parser.add_argument('--hdf5_path', type=str, default='time_vs_error.hdf5')
    args = parser.parse_args()

    path = args.hdf5_path

    if not os.path.exists(path):
        print("%s doesn't exist... creating" % path)
        with h5py.File(path, 'w') as f:
            ns = get_ns(args)
            prod = product(marchers, speed_funcs())
            for M, s in prod:
                name = get_dataset_name(M, s)
                print(name)
                f.create_dataset(name + '/rms', data=get_rms(M, s, ns))
                f.create_dataset(name + '/t', data=get_times(M, s, ns))

    global hdf5_file
    hdf5_file = h5py.File(path, 'r')

    app = PlotApp(size=INITIAL_SIZE)
    app.MainLoop()
