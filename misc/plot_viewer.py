#!/usr/bin/env python3

import sys

sys.path.insert(0, '../build/Release')

import argparse
import glob
import h5py
import numpy as np
import os.path

import wx
import wx.lib.agw.aui as aui

from wx.lib.mixins.listctrl import \
    CheckListCtrlMixin, \
    ColumnSorterMixin, \
    ListCtrlAutoWidthMixin

import matplotlib as mpl
import matplotlib.figure
import matplotlib.backends.backend_wxagg as mpl_backend

from itertools import product

plotdata = dict()

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

        d = dict()
        d['fig'] = mpl.figure.Figure()
        d['ax'] = d['fig'].add_subplot(111)
        d['ax'].set_xlabel(self._xlabel)
        d['ax'].set_ylabel(self._ylabel)
        d['legend'] = d['ax'].legend()
        assert d['legend'] is not None
        d['legend'].set_visible(False)
        d['canvas'] = mpl_backend.FigureCanvas(self, wx.ID_ANY, d['fig'])

        plotdata[self._plotdata_key] = d

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(d['canvas'], flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

class LogLogTimeVsRmsErrorPanel(LogLogPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_time_vs_rms_error'
        self._xlabel = 'Time'
        self._ylabel = 'RMS Error'
        LogLogPanel.__init__(self, parent, *args, **kwargs)

class LogLogTimeVsMaxErrorPanel(LogLogPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_time_vs_max_error'
        self._xlabel = 'Time'
        self._ylabel = 'Max Error'
        LogLogPanel.__init__(self, parent, *args, **kwargs)

class LogLogSizeVsRmsErrorPanel(LogLogPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_size_vs_rms_error'
        self._xlabel = 'Problem Size (n)'
        self._ylabel = 'RMS Error'
        LogLogPanel.__init__(self, parent, *args, **kwargs)

class LogLogSizeVsMaxErrorPanel(LogLogPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_size_vs_max_error'
        self._xlabel = 'Problem Size (n)'
        self._ylabel = 'Max Error'
        LogLogPanel.__init__(self, parent, *args, **kwargs)

class LogLogSizeVsTimePanel(LogLogPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_size_vs_time'
        self._xlabel = 'Problem Size (n)'
        self._ylabel = 'Time'
        LogLogPanel.__init__(self, parent, *args, **kwargs)

class ErrorVolPanel(wx.SplitterWindow):
    def __init__(self, parent, *args, **kwargs):
        wx.SplitterWindow.__init__(
            self, parent, style=wx.SP_THIN_SASH, *args, **kwargs)

        self._top_panel = wx.Panel(self, style=wx.SUNKEN_BORDER)
        self.Initialize(self._top_panel)

        self._bottom_panel = wx.Panel(self, style=wx.SUNKEN_BORDER)
        self.Initialize(self._bottom_panel)

        self.SplitHorizontally(self._top_panel, self._bottom_panel, -60)

        # Top panel widgets

        d = dict()
        d['fig'] = mpl.figure.Figure()
        d['ax'] = d['fig'].add_subplot(111)
        d['legend'] = d['ax'].legend()
        assert d['legend'] is not None
        d['legend'].set_visible(False)
        d['canvas'] = mpl_backend.FigureCanvas(
            self._top_panel, wx.ID_ANY, d['fig'])
        d['im'] = d['ax'].imshow([[1]])
        d['colorbar'] = d['fig'].colorbar(d['im'])
        plotdata['errorvol'] = d

        self._top_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._top_panel_sizer.Add(d['canvas'], flag=wx.EXPAND, proportion=1)
        self._top_panel.SetSizer(self._top_panel_sizer)

        # Bottom panel widgets

        self._axis_names = ['xy', 'yz', 'xz']
        self._axis_radio_box = wx.RadioBox(
            self._bottom_panel, wx.ID_ANY, 'Slice Plane', wx.DefaultPosition,
            wx.DefaultSize, self._axis_names)

        self.Bind(wx.EVT_RADIOBOX, self.HandleRadioEvent, self._axis_radio_box)

        self._slice_slider_style = \
            wx.SL_HORIZONTAL | wx.SL_BOTTOM | wx.SL_LABELS
        self._slice_slider = wx.Slider(
            self._bottom_panel, wx.ID_ANY, 128, 0, 256, wx.DefaultPosition,
            wx.DefaultSize, self._slice_slider_style)

        self.Bind(wx.EVT_SLIDER, self.HandleSliderEvent, self._slice_slider)

        self._size_text = wx.StaticText(
            self._bottom_panel, wx.ID_ANY, "Size (n):", wx.DefaultPosition)

        self._sizes = list(map(str, [5, 9, 17, 33]))
        self._size_chooser = wx.Choice(
            self._bottom_panel, wx.ID_ANY, wx.DefaultPosition,
            choices=self._sizes)

        self.Bind(wx.EVT_CHOICE, self.HandleChoiceEvent, self._size_chooser)

        self._bottom_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self._bottom_panel_sizer.Add(self._axis_radio_box)
        self._bottom_panel_sizer.Add(self._slice_slider)
        self._bottom_panel_sizer.Add(self._size_text)
        self._bottom_panel_sizer.Add(self._size_chooser)

        self._bottom_panel.SetSizer(self._bottom_panel_sizer)

    def HandleRadioEvent(self, event):
        index = event.GetInt()
        d = plotdata['errorvol']
        d['slice_plane'] = self._axis_names[index]
        d['im'].set_data(d['data'][[
            np.index_exp[:, :, d['slice']['xy']],
            np.index_exp[d['slice']['yz'], :, :],
            np.index_exp[:, d['slice']['xz'], :]
        ][index]])
        d['canvas'].draw()
        self._slice_slider.SetValue(d['slice'][d['slice_plane']])

    def HandleSliderEvent(self, event):
        index = event.GetInt()
        d = plotdata['errorvol']
        d['slice'][d['slice_plane']] = index
        d['im'].set_data(d['data'][[
            np.index_exp[:, :, index],
            np.index_exp[index, :, :],
            np.index_exp[:, index, :]
        ][{'xy': 0, 'yz': 1, 'xz': 2}[d['slice_plane']]]])
        d['canvas'].draw()

    def HandleChoiceEvent(self, event):
        index = event.GetInt()
        size = self._sizes[index]
        self.SetSelectedSize(size)

    def SetSizes(self, sizes):
        self._sizes = sizes
        self._size_chooser.Clear()
        self._size_chooser.SetItems([str(s) for s in sizes])

    def SetSelectedSize(self, n):
        assert(n in self._sizes)

        self._slice_slider.SetMax(n - 1)

        d = plotdata['errorvol']
        s = d['slice'][d['slice_plane']]
        self._slice_slider.SetValue(s)

        u = np.array(hdf5_file[d['dsetname'] + '/u%d' % n])
        U = np.array(hdf5_file[d['dsetname'] + '/U%d' % n])
        d['data'] = u - U

        d['slice']['xy'] = int(n/2)
        d['slice']['yz'] = int(n/2)
        d['slice']['xz'] = int(n/2)

        d['im'].set_data(d['data'][:, :, d['slice'][d['slice_plane']]])

        vmin = d['data'].min()
        vmax = d['data'].max()
        d['colorbar'].set_clim(vmin, vmax)
        d['colorbar'].draw_all()

        d['canvas'].draw()

class LogLogControlPanel(wx.Panel, ColumnSorterMixin):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self._list = CheckListCtrl(self)
        self._list.InsertColumn(0, 'OLIM')
        self._list.InsertColumn(1, 'Type')
        self._list.InsertColumn(2, 'Function')

        # Initialize dictionary used to keep track of plots
        self._plots = dict()

        self.itemDataMap = dict()
        if hdf5_file is not None:
            for dataset_name in get_dataset_names(hdf5_file):
                mname, sname = dataset_name.split('/')
                if mname == 'basic_3d':
                    olim, type_ = 'basic_3d', 'finite_diff'
                elif 'hu' in mname:
                    olim, type_ = mname[:9], mname[10:]
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
        if hdf5_file is not None:
            dsetname = get_dataset_names(hdf5_file)[index]
            d = plotdata[self._plotdata_key]
            if not flag:
                self._plots[dsetname].remove()
                del self._plots[dsetname]
                if len(self._plots) == 0:
                    assert d['legend'] is not None
                    assert d['legend'].get_visible()
                    d['legend'].set_visible(False)
            else:
                X = hdf5_file[dsetname + '/' + self._x_axis_key]
                Y = hdf5_file[dsetname + '/' + self._y_axis_key]
                line, = d['ax'].loglog(X, Y, 'o-', label=dsetname)
                self._plots[dsetname] = line
                if len(self._plots) == 1:
                    assert d['legend'] is not None
                    assert not d['legend'].get_visible()
                    d['legend'].set_visible(True)
            d['ax'].legend()
            d['canvas'].draw()

class LogLogTimeVsRmsErrorControlPanel(LogLogControlPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_time_vs_rms_error'
        self._x_axis_key = 't'
        self._y_axis_key = 'rms'
        LogLogControlPanel.__init__(self, parent, *args, **kwargs)

class LogLogTimeVsMaxErrorControlPanel(LogLogControlPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_time_vs_max_error'
        self._x_axis_key = 't'
        self._y_axis_key = 'max'
        LogLogControlPanel.__init__(self, parent, *args, **kwargs)

class LogLogSizeVsRmsErrorControlPanel(LogLogControlPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_size_vs_rms_error'
        self._x_axis_key = 'n'
        self._y_axis_key = 'rms'
        LogLogControlPanel.__init__(self, parent, *args, **kwargs)

class LogLogSizeVsMaxErrorControlPanel(LogLogControlPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_size_vs_max_error'
        self._x_axis_key = 'n'
        self._y_axis_key = 'max'
        LogLogControlPanel.__init__(self, parent, *args, **kwargs)

class LogLogSizeVsTimeControlPanel(LogLogControlPanel):
    def __init__(self, parent, *args, **kwargs):
        self._plotdata_key = 'loglog_size_vs_time'
        self._x_axis_key = 'n'
        self._y_axis_key = 't'
        LogLogControlPanel.__init__(self, parent, *args, **kwargs)

class ErrorVolControlPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self._list = wx.ListCtrl(
            self, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)

        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnItemSelected,
                  self._list)

        self._list.InsertColumn(0, 'OLIM')
        self._list.InsertColumn(1, 'Type')
        self._list.InsertColumn(2, 'Function')

        if hdf5_file is not None:
            for dsetname in get_dataset_names(hdf5_file):
                mname, sname = dsetname.split('/')
                if mname == 'basic_3d':
                    olim, type_ = 'basic_3d', 'finite_diff'
                elif 'hu' in mname:
                    olim, type_ = mname[:9], mname[10:]
                else:
                    olim, type_ = mname.split('_')
                index = self._list.InsertItem(sys.maxsize, olim)
                self._list.SetItem(index, 1, type_)
                self._list.SetItem(index, 2, sname)

        self._list.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self._list.SetColumnWidth(1, wx.LIST_AUTOSIZE)
        self._list.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(self._list, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

    def SetErrorVolPanel(self, error_vol_panel):
        self._error_vol_panel = error_vol_panel

    def OnItemSelected(self, event):
        index = event.GetIndex()
        dsetname = get_dataset_names(hdf5_file)[index]
        sizes = sorted({
            int(k[1:]) for k in hdf5_file[dsetname].keys() if
            k[0] in {'u', 'U'}})

        self._error_vol_panel.SetSizes(sizes)

        d = plotdata['errorvol']
        d['dsetname'] = dsetname
        d['sizes'] = sizes

        n = sizes[int(len(sizes)/2)] # middle element... why not
        u = np.array(hdf5_file[dsetname + '/u%d' % n])
        U = np.array(hdf5_file[dsetname + '/U%d' % n])
        d['data'] = u - U

        d['slice_plane'] = 'xy'

        d['slice'] = dict()

        d['slice']['xy'] = int(n/2)
        d['slice']['yz'] = int(n/2)
        d['slice']['xz'] = int(n/2)

        d['im'] = d['ax'].imshow(d['data'][:, :, d['slice']['xy']])

        vmin = d['data'].min()
        vmax = d['data'].max()
        d['colorbar'].set_clim(vmin, vmax)

        d['canvas'].draw()

        self._error_vol_panel.SetSelectedSize(n)

class DatasetSelectPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)

        self._panels = dict()
        self._panels['loglog_time_vs_rms_error'] = LogLogTimeVsRmsErrorControlPanel(self)
        self._panels['loglog_time_vs_max_error'] = LogLogTimeVsMaxErrorControlPanel(self)
        self._panels['loglog_size_vs_rms_error'] = LogLogSizeVsRmsErrorControlPanel(self)
        self._panels['loglog_size_vs_max_error'] = LogLogSizeVsMaxErrorControlPanel(self)
        self._panels['loglog_size_vs_time'] = LogLogSizeVsTimeControlPanel(self)
        self._panels['error_vol_control'] = ErrorVolControlPanel(self)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        for panel in self._panels.values():
            panel.Hide()
            self._sizer.Add(panel, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

        self._panels['loglog_time_vs_rms_error'].Show()

    def SetErrorVolPanel(self, error_vol_panel):
        self._panels['error_vol_control'].SetErrorVolPanel(error_vol_panel)

    def switch_control_panels(self, key):
        k = next((k for k in self._panels if self._panels[k].IsShown()), None)
        if k == key:
            return
        self._panels[k].Hide()
        self._panels[key].Show()
        self.Layout()

class NotebookPanel(wx.Panel):
    def __init__(self, parent, dataset_select_panel, *args, **kwargs):
        self._parent = parent
        self._dataset_select_panel = dataset_select_panel

        wx.Panel.__init__(self, parent, *args, **kwargs)

        self._notebook = aui.AuiNotebook(self)
        self._notebook.Bind(aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.page_changed)

        error_vol_panel = ErrorVolPanel(self._notebook)
        self._dataset_select_panel.SetErrorVolPanel(error_vol_panel)

        self._pages_and_labels = [
            (LogLogTimeVsRmsErrorPanel(self._notebook), 'Time vs. RMS Error'),
            (LogLogTimeVsMaxErrorPanel(self._notebook), 'Time vs. Max Error'),
            (LogLogSizeVsRmsErrorPanel(self._notebook), 'Size vs. RMS Error'),
            (LogLogSizeVsMaxErrorPanel(self._notebook), 'Size vs. Max Error'),
            (LogLogSizeVsTimePanel(self._notebook), 'Size vs. Time'),
            (error_vol_panel, "3D Error")
        ]
        for page, label in self._pages_and_labels:
            self._notebook.AddPage(page, label)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(self._notebook, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

    def page_changed(self, event):
        keys = ['loglog_time_vs_rms_error', 'loglog_time_vs_max_error',
                'loglog_size_vs_rms_error', 'loglog_size_vs_max_error',
                'loglog_size_vs_time', 'error_vol_control']
        old_sel = event.GetOldSelection()
        if old_sel == -1:
            return
        sel = event.GetSelection()
        self._dataset_select_panel.switch_control_panels(keys[sel])

class PlotFrame(wx.Frame):
    def __init__(self, parent, id=-1, title='Plot Viewer',
                 pos=wx.DefaultPosition, size=(800, 600),
                 style=wx.DEFAULT_FRAME_STYLE):
        wx.Frame.__init__(self, parent, id, title, pos, size, style)

        self._mgr = aui.AuiManager()
        self._mgr.SetManagedWindow(self)

        self._panels = dict()
        self._panels['dataset_select'] = DatasetSelectPanel(self)
        self._panels['notebook'] = NotebookPanel(
            self, self._panels['dataset_select'])
        self._mgr.AddPane(
            self._panels['notebook'],
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
        
INITIAL_SIZE = (1000, 800)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    args = parser.parse_args()

    global hdf5_file
    if args.path is not None:
        hdf5_file = h5py.File(args.path, 'r')
    else:
        hdf5_file = None

    app = PlotApp(size=INITIAL_SIZE)
    app.MainLoop()
