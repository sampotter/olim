import argparse
import enum
import os.path
import subprocess
import sys

import git

import wx
import wx.lib.agw.aui as aui

import matplotlib as mpl
import matplotlib.figure
import matplotlib.backends.backend_wxagg as mpl_backend
import numpy as np

class CMakeBuildType(enum.Enum):
    Debug = 1
    RelWithDebInfo = 2
    Release = 3

    def __str__(self, value):
        return {
            1: 'Debug',
            2: 'RelWithDebInfo',
            3: 'Release'
        }[value]

class CMakeBackend(enum.Enum):
    Make = 1
    Ninja = 2

    def __str__(self, value):
        return {
            1: 'Make',
            2: 'Ninja',
        }[value]

class CMakeBuild(object):
    def __init__(self, path=None, build_type=CMakeBuildType.Release,
                 definitions=None, backend=CMakeBackend.Make):
        if not path:
            raise 'missing path'
        self._path = path
        self._build_type = build_type
        self._definitions = definitions
        self._backend = backend

    def get_build_cmd(self):
        args = ['cmake']
        args += ['-D%s=%s' % (k, v) for k, v in self._definitions.items()]
        args =+ ['-DCMAKE_BUILD_TYPE=%s' % str(self._build_type)]
        args += ['-G%s' % str(self._backend)]
        args += [path]
        return args

    def build(self):
        curdir = os.curdir()
        os.chdir(self._path)
        subprocess.run(self.get_build_cmd())
        os.chdir(curdir)

class CMakeProject(object):
    def __init__(self, git_url=None, repo_path=None):
        self._repo = git.Repo.clone_from(git_url, repo_path)
        self._repo_path = repo_path

    def build_release(self, path='build/Release'):
        build_path = os.path.join(self._repo_path, path)
        os.makedirs(build_path)

class ErrorSlicePanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self._figure = mpl.figure.Figure()
        self._axes = self._figure.add_subplot(111)
        self._axes.plot(np.arange(10), np.arange(10))
        self._figure_canvas = mpl_backend.FigureCanvas(
            self, wx.ID_ANY, self._figure)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(self._figure_canvas, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

class CommitListPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self._list = wx.ListCtrl(self, wx.ID_ANY, style=wx.LC_REPORT)
        self._list.InsertColumn(0, "Git SHA-1 Hash")
        self._list.InsertColumn(1, "Message")
        for i in range(10):
            index = self._list.InsertStringItem(sys.maxsize, str(i))
            self._list.SetStringItem(index, 1, "blah")
        self._list.SetColumnWidth(0, wx.LIST_AUTOSIZE_USEHEADER)
        self._list.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(self._list, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

class TextOutputPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self._text = wx.TextCtrl(self, wx.ID_ANY, style=wx.TE_MULTILINE)

        self._sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._sizer.Add(self._text, flag=wx.EXPAND, proportion=1)
        self.SetSizer(self._sizer)

class VolPlotFrame(wx.Frame):
    def __init__(self, parent, id=-1, title="AUI Test", pos=wx.DefaultPosition,
                 size=(800, 600), style=wx.DEFAULT_FRAME_STYLE):
        wx.Frame.__init__(self, parent, id, title, pos, size, style)
        
        # Setup AUI
        self._mgr = aui.AuiManager()
        self._mgr.SetManagedWindow(self)

        # Create panels
        self._panels = dict()
        self._panels['error_slice'] = ErrorSlicePanel(self)
        self._panels['commit_list'] = CommitListPanel(self)
        self._panels['output'] = TextOutputPanel(self)

        # Add the panes to the manager
        self._mgr.AddPane(
            self._panels['error_slice'],
            aui.AuiPaneInfo().CenterPane())
        self._mgr.AddPane(
            self._panels['commit_list'],
            aui.AuiPaneInfo().Left().Caption("Commit List"))
        self._mgr.AddPane(
            self._panels['output'],
            aui.AuiPaneInfo().Bottom().Caption("Output"))

        # Commit all changes just made
        self._mgr.Update()

        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, event):
        self._mgr.UnInit()
        event.Skip()

class VolPlotApp(wx.App):
    def __init__(self, size=(800, 600), redirect=False, filename=None,
                 useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self, redirect, filename, useBestVisual,
                        clearSigInt)
        self._frame = VolPlotFrame(None, size=size)
        self.SetTopWindow(self._frame)
        self._frame.Show()

    def OnInit(self):
        return True
        
INITIAL_SIZE = (1000, 800)

if __name__ == '__main__':
    app = VolPlotApp(size=INITIAL_SIZE)
    app.MainLoop()
