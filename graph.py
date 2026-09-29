# -*- coding: utf-8 -*-
# matplotlib.org/stable/gallery/user_interfaces/embedding_in_qt_sgskip.html
import sys
import time
import re
import os
import functools
import collections
import numpy as np
import pandas as pd

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.figure import Figure

def get_files(path, suffix):
    '''
    for a given path and suffix, pull all the matching filenames
    and collect them into a dict grouped by run number
    '''
    pattern = r'.*run_([0-9]+)_proc_([0-9]+)_salt_([0-9]+)_' + suffix
    matches = [re.search(pattern, f) for f in os.listdir(path)
            if re.search(pattern, f)] # re.search() returns None if no match
    d = collections.defaultdict(list)
    for m, r in [(m.group(0), m.group(1)) for m in matches]:
        d[int(r)].append(m)
    return d

def reduce_files(path, suffix, run, xcol, **pd_kwargs):
    '''
    for a given run, read in all the files and reduce all the columns.
    do not reduce the column named "xcol" in the returned dataframe
    (the time or the rep number shouldn't be summed over processes)
    NB: for population simply do (path, "population.csv", run, "Rep")
        for a hist: (path, "rep_N.csv", run, "Time(s)", skiprows=[1])
        where N is the rep number you want, or "final.csv" for the final ones.
    '''
    filedict = get_files(path, suffix)
    ll = []
    files = filedict[run]
    for f in files:
        ll.append(pd.read_csv(os.path.join(path, f), sep=r'\s+', **pd_kwargs))
    xdata = ll[0][xcol].to_numpy()
    df = functools.reduce(lambda x, y: x.add(y, fill_value=np.nan), ll)
    return df.div(len(files))

class GraphWindow(QtWidgets.QMainWindow):
    def __init__(self, folder, active):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.folder = folder
        self.active = active
        '''
        should have: a load folder option, or to be able to pass the
        active folder to the constructor if a simulation is running.
        maybe watchdog to watch the folder and regexp matches to
        find the per-run files as needed? i think one tab per run
        and then each tab should have a top panel for the population
        of each state as a function of rep (reduced across processes),
        with a vline at burn_reps. then below that a larger panel with
        the current histogram of all processes, updated when the fortran
        updates it? hmmmm
        '''
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)

        pop_canvas = FigureCanvas(Figure(figsize=(9, 3)))
        # Ideally one would use self.addToolBar here, but it is slightly
        # incompatible between PyQt6 and other bindings, so we just add the
        # toolbar as a plain widget instead.
        layout.addWidget(NavigationToolbar(pop_canvas, self))
        layout.addWidget(pop_canvas)

        dynamic_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(dynamic_canvas)
        layout.addWidget(NavigationToolbar(dynamic_canvas, self))

        self._pop_ax = pop_canvas.figure.subplots()
        t = np.linspace(0, 10, 501)
        self._pop_ax.plot(t, np.tan(t), ".")

        self._dynamic_ax = dynamic_canvas.figure.subplots()
        # Set up a Line2D.
        self.xdata = np.linspace(0, 10, 101)
        self._update_ydata()
        self._line, = self._dynamic_ax.plot(self.xdata, self.ydata)
        # The below two timers must be attributes of self, so that the garbage
        # collector won't clean them after we finish with __init__...

        # The data retrieval may be fast as possible (Using QRunnable could be
        # even faster).
        self.data_timer = dynamic_canvas.new_timer(1)
        self.data_timer.add_callback(self._update_ydata)
        self.data_timer.start()
        # Drawing at 50Hz should be fast enough for the GUI to feel smooth, and
        # not too fast for the GUI to be overloaded with events that need to be
        # processed while the GUI element is changed.
        self.drawing_timer = dynamic_canvas.new_timer(20)
        self.drawing_timer.add_callback(self._update_canvas)
        self.drawing_timer.start()

    def plot_pop_df(df, ax):
        for c in df.columns:
            if c not in ['Rep', 'rep_end_time']:
                ax.plot(df['Rep'], df[c], label=c)

    def _update_ydata(self):
        # Shift the sinusoid as a function of time.
        self.ydata = np.sin(self.xdata + time.time())

    def _update_canvas(self):
        self._line.set_data(self.xdata, self.ydata)
        # It should be safe to use the synchronous draw() method for most drawing
        # frequencies, but it is safer to use draw_idle().
        self._line.figure.canvas.draw_idle()
