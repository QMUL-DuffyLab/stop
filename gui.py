# callum - 30/07/2026
# -*- coding: utf-8 -*-
# JsonModel taken from
# doc.qt.io/qtforpython-6/examples/example_widgets_itemviews_jsonmodel.html

import sys
import os
import random
import numpy as np
from PyQt6 import QtCore, QtWidgets, QtGui
from PyQt6.QtCore import Qt
import json
from PyQt6.QtWidgets import *
from typing import Any

'''

ProteinBuilder here is a wizard that takes the user through
the various quantities they need to define in order for the
TCSPC simulation to run. The problem is that
a.) many of the quantities depend on the values of previous quantities,
which would naturally suggest a QWizard, but
b.) the data is too annoying and heterogenous to really fit naturally
in the field mechanism of QWizards. several quantities are matrices
whose size is given at runtime by the user, and so on. so:
i'm instead overriding the validatePage() method to add various
quantities to the parent QWizard when the user presses next on each screen

'''

class loadExisting(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Load existing protein data.")
        self.setSubTitle('''
If you'd like to import existing protein data, you can do that here.
Otherwise, click Next to start specifying the parameters of your protein.
                         ''')

    def initializePage(self):
        layout = QVBoxLayout()
        gl = QGridLayout()
        self.filename = QLineEdit(os.getcwd())
        self.load_success = False
        self.browseButton = QPushButton("Browse") 
        self.browseButton.setToolTip("Search for a JSON file")
        self.browseButton.clicked.connect(self.onBrowseButton)
        self.loadButton = QPushButton("Load") 
        self.loadButton.setToolTip("Click to import JSON data from this file")
        self.loadButton.clicked.connect(self.onLoadButton)
        self.proteinChooser = QComboBox()
        self.resetButton = QPushButton("Reset")
        self.resetButton.setToolTip("Delete loaded JSON data and start again")
        self.resetButton.clicked.connect(self.onResetButton)
        gl.addWidget(QLabel("Filename:"), 0, 0)
        gl.addWidget(self.filename, 0, 1)
        gl.addWidget(self.browseButton, 0, 2)
        gl.addWidget(self.loadButton, 0, 3)
        self.plabel = QLabel("Protein name:")
        self.plabel.setToolTip('''Give your protein a short, descriptive name.
Will be used to generate output directory structure.''')
        gl.addWidget(self.plabel, 1, 0)
        gl.addWidget(self.proteinChooser, 1, 1)
        gl.addWidget(self.resetButton, 1, 2)
        layout.addLayout(gl)
        self.setLayout(layout)

    def load_from_file(self):
        print(self.filename.text())
        with open(self.filename.text()) as f:
            try:
                self.data = json.load(f)
                print(self.data)
                success = True
            except:
                self.data = {}
                print("JSON load failed.")
                # TODO: make a QMessageBox for this
                success = False
        return success

    def onBrowseButton(self):
        self.fn, _ = QFileDialog.getOpenFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.filename.setText(self.fn)

    def onLoadButton(self):
        self.load_success = self.load_from_file()
        if self.load_success:
            protein_names = self.data.keys()
            for name in protein_names:
                self.proteinChooser.addItem(name)
            # make a QMessageBox here explaining if it fails

    def onResetButton(self):
        self.load_success = False
        self.filename.setText("")
        self.proteinChooser.clear()
        self.data = {}
        self.parent.data = {}

    def updateData(self):
        if self.load_success:
            name = self.proteinChooser.currentText()
            self.parent.data = self.data[name]
            self.parent.data['name'] = name
            self.parent.data['filename'] = self.filename.text()
        else:
            self.parent.data = {}

    def validatePage(self):
        self.updateData()
        print(f"LE exit: data = {self.parent.data}")
        return True

class nameNumber(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        layout = QVBoxLayout()
        gl = QGridLayout()
        self.protein_name = QLineEdit()
        self.registerField('protein_name', self.protein_name, "text")
        self.n_p = QSpinBox()
        self.registerField('n_p', self.n_p)
        self.n_s = QSpinBox()
        self.registerField('n_s', self.n_s)
        gl.addWidget(QLabel("Protein name:"), 0, 0)
        gl.addWidget(self.protein_name, 0, 1)
        gl.addWidget(QLabel("Number of pigments:"), 1, 0)
        gl.addWidget(self.n_p, 1, 1)
        gl.addWidget(QLabel("Number of states:"), 2, 0)
        gl.addWidget(self.n_s, 2, 1)
        layout.addLayout(gl)
        self.setLayout(layout)

    def initializePage(self):
        if 'name' in self.parent.data.keys():
            self.protein_name.setText(self.parent.data['name'])
        else:
            self.protein_name.setText("")
        if 'n_p' in self.parent.data.keys():
            self.n_p.setValue(self.parent.data['n_p'])
        else:
            self.n_p.setValue(0)
        if 'n_s' in self.parent.data.keys():
            self.n_s.setValue(self.parent.data['n_s'])
        else:
            self.n_s.setValue(0)
        self.n_s.setRange(self.n_p.value(), 20)

    def cleanupPage(self):
        self.protein_name.setText("")
        self.n_p.setValue(0)
        self.n_s.setValue(0)
        if getattr(self, "fields", False):
            for field in self.fields:
                del self.parent.data[field]

    def updateData(self):
        '''
        update the parent QWizard's data struct with the
        data that's been entered here. also keep track of
        the names of the fields so that they can be cleaned
        up by cleanupPage()
        '''
        self.parent.data["name"] = self.field('protein_name')
        self.parent.data["n_p"]  = self.field('n_p')
        self.parent.data["n_s"]  = self.field('n_s')
        self.fields = ["name", "n_p", "n_s"]

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []
        if len(self.parent.data["name"]) == 0:
            msgs.append("Length of protein name must be > 0.")
            validated = False
        if self.parent.data["n_p"] <= 0:
            msgs.append("Number of pigments must be > 0.")
            validated = False
        if self.parent.data["n_s"] <= 0:
            msgs.append("Number of states must be > 0.")
            validated = False
        return validated, msgs

    def validatePage(self):
        '''
        update the parent's data, check it, print the current
        dict for my benefit, then carry on if all is well
        '''
        self.updateData()
        validated, msgs = self.checkData()
        if not validated:
            self.errors = QMessageBox.critical(self,
            "whoospy daisy", ('\n').join(msgs))
        print(f"NN exit: data = {self.parent.data}")
        return validated

class namePigmentsStates(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.layout = QVBoxLayout()
        self.pl = QGridLayout()
        self.sl = QGridLayout()
        self.nl = QGridLayout()
        self.layout.addLayout(self.pl)
        self.layout.addLayout(self.sl)
        self.layout.addLayout(self.nl)
        self.setLayout(self.layout)

    def initializePage(self):
        self.n_p = self.field('n_p')
        self.n_s = self.field('n_s')
        self.pigment_names = []
        self.state_names = []
        self.n_tot = []
        self.n_thermal = []
        for i in range(self.n_p):
            self.pl.addWidget(QLabel(f"Name of pigment {i + 1:d}:"), i, 0)
            current_name = QLineEdit()
            self.pigment_names.append(current_name)
            self.pl.addWidget(current_name, i, 1)
        for i in range(self.n_s):
            self.sl.addWidget(QLabel(f"Name of state {i + 1:d}:"), i, 0)
            current_state = QLineEdit()
            self.state_names.append(current_state)
            self.sl.addWidget(current_state, i, 1)
        self.nl.addWidget(QLabel("Total number of pigments"), 0, 1)
        self.nl.addWidget(QLabel("Thermally accessible pigments"), 0, 2)
        for i in range(self.n_p):
            self.nl.addWidget(QLabel(f"Pigment {i + 1:d}:"), i + 1, 0)
            self.n_tot.append(QSpinBox())
            self.n_thermal.append(QSpinBox())
            self.nl.addWidget(self.n_tot[i], i + 1, 1)
            self.nl.addWidget(self.n_thermal[i], i + 1, 2)
        '''
        check if there's data loaded and fill values if so
        '''
        keys = ["pigment_names", "state_names"]
        boxlists = [self.pigment_names, self.state_names]
        for k, b in zip(keys, boxlists):
            if k in self.parent.data:
                names = self.parent.data[k]
            else:
                names = ["" for _ in range(self.n_p)]
            for i, n in enumerate(names):
                b[i].setText(n)
        keys = ["n_tot", "n_thermal"]
        boxlists = [self.n_tot, self.n_thermal]
        for k, b in zip(keys, boxlists):
            if k in self.parent.data:
                vals = self.parent.data[k]
            else:
                vals = [0 for _ in range(self.n_p)]
            for i, n in enumerate(vals):
                b[i].setValue(n)

    def cleanupPage(self):
        self.n_p = 0
        self.n_s = 0
        self.pigment_names = []
        self.state_names = []
        self.n_tot = []
        self.n_thermal = []
        for item in self.pigment_names:
            item.setText("")
        for item in self.state_names:
            item.setText("")
        for layout in self.pl, self.sl, self.nl:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget:
                    child.widget().deleteLater()
        if getattr(self, "fields", False):
            for field in self.fields:
                del self.parent.data[field]

    def updateData(self):
        self.parent.data["pigment_names"] = [p.text()
                                        for p in self.pigment_names]
        self.parent.data["state_names"]   = [s.text()
                                        for s in self.state_names]
        self.parent.data["n_tot"] = [int(p.value()) for p in self.n_tot]
        self.parent.data["n_thermal"] = [int(p.value())
                                         for p in self.n_thermal]
        self.fields = ["pigment_names", "state_names", "n_tot", "n_thermal"]

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []
        data = self.parent.data

        for name_string, n_string in zip(
            ['pigment_names', 'state_names'], ['n_p', 'n_s']):
            nsmsg = name_string.replace("_", " ")
            if any([len(p) == 0 for p in data[name_string]]):
                msgs.append(f"Length of {nsmsg} must be > 0.")
                validated = False
            if len(data[name_string]) != data[n_string]:
                msgs.append(f"Number of {nsmsg} doesn't match {n_string}.")
                validated = False
        for i in range(data['n_p']):
            if self.n_tot[i].value() <= 0:
                msgs.append(f"Number of pigments for pigment {i} must be > 0.")
                validated = False
            if self.n_thermal[i].value() <= 0:
                msgs.append(f"Number of thermally accessible pigments for pigment {i} must be > 0.")
                validated = False
            if self.n_tot[i].value() < self.n_thermal[i].value():
                msgs.append(f"Number of thermally accessible states for pigment {i} is larger than total.")
                validated = False
        return validated, msgs

    def validatePage(self):
        '''
        update the parent's data, check it, print the current
        dict for my benefit, then carry on if all is well
        '''
        self.updateData()
        validated, msgs = self.checkData()
        if not validated:
            self.errors = QMessageBox.critical(self,
            "whoospy daisy", ('\n').join(msgs))
        print(f"NPS exit: data = {self.parent.data}")
        return validated

class pigmentProperties(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.layout = QVBoxLayout()
        self.pl = QGridLayout()
        self.layout.addLayout(self.pl)
        self.setLayout(self.layout)
        self.n_s = 0
        # column headers
        self.hopLabel = QLabel("Hopping time (s)")
        self.hopLabel.setToolTip("The hopping time for each state from one protein to its neighbours, in seconds. e.g. for 1ps, enter 1e-12.")
        self.pl.addWidget(self.hopLabel, 0, 1)
        self.pl.addWidget(QLabel("Decay time (s)"), 0, 2)
        self.pl.addWidget(QLabel("Cross-section (cm^{-1})"), 0, 3)
        self.emissiveLabel = QLabel("Emissive decay?")
        self.emissiveLabel.setToolTip("At least one decay must be emissive; that is, visible to the detector. Multiple boxes can be checked here if there are multiple decay pathways.")
        self.pl.addWidget(self.emissiveLabel, 0, 4)
        self.pigmentLabel = QLabel("Pigment") 
        self.pigmentLabel.setToolTip("Which pigment does each state belong to?")
        self.pl.addWidget(self.pigmentLabel, 0, 5)
        self.abundanceLabel = QLabel("Abundance") 
        self.pigmentLabel.setToolTip("What fraction of sites have this state present?")
        self.pl.addWidget(self.abundanceLabel, 0, 6)

    def initializePage(self):
        self.n_s = self.field("n_s")
        self.hop       = []
        self.decay     = []
        self.xsec      = []
        self.emissive  = []
        self.which_p   = []
        self.abundance = []
        for i in range(self.n_s):
            row = i + 1
            state_name = self.parent.data["state_names"][i]
            self.hop.append(QLineEdit("0.0"))
            self.decay.append(QLineEdit("0.0"))
            self.xsec.append(QLineEdit("0.0"))
            self.emissive.append(QCheckBox())
            self.which_p.append(QComboBox())
            self.abundance.append(QLineEdit("0.0"))
            self.pl.addWidget(QLabel(state_name), row, 0)
            self.pl.addWidget(self.hop[i], row, 1)
            self.pl.addWidget(self.decay[i], row, 2)
            self.pl.addWidget(self.xsec[i], row, 3)
            self.pl.addWidget(self.emissive[i], row, 4)
            # need to find out how to centre the checkboxes. it's annoying
            self.pl.setAlignment(self.emissive[i],
                                 Qt.AlignmentFlag.AlignHCenter)
            for j in range(self.field("n_p")):
                name = self.parent.data["pigment_names"][j]
                self.which_p[i].addItem(name)
            self.pl.addWidget(self.which_p[i], row, 5)
            self.pl.addWidget(self.abundance[i], row, 6)
        '''
        check if there's data loaded and fill values if so
        '''
        keys = ["hop", "xsec", "abundance"]
        boxlists = [self.hop, self.xsec, self.abundance]
        for k, b in zip(keys, boxlists):
            if k in self.parent.data:
                names = self.parent.data[k]
            else:
                if k == 'abundance':
                    names = [1.0 for _ in range(self.n_s)]
                else:
                    names = [0.0 for _ in range(self.n_s)]
            for i, n in enumerate(names):
                b[i].setText(str(n))
        if "intra" in self.parent.data:
            for i in range(self.n_s):
                self.decay[i].setText(str(self.parent.data["intra"][i][i]))
        if "emissive" in self.parent.data:
            ea = self.parent.data["emissive"]
        else:
            ea = [False for _ in range(self.n_s)]
        if "which_pigment" in self.parent.data:
            which = self.parent.data["which_pigment"]
        else:
            # python's 0-based and fortran is 1-based
            # so the conversion has to be done somewhere; i do it
            # on the python side and put the 1-based arrays in JSON
            which = [1 for _ in range(self.n_s)]
        for i in range(self.n_s):
            self.emissive[i].setChecked(ea[i])
            self.which_p[i].setCurrentIndex(which[i] - 1)

    def cleanupPage(self):
        while self.pl.count():
            child = self.pl.takeAt(0)
            if child.widget:
                child.widget().deleteLater()
        self.hop       = []
        self.decay     = []
        self.xsec      = []
        self.emissive  = []
        self.which_p   = []
        self.abundance = []
        self.n_s = 0
        self.n_p = 0
        if getattr(self, "fields", False):
            for field in self.fields:
                del self.parent.data[field]

    def updateData(self):
        self.parent.data["hop"] = [float(p.text())
                            if p.text() != '' else 0.0 for p in self.hop]
        self.parent.data["decay"] = [float(p.text())
                            if p.text() != '' else 0.0 for p in self.decay]
        self.parent.data["xsec"]  = [float(p.text())
                            if p.text() != '' else 0.0 for p in self.xsec]
        self.parent.data["emissive"] = [p.isChecked() for p in self.emissive]
        self.parent.data["which_pigment"] = [p.currentIndex() + 1
                                      for p in self.which_p]
        dist = []
        which = self.parent.data["which_pigment"]
        for i in range(self.parent.data['n_s']):
            row = []
            for j in range(self.parent.data['n_s']):
               row.append(False if which[i] == which[j] else True) 
            dist.append(row)
        self.parent.data["dist"] = dist
        self.parent.data["abundance"]  = [float(p.text())
                        if p.text() != '' else 0.0 for p in self.abundance]
        self.fields = ["hop", "decay", "xsec", "emissive",
                       "which_pigment", "abundance"]

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []
        data = self.parent.data
        if not(any(data["emissive"])):
            msgs.append("At least one decay should be emissive!")
            validated = False
        if all([d == 0.0 for d in data['xsec']]):
            msgs.append("At least one state should have a non-zero cross section.")
            validated = False
        if all([d == 0.0 for d in data['abundance']]):
            msgs.append("At least one state should have a non-zero abundance.")
            validated = False
        return validated, msgs
        
    def validatePage(self):
        '''
        update the parent's data, check it, print the current
        dict for my benefit, then carry on if all is well
        '''
        self.updateData()
        validated, msgs = self.checkData()
        if not validated:
            self.errors = QMessageBox.critical(self,
            "whoospy daisy", ('\n').join(msgs))
        print(f"PP exit: data = {self.parent.data}")
        return validated

class MatrixTable(QTableWidget):
    def __init__(self, data, key, symmetric=False, block_diagonal=False):
        QTableWidget.__init__(self)
        self.state_names = data["state_names"]
        self.n_states = len(self.state_names)
        self.data = data
        self.setRowCount(self.n_states)
        self.setColumnCount(self.n_states)
        self.setHorizontalHeaderLabels(self.state_names)
        self.setVerticalHeaderLabels(self.state_names)
        for i in range(self.n_states):
            for j in range(self.n_states):
                item = QTableWidgetItem()
                if key in self.data:
                    item.setText(str(self.data[key][i][j]))
                self.setItem(i, j, item)
        if symmetric:
            for i in range(self.n_states):
                for j in range(self.n_states):
                    if i > j:
                        item = QtWidgets.QTableWidgetItem()
                        item.setBackground(QtGui.QColor("darkGray"))
                        item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
                        self.setItem(i, j, item)
        if block_diagonal:
            for i in range(self.n_states):
                item = QtWidgets.QTableWidgetItem()
                item.setBackground(QtGui.QColor("darkGray"))
                item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
                self.setItem(i, i, item)

class RemainderTable(QGridLayout):
    def __init__(self, data):
        QGridLayout.__init__(self)
        self.state_names = data["state_names"]
        self.n_states = len(self.state_names)
        self.data = data
        for i in range(self.n_states):
            self.addWidget(QLabel(self.state_names[i]), 0, i + 1)
            self.addWidget(QLabel(self.state_names[i]), i + 1, 0)
            for j in range(self.n_states):
                current = QtWidgets.QComboBox()
                current.addItem("None")
                for k in range(self.n_states):
                    current.addItem(self.state_names[k])
                self.addWidget(current, i + 1, j + 1)
                if "ann_remainder" in data:
                    # 0 is the None index is the item list
                    # so we don't need to mess about here
                    ci = self.data["ann_remainder"][i][j]
                    current.setCurrentIndex(ci)

    def cleanup(self):
        while self.count():
            child = self.takeAt(0)
            if child.widget:
                child.widget().deleteLater()

class matrixTables(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.layout = QVBoxLayout()

    def initializePage(self):
        self.intra = MatrixTable(self.parent.data, "intra", block_diagonal=True)
        self.ann   = MatrixTable(self.parent.data, "ann", symmetric=True)
        self.ann_rem_layout = RemainderTable(self.parent.data)
        self.intra_label = QLabel("Transfer times between states (s)")
        self.layout.addWidget(self.intra_label)
        self.layout.addWidget(self.intra)
        self.ann_label = QLabel("Annihilation times between states (s)")
        self.layout.addWidget(self.ann_label)
        self.layout.addWidget(self.ann)
        self.ann_rem_label = QLabel("Remaining state after annihilation event")
        self.layout.addWidget(self.ann_rem_label)
        self.layout.addLayout(self.ann_rem_layout)
        self.setLayout(self.layout)

    def cleanupPage(self):
        self.ann_rem_layout.cleanup()
        self.intra  = []
        self.ann    = []
        self.intra_label = None
        self.ann_label = None
        if getattr(self, "fields", False):
            for field in self.fields:
                del self.parent.data[field]

    def updateData(self):
        n_s = self.parent.data["n_s"]
        intra = np.zeros((n_s, n_s), dtype=float)
        ann = np.zeros_like(intra)
        ann_rem = np.zeros((n_s, n_s), dtype=int)
        for i in range(n_s):
            for j in range(n_s):
                if i == j:
                    intra[i, i] = self.parent.data['decay'][i]
                else:
                    intra[i, j] = float(self.intra.item(i, j).text())
                if j >= i:
                    # annihilation matrix must be symmetric
                    ann[i, j] = float(self.ann.item(i, j).text())
                    ann[j, i] = float(self.ann.item(i, j).text())
                current = self.ann_rem_layout.itemAtPosition(i + 1, j + 1)
                # None is always added as the first (0) index
                ann_rem[i, j] = current.widget().currentIndex()
        self.parent.data["intra"] = intra.tolist()
        self.parent.data["ann"] = ann.tolist()
        self.parent.data["ann_remainder"] = ann_rem.tolist()
        self.fields = ["intra", "ann", "ann_remainder"]
        print(f"MT exit: data = {self.parent.data}")

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []
        data = self.parent.data
        intra = data["intra"]
        ann = data["ann"]
        ann_rem = data["ann_remainder"]
        for i in range(data["n_s"]):
            for j in range(data["n_s"]):
                s1 = data["state_names"][i]
                s2 = data["state_names"][j]
                if intra[i][j] < 0.0:
                    msgs.append(f"Transfer times from {s1} to {s2} cannot be < 0.")
                    validated = False
                if ann[i][j] < 0.0:
                    msgs.append(f"Annihilation time for {s1} and {s2} cannot be < 0.")
                    validated = False
                if ann[i][j] != ann[j][i]:
                    msgs.append("Annihilation matrix is not symmetric.")
                    validated = False
                if ann_rem[i][j] < 0:
                    msgs.append(f"Annihilation remainder for {s1} and {s2} invalid.")
                    validated = False
                if ann[i][j] == 0.0 and ann_rem[i][j] > 0:
                    msgs.append(f"States {s1} and {s2} have an annihilation remainder set but a zero annihilation rate.")
                    validated = False
        return validated, msgs

    def validatePage(self):
        '''
        update the parent's data, check it, print the current
        dict for my benefit, then carry on if all is well
        '''
        self.updateData()
        validated, msgs = self.checkData()
        if not validated:
            self.errors = QMessageBox.critical(self,
            "whoospy daisy", ('\n').join(msgs))
        print(f"PP exit: data = {self.parent.data}")
        return validated

class savePage(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Save")
        self.setSubTitle("Save protein data to file.")

    def initializePage(self):
        layout = QVBoxLayout()
        gl = QGridLayout()
        if 'filename' in self.parent.data:
            self.filename = QLineEdit(self.parent.data['filename'])
        else:
            self.filename = QLineEdit(os.getcwd())
        self.save_success = False
        self.load_success = False
        self.existing_data = {}
        self.browseButton = QPushButton("Browse") 
        self.browseButton.setToolTip("Search for a JSON file")
        self.browseButton.clicked.connect(self.onBrowseButton)
        self.saveButton = QPushButton("Save") 
        self.saveButton.setToolTip("Click to save JSON data to this file")
        self.saveButton.clicked.connect(self.onSaveButton)
        gl.addWidget(QLabel("Filename:"), 0, 0)
        gl.addWidget(self.filename, 0, 1,)
        gl.addWidget(self.browseButton, 0, 2)
        gl.addWidget(self.saveButton, 0, 3)
        layout.addLayout(gl)
        self.setLayout(layout)

    def save_to_file(self):
        dd = self.parent.data
        name = dd.pop('name')
        final_data = {name: dd}
        print(final_data)
        dd['name'] = name
        # don't need these in the JSON
        print(final_data)
        if 'filename' in final_data:
            # filename is only a key if a file was loaded at the start
            del final_data[name]['filename']
        del final_data[name]['decay']
        success = True
        self.existing_data = {}
        # if the filename exists, try to open it and parse the JSON
        if os.path.isfile(self.filename.text()):
            with open(self.filename.text(), "r+", encoding='utf-8') as f:
                try:
                    self.existing_data = json.load(f)
                    print(self.existing_data)
                    self.load_success = True
                except:
                    print("Failed to load existing protein data from JSON.")
                    self.load_success = False
        # if the protein name matches one that's already there and we just
        # merge the dicts, the original will be overwritten, so check
        if self.parent.data["name"] in self.existing_data.keys():
            overwrite = True
            self.overwriteCheck = QMessageBox.question(self,
                "", "Protein name already exists in data file. Overwrite?")

            if self.overwriteCheck == QMessageBox.StandardButton.NoButton:
                overwrite = False

            if overwrite:
                self.total_data = self.existing_data | final_data
            else:
                success = False
        else:
            final_data = self.existing_data | final_data
            with open(self.filename.text(), "w") as f:
                json.dump(final_data, f)
        return success

    def onBrowseButton(self):
        self.fn, _ = QFileDialog.getSaveFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.filename.setText(self.fn)

    def onSaveButton(self):
        self.save_success = self.save_to_file()

    def validatePage(self):
        return self.save_success

class ProteinBuilder(QWizard):
    def __init__(self):
        super().__init__()
        self.data = {}
        self.addPage(loadExisting(self))
        self.addPage(nameNumber(self))
        self.addPage(namePigmentsStates(self))
        self.addPage(pigmentProperties(self))
        self.addPage(matrixTables(self))
        self.addPage(savePage(self))
        self.setWindowTitle("Protein builder for STOP")

if __name__ == "__main__":
    app =QtWidgets.QApplication([])

    widget = ProteinBuilder()
    widget.show()

    sys.exit(app.exec())
