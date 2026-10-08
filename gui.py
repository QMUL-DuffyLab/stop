# callum - 30/07/2026
# -*- coding: utf-8 -*-
# shout out to
# pythonguis.com/faq/constantly-print-subprocess-output-while-process-is-running/
# for the implementation of runPage.run() at the bottom here

import sys
import os
import random
import numpy as np
from PyQt6 import QtCore, QtWidgets, QtGui, uic
from PyQt6.QtCore import Qt
import json
from PyQt6.QtWidgets import *
from typing import Any
import parse

'''

STOPSetup is a wizard that takes the user through
the various quantities they need to define in order for the
TCSPC simulation to run. The problem is that
a.) many of the quantities depend on the values of previous quantities,
which would naturally suggest a QWizard, but
b.) the data is too annoying and heterogenous to really fit naturally
in the field mechanism of QWizards. several quantities are matrices
whose size is given at runtime by the user, and so on. so:
i've added a data member to the parent QWizard class and each page adds
relevant quantities to that as the user goes; I've added updateData()
and checkData() methods to each page to validate that data as necessary,
checkData() returns a bool along with any error messages,
and validatePage() is overridden to check that bool and print the error
messages in a QMessageBox if there are any.

TODO:
    - maybe (MAYBE) get the fortran to print out intermediate histograms
      as it goes, and then plot them in a separate window, possibly along
      with the population per rep (checked at intervals; could do this
      by checking when runPage.output_area is updated, since the fortran
      prints to stdout every 100 reps).
    - test suite. write a bunch of toy protein and simulation JSON files
      which either should or shouldn't parse, put them in a tests dir,
      glob them and run them through the parsers one by one. GUI testing
      for things like going back and forth through the wizard is harder
      to standardise, I've been trying to test as i go
'''

def load_from_file(widget, filename):
    data = {}
    success = False
    if os.path.isfile(filename):
        with open(filename) as f:
            try:
                data = json.load(f)
                success = True
            except:
                box = QMessageBox.critical(widget,
                "JSON import failed", f"Failed to load data from {filename}.")
    else:
        box = QMessageBox.critical(widget,
        "Load failed", f"Filename {filename} does not exist.")
    return success, data

main_ui_class, main_ui_widget = uic.loadUiType("main.ui")

class loadExisting(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Load existing protein data.")
        self.setSubTitle("If you'd like to import existing protein data "
                "to edit, you can do that here. Otherwise, click Next "
                "to start specifying the parameters of your protein.")
        self.data = self.parent.protein_data

    def initializePage(self):
        layout = QVBoxLayout()
        gl = QGridLayout()
        fstr = os.path.join(os.getcwd(), "protein.json")
        self.filename = QLineEdit(fstr)
        self.load_success = False
        self.browseButton = QPushButton("Browse") 
        self.browseButton.setToolTip("Search for a JSON file")
        self.browseButton.clicked.connect(self.onBrowseButton)
        self.proteinChooser = QComboBox()
        self.resetButton = QPushButton("Reset")
        self.resetButton.setToolTip("Delete loaded JSON data and start again")
        self.resetButton.clicked.connect(self.onResetButton)
        gl.addWidget(QLabel("Filename:"), 0, 0)
        gl.addWidget(self.filename, 0, 1)
        gl.addWidget(self.browseButton, 0, 2)
        self.plabel = QLabel("Protein name:")
        self.plabel.setToolTip('''Give your protein a short, descriptive name.
Will be used to generate output directory structure.''')
        gl.addWidget(self.plabel, 1, 0)
        gl.addWidget(self.proteinChooser, 1, 1)
        gl.addWidget(self.resetButton, 1, 2)
        layout.addLayout(gl)
        self.setLayout(layout)

    def onBrowseButton(self):
        self.proteinChooser.clear()
        self.fn, _ = QFileDialog.getOpenFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.filename.setText(self.fn)
        self.load_success, self.all_json_data = load_from_file(self,
                self.filename.text())
        if self.load_success:
            protein_names = self.all_json_data.keys()
            for name in protein_names:
                self.proteinChooser.addItem(name)

    def onResetButton(self):
        self.load_success = False
        self.filename.setText("")
        self.proteinChooser.clear()
        self.all_json_data = {}
        self.parent.protein_data = {}

    def updateData(self):
        if self.load_success:
            name = self.proteinChooser.currentText()
            self.data = self.all_json_data[name]
            self.parent.protein_data = self.data
            self.parent.protein_name = name
            self.parent.protein_file = self.filename.text()
        else:
            self.data = {}

    def validatePage(self):
        self.updateData()
        return True

class nameNumber(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.setTitle("Name and basic protein properties")
        self.setSubTitle("Enter the name of the protein "
                "if not yet given, the number of different pigments, "
                "and the number of states in total across them.")
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
        self.data = self.parent.protein_data
        self.protein_name.setText(self.parent.protein_name)
        if 'n_p' in self.data.keys():
            self.n_p.setValue(self.data['n_p'])
            self.n_s.setRange(self.n_p.value(), 20)
        if 'n_s' in self.data.keys():
            self.n_s.setValue(self.data['n_s'])

    def cleanupPage(self):
        self.protein_name.setText("")
        self.n_p.setValue(0)
        self.n_s.setRange(0, 20)
        self.n_s.setValue(0)
        for key in self.updated_keys:
            if key in self.data:
                del self.data[key]
        self.parent.protein_data = self.data
        self.updated_keys = []

    def updateData(self):
        '''
        update the parent QWizard's data struct with the
        data that's been entered here. also keep track of
        the names of the fields so that they can be cleaned
        up by cleanupPage()
        '''
        self.parent.protein_name = self.field('protein_name')
        self.data["n_p"]  = self.field('n_p')
        self.data["n_s"]  = self.field('n_s')
        self.parent.protein_data = self.data
        self.updated_keys = ["n_p", "n_s"]

    def checkData(self):
        return parse.parse_protein(self.parent.protein_data,
                keys=self.updated_keys)

    def validatePage(self):
        self.updateData()
        if self.parent.protein_name == "":
            self.errors = QMessageBox.critical(self,
            self.title(), "Protein name cannot be blank.")
            return False
        valid, msgs = self.checkData()
        if not valid:
            self.errors = QMessageBox.critical(self,
            self.title(), ('\n').join(msgs))
        return valid

class namePigmentsStates(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.setTitle("Names of pigments and states")
        self.setSubTitle("These names will be used in the "
            "output files from the fortran. Short names with "
            "no spaces recommended, otherwise the fitting code will "
            "get confused about which histogram column is which."
                )
        self.parent = parent
        self.updated_keys = []
        self.layout = QVBoxLayout()
        self.pl = QGridLayout()
        self.sl = QGridLayout()
        self.nl = QGridLayout()
        self.layout.addLayout(self.pl)
        self.layout.addLayout(self.sl)
        self.layout.addLayout(self.nl)
        self.setLayout(self.layout)

    def initializePage(self):
        self.data = self.parent.protein_data
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
            if k in self.data:
                names = self.data[k]
            else:
                names = ["" for _ in range(self.n_p)]
            for i, n in enumerate(names):
                b[i].setText(n)
        keys = ["n_tot", "n_thermal"]
        boxlists = [self.n_tot, self.n_thermal]
        for k, b in zip(keys, boxlists):
            if k in self.data:
                vals = self.data[k]
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
        for key in self.updated_keys:
            if key in self.data:
                del self.data[key]
        self.parent.protein_data = self.data
        self.updated_keys = []

    def updateData(self):
        self.data["pigment_names"] = [p.text() for p in self.pigment_names]
        self.data["state_names"]   = [s.text() for s in self.state_names]
        self.data["n_tot"]         = [int(p.value()) for p in self.n_tot]
        self.data["n_thermal"]     = [int(p.value()) for p in self.n_thermal]
        self.parent.protein_data = self.data
        self.updated_keys = ["pigment_names", "state_names",
                "n_tot", "n_thermal"]

    def checkData(self):
        return parse.parse_protein(self.parent.protein_data,
                keys=self.updated_keys)

    def validatePage(self):
        '''
        update the parent's data, check it, print the current
        dict for my benefit, then carry on if all is well
        '''
        self.updateData()
        valid, msgs = self.checkData()
        if not valid:
            self.errors = QMessageBox.critical(self,
            self.title(), ('\n').join(msgs))
        return valid

class stateProperties(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.setTitle("State properties.")
        self.setSubTitle("Per-state properties like decay time, "
        "cross-section, etc.")
        self.parent = parent
        self.updated_keys = []
        self.layout = QVBoxLayout()
        self.pl = QGridLayout()
        self.layout.addLayout(self.pl)
        self.setLayout(self.layout)
        self.n_s = 0
        # column headers
        self.hopLabel = QLabel("Hopping time (s)")
        self.hopLabel.setToolTip(
'''The hopping time for each state from one protein to its neighbours,
in seconds. e.g. for 1ps, enter 1e-12.''')
        self.pl.addWidget(self.hopLabel, 0, 1)
        self.pl.addWidget(QLabel("Decay time (s)"), 0, 2)
        self.pl.addWidget(QLabel("Cross-section (cm^{-1})"), 0, 3)
        self.emissiveLabel = QLabel("Emissive decay?")
        self.emissiveLabel.setToolTip(
'''At least one decay must be emissive; that is, visible to the detector.
Multiple boxes can be checked here if there are multiple decay pathways.''')
        self.pl.addWidget(self.emissiveLabel, 0, 4)
        self.pigmentLabel = QLabel("Pigment") 
        self.pigmentLabel.setToolTip(
                "Which pigment does each state belong to?")
        self.pl.addWidget(self.pigmentLabel, 0, 5)
        self.abundanceLabel = QLabel("Abundance") 
        self.pigmentLabel.setToolTip(
                "What fraction of sites have this state present?")
        self.pl.addWidget(self.abundanceLabel, 0, 6)

    def initializePage(self):
        self.data = self.parent.protein_data
        self.n_s = self.field("n_s")
        self.hop       = []
        self.decay     = []
        self.xsec      = []
        self.emissive  = []
        self.which_p   = []
        self.abundance = []
        for i in range(self.n_s):
            row = i + 1
            state_name = self.data["state_names"][i]
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
                name = self.data["pigment_names"][j]
                self.which_p[i].addItem(name)
            self.pl.addWidget(self.which_p[i], row, 5)
            self.pl.addWidget(self.abundance[i], row, 6)
        '''
        check if there's data loaded and fill values if so
        '''
        keys = ["hop", "xsec", "abundance"]
        boxlists = [self.hop, self.xsec, self.abundance]
        for k, b in zip(keys, boxlists):
            if k in self.data:
                names = self.data[k]
            else:
                if k == 'abundance':
                    names = [1.0 for _ in range(self.n_s)]
                else:
                    names = [0.0 for _ in range(self.n_s)]
            for i, n in enumerate(names):
                b[i].setText(str(n))
        if "intra" in self.data:
            for i in range(self.n_s):
                self.decay[i].setText(str(self.data["intra"][i][i]))
        if "emissive" in self.data:
            ea = self.data["emissive"]
        else:
            ea = [False for _ in range(self.n_s)]
        if "which_pigment" in self.data:
            which = self.data["which_pigment"]
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
        for key in self.updated_keys:
            del self.data[key]
        self.parent.protein_data = self.data
        self.updated_keys = []

    def updateData(self):
        self.data["hop"] = [float(p.text())
                            if p.text() != '' else 0.0 for p in self.hop]
        self.data["decay"] = [float(p.text())
                            if p.text() != '' else 0.0 for p in self.decay]
        self.data["xsec"]  = [float(p.text())
                            if p.text() != '' else 0.0 for p in self.xsec]
        self.data["emissive"] = [p.isChecked() for p in self.emissive]
        self.data["which_pigment"] = [p.currentIndex() + 1
                                      for p in self.which_p]
        dist = []
        which = self.data["which_pigment"]
        for i in range(self.data['n_s']):
            row = []
            for j in range(self.data['n_s']):
               row.append(False if which[i] == which[j] else True) 
            dist.append(row)
        self.data["dist"] = dist
        self.data["abundance"]  = [float(p.text())
                        if p.text() != '' else 0.0 for p in self.abundance]
        self.parent.protein_data = self.data
        self.updated_keys = ["hop", "decay", "xsec", "emissive",
                       "which_pigment", "abundance"]

    def checkData(self):
        return parse.parse_protein(self.parent.protein_data,
                keys=self.updated_keys)
        
    def validatePage(self):
        self.updateData()
        valid, msgs = self.checkData()
        if not valid:
            self.errors = QMessageBox.critical(self,
            self.title(), ('\n').join(msgs))
        return valid

class matrixTables(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.setTitle("Intra-state rates")
        self.setSubTitle("Details of transfer and annihilation rates.")
        self.parent = parent
        self.updated_keys = []
        self.layout = QGridLayout()
        self.intra_label = QLabel("Transfer times between states (s)")
        self.ann_label = QLabel("Annihilation times between states (s)")
        self.ann_rem_label = QLabel("Remaining state after annihilation event")
        self.layout.addWidget(self.intra_label, 0, 0)
        self.layout.addWidget(self.ann_label, 1, 0)
        self.layout.addWidget(self.ann_rem_label, 2, 0)
        self.intra = QTableWidget()
        self.ann = QTableWidget()
        self.ann_rem = QGridLayout()
        self.layout.addWidget(self.intra, 0, 1)
        self.layout.addWidget(self.ann, 1, 1)
        self.layout.addLayout(self.ann_rem, 2, 1)
        self.setLayout(self.layout)

    def initializePage(self):
        self.data = self.parent.protein_data
        state_names = self.data["state_names"]
        n_s = self.data["n_s"]
        for widget, key in zip([self.intra, self.ann], ["intra", "ann"]):
            widget.setRowCount(n_s)
            widget.setColumnCount(n_s)
            widget.setHorizontalHeaderLabels(state_names)
            widget.setVerticalHeaderLabels(state_names)
            for i in range(n_s):
                for j in range(n_s):
                    item = QTableWidgetItem()
                    if key in self.data:
                        item.setText(str(self.data[key][i][j]))
                    widget.setItem(i, j, item)
            for i in range(n_s):
                # ann is symmetric - block below the diagonal
                if key == "ann":
                    for j in range(n_s):
                        if i > j:
                            item = QTableWidgetItem()
                            item.setBackground(QtGui.QColor("darkGray"))
                            item.setFlags(Qt.ItemFlag.ItemIsSelectable | 
                                          Qt.ItemFlag.ItemIsEditable)
                            widget.setItem(i, j, item)
                # intra should have diagonal blocked - decays already given
                if key == "intra":
                    item = QTableWidgetItem()
                    item.setBackground(QtGui.QColor("darkGray"))
                    item.setFlags(
                    Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
                    widget.setItem(i, i, item)
        # ann_rem
        for i in range(n_s):
            self.ann_rem.addWidget(QLabel(state_names[i]), 0, i + 1)
            self.ann_rem.addWidget(QLabel(state_names[i]), i + 1, 0)
            for j in range(n_s):
                current = QComboBox()
                current.addItem("None")
                for k in range(n_s):
                    current.addItem(state_names[k])
                self.ann_rem.addWidget(current, i + 1, j + 1)
                if "ann_remainder" in self.data:
                    # 0 is the None index in the item list
                    # so we don't need to mess about here
                    ci = self.data["ann_remainder"][i][j]
                    current.setCurrentIndex(ci)

    def cleanupPage(self):
        self.intra.clear()
        self.ann.clear()
        while self.ann_rem.count():
            child = self.ann_rem.takeAt(0)
            if child.widget:
                child.widget().deleteLater()
        for key in self.updated_keys:
            if key in self.data:
                del self.data[key]
        self.parent.protein_data = self.data
        self.updated_keys = []

    def updateData(self):
        n_s = self.data["n_s"]
        intra = np.zeros((n_s, n_s), dtype=float)
        ann = np.zeros_like(intra)
        ann_rem = np.zeros((n_s, n_s), dtype=int)
        for i in range(n_s):
            for j in range(n_s):
                if i == j:
                    intra[i, i] = self.data['decay'][i]
                else:
                    intra[i, j] = float(self.intra.item(i, j).text())
                if j >= i:
                    # annihilation matrix must be symmetric
                    ann[i, j] = float(self.ann.item(i, j).text())
                    ann[j, i] = float(self.ann.item(i, j).text())
                current = self.ann_rem.itemAtPosition(i + 1, j + 1)
                # None is always added as the first (0) index
                ann_rem[i, j] = current.widget().currentIndex()
        self.data["intra"] = intra.tolist()
        self.data["ann"] = ann.tolist()
        self.data["ann_remainder"] = ann_rem.tolist()
        self.parent.protein_data = self.data
        self.updated_keys = ["intra", "ann", "ann_remainder"]

    def checkData(self):
        '''
        special extra check here that the remainders and rates match
        (there are no nonzero rates set with zero remainder or vice-versa)
        '''
        msgs = []
        ann = self.data["ann"]
        ann_rem = self.data["ann_remainder"]
        for i in range(self.data["n_s"]):
            s1 = self.data["state_names"][i]
            for j in range(self.data["n_s"]):
                s2 = self.data["state_names"][j]
                if ann[i][j] == 0.0 and ann_rem[i][j] > 0:
                    msgs.append(f"States {s1} and {s2} have an annihilation "
                    "remainder set but a zero annihilation rate. "
                    "Either set the remainder to None in the bottom right, "
                    "or provide a non-zero annihilation rate between them.")
                    return False, msgs
        return parse.parse_protein(self.parent.protein_data, 
                keys=self.updated_keys)

    def validatePage(self):
        self.updateData()
        valid, msgs = self.checkData()
        if not valid:
            self.errors = QMessageBox.critical(self,
            self.title(), ('\n').join(msgs))
        return valid

class saveProteinPage(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Save")
        self.setSubTitle("Save protein data to file.")

    def initializePage(self):
        self.data = self.parent.protein_data
        layout = QVBoxLayout()
        gl = QGridLayout()
        pf = self.parent.protein_file
        default = os.path.join(os.getcwd(), f"{self.parent.protein_name}.json")
        fstr = default if pf == "" else pf
        self.filename = QLineEdit(fstr)
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
        name = self.parent.protein_name
        final_data = {name: self.data}
        success = True
        self.existing_data = {}
        if os.path.isfile(self.filename.text()):
            with open(self.filename.text(), "r+", encoding='utf-8') as f:
                try:
                    self.existing_data = json.load(f)
                    self.load_success = True
                except:
                    box = QMessageBox.critical(self,
                    "JSON save failed", 
                    "Failed to load existing protein data from JSON to merge.")
                    self.load_success = False
        # if the protein name matches one that's already there and we just
        # merge the dicts, the original will be overwritten, so check
        if self.parent.protein_name in self.existing_data.keys():
            overwrite = True
            self.overwriteCheck = QMessageBox.question(self,
                "", "Protein name already exists in data file. Overwrite?")

            if self.overwriteCheck == QMessageBox.StandardButton.NoButton:
                overwrite = False
            if overwrite:
                total_data = self.existing_data | final_data
                with open(self.filename.text(), "w") as f:
                    json.dump(total_data, f)
            else:
                success = False
        else:
            final_data = self.existing_data | final_data
            with open(self.filename.text(), "w") as f:
                json.dump(final_data, f)
        if success:
            # these will be needed on runPage later on
            self.parent.protein_file = self.filename.text()
        return success

    def onBrowseButton(self):
        self.fn, _ = QFileDialog.getSaveFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.filename.setText(self.fn)

    def onSaveButton(self):
        self.save_success = self.save_to_file()

    def checkData(self):
        return parse.parse_protein(self.data)

    def validatePage(self):
        if not self.save_success:
            self.save_success = self.save_to_file()
            if not self.save_success:
                box = QMessageBox.critical(self,
                "Protein not saved", 
                "Protein data has not been saved successfully.")
                return False
        valid, msgs = self.checkData()
        if not valid:
            self.errors = QMessageBox.critical(self,
            self.title(), ('\n').join(msgs))
        return valid

class loadSimulation(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Load existing simulation parameters.")
        self.data = {}

    def initializePage(self):
        layout = QVBoxLayout()
        gl = QGridLayout()
        fstr = os.path.join(os.getcwd(), "simulation.json")
        self.filename = QLineEdit(fstr)
        self.load_success = False
        self.browseButton = QPushButton("Browse") 
        self.browseButton.setToolTip("Search for a JSON file")
        self.browseButton.clicked.connect(self.onBrowseButton)
        self.resetButton = QPushButton("Reset")
        self.resetButton.setToolTip("Delete loaded JSON data and start again")
        self.resetButton.clicked.connect(self.onResetButton)
        gl.addWidget(QLabel("Filename:"), 0, 0)
        gl.addWidget(self.filename, 0, 1)
        gl.addWidget(self.browseButton, 0, 2)
        gl.addWidget(self.resetButton, 0, 3)
        layout.addLayout(gl)
        self.setLayout(layout)

    def onBrowseButton(self):
        self.fn, _ = QFileDialog.getOpenFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.filename.setText(self.fn)
        self.load_success, self.data = load_from_file(self,
                                                      self.filename.text())

    def onResetButton(self):
        self.load_success = False
        self.filename.setText("")
        self.data = {}

    def updateData(self):
        if self.load_success:
            self.parent.sim_file = self.filename.text()
            self.parent.sim_data = self.data

    def validatePage(self):
        self.updateData()
        return True

class simulationParameters(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Enter simulation parameters.")
        layout = QVBoxLayout()
        gl = QGridLayout()
        self.fwhmLabel     = QLabel("Pulse FWHM(s):")
        self.fluenceLabel  = QLabel("Fluence (photons per pulse):")
        self.nSitesLabel   = QLabel("Number of sites in lattice:")
        self.latticeLabel  = QLabel("Lattice type:")
        self.repRateLabel  = QLabel("Rep rate (Hz):")
        self.burnRepsLabel = QLabel("Number of burn reps:")
        self.tMaxLabel     = QLabel("Maximum binning time (s):")
        self.dt1Label      = QLabel("Binning time step (s):")
        self.dt2Label      = QLabel("Dark time step (s):")
        self.binwidthLabel = QLabel("Binwidth (s):")
        self.nCountsLabel  = QLabel("Number of counts:")
        self.nRepeatsLabel = QLabel("Number of repeats:")
        self.debugLabel    = QLabel("Debug mode:")

        self.fwhmLabel.setToolTip("The FWHM of the pulse, in seconds.")
        self.fluenceLabel.setToolTip("The fluence in photons per pulse.")
        self.nSitesLabel.setToolTip(
'''Number of sites in the lattice. Generally unless you are working with
aggregates of a known, specific size, it's best to leave this on the order
of 100 (especially if hopping is allowed, and if some states are not 
present on every site, for statistical reasons.)''')
        self.latticeLabel.setToolTip(
'''Sets the connectivity of the lattice. Unless you have good reason to
think your aggregate is a line, you can probably ignore this; changing 
from honeycomb to square to hex generally does not make a qualitative
difference to the results.''')
        self.repRateLabel.setToolTip("Laser repetition rate in Hz.")
        self.burnRepsLabel.setToolTip(
'''"Burn reps" are laser reps performed at the start of the simulation
*without binning decays*, in order to allow steady state populations to 
develop and simulate the experimental situation, where the system's probably
already in the steady state before measurement begins. The actual number 
you need is dependent on the system to an extent; see the output files
labelled "_population.csv" for the total population of each state at the
end of each rep to get a better idea.''')
        self.tMaxLabel.setToolTip("Time over which decays will be binned.")
        self.dt1Label.setToolTip(
'''Time step to use during the measurement time, i.e. from t = 0 
to whatever value you give t_max directly above.''')
        self.dt2Label.setToolTip(
'''Time step to use in between measurement phases, i.e. from t = t_max 
up to t = 1 / rep_rate, when the next pulse begins. This time step 
should generally be longer since we're not binning so knowledge
of the precise timing of events is unimportant, and also the shorter this
time step is the longer the simulation will take.''')
        self.binwidthLabel.setToolTip("Binwidth in seconds.")
        self.nCountsLabel.setToolTip(
'''Number of counts to simulate to. Note that this is not a strict limit;
the code will periodically check what the total number of emissive decays
(those visible to the detector) is per bin, and will stop once any
bin reaches n_counts or greater.''')
        self.nRepeatsLabel.setToolTip(
'''Number of repeats to do. The way this works is that each repeat is 
set up with a new RNG seed; additionally, if there are some 
abundances < 1.0, the code will randomise the locations of
proteins with and without the given states for each repeat.''')
        self.debugLabel.setToolTip(
'''Set debug mode on or off. Generally you can probably leave it off;
it's mostly for my benefit, or for if you want to modify the guts of 
the fortran code in some way.''')

        gl.addWidget(self.fwhmLabel,     0, 0)
        gl.addWidget(self.fluenceLabel,  1, 0)
        gl.addWidget(self.nSitesLabel,   2, 0)
        gl.addWidget(self.latticeLabel,  3, 0)
        gl.addWidget(self.repRateLabel,  4, 0)
        gl.addWidget(self.burnRepsLabel, 5, 0)
        gl.addWidget(self.tMaxLabel,     6, 0)
        gl.addWidget(self.dt1Label,      7, 0)
        gl.addWidget(self.dt2Label,      8, 0)
        gl.addWidget(self.binwidthLabel, 9, 0)
        gl.addWidget(self.nCountsLabel,  10, 0)
        gl.addWidget(self.nRepeatsLabel, 11, 0)
        gl.addWidget(self.debugLabel,    12, 0)

        self.fwhmBox     = QLineEdit()
        self.fluenceBox  = QLineEdit()
        self.nSitesBox   = QLineEdit()
        self.latticeBox  = QComboBox()
        self.repRateBox  = QLineEdit()
        self.burnRepsBox = QLineEdit()
        self.tMaxBox     = QLineEdit()
        self.dt1Box      = QLineEdit()
        self.dt2Box      = QLineEdit()
        self.binwidthBox = QLineEdit()
        self.nCountsBox  = QLineEdit()
        self.nRepeatsBox = QLineEdit()
        self.debugBox    = QCheckBox()
        self.latticeBox.addItem("hex")
        self.latticeBox.addItem("square")
        self.latticeBox.addItem("honeycomb")
        self.latticeBox.addItem("line")
        self.boxes = [self.fwhmBox, self.fluenceBox, self.nSitesBox,
                 self.latticeBox, self.repRateBox, self.burnRepsBox,
                 self.tMaxBox, self.dt1Box, self.dt2Box,
                 self.binwidthBox, self.nCountsBox,
                 self.nRepeatsBox, self.debugBox]
        self.sim_keys = ["fwhm",  "fluence",  "n_sites",  "lattice", 
                    "rep_rate",  "burn_reps",  "tmax",  "dt1",  "dt2",
                    "binwidth",  "n_counts",  "n_repeats", "debug"]

        gl.addWidget(self.fwhmBox,     0, 1)
        gl.addWidget(self.fluenceBox,  1, 1)
        gl.addWidget(self.nSitesBox,   2, 1)
        gl.addWidget(self.latticeBox,  3, 1)
        gl.addWidget(self.repRateBox,  4, 1)
        gl.addWidget(self.burnRepsBox, 5, 1)
        gl.addWidget(self.tMaxBox,     6, 1)
        gl.addWidget(self.dt1Box,      7, 1)
        gl.addWidget(self.dt2Box,      8, 1)
        gl.addWidget(self.binwidthBox, 9, 1)
        gl.addWidget(self.nCountsBox,  10, 1)
        gl.addWidget(self.nRepeatsBox, 11, 1)
        gl.addWidget(self.debugBox,    12, 1)

        layout.addLayout(gl)
        self.setLayout(layout)

    def initializePage(self):
        if len(self.parent.sim_data) > 0:
            parent_keys = self.parent.sim_data.keys()
            for sk, b in zip(self.sim_keys, self.boxes):
                if sk in parent_keys:
                    if sk == "lattice":
                        for i in range(b.count()):
                            if b.itemText(i) == self.parent.sim_data[sk]:
                                b.setCurrentIndex(i)
                    elif sk == 'debug':
                        if type(self.parent.sim_data[sk]) == bool:
                            b.setChecked(self.parent.sim_data[sk])
                    else:
                        b.setText(str(self.parent.sim_data[sk]))

    def cleanupPage(self):
        for sk, b in zip(self.sim_keys, self.boxes):
            if sk == 'debug':
                b.setChecked(False)
            if sk == 'lattice':
                b.setCurrentIndex(0)
            else:
                b.setText("")

    def updateData(self):
        for sk, b in zip(self.sim_keys, self.boxes):
            if sk == 'debug':
                v = b.isChecked()
            elif sk == 'lattice':
                v = b.currentText()
            elif sk in ['n_sites', 'burn_reps', 'n_counts', 'n_repeats']:
                v = int(b.text() if b.text() != '' else 0)
            else:
                v = float(b.text() if b.text() != '' else 0.0)
            self.parent.sim_data[sk] = v

    def checkData(self):
        msgs = []
        if self.parent.sim_data['dt1'] >= self.parent.sim_data['binwidth']:
            msgs.append("Measurement time step dt1 should be smaller "
            "than binwidth.")
            return False, msgs
        if self.parent.sim_data['dt1'] > self.parent.sim_data['dt2']:
            msgs.append("Measurement time step dt1 should be smaller "
            "than dark time step dt2.")
            return False, msgs
        return parse.parse_simulation(self.parent.sim_data)

    def validatePage(self):
        self.updateData()
        valid, msgs = self.checkData()
        if not valid:
            self.errors = QMessageBox.critical(self,
            self.title(), ('\n').join(msgs))
        return valid

class saveSimPage(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.setTitle("Save")
        self.setSubTitle("Save simulation data to file.")

    def initializePage(self):
        layout = QVBoxLayout()
        gl = QGridLayout()
        sf = self.parent.sim_file
        fstr = os.getcwd() if sf == "" else sf
        self.filename = QLineEdit(fstr)
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
        dd = self.parent.sim_data
        # don't need these in the JSON
        if 'filename' in dd:
            # filename is only a key if a file was loaded at the start
            del dd['filename']
        success = True
        overwrite = True
        # if the filename exists, warn user
        if os.path.isfile(self.filename.text()):
            self.overwriteCheck = QMessageBox.question(self,
                "", "File already exists. Overwrite?")
            if self.overwriteCheck == QMessageBox.StandardButton.NoButton:
                overwrite = False
            if overwrite:
                try:
                    with open(self.filename.text(), "w") as f:
                        json.dump(dd, f)
                except:
                    box = QMessageBox.critical(self,
                    "JSON overwrite failed",
                    "Failed to overwrite simulation JSON")
                    success = False
            else:
                success = False
        else:
            try:
                with open(self.filename.text(), "w") as f:
                    json.dump(dd, f)
            except:
                box = QMessageBox.critical(self,
                "JSON save failed", "Failed to save simulation JSON")
                success = False
        if success:
            self.parent.sim_file = self.filename.text()
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

class ProteinDataBuilder(QWizard):
    def __init__(self):
        super().__init__()
        self.resize(QtCore.QSize(800, 600))
        self.protein_data = {}
        self.protein_name = ""
        self.protein_file = ""
        self.addPage(loadExisting(self))
        self.addPage(nameNumber(self))
        self.addPage(namePigmentsStates(self))
        self.addPage(stateProperties(self))
        self.addPage(matrixTables(self))
        self.addPage(saveProteinPage(self))
        self.setWindowTitle("Protein parameter wizard")

class SimulationDataBuilder(QWizard):
    def __init__(self):
        super().__init__()
        self.resize(QtCore.QSize(800, 600))
        self.sim_data = {}
        self.sim_file = ""
        self.addPage(loadSimulation(self))
        self.addPage(simulationParameters(self))
        self.addPage(saveSimPage(self))

class main_window(main_ui_class, main_ui_widget):
    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.browseProteinButton.clicked.connect(self.browseProtein)
        self.createProteinButton.clicked.connect(self.createProtein)
        self.browseSimulationButton.clicked.connect(self.browseSimulation)
        self.createSimulationButton.clicked.connect(self.createSimulation)
        self.proteinChoice.currentTextChanged.connect(self.updateConnected)
        self.runButton.clicked.connect(self.run)
        self.stopButton.clicked.connect(self.kill)
        self.quitButton.clicked.connect(self.quit)
        self.protein_window = None
        self.simulation_window = None
        self.process = None
        self.all_pd = {}
        self.pd = {}
        self.sd = {}
        cs = f"Recommended max cores (os.cpu_count()): {os.cpu_count()}"
        self.coresGuideLabel.setText(cs)

    def updateConnected(self):
        p = self.proteinChoice.currentText()
        if p in self.all_pd:
            hh = self.all_pd[p]['hop']
            b = any([h > 0.0 for h in hh])
            self.connectedBox.setChecked(b)

    def browseProtein(self):
        self.proteinChoice.clear()
        self.fn, _ = QFileDialog.getOpenFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.proteinFileBox.setText(self.fn)
        self.load_success, self.all_pd = load_from_file(self,
                        self.proteinFileBox.text())
        if self.load_success:
            for k in self.all_pd.keys():
                self.proteinChoice.addItem(k)

    def createProtein(self):
        if self.protein_window is None:
            self.protein_window = ProteinDataBuilder()
            protein_filename = self.protein_window.protein_file
            self.protein_window.show()
            self.proteinFileBox.setText(protein_filename)
        else:
            self.protein_window.close()
            self.protein_window = None

    def browseSimulation(self):
        self.fn, _ = QFileDialog.getOpenFileName(self, "Select JSON file",
                                              os.getcwd(),
                                              "JSON file (*.json)")
        self.simulationFileBox.setText(self.fn)

    def createSimulation(self):
        if self.simulation_window is None:
            self.simulation_window = SimulationDataBuilder()
            self.simulation_window.show()
        else:
            self.simulation_window.close()
            self.simulation_window = None

    def run(self):
        '''
        run the TCSPC simulation with the options below.
        '''
        if self.process is not None:
            return
        pls, self.all_pd = load_from_file(self, self.proteinFileBox.text())
        if not pls:
            return
        sls, self.sd = load_from_file(self, self.simulationFileBox.text())
        if not sls:
            return

        p = f"{self.proteinChoice.currentText()}"
        self.pd = self.all_pd[p]
        valid, msgs = parse.parse_protein(self.pd)
        if not valid:
            msgs.append("Try clicking Create and fixing JSON data.")
            self.simulationErrors = QMessageBox.critical(self,
            "Protein load errors", ('\n').join(msgs))
            return

        valid, msgs = parse.parse_simulation(self.sd)
        if not valid:
            msgs.append("Try clicking Create and fixing JSON data.")
            self.simulationErrors = QMessageBox.critical(self,
            "Simulation errors", ('\n').join(msgs))
            return

        self.process = QtCore.QProcess(self)
        self.outputArea.clear()
        self.runButton.setEnabled(False)
        self.process.readyReadStandardOutput.connect(self.read_stdout)
        self.process.readyReadStandardError.connect(self.read_stderr)
        self.process.finished.connect(self.run_finished)

        pf = f"{self.proteinFileBox.text()}"
        sf = f"{self.simulationFileBox.text()}"
        n = f"{self.numCores.value()}"

        if self.connectedBox.isChecked():
            c = "--connection"
            if all([h == 0.0 for h in self.pd['hop']]):
                box = QMessageBox.warning(self,
                "Simulation errors","All hopping rates are set to zero "
                "but connection box is checked. Proceed?", 
                QMessageBox.StandardButton.Yes | 
                QMessageBox.StandardButton.No)
                if box == QMessageBox.StandardButton.No:
                    valid = False
        else:
            c = "--no-connection"
            if any([h > 0.0 for h in self.pd['hop']]):
                box = QMessageBox.warning(self,
                "Simulation errors","There are non-zero hopping rates, "
                "but connection box is set to unconnected (detergent)"
                ". Proceed?", QMessageBox.StandardButton.Yes |
                QMessageBox.StandardButton.No)
                if box == QMessageBox.StandardButton.No:
                    valid = False

        args = ["main.py", "-pf", pf, "-sf", sf, "-p", p, c, "-n", n]
        if self.outputDirBox.text() != "":
            args.append("-o")
            args.append(self.outputDirBox.text())

        if valid:
            self.process.start("python", args)
        else:
            self.process.kill()
            self.process = None
            self.runButton.setEnabled(False)

    def read_stdout(self):
        data = self.process.readAllStandardOutput()
        text = bytes(data).decode("utf-8")
        self.outputArea.appendPlainText(text.rstrip())

    def read_stderr(self):
        data = self.process.readAllStandardError()
        text = bytes(data).decode("utf-8")
        self.outputArea.appendPlainText(text.rstrip())

    def run_finished(self):
        output = self.outputArea.toPlainText()
        args = [self.pd,
                self.sd, self.proteinChoice.currentText(),
                self.connectedBox.isChecked()]
        if self.outputDirBox.text() != "":
            args.append(self.outputDirBox.text())
        outdir = parse.generate_output_dirs(*args)
        outfile = os.path.join(outdir, "stdout.log")
        try:
            with open(outfile, "w") as f:
                f.write(output)
        except:
            box = QMessageBox(self, "Logfile save error",
                "Unable to save log file.")
        self.process = None
        self.runButton.setEnabled(True)

    def kill(self):
        '''
        stop the TCSPC simulation
        '''
        if self.process is None:
            return
        self.process.kill()
        self.outputArea.appendPlainText("Simulation stopped.")
        self.process = None

    def quit(self):
        if self.process is not None:
            box = QMessageBox.critical(self,
            "Really quit?", "Simulation is running. Really quit?")
            box.setStandardButtons(QtMessageBox.Yes | QtMessageBox.No)
            box.setDefaultButton(QtMessageBox.StandardButton.No)
            button = box.exec()
            if button == QMessageBox.Yes:
                self.kill()
                self.close()
        else:
            self.kill()
            self.close()

def start():
    app =QtWidgets.QApplication([])
    widget = main_window(None)
    widget.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    start()
