# callum - 30/07/2026
# -*- coding: utf-8 -*-
# JsonModel taken from
# doc.qt.io/qtforpython-6/examples/example_widgets_itemviews_jsonmodel.html

import sys
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

class nameNumber(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        layout = QVBoxLayout()
        gl = QGridLayout()
        protein_name = QLineEdit()
        self.registerField('protein_name', protein_name, "text")
        n_p = QSpinBox()
        n_p.setRange(0, 20)
        self.registerField('n_p', n_p)
        n_s = QSpinBox()
        n_s.setRange(n_p.value(), 20)
        self.registerField('n_s', n_s)
        gl.addWidget(QLabel("Protein name:"), 0, 0)
        gl.addWidget(protein_name, 0, 1)
        gl.addWidget(QLabel("Number of pigments:"), 1, 0)
        gl.addWidget(n_p, 1, 1)
        gl.addWidget(QLabel("Number of states:"), 2, 0)
        gl.addWidget(n_s, 2, 1)
        layout.addLayout(gl)
        self.setLayout(layout)

    def updateData(self):
        '''
        update the parent QWizard's data struct with the
        data that's been entered here. also keep track of
        the names of the fields so that they can be cleaned
        up by cleanupPage()
        '''
        self.parent.data["name"]       = self.field('protein_name')
        self.parent.data["n_pigments"] = self.field('n_p')
        self.parent.data["n_states"]   = self.field('n_s')
        self.fields = ["name", "n_pigments", "n_states"]

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []

        if len(self.parent.data["name"]) == 0:
            msgs.append("Length of protein name must be > 0.")
            validated = False
            # NB: checking string legality??

        if self.parent.data["n_pigments"] <= 0:
            msgs.append("Number of pigments must be > 0.")
            validated = False

        if self.parent.data["n_states"] <= 0:
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
            print(msgs)
        print(f"NN exit: data = {self.parent.data}")
        return validated

class namePigmentsStates(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent

    def initializePage(self):
        n_p = self.field('n_p')
        n_s = self.field('n_s')
        layout = QVBoxLayout()
        pl = QGridLayout()
        sl = QGridLayout()

        self.pigment_names = []
        self.state_names = []
        for i in range(n_p):
            pl.addWidget(QLabel(f"Name of pigment {i + 1:d}:"), i, 0)
            current_name = QLineEdit()
            self.pigment_names.append(current_name)
            self.registerField(f"pigment_name{i + 1:d}", current_name)
            pl.addWidget(current_name, i, 1)
        for i in range(n_s):
            sl.addWidget(QLabel(f"Name of state {i + 1:d}:"), i, 0)
            current_state = QLineEdit()
            self.state_names.append(current_state)
            sl.addWidget(current_state, i, 1)
            self.registerField(f"state_name{i + 1:d}", current_state)
        layout.addLayout(pl)
        layout.addLayout(sl)
        self.setLayout(layout)

    def updateData(self):
        print(self.pigment_names)
        self.parent.data["pigment_names"] = [p.text()
                                        for p in self.pigment_names]
        self.parent.data["state_names"]   = [s.text()
                                        for s in self.state_names]

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []
        data = self.parent.data

        for name_string, n_string in zip(
            ['pigment_names', 'state_names'], ['n_pigments', 'n_states']):
            nsmsg = name_string.replace("_", " ")

            if any([len(p) == 0 for p in data[name_string]]):
                msgs.append(f"Length of {nsmsg} must be > 0.")
                validated = False
                # NB: checking string legality??

            if len(data[name_string]) != data[n_string]:
                msgs.append(f"Number of {nsmsg} doesn't match {n_string}.")
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
            print(msgs)
        print(f"NPS exit: data = {self.parent.data}")
        return validated

class pigmentProperties(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)

    def initializePage(self):
        n_s = self.field("n_s")
        layout = QVBoxLayout()
        pl = QGridLayout()
        self.hop      = []
        self.decay    = []
        self.xsec     = []
        self.emissive = []
        self.which_p  = []

        # column headers
        pl.addWidget(QLabel("Hopping time (s)"), 0, 1)
        pl.addWidget(QLabel("Decay time (s)"), 0, 2)
        pl.addWidget(QLabel("Cross-section (cm^{-1})"), 0, 3)
        pl.addWidget(QLabel("Emissive decay?"), 0, 4)
        pl.addWidget(QLabel("Pigment"), 0, 5)
        for i in range(n_s):
            row = i + 1
            state_name = self.field(f"state_name{i + 1:d}")
            pl.addWidget(QLabel(state_name), row, 0)
            self.hop.append(QLineEdit("0.0"))
            pl.addWidget(self.hop[i], row, 1)
            self.decay.append(QLineEdit("0.0"))
            pl.addWidget(self.decay[i], row, 2)
            self.xsec.append(QLineEdit("0.0"))
            pl.addWidget(self.xsec[i], row, 3)
            self.emissive.append(QCheckBox())
            pl.addWidget(self.emissive[i], row, 4)
            self.which_p.append(QComboBox())
            for j in range(self.field("n_p")):
                name = self.field(f"pigment_name{j + 1:d}")
                self.which_p[i].addItem(name)
            pl.addWidget(self.which_p[i], row, 5)
        layout.addLayout(pl)
        self.setLayout(layout)

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
        self.fields = ["hop", "decay", "xsec", "emissive", "which_pigment"]

    def checkData(self):
        '''
        check the entered data to make sure it's all tickety boo
        '''
        validated = True
        msgs = []
        data = self.parent.data

        return validated, msgs
        
    def validatePage(self):
        '''
        update the parent's data, check it, print the current
        dict for my benefit, then carry on if all is well
        '''
        self.updateData()
        validated, msgs = self.checkData()
        if not validated:
            print(msgs)
        print(f"PP exit: data = {self.parent.data}")
        return validated

    def cleanupPage(self):
        '''
        if the user goes back, delete the fields we've
        added here so that there's no problem with duplicates
        '''
        if getattr(self, "fields", False):
            for field in self.fields:
                del self.parent.data[field]

class MatrixTable(QTableWidget):
    def __init__(self, state_names, block_diagonal=False):
        QTableWidget.__init__(self)
        self.state_names = state_names
        self.n_states = len(state_names)
        self.setRowCount(self.n_states)
        self.setColumnCount(self.n_states)
        self.setHorizontalHeaderLabels(self.state_names)
        self.setVerticalHeaderLabels(self.state_names)
        if block_diagonal:
            for i in range(self.n_states):
                item = QtWidgets.QTableWidgetItem()
                item.setBackground(QtGui.QColor("darkGray"))
                item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
                self.setItem(i, i, item)

class RemainderTable(QTableWidget):
    def __init__(self, state_names):
        QTableWidget.__init__(self)
        self.state_names = state_names
        self.n_states = len(state_names)
        self.setRowCount(self.n_states)
        self.setColumnCount(self.n_states)
        self.setHorizontalHeaderLabels(self.state_names)
        self.setVerticalHeaderLabels(self.state_names)
        for i in range(self.n_states):
            for j in range(self.n_states):
                item = QtWidgets.QComboBox()
                for k in range(self.n_states):
                    item.addItem(self.state_names[i])
                self.setItem(i, j, item)

class matrixTables(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)

    def initializePage(self):
        layout = QVBoxLayout()
        state_names = [self.field(f"state_name{i + 1:d}")
                       for i in range(self.field("n_s"))]
        intra = MatrixTable(state_names, block_diagonal=True)
        ann   = MatrixTable(state_names)
        ann_rem = RemainderTable(state_names)
        layout.addWidget(intra)
        layout.addWidget(ann)
        layout.addWidget(ann_rem)
        self.setLayout(layout)

    def validatePage(self):
        intra = np.zeros((self.data["n_states"], self.data["n_states"]), dtype=float)
        ann = np.zeros_like(intra)
        ann_rem = np.zeros((self.data["n_states"], self.data["n_states"]), dtype=int)
        for i in range(self.data["n_states"]):
            for j in range(self.data["n_states"]):
                intra[i, j] = float(self.intra.item(i, j).text())
                ann[i, j] = float(self.ann.item(i, j).text())
                # ann_rem[i, j] = self.ann_rem.item(i, j).currentIndex()
        self.parent.data["intra"] = intra
        self.parent.data["ann"] = ann
        self.parent.data["ann_remainder"] = ann_rem
        print(f"MT exit: data = {self.parent.data}")

class ProteinBuilder(QWizard):
    def __init__(self):
        super().__init__()
        self.data = {}
        self.addPage(nameNumber(self))
        self.addPage(namePigmentsStates(self))
        self.addPage(pigmentProperties(self))
        self.addPage(matrixTables(self))


if __name__ == "__main__":
    app =QtWidgets.QApplication([])

    widget = ProteinBuilder()
    widget.show()

    sys.exit(app.exec())
