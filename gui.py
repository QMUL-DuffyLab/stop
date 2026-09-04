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
        self.setTitle("Protein Builder")
        self.setSubTitle('''
If you'd like to import existing protein data, you can do that here.
Otherwise, click Next to start specifying the parameters of your protein.
                         ''')

    def initializePage(self):
        layout = QVBoxLayout()
        gl = QGridLayout()
        self.filename = QLineEdit()
        self.load_success = False
        self.browseButton = QPushButton("Browse") 
        self.browseButton.clicked.connect(self.onBrowseButton)
        self.loadButton = QPushButton("Load") 
        self.loadButton.clicked.connect(self.onLoadButton)
        self.proteinChooser = QComboBox()
        self.resetButton = QPushButton("Reset")
        self.resetButton.clicked.connect(self.onResetButton)
        gl.addWidget(QLabel("Filename:"), 0, 0)
        gl.addWidget(self.filename, 0, 1)
        gl.addWidget(self.browseButton, 0, 2)
        gl.addWidget(self.loadButton, 0, 3)
        gl.addWidget(QLabel("Protein name:"), 1, 0)
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
            # NB: checking string legality??

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
            print(msgs)
        print(f"NN exit: data = {self.parent.data}")
        return validated

class namePigmentsStates(QWizardPage):
    def __init__(self, parent):
        QWizardPage.__init__(self, parent)
        self.parent = parent
        self.layout = QVBoxLayout()
        self.pl = QGridLayout()
        self.sl = QGridLayout()
        self.layout.addLayout(self.pl)
        self.layout.addLayout(self.sl)
        self.setLayout(self.layout)

    def initializePage(self):
        self.n_p = self.field('n_p')
        self.n_s = self.field('n_s')
        self.pigment_names = []
        self.state_names = []
        for i in range(self.n_p):
            self.pl.addWidget(QLabel(f"Name of pigment {i + 1:d}:"), i, 0)
            current_name = QLineEdit()
            self.pigment_names.append(current_name)
            #self.registerField(f"pigment_name{i + 1:d}", current_name)
            self.pl.addWidget(current_name, i, 1)
        for i in range(self.n_s):
            self.sl.addWidget(QLabel(f"Name of state {i + 1:d}:"), i, 0)
            current_state = QLineEdit()
            self.state_names.append(current_state)
            #self.registerField(f"state_name{i + 1:d}", current_state)
            self.sl.addWidget(current_state, i, 1)
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
    
    def cleanupPage(self):
        self.pigment_names = []
        self.state_names = []
        self.n_p = 0
        self.n_s = 0
        for item in self.pigment_names:
            item.setText("")
        for item in self.state_names:
            item.setText("")
        for layout in self.pl, self.sl:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget:
                    child.widget().deleteLater()


    def updateData(self):
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
            ['pigment_names', 'state_names'], ['n_p', 'n_s']):
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
        self.parent = parent
        self.layout = QVBoxLayout()
        self.pl = QGridLayout()
        self.layout.addLayout(self.pl)
        self.setLayout(self.layout)
        self.n_s = 0
        # column headers
        self.pl.addWidget(QLabel("Hopping time (s)"), 0, 1)
        self.pl.addWidget(QLabel("Decay time (s)"), 0, 2)
        self.pl.addWidget(QLabel("Cross-section (cm^{-1})"), 0, 3)
        self.pl.addWidget(QLabel("Emissive decay?"), 0, 4)
        self.pl.addWidget(QLabel("Pigment"), 0, 5)

    def initializePage(self):
        self.n_s = self.field("n_s")
        self.hop      = []
        self.decay    = []
        self.xsec     = []
        self.emissive = []
        self.which_p  = []
        for i in range(self.n_s):
            row = i + 1
            state_name = self.parent.data["state_names"][i]
            self.hop.append(QLineEdit("0.0"))
            self.decay.append(QLineEdit("0.0"))
            self.xsec.append(QLineEdit("0.0"))
            self.emissive.append(QCheckBox())
            self.which_p.append(QComboBox())
            self.pl.addWidget(QLabel(state_name), row, 0)
            self.pl.addWidget(self.hop[i], row, 1)
            self.pl.addWidget(self.decay[i], row, 2)
            self.pl.addWidget(self.xsec[i], row, 3)
            self.pl.addWidget(self.emissive[i], row, 4)
            # need to find out how to centre the checkboxes. it's annoying
            # self.pl.setAlignment(self.emissive[i], Qt.AlignHCenter)
            for j in range(self.field("n_p")):
                name = self.parent.data["pigment_names"][j]
                self.which_p[i].addItem(name)
            self.pl.addWidget(self.which_p[i], row, 5)

    def cleanupPage(self):
        while self.pl.count():
            child = self.pl.takeAt(0)
            if child.widget:
                child.widget().deleteLater()
        self.hop      = []
        self.decay    = []
        self.xsec     = []
        self.emissive = []
        self.which_p  = []
        self.n_s = 0

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
        self.addPage(loadExisting(self))
        self.addPage(nameNumber(self))
        self.addPage(namePigmentsStates(self))
        self.addPage(pigmentProperties(self))
        self.addPage(matrixTables(self))


if __name__ == "__main__":
    app =QtWidgets.QApplication([])

    widget = ProteinBuilder()
    widget.show()

    sys.exit(app.exec())
