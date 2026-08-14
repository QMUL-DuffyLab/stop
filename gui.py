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

class MatrixTable(QTableWidget):
    def __init__(self, data, block_diagonal=False):
        QTableWidget.__init__(self)
        self.data = data
        self.n_states = self.data["n_states"]
        self.setRowCount(self.n_states)
        self.setColumnCount(self.n_states)
        self.setHorizontalHeaderLabels(self.data["state_names"])
        self.setVerticalHeaderLabels(self.data["state_names"])
        if block_diagonal:
            for i in range(self.n_states):
                item = QtWidgets.QTableWidgetItem()
                item.setBackground(QtGui.QColor("darkGray"))
                item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
                self.setItem(i, i, item)

class RemainderTable(QTableWidget):
    def __init__(self, data):
        QTableWidget.__init__(self)
        self.data = data
        self.n_states = self.data["n_states"]
        self.setRowCount(self.n_states)
        self.setColumnCount(self.n_states)
        self.setHorizontalHeaderLabels(self.data["state_names"])
        self.setVerticalHeaderLabels(self.data["state_names"])
        for i in range(self.n_states):
            for j in range(self.n_states):
                item = QtWidgets.QComboBox()
                for k in range(self.n_states):
                    item.addItem(self.data["state_names"][i])
                self.setItem(i, j, item)

class MatrixTables(QtWidgets.QWidget):
    def __init__(self, data, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        layout = QVBoxLayout(self)
        self.data = data
        print(f"MT init: self.data = {self.data}")
        self.intra = MatrixTable(self.data, block_diagonal=True)
        self.ann   = MatrixTable(self.data)
        # self.ann_rem = RemainderTable(self.data)
        layout.addWidget(self.intra)
        layout.addWidget(self.ann)
        # layout.addWidget(self.ann_rem)

        btnNext = QtWidgets.QPushButton("Next")
        btnNext.clicked.connect(self.showNext)
        btnLayout = QHBoxLayout()
        btnLayout.addWidget(btnNext)

        layout.addLayout(btnLayout)
        self.setLayout(layout)

    def updateData(self):
        intra = np.zeros((self.data["n_states"], self.data["n_states"]), dtype=float)
        ann = np.zeros_like(intra)
        ann_rem = np.zeros((self.data["n_states"], self.data["n_states"]), dtype=int)
        for i in range(self.data["n_states"]):
            for j in range(self.data["n_states"]):
                intra[i, j] = float(self.intra.item(i, j).text())
                ann[i, j] = float(self.ann.item(i, j).text())
                # ann_rem[i, j] = self.ann_rem.item(i, j).currentIndex()
        self.data["intra"] = intra
        self.data["ann"] = ann
        self.data["ann_remainder"] = ann_rem
        print(f"MT exit: self.data = {self.data}")

    def showNext(self):
        self.updateData()

class NameNumber(QtWidgets.QWidget):
    def __init__(self, data, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        layout = QGridLayout(self)
        self.protein_name = QLineEdit()
        self.n_p = QSpinBox()
        self.n_p.setRange(1, 20)
        self.n_p.setValue(1)
        self.n_s = QSpinBox()
        self.n_s.setRange(self.n_p.value(), 20)
        self.n_s.setValue(2)
        layout.addWidget(QLabel("Protein name:"), 0, 0)
        layout.addWidget(self.protein_name, 0, 1)
        layout.addWidget(QLabel("Number of pigments:"), 1, 0)
        layout.addWidget(self.n_p, 1, 1)
        layout.addWidget(QLabel("Number of states:"), 2, 0)
        layout.addWidget(self.n_s, 2, 1)

    def updateData(self, data):
        data["name"] = self.protein_name.text()
        data["n_pigments"] = self.n_p.value()
        data["n_states"] = self.n_s.value()

class NamePigmentsStates(QtWidgets.QWidget):
    def __init__(self, data, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        layout = QVBoxLayout()
        pl = QGridLayout()
        sl = QGridLayout()
        self.pigmentProperties = None
        self.pigment_names = []
        self.state_names = []
        if "n_pigments" in data:
            self.n_p = data["n_pigments"]
        else:
            self.n_p = 0
        if "n_states" in data:
            self.n_s = data["n_states"]
        else:
            self.n_s = 0
        print(f"NPS init: data = {data}")
        self.data = data
        print(f"NPS init: self.data = {self.data}")

        for i in range(self.n_p):
            pl.addWidget(QLabel(f"Name of pigment {i + 1:d}:"), i, 0)
            self.pigment_names.append(QLineEdit())
            pl.addWidget(self.pigment_names[i], i, 1)
        for i in range(self.n_s):
            sl.addWidget(QLabel(f"Name of state {i + 1:d}:"), i, 0)
            self.state_names.append(QLineEdit())
            sl.addWidget(self.state_names[i], i, 1)
        layout.addLayout(pl)
        layout.addLayout(sl)

        btnNext = QtWidgets.QPushButton("Next")
        btnNext.clicked.connect(self.showPigmentProperties)
        btnLayout = QHBoxLayout()
        btnLayout.addWidget(btnNext)

        layout.addLayout(btnLayout)
        self.setLayout(layout)

    def updateData(self):
        self.data["pigment_names"] = [p.text() for p in self.pigment_names]
        self.data["state_names"]   = [s.text() for s in self.state_names]

    def showPigmentProperties(self):
        self.updateData()
        if self.pigmentProperties is None:
            self.pigmentProperties = PigmentProperties(self.data)
        self.pigmentProperties.show()
        print(f"NPS exit: self.data = {self.data}")

class PigmentProperties(QtWidgets.QWidget):
    def __init__(self, data, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        layout = QVBoxLayout()
        pl = QGridLayout()
        self.hop = []
        self.decay = []
        self.xsec = []
        self.emissive = []
        self.which_p = []
        print(f"PP init: data = {data}")
        self.data = data
        print(f"PP init: self.data = {self.data}")
        self.matrixTables = None

        # column headers
        pl.addWidget(QLabel("Hopping time (s)"), 0, 1)
        pl.addWidget(QLabel("Decay time (s)"), 0, 2)
        pl.addWidget(QLabel("Cross-section (cm^{-1})"), 0, 3)
        pl.addWidget(QLabel("Emissive decay?"), 0, 4)
        pl.addWidget(QLabel("Pigment"), 0, 5)
        for i in range(self.data["n_states"]):
            row = i + 1
            pl.addWidget(QLabel(f"{self.data['state_names'][i]}:"), row, 0)
            self.hop.append(QLineEdit("0.0"))
            pl.addWidget(self.hop[i], row, 1)
            self.decay.append(QLineEdit("0.0"))
            pl.addWidget(self.decay[i], row, 2)
            self.xsec.append(QLineEdit("0.0"))
            pl.addWidget(self.xsec[i], row, 3)
            self.emissive.append(QCheckBox())
            pl.addWidget(self.emissive[i], row, 4)
            self.which_p.append(QComboBox())
            for name in self.data['pigment_names']:
                self.which_p[i].addItem(name)
            pl.addWidget(self.which_p[i], row, 5)

        layout.addLayout(pl)

        btnNext = QtWidgets.QPushButton("Next")
        btnNext.clicked.connect(self.showMatrixTables)
        btnLayout = QHBoxLayout()
        btnLayout.addWidget(btnNext)

        layout.addLayout(btnLayout)
        self.setLayout(layout)

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

    def showMatrixTables(self):
        self.updateData()
        if self.matrixTables is None:
            self.matrixTables = MatrixTables(self.data)
        self.matrixTables.show()
        print(f"PP exit: self.data = {self.data}")

class ProteinBuilder(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.data = {}
        layout = QVBoxLayout(self)

        self.nameNumber = NameNumber(self.data)
        self.namePigmentStates = None
        self.pigmentProperties = None
        layout.addWidget(self.nameNumber)

        btnNext = QtWidgets.QPushButton("Next")
        btnNext.clicked.connect(self.showNamePigmentStates)
        btnLayout = QHBoxLayout()
        btnLayout.addWidget(btnNext)

        layout.addLayout(btnLayout)
        self.setLayout(layout)

    def showNamePigmentStates(self):
        self.nameNumber.updateData(self.data)
        if self.namePigmentStates is None:
            self.namePigmentStates = NamePigmentsStates(self.data)
        self.namePigmentStates.show()

if __name__ == "__main__":
    app =QtWidgets.QApplication([])

    widget = ProteinBuilder()
    widget.show()

    sys.exit(app.exec())
