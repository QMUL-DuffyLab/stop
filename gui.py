# callum - 30/07/2026
# -*- coding: utf-8 -*-
import sys
import random
from PyQt6 import QtCore, QtWidgets, QtGui
from PyQt6.QtCore import Qt
import json
from PyQt6.QtWidgets import (QLineEdit, QPushButton, QApplication,
    QVBoxLayout, QHBoxLayout, QSpinBox, QTableWidget, QDialog)


class Widge(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.hello = ["allo m8", "fook off m8", "wassup lad", "piss off idiot"]
        
        self.button = QtWidgets.QPushButton("Click me!")
        self.text = QtWidgets.QLabel("Hello World")

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button)

        self.button.clicked.connect(self.magic)

    def magic(self):
        self.text.setText(random.choice(self.hello))

class MatrixTable(QTableWidget):
    def __init__(self, table_name, state_names, block_diagonal=False):
        QTableWidget.__init__(self)
        if type(table_name) is not str:
            raise TypeError("Table name passed to MatrixTable must be str.")
        if type(state_names) is not list:
            raise TypeError("State names passed to MatrixTable must be list.")
        self.table_name = table_name
        self.state_names = state_names
        self.n_states = len(self.state_names)
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

class NameNumber(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        layout = QVBoxLayout(self)
        self.protein_name = QLineEdit("Enter protein name")
        self.n_p = QSpinBox()
        self.n_p.setRange(1, 20)
        self.n_p.setValue(1)
        self.n_s = QSpinBox()
        self.n_s.setRange(self.n_p.value(), 20)
        self.n_s.setValue(2)
        layout.addWidget(self.protein_name)
        layout.addWidget(self.n_p)
        layout.addWidget(self.n_s)

class ProteinBuilder(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        layout = QVBoxLayout()

        nameNumber = NameNumber()
        self.protein_name = nameNumber.protein_name.value()
        self.n_p = nameNumber.n_p.value()
        self.n_s = nameNumber.n_s.value()

        self.Stack = QtWidgets.QStackedWidget()
        self.Stack.addWidget(NameNumber())
        self.Stack.addWidget(MatrixTable("LHCII",
                            ["Chl_S", "Chl_T", "Car_S", "Car_T"], True))

        btnNext = QtWidgets.QPushButton("Next")
        btnNext.clicked.connect(self.onNext)
        btnPrevious = QtWidgets.QPushButton("Previous")
        btnPrevious.clicked.connect(self.onPrevious)
        btnLayout = QHBoxLayout()
        btnLayout.addWidget(btnPrevious)
        btnLayout.addWidget(btnNext)

        layout.addWidget(self.Stack)
        layout.addLayout(btnLayout)
        self.setLayout(layout)

    def onNext(self):
        self.Stack.setCurrentIndex(self.Stack.currentIndex() + 1)

    def onPrevious(self):
        self.Stack.setCurrentIndex(self.Stack.currentIndex() - 1)

if __name__ == "__main__":
    app =QtWidgets.QApplication([])

    widget = ProteinBuilder()
    widget.show()

    sys.exit(app.exec())
