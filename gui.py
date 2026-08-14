# callum - 30/07/2026
# -*- coding: utf-8 -*-
# JsonModel taken from
# doc.qt.io/qtforpython-6/examples/example_widgets_itemviews_jsonmodel.html

import sys
import random
from PyQt6 import QtCore, QtWidgets, QtGui
from PyQt6.QtCore import Qt
import json
from PyQt6.QtWidgets import *
from typing import Any

class TreeItem:
    """A Json item corresponding to a line in QTreeView"""

    def __init__(self, parent: "TreeItem" = None):
        self._parent = parent
        self._key = ""
        self._value = ""
        self._value_type = None
        self._children = []

    def appendChild(self, item: "TreeItem"):
        """Add item as a child"""
        self._children.append(item)

    def child(self, row: int) -> "TreeItem":
        """Return the child of the current item from the given row"""
        return self._children[row]

    def parent(self) -> "TreeItem":
        """Return the parent of the current item"""
        return self._parent

    def childCount(self) -> int:
        """Return the number of children of the current item"""
        return len(self._children)

    def row(self) -> int:
        """Return the row where the current item occupies in the parent"""
        return self._parent._children.index(self) if self._parent else 0

    @property
    def key(self) -> str:
        """Return the key name"""
        return self._key

    @key.setter
    def key(self, key: str):
        """Set key name of the current item"""
        self._key = key

    @property
    def value(self) -> str:
        """Return the value name of the current item"""
        return self._value

    @value.setter
    def value(self, value: str):
        """Set value name of the current item"""
        self._value = value

    @property
    def value_type(self):
        """Return the python type of the item's value."""
        return self._value_type

    @value_type.setter
    def value_type(self, value):
        """Set the python type of the item's value."""
        self._value_type = value

    @classmethod
    def load(
        cls, value: list | dict, parent: "TreeItem" = None, sort=True
    ) -> "TreeItem":
        """Create a 'root' TreeItem from a nested list or a nested dictonary

        Examples:
            with open("file.json") as file:
                data = json.dump(file)
                root = TreeItem.load(data)

        This method is a recursive function that calls itself.

        Returns:
            TreeItem: TreeItem
        """
        rootItem = TreeItem(parent)
        rootItem.key = "root"

        if isinstance(value, dict):
            items = sorted(value.items()) if sort else value.items()

            for key, value in items:
                child = cls.load(value, rootItem)
                child.key = key
                child.value_type = type(value)
                rootItem.appendChild(child)

        elif isinstance(value, list):
            for index, value in enumerate(value):
                child = cls.load(value, rootItem)
                child.key = index
                child.value_type = type(value)
                rootItem.appendChild(child)

        else:
            rootItem.value = value
            rootItem.value_type = type(value)

        return rootItem

class JsonModel(QAbstractItemModel):
    """ An editable model of Json data """

    def __init__(self, parent: QObject = None):
        super().__init__(parent)

        self._rootItem = TreeItem()
        self._headers = ("key", "value")

    def clear(self):
        """ Clear data from the model """
        self.load({})

    def load(self, document: dict):
        """Load model from a nested dictionary returned by json.loads()

        Arguments:
            document (dict): JSON-compatible dictionary
        """

        assert isinstance(
            document, (dict, list, tuple)
        ), "`document` must be of dict, list or tuple, " f"not {type(document)}"

        self.beginResetModel()

        self._rootItem = TreeItem.load(document)
        self._rootItem.value_type = type(document)

        self.endResetModel()

        return True

    def data(self, index: QModelIndex, role: Qt.ItemDataRole) -> Any:
        """Override from QAbstractItemModel

        Return data from a json item according index and role

        """
        ret = None
        if index.isValid():
            item = index.internalPointer()
            match role:
                case Qt.ItemDataRole.DisplayRole:
                    match index.column():
                        case 0:
                            ret = item.key
                        case 1:
                            ret = item.value
                case Qt.ItemDataRole.EditRole:
                    if index.column() == 1:
                        ret = item.value
        return ret

    def setData(self, index: QModelIndex, value: Any, role: Qt.ItemDataRole):
        """Override from QAbstractItemModel

        Set json item according index and role

        Args:
            index (QModelIndex)
            value (Any)
            role (Qt.ItemDataRole)

        """
        if role == Qt.ItemDataRole.EditRole:
            if index.column() == 1:
                item = index.internalPointer()
                item.value = str(value)

                self.dataChanged.emit(index, index, [Qt.ItemDataRole.EditRole])

                return True

        return False

    def headerData(
        self, section: int, orientation: Qt.Orientation, role: Qt.ItemDataRole
    ):
        """Override from QAbstractItemModel

        For the JsonModel, it returns only data for columns (orientation = Horizontal)

        """
        if role != Qt.ItemDataRole.DisplayRole:
            return None

        if orientation == Qt.Orientation.Horizontal:
            return self._headers[section]

    def index(self, row: int, column: int, parent=QModelIndex()) -> QModelIndex:
        """Override from QAbstractItemModel

        Return index according row, column and parent

        """
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self._rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()

    def parent(self, index: QModelIndex) -> QModelIndex:
        """Override from QAbstractItemModel

        Return parent index of index

        """

        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parent()

        if parentItem == self._rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent=QModelIndex()):
        """Override from QAbstractItemModel

        Return row count from parent index
        """
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            parentItem = self._rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount()

    def columnCount(self, parent=QModelIndex()):
        """Override from QAbstractItemModel

        Return column number. For the model, it always return 2 columns
        """
        return 2

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:
        """Override from QAbstractItemModel

        Return flags of index
        """
        flags = super(JsonModel, self).flags(index)

        if index.column() == 1:
            return Qt.ItemFlag.ItemIsEditable | flags
        else:
            return flags

    def to_json(self, item=None):

        if item is None:
            item = self._rootItem

        nchild = item.childCount()

        if item.value_type is dict:
            document = {}
            for i in range(nchild):
                ch = item.child(i)
                document[ch.key] = self.to_json(ch)
            return document

        elif item.value_type == list:
            document = []
            for i in range(nchild):
                ch = item.child(i)
                document.append(self.to_json(ch))
            return document

        else:
            return item.value

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
        layout = QVBoxLayout(self)
        pl = QGridLayout(self)
        sl = QGridLayout(self)
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
        print(f"data = {data}")

        for i in range(self.n_p):
            pl.addWidget(QLabel(f"Name of pigment {i:d}:"), i, 0)
            self.pigment_names.append(QLineEdit())
            pl.addWidget(self.pigment_names[i], i, 1)
        for i in range(self.n_s):
            sl.addWidget(QLabel(f"Name of state {i:d}:"), i, 0)
            self.state_names.append(QLineEdit())
            sl.addWidget(self.state_names[i], i, 1)
            sl.addWidget(QComboBox(), i, 2)
        layout.addLayout(pl)
        layout.addLayout(sl)

    def updateData(self, data):
        data["pigment_names"] = [p.text() for p in self.pigment_names]
        data["state_names"] = [s.text() for s in self.state_names]


class ProteinBuilder(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.data = {}
        layout = QVBoxLayout(self)

        self.Stack = QtWidgets.QStackedWidget()
        self.Stack.addWidget(NameNumber(self.data))
        self.Stack.addWidget(NamePigmentsStates(self.data))
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
        self.Stack.currentWidget().updateData(self.data)
        nxi = self.Stack.currentIndex() + 1
        self.Stack.widget(nxi).update()
        self.Stack.setCurrentIndex(nxi)
        print(f"data = {self.data}")

    def onPrevious(self):
        self.Stack.currentWidget().updateData(self.data)
        pxi = self.Stack.currentIndex() - 1
        self.Stack.widget(pxi).update()
        self.Stack.setCurrentIndex(pxi)
        print(f"data = {self.data}")

if __name__ == "__main__":
    app =QtWidgets.QApplication([])

    widget = ProteinBuilder()
    widget.show()

    sys.exit(app.exec())
