#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to plot output results in the interface. 

@author: maxime - maxime.bernard@uni-potsdam.de

"""


import os
from PyQt5.QtWidgets import (QWidget, QScrollArea,QGridLayout,
                             QVBoxLayout,QHBoxLayout,QGroupBox,QLabel,QComboBox,
                             QCheckBox,QLineEdit,QRadioButton,QFrame,QSpacerItem,
                             QSizePolicy)
from PyQt5.QtGui import (QIntValidator)
from PyQt5.QtCore import Qt

import configs as conf

#-----------------------------------------------------------------------------
#-------------------------------- Classes ------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

class CoolingRates_Parameters(QWidget):
     """ 
     This class defines the parameters to plot 2D map of cooling rates.
     It asks the user to provide range of temperature or time values from
     which the cooling rates will be calculated.
     
     #see class GraphWin, def "loadData"
     """
     def __init__(self,parent):
         super().__init__()   
         
         # Parameters
         self.data = 'temperature'
         self.TitleLabel = QLabel("Parameters for 2D map of cooling rates")
         self.TitleLabel.setFont(conf.fontBold12)
         self.minTempLabel = QLabel("min:")
         self.minTempLabel.setFont(conf.font10)
         self.minTempLabel.setAlignment(Qt.AlignCenter)
         self.minTempLabel.setToolTip("Temperature of the lower limit.")
         self.minTemp = QLineEdit()
         self.maxTempLabel = QLabel("max:")
         self.maxTempLabel.setAlignment(Qt.AlignCenter)
         self.maxTempLabel.setFont(conf.font10)
         self.maxTemp = QLineEdit()
         self.maxTempLabel.setToolTip("Temperature of the higher limit.")
         self.ThemeSelectionLabel = QLabel("Calculate cooling rates from:")
         self.ThemeSelectionLabel.setFont(conf.font11)
         self.ThemeSelectionLabel.setAlignment(Qt.AlignLeft)
         self.TempSelection = QRadioButton("temperature")
         self.TempSelection.setChecked(True)
         self.timeSelection = QRadioButton("time")
         self.timeSelection.setChecked(False)
         self.rangeLabel = QLabel("range (°C):")
         self.rangeLabel.setFont(conf.font10)
         self.rangeLabel.setAlignment(Qt.AlignRight)
         
         # Interactive slots
         self.TempSelection.toggled.connect(lambda: self.updateLabel(self.TempSelection))
         self.timeSelection.toggled.connect(lambda: self.updateLabel(self.timeSelection))

         # Add to layout
         self.Box1 = QGroupBox()
         layout = QGridLayout()
         self.vBox = QVBoxLayout()
         
         layout.addWidget(self.TempSelection,1,1)
         layout.addWidget(self.timeSelection,1,2)
         layout.addWidget(self.minTempLabel,2,1)
         layout.addWidget(self.maxTempLabel,2,2)
         layout.addWidget(self.rangeLabel,3,0)
         layout.addWidget(self.minTemp,3,1)
         layout.addWidget(self.maxTemp,3,2)
         layout.setSpacing(10)
         layout.setColumnStretch(0,1)
         layout.setColumnStretch(4,1)
         self.vBox.addWidget(self.TitleLabel)
         self.vBox.addStretch(2)
         self.vBox.addWidget(self.ThemeSelectionLabel)
         self.vBox.addStretch(1)
         self.vBox.addLayout(layout)
         self.vBox.addSpacerItem(QSpacerItem(150, 20, QSizePolicy.Expanding))

     #-----------------------------------------------------
     def updateLabel(self, item):
        """ To change some label """
        if item.text() == "temperature":
            self.rangeLabel.setText("range (°C):")
            self.data = 'temperature'
        else:
            self.rangeLabel.setText("range (Myr):")
            self.data = 'time' 




##################################################################################           
class plot2DVTK(QWidget):
    """
    This class defines parameters to plot 2D maps of data contained in a vtk file
    #see class GraphWin, def "plotData"
    
    """
    def __init__(self,parent,datalist):
        super().__init__()
        
        # Signals
        self.interpolate = 0
        self.data = 'Depth'
        self.interpolation = 'linear'
        self.normalized = 0
        self.datalist = datalist
        
        # Layouts
        self.ScrollArea = QScrollArea()
        self.layout = QGridLayout()
        layout2 = QVBoxLayout()
        layout2.setSpacing(20)
        Grid1 = QGridLayout()
        HBox = QHBoxLayout()
        HBox2 = QHBoxLayout()
        HBox3 = QHBoxLayout()
        self.Box1 = QGroupBox("Temperature")
        self.Box1.setCheckable(True)
        self.Box1.setChecked(False)
        self.Box1.setEnabled(False) 
        if parent.dataVTK_name == 'Temperature':
            self.Box1.setEnabled(True)   
        
        #### Parameters ####
        # Data
        self.Label = QLabel("Please choose the data to plot:")
        self.Label.setFont(conf.fontBold12)
        self.DataSelection = QComboBox()
        self.DataSelection.addItems(self.datalist)
        self.normalizedCheckBox = QCheckBox("use normalized values")
        # Contour plot
        self.PlotDetailsLabel = QLabel('Set the plot characteristics:')
        self.PlotDetailsLabel.setFont(conf.fontBold12)
        self.ContourLabel = QLabel("Draw contour plot?")
        self.ContourLabel.setFont(conf.font11)
        self.ContourCheck = QCheckBox()
        self.ContourCheck.setChecked(False)
        self.nconLabel = QLabel("Number of contours:")
        self.nconLabel.setFont(conf.font10)
        self.nconLabel.setAlignment(Qt.AlignRight)
        self.ncon = QLineEdit()
        self.ncon.setValidator(QIntValidator())
        self.ncon.setEnabled(False)
        # Temperature field
        self.ThemeSelectionLabel = QLabel("data to interpolate:")
        self.ThemeSelectionLabel.setFont(conf.font10)
        self.ThemeSelectionLabel.setAlignment(Qt.AlignRight)
        self.depthSelection = QRadioButton("depth")
        self.depthSelection.setChecked(True)
        self.isothermSelection = QRadioButton("isotherm")
        self.isothermSelection.setChecked(False)
        self.isothermLabel = QLabel("set depth (km):")
        self.isothermLabel.setFont(conf.font10)
        self.isothermLabel.setAlignment(Qt.AlignRight)
        self.isotherm = QLineEdit()
        self.isotherm.setValidator(conf.DoubleValidator)
        self.interpMethodLabel = QLabel("interpolation:")
        self.interpMethodLabel.setFont(conf.font10)
        self.interpMethodLabel.setAlignment(Qt.AlignRight)
        self.interpMethodLabel.setToolTip("for more details, see Scipy.interp1d")
        self.interpMethod = QComboBox()
        items = ['linear','nearest','slinear','quadratic','cubic']
        self.interpMethod.addItems(items)
        
        # Interactive slots
        self.contours = 0
        self.ContourCheck.toggled.connect(lambda: self.updateCheckBox(self.ContourCheck))
        self.DataSelection.currentIndexChanged.connect(lambda:self.updateBoxes(self.DataSelection))
        self.Box1.toggled.connect(lambda: self.activateInterp(self.Box1))
        self.isothermSelection.toggled.connect(lambda: self.updateLabel(self.isothermSelection))
        self.depthSelection.toggled.connect(lambda: self.updateLabel(self.depthSelection))
        self.interpMethod.currentIndexChanged.connect(lambda: self.updateInterp(self.interpMethod))
        self.normalizedCheckBox.stateChanged.connect(lambda: self.do_normalization())
        
        # Add to layout
        # data
        HBox.addWidget(self.DataSelection)
        HBox.addWidget(self.normalizedCheckBox)
        HBox.addStretch(1)
        # Contour
        HBox2.addWidget(self.ContourLabel)
        HBox2.addWidget(self.ContourCheck,0)
        HBox2.addStretch(1)
        HBox3.addStretch(1)
        HBox3.addWidget(self.nconLabel)
        HBox3.addWidget(self.ncon,)
        HBox3.addStretch(1)
        # Temperature field
        Grid1.addWidget(self.ThemeSelectionLabel, 0, 0)
        Grid1.addWidget(self.depthSelection, 0, 1)
        Grid1.addWidget(self.isothermSelection, 0, 2)
        Grid1.addWidget(self.isothermLabel, 2, 0)
        Grid1.addWidget(self.isotherm, 2, 1)
        Grid1.addWidget(self.interpMethodLabel, 3, 0)
        Grid1.addWidget(self.interpMethod, 3, 1)
        Grid1.setSpacing(10)
        Grid1.setColumnStretch(3,1)
        self.Box1.setLayout(Grid1)
        
        # Separator
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        
        # then the main layout
        Q = QWidget()
        layout2.addWidget(self.Label)
        layout2.addLayout(HBox)
        layout2.setStretch(2,1)
        layout2.addWidget(Separator)
        layout2.addWidget(self.PlotDetailsLabel)
        layout2.addLayout(HBox2)
        layout2.addLayout(HBox3)
        layout2.setStretch(7,2)
        layout2.addWidget(self.Box1)
        Q.setLayout(layout2)
        self.ScrollArea.setWidget(Q)
        self.layout.addWidget(self.ScrollArea,0,0,6,1)
    
    #-----------------------------------------------------
    def do_normalization(self):
        """to normalize the data to the range of values"""
        if self.normalizedCheckBox.isChecked():
            self.normalized = 1
        else:
            self.normalized = 0
            
    #-----------------------------------------------------
    def updateInterp(self, selection):
        """ To set the interpolation method """
        self.interpolation = selection.currentText()
        
    #-----------------------------------------------------
    def activateInterp(self, box):
        """ To signal the user wants to interpolate data """
        if box.isChecked() == True:
            self.interpolate = 1
            self.depthSelection.setEnabled(True)
            self.isothermSelection.setEnabled(True)
        else:
            self.interpolate = 0
            self.depthSelection.setEnabled(False)
            self.isothermSelection.setEnabled(False)
        
    #-----------------------------------------------------
    def updateLabel(self, item):
        """ To change some label """
        if item.text() == "depth":
            self.isothermLabel.setText("set depth (km):")
            self.data = 'Depth'
        else:
            self.isothermLabel.setText("set isotherm (°C):")
            self.data = 'Isotherm'
            
    #-----------------------------------------------------
    def updateBoxes(self, item):
        """ To enable/disable some boxes """
        if item.currentText() == "Temperature":
            self.Box1.setEnabled(True)
        else:
            self.Box1.setEnabled(False)
            
    #-----------------------------------------------------
    def updateCheckBox(self,checkBox):
        if checkBox.isChecked() == True:
            self.contours = 1
            self.ncon.setEnabled(True)
        else:
            self.contours = 0
            self.ncon.setEnabled(False)