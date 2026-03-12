#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This module contains functions and classes used in the interface to plot model outputs.

@author: Maxime Bernard

"""


import os
import configs as conf
from PyQt5.QtWidgets import (QWidget, QPushButton,QVBoxLayout, QHBoxLayout, QGridLayout,
                             QTabWidget, QLabel,QMessageBox,QComboBox,QFileDialog,
                             QDockWidget,QTreeView,QListView,QGroupBox,QErrorMessage)

from PyQt5.QtCore import Qt, QSize
from PyQt5.Qt import QStandardItemModel, QStandardItem

from math import *
import pandas as pd

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.collections import LineCollection
import traceback
import pyvista as pv
import pyvistaqt as pvqt
import xarray as xr
import Utils.PGUI_utils as pgu
from Utils import misfits
import Utils.interface as U_interface
import Utils.graphics.set_Thermochronometers as plot_thermo
import Utils.graphics.batch_results as batch
import Utils.graphics.plot_NA as plot_NA
import Utils.graphics.load_results as load_results
import Utils.graphics.plot_2Dmaps as plot_2D
import Utils.interface_inputs as GUI_inputs


import Thermochronology.thermochronometers as th
import Thermochronology.settings as Thermo_settings


# Path of executable
FolderPath = conf.FolderPath


class GraphWin(pvqt.MainWindow):
    """
    This class handle the Window hosting the plots.
    parent = MainWindow
    
    """
    def __init__(self, parent):
        super().__init__()
        
        # -------------- Define plot area ----------------------
        # Define new plot space
        self.MainWin = parent
        self.UsePredictedElevation = self.MainWin.UsePredictedElevation
        self.PecubePath = conf.PecubeFolderPath
        self.FolderPath = FolderPath
        self.Object3D = 0.
        self.data = 0
        self.dataFiles = []
        self.CountPlotsRow = 0
        self.CountPlotsColumn = 0
        self.SavePlots = []
        self.SaveToolbar = []
        self.List3DProp = []
        self.Canvas3D = []
        self.He4He3 = 0
        self.FitFunct = 1
        self.poly = 0
        self.CurrentFolder = 'Nil'
        self.TabPlot = QTabWidget()
        self.TabPlot.setTabBar(U_interface.HorizontalTabWidget())
        self.TabPlot.setTabsClosable(True)
        self.TabPlot.tabCloseRequested.connect(self.close_current_tab)
        self.tab2D = QWidget()  # for 2D plots
        self.TabPlot.addTab(self.tab2D, "Plots")
        self.TabPlot.currentChanged.connect(lambda:self.show_Properties(self.TabPlot))
        self.MainMesh = None #set Main Mesh to None
        self.Samples_loc = None #set the samples location to none
        self.indexSplitter = 0
        #Colorbar properties
        self.sargs = dict(height=0.1, width=0.2, vertical=True, position_x=0.95, position_y=0.05,
                      label_font_size=14,interactive=True,color=conf.Colors3Dplot['colorbarLabels'])
        #Config state: 2D (=0) or 3D plot (=1)
        self.activePlot = 0
        # #Initiate list of items ID (2D or 3D plot?) that will be added to plot
        self.itemID = []
        
        # -------- Define side bar -------------
        # Define Data Tab
        # Add data button
        AddDataButton = QPushButton("Load project...")
        AddDataButton.resize(QSize(30, 20))
        AddDataButton.adjustSize()
        AddDataButton.clicked.connect(lambda: self.loadData(1))
        ModelButton = QPushButton("Add 3D model...")
        ModelButton.resize(QSize(30, 20))
        ModelButton.adjustSize()
        ModelButton.clicked.connect(lambda: self.load3DData())
        RemoveButton = QPushButton("Remove data...")
        RemoveButton.resize(QSize(30, 20))
        RemoveButton.clicked.connect(lambda: self.RemoveData())
        self.DataCombo = QComboBox()
        self.DataCombo.addItem("Further data...\t")
        if self.MainWin.ParamTable:
            self.folderName = self.MainWin.ParamTable.PNAME
            self.FolderPath = FolderPath
            self.get_dataList()         
                
        self.DataCombo.activated.connect(lambda:self.Show_data(self.DataCombo))
        HLay = QHBoxLayout()
        HLay.addWidget(AddDataButton)
        HLay.addWidget(ModelButton)
        HLay.addStretch(1)
        # Treeview with data loaded
        self.treeView = QTreeView()
        self.treeView.setHeaderHidden(True)
        #for 3d Models
        self.treeView3D = QTreeView()
        self.treeView3D.setHeaderHidden(True)
        self.rootNode = QListView()
        self.rootNode3D = QListView()
        self.treeModel = QStandardItemModel(self.rootNode)
        self.treeModel.itemChanged.connect(lambda: self.ShowData())
        self.treeModel3D = QStandardItemModel(self.rootNode3D)
        self.treeView.setModel(self.treeModel)
        self.treeView.expandAll()
        
        #Labels for TreeViews
        Label2D = QLabel("data:")
        Label2D.setFont(conf.fontBold10)
        Label3D = QLabel("3D models:")
        Label3D.setFont(conf.fontBold10)

        self.VLay = QVBoxLayout()
        self.VLay.addLayout(HLay)
        HLay2 = QHBoxLayout()
        HLay2.addWidget(RemoveButton)
        HLay2.addStretch(2)
        self.VLay.addLayout(HLay2)
        HLay3 = QHBoxLayout()
        HLay3.addWidget(self.DataCombo)
        HLay3.addStretch(1)
        self.VLay.addLayout(HLay3)
        self.VLay.addWidget(Label2D)
        self.VLay.addWidget(self.treeView)
        self.VLay.addWidget(Label3D)
        self.VLay.addWidget(self.treeView3D)

        widgetDock = QWidget()
        widgetDock.setLayout(self.VLay)
        self.DockData = QDockWidget("Data", self)
        self.DockData.setWidget(widgetDock)
        self.DockData.setMinimumWidth(300)

        # Plot Area
        # Define Properties Tab
        ########## Create plot properties ##########
        ####### 2D plot properties  ######
        self.Dock2dProperties = QDockWidget("Properties", self)
        self.Dock2dProperties.setMinimumWidth(300)
        
        #Data Properties
        DataProp = QLabel("Data properties")
        DataProp.setFont(conf.fontBold12)
        DataProp.setAlignment(Qt.AlignLeft)
        #Plot properties
        PlotProp = QLabel("Plot properties")
        PlotProp.setFont(conf.fontBold12)
        PlotProp.setAlignment(Qt.AlignLeft)
        #Colorbar range
        CbRange = QLabel("data range:")
        CbRange.setFont(conf.fontBold8)
        CbRange.setAlignment(Qt.AlignLeft)
        
        ####### Build the grid #######
        VBoxMain = QVBoxLayout()
        #Data Properties
        BoxData = QGroupBox("Data properties")
        GboxDataProp = QGridLayout()
        GboxDataProp.setRowStretch(7,1)
        BoxData.setLayout(GboxDataProp)
        VBoxMain.addWidget(BoxData)
        #Plot properties
        BoxPlot = QGroupBox("Plot properties")
        GboxPlotProp = QGridLayout()
        GboxPlotProp.setRowStretch(5,1)
        GboxPlotProp.setSpacing(10)
        BoxPlot.setLayout(GboxPlotProp)
        # VBoxMain.addStretch(1)
        VBoxMain.addWidget(BoxPlot)
        widgetDockProp = QWidget()
        widgetDockProp.setLayout(VBoxMain)
        self.Dock2dProperties.setWidget(widgetDockProp)
        
        # ------------- Combine layouts ---------------------
        self.VBox = QVBoxLayout()
        # self.MainLay = QGridLayout()
        # Create the 6 spots for 2D plots
        self.MainLay = U_interface.QCustomSplitter(Qt.Vertical)
        self.Splitter1 = U_interface.QCustomSplitter(Qt.Horizontal)
        self.Splitter2 = U_interface.QCustomSplitter(Qt.Horizontal)
        FakeLabel = QWidget()
        FakeLabel2 = QWidget()
        FakeLabel3 = QWidget()
        FakeLabel4 = QWidget()
        FakeLabel5 = QWidget()
        FakeLabel6 = QWidget()
        self.Splitter1.addWidget(FakeLabel)
        self.Splitter1.addWidget(FakeLabel2)
        self.Splitter1.addWidget(FakeLabel3)
        self.Splitter2.addWidget(FakeLabel4)
        self.Splitter2.addWidget(FakeLabel5)
        self.Splitter2.addWidget(FakeLabel6)
        # Hide them
        self.Splitter1.widget(0).setVisible(False)
        self.Splitter1.widget(1).setVisible(False)
        self.Splitter1.widget(2).setVisible(False)
        self.Splitter2.widget(0).setVisible(False)
        self.Splitter2.widget(1).setVisible(False)
        self.Splitter2.widget(2).setVisible(False)
        self.Splitter2.setVisible(False)
        # Add to main splitter
        self.MainLay.addWidget(self.Splitter1)
        self.MainLay.addWidget(self.Splitter2)
        VLayout = QVBoxLayout()
        VLayout.addWidget(self.MainLay)
        # MainLay.addLayout(SideBar.VLay,1)
        self.MainWin.addDockWidget(Qt.LeftDockWidgetArea, self.DockData)
        self.MainWin.addDockWidget(Qt.LeftDockWidgetArea, self.Dock2dProperties)
        self.MainWin.tabifyDockWidget(self.DockData, self.Dock2dProperties)
        self.MainWin.setTabPosition(Qt.AllDockWidgetAreas, QTabWidget.North)
        self.MainWin.setDockOptions(
           self. MainWin.GroupedDragging | self.MainWin.AllowTabbedDocks | self.MainWin.AllowNestedDocks)
        self.widget = QWidget()
        self.tab2D.setLayout(VLayout)
        self.VBox.addWidget(self.TabPlot)
        #Set the Layout to the Output layout of the main Window
        self.MainWin.OutputWidget.setLayout(self.VBox)
    
    #----------------------------------------------------------
    def get_dataList(self):
        """
        To get the list of available data within the current project
        data can be: Age-elevation, Tt paths, Map cooling rates, Map Temperature, Map ages, date-eU
        
        """
        # Load inputs parameters from files and assign useful parameters
        try:
            oldParameters = pgu.readOldInput(self.folderName)
        except Exception as E:
            tb = traceback.format_exc()
            QErrorMessage(self).showMessage('Error with loading data:' + str(tb) + str(E))
            return
        self.ParametersInput = oldParameters.ParametersInput
        self.MainWin.ParamTable = GUI_inputs.ParamWin(self.MainWin,self.folderName, self.ParametersInput, oldParameters, 2)
        # The number of samples
        self.nbSamples = int(self.ParametersInput.DParameters['Nb_Samples'])
        self.nbGrain_per_sample = self.ParametersInput.TableParameters['nb_grains']
        
        self.DataCombo.clear()
        self.dataFiles = ["CompareAGE.csv"]
        path = os.path.join(self.PecubePath,self.folderName,"output","CompareAGE.csv")

        if os.path.exists(path):
            self.DataCombo.addItem("Age comparison")
            self.DataCombo.addItem("Elevation comparison")
            self.DataCombo.addItem("Date-eU")
            self.DataCombo.addItem("Age transect")
            if int(self.ParametersInput.DParameters[conf.Variable_names['Use Fission Track']]) > 0:
                self.DataCombo.addItem("MFTL comparison")
                self.DataCombo.addItem("MFTL vs elevation")
        # Is there any trapped-charge predictions ?
        path = os.path.join(self.PecubePath,self.folderName,"output",conf.Variable_names["ThL_file"])
        path2 = os.path.join(self.PecubePath,self.folderName,"output",conf.Variable_names["OSL_file"])
        path3 = os.path.join(self.PecubePath,self.folderName,"output",conf.Variable_names["ESR_file"])
        if os.path.exists(path) or os.path.exists(path2) or os.path.exists(path3):
            self.DataCombo.addItem("Trapped-charge comparison")
            self.DataCombo.addItem("nN_elevation")
        if os.path.exists(path):
            self.dataFiles.append(conf.Variable_names["ThL_file"])
            self.dataFiles.append(conf.Variable_names["ThL_file_input"])
        if os.path.exists(path2):
            self.dataFiles.append(conf.Variable_names["OSL_file"])
            self.dataFiles.append(conf.Variable_names["OSL_file_input"])
        if os.path.exists(path3):
            self.dataFiles.append(conf.Variable_names["ESR_file"])
            self.dataFiles.append(conf.Variable_names["ESR_file_input"])
        # 
        path = os.path.join(self.PecubePath,self.folderName,"output")
        if os.path.exists(path):
            filenames = os.listdir(path)
            target1 = [file for file in filenames if file.startswith('TimeAge')]
            target2 = [file for file in filenames if file.startswith('Ages0')]
            if target1 or target2:
                self.DataCombo.addItem("Age-elevation")
                self.dataFiles.extend(target1)
                self.dataFiles.extend(target2)
                
        # Is there any CoolingRates.csv ?
        if os.path.exists(os.path.join(self.PecubePath,self.folderName,"output","CoolingRates.csv")):
            self.DataCombo.addItem("Map cooling rates")
            self.dataFiles.append("CoolingRates.csv")
            
        # Is there any TimeTemperaturePaths.csv ?
        if os.path.exists(os.path.join(self.PecubePath,self.folderName,"output","TimeTemperaturePaths.csv")):
            self.DataCombo.addItem("Tt paths") 
            self.dataFiles.append("TimeTemperaturePaths.csv")
            
        # Is there any luminescence prediction ?
        path = os.path.join(self.PecubePath,self.folderName,"output")
        if os.path.exists(path):
            filenames = os.listdir(path)
            target1 = [file for file in filenames if file.startswith('ESRTime') or file.startswith('OSLTime') or file.startswith('TLTime')]
            if target1:
                self.DataCombo.addItem("nN vs time")
                self.dataFiles.extend(target1)
                
        # Is there any PecubeXXX.vtk ?
        path = os.path.join(self.PecubePath,self.folderName,"VTK")
        if os.path.exists(path):
            filenames = os.listdir(path)
            target = [file for file in filenames if file.startswith('Pecube')]
            if target:
                self.DataCombo.addItem("Map temperatures") 
            # Is there any AgeXXX.vtk ?
            target = [file for file in filenames if file.startswith('Age')]
            if target:
                self.DataCombo.addItem("Map ages") 
        self.CurrentFolder = self.folderName #store the current project
        if int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2:
            self.dataFiles.append(self.ParametersInput.DParameters['data_folder']+".csv")
 
        if os.path.exists(os.path.join(self.PecubePath,self.folderName,"output","Compare43He.csv")):
            self.DataCombo.addItem("4He/3He data")
            self.dataFiles.append("Compare43He.csv")
            self.dataFiles.append(conf.Variable_names["43_inputFile"])
            
        if not self.DataCombo.currentText():
            self.DataCombo.addItem("Further data...\t")
        
        # Create dataset - only for sample-specific predictions
        if int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2:
            self.dataset = load_results.dataSets(self.dataFiles,self.folderName,self.ParametersInput.DParameters['data_folder'])
        
        # Is there any inversion file ?
        path = os.path.join(self.PecubePath,self.folderName,"NA","NA_int_res.csv")
        path2 = os.path.join(self.PecubePath,self.folderName,"NA","nab.out")
        if int(self.ParametersInput.DParameters[conf.Variable_names['use inversion or batch']]) < 2: # inversion mode
            if os.path.exists(path) or os.path.exists(path2):
                self.DataCombo.addItem("Inversion results")
        
        elif int(self.ParametersInput.DParameters[conf.Variable_names['use inversion or batch']]) == 2: # batch mode
            if os.path.exists(os.path.join(self.PecubePath,self.folderName,"NA","models.in")):
                self.DataCombo.addItem("Batch results")
        
        
    
    #----------------------------------------------------------
    def Show_data(self,d):
        """
        To plot data from the combo box. The data are chosen by the user.
        
        """
        # For sample-specific predictions
        if d.currentText() == 'Age-elevation' and int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2:
            # Update data Combo box to plot Age-elevation profile
            filename = os.path.join(self.PecubePath,self.FolderPath,"output","CompareAGE.csv")
            # Get AHe ages from TimeAge***.csv file
            Ages, AgeError = self.get_specific_ages()
            # Get observed elevations of samples
            samples_elevations = self.dataset.MainDataset.HEIGHTPRED.values
            # Get observed ages and errors
            AgesObs, AgeErrorObs = self.getObsAges()
            # Get dic Pecube
            agecol, errname, agename, predname, colores, profdict, He43dict,coloramp, ageMarker  =  Thermo_settings.dict_pecube()
            CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
            inputData = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
            # Plot Age-elevation profile
            try:
                self.plot_AgeData('Age elevation',self.thermochronometers,CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                        predname=predname,ageMarker=ageMarker)
            except Exception as E:
                QErrorMessage(self).showMessage('Something wrong happen with plot Age elevation:' + str(E) + '. Please, make sure your input data location are all included in your DEM.')
                return
        
        # For all model nodes predictions
        elif d.currentText() == 'Age-elevation' and int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 1:
            # Load data from Ages_n.csv
            # Check number of Ages step in the folder
            files = os.listdir(os.path.join(self.PecubePath,self.FolderPath,"output"))
            files_list = []
            for i in range(len(files)):
                if files[i][0:3] == 'Age':
                    files_list.append(files[i])
            
            # Ask the user to choose the file
            # Select file to open
            msg = QMessageBox.question(self, 'Pecube Message',
                                        'Several files have been found... Please select one or more.',
                                        QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Ok)
            if msg == QMessageBox.Ok:
                dialog = QFileDialog()
                dialog.setDirectory(os.path.join(self.PecubePath,self.folderName,"Output"))
                options = QFileDialog.Options()
                options |= QFileDialog.DontUseNativeDialog
                self.filePath, _ = dialog.getOpenFileNames(
                    self, "QFileDialog.getOpenFileName()", "", "CSV Files (Ages*.csv)", options=options)
                if self.filePath:
                    # Get thermochronometers chose by the user
                    self.System_Thermo = th.ThermoAges()
                    self.wThermo = U_interface.WinExtraParameters(500, 400, text='Thermochronological systems')
                    self.wThermo.OkButton.clicked.connect(lambda: self.get_ThermoSys(self.filePath))
                    Label = QLabel("Which thermochronological ages, do you want to plot?")
                    Label.setFont(conf.fontBold12)
                    self.wThermo.layout.addWidget(Label, 0, 0)
                    self.wThermo.layout.addLayout(self.System_Thermo.VMainBox, 1, 0)
                    Q2 = QWidget()
                    Q2.setLayout(self.wThermo.layout)
                    self.wThermo.setCentralWidget(Q2)
                    conf.WindowsOpen.append(self.wThermo)
                    self.wThermo.show()
            else:
                return
        
        # Plot observations vs predictions
        elif d.currentText() == "Age comparison":
            #Get ages from TimeAge.csv file
            Ages,AgeError = self.get_specific_ages()
            # Get observed ages
            AgesObs, AgeErrorObs = self.getObsAges()
            # Get dic Pecube
            agecol, errname, agename, predname, colores, profdict, He43dict,coloramp,ageMarker =  Thermo_settings.dict_pecube()
            CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
            inputData = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
            try:
                self.plotComparisons('Age comparison',self.thermochronometers,CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                     ageMarker=ageMarker)
            except Exception as E:
                QErrorMessage(self).showMessage('Something wrong happen with pplot obs vs pred: '+ str(E)+ '. Please, make sure your input data location are all included in your DEM.')
                return
        
        # Plot observations vs predictions for fission track
        elif d.currentText() == "MFTL comparison":
            self.windowParam = U_interface.WinExtraParameters(400, 200, text='MFTL data',fixed=False)
            self.parameters = plot_thermo.ChooseMFTL()
            self.windowParam.layout.addLayout(self.parameters.Layout, 0, 0)
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_MFTL("MFTL comparison",self.parameters.Data2Plot)) 
            
        
        # Plot MFTL vs elevation for fission track
        elif d.currentText() == "MFTL vs elevation" and int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2:
            self.windowParam = U_interface.WinExtraParameters(400, 200, text='MFTL data',fixed=False)
            self.parameters = plot_thermo.ChooseMFTL()
            self.windowParam.layout.addLayout(self.parameters.Layout, 0, 0)
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_MFTL("MFTL elevation",self.parameters.Data2Plot))
            
        
        # Plot observations vs predictions
        elif d.currentText() == "Trapped-charge comparison":
            #Get data from ThL_data.csv and/or OSL_data and/or ESR_data.csv file(s)
            # Ages,AgeError = self.get_Trapped_charge()
            # Get dic Pecube
            agecol, errname, agename, predname, colores, profdict, He43dict,coloramp,ageMarker =  Thermo_settings.dict_pecube()
            CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
            #Empyt dataset
            Pred = 0
            Obs = 0
            TrapCharSys = []
            ThL = 0
            try:
                ThLData = self.dataset.MainDataset.filter_by_attrs(description="ThL_data_input")
                Obs = xr.merge([Obs,ThLData])
                ThLPred = self.dataset.MainDataset.filter_by_attrs(description="ThL_data")
                Pred = xr.merge([Pred,ThLPred])
                if ThLData:
                    TrapCharSys.append('THL')
                ThL = 1
            except Exception:
                pass
            OSL = 0
            try:
                OSLData = self.dataset.MainDataset.filter_by_attrs(description="OSL_data_input")
                OSLPred = self.dataset.MainDataset.filter_by_attrs(description="OSL_data")
                if not Obs:
                    Obs = OSLData
                    Pred = OSLPred
                else:
                    Obs = xr.merge([Obs,OSLData])
                    Pred = xr.merge([Pred,OSLPred])
                if OSLData:
                    TrapCharSys.append('OSL')
                OSL = 1
            except:
                pass
            ESR = 0
            try:
                ESRData = self.dataset.MainDataset.filter_by_attrs(description="ESR_data_input")
                ESRPred = self.dataset.MainDataset.filter_by_attrs(description="ESR_data")
                if not Obs:
                    Obs = ESRData
                    Pred = ESRPred
                else:
                    Obs = xr.merge([Obs,ESRData])
                    Pred = xr.merge([Pred,ESRPred])
                if ESRData:
                    TrapCharSys.append('ESR')
                ESR = 1
            except:
                pass
            
            try:
                self.plotComparisons('Trapped-charge comparison',TrapCharSys,Pred,inputc=Obs,agecol=agecol,errname=errname,colores=colores,ageMarker=ageMarker)
            except Exception as E:
                QErrorMessage(self).showMessage('Something wrong happen with pplot obs vs pred: ' + str(E) + '. Please, make sure your input data location are all included in your DEM.')
                return
        
        # Plot nN_elevation
        elif d.currentText() == "nN_elevation":
            #Get data from ThL_data.csv and/or OSL_data and/or ESR_data.csv file(s)
            # Get dic Pecube
            agecol, errname, agename, predname, colores, profdict, He43dict,coloramp,ageMarker =  Thermo_settings.dict_pecube()
            CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
            #Empyt dataset
            Pred = 0
            Obs = 0
            TrapCharSys = []
            ThL = 0
            try:
                ThLData = self.dataset.MainDataset.filter_by_attrs(description="ThL_data_input")
                Obs = xr.merge([Obs,ThLData])
                ThLPred = self.dataset.MainDataset.filter_by_attrs(description="ThL_data")
                Pred = xr.merge([Pred,ThLPred])
                if ThLData:
                    TrapCharSys.append('THL')
                ThL = 1
            except Exception:
                pass
            OSL = 0
            try:
                OSLData = self.dataset.MainDataset.filter_by_attrs(description="OSL_data_input")
                OSLPred = self.dataset.MainDataset.filter_by_attrs(description="OSL_data")
                if not Obs:
                    Obs = OSLData
                    Pred = OSLPred
                else:
                    Obs = xr.merge([Obs,OSLData])
                    Pred = xr.merge([Pred,OSLPred])
                if OSLData:
                    TrapCharSys.append('OSL')
                OSL = 1
            except:
                pass
            ESR = 0
            try:
                ESRData = self.dataset.MainDataset.filter_by_attrs(description="ESR_data_input")
                ESRPred = self.dataset.MainDataset.filter_by_attrs(description="ESR_data")
                if not Obs:
                    Obs = ESRData
                    Pred = ESRPred
                else:
                    Obs = xr.merge([Obs,ESRData])
                    Pred = xr.merge([Pred,ESRPred])
                if ESRData:
                    TrapCharSys.append('ESR')
                ESR = 1
            except:
                pass
            
            try:
                self.plot_AgeData('Trapped charge',TrapCharSys,Pred,inputc=Obs,agecol=agecol,errname=errname,colores=colores,
                                        predname=predname,ageMarker=ageMarker)
            except:
                QErrorMessage(self).showMessage('Something wrong happen with pplot obs vs pred. Please, make sure your input data location are all included in your DEM.')
                return
            
        # Plot observations vs predictions
        elif d.currentText() == "Elevation comparison":
             #Get ages from TimeAge.csv file
             Ages,AgeError = self.get_specific_ages()
             # Get observed ages
             AgesObs, AgeErrorObs = self.getObsAges()
             # Get dic Pecube
             agecol, errname, agename, predname, colores, profdict,He43dict,coloramp,ageMarker  =  Thermo_settings.dict_pecube()
             CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
             inputData = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
             try:
                 self.plotComparisons('Elevation comparison',self.thermochronometers,CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,ageMarker=ageMarker)
             except Exception as E:
                 QErrorMessage(self).showMessage('Something wrong happen with pplot obs vs pred: ' + str(E))
                 return
            
        # Plot Date vs eU, Dates vs distance 
        elif d.currentText() == "Date-eU":
            if int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2:
                #Get ages from TimeAge.csv file
                Ages,AgeError = self.get_specific_ages()
                # Get observed ages
                AgesObs, AgeErrorObs = self.getObsAges()
                #Get eU from Samples_settings file
                eU_array = self.get_eU_array()

                # Get dic Pecube
                agecol, errname, agename, predname, colores, profdict,He43dict,coloramp,ageMarker =  Thermo_settings.dict_pecube()
                CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
                inputData = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
                # Plot Age-elevation profile
                try:
                    self.plot_AgeData('Date vs eU',self.thermochronometers,CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                            predname=predname,ageMarker=ageMarker)
                except Exception as E:
                    QErrorMessage(self).showMessage('Something wrong happen with plot Age elevation:' + str(E) + 
                                                    '. Please, make sure your input data location are all included in your DEM.')
                    return
                # try:
                #     self.plot_AgeElevation('Date-eU',eU_array,AgeError,Ages,AgesObs,AgeErrorObs, 'eU (ppm)','Age (Ma)')
                # except:
                #     QErrorMessage(self).showMessage('Something wrong happen with plot_AgeElevation')
                #     return
            else:
                return
            
        elif d.currentText() == "Age transect":
            # Get ages accross distance
            Ages, AgeError = self.get_specific_ages()
            # Get distances
            # Distance_array = self.get_sample_location()
            # Get observed ages
            AgesObs, AgeErrorObs = self.getObsAges()
            # Get dic Pecube
            agecol, errname, agename, predname, colores, profdict,He43dict,coloramp,ageMarker  =  Thermo_settings.dict_pecube()
            CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
            inputc = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
            # Create extra window to ask information from the user
            self.windowParam = U_interface.WinExtraParameters(700, 600, text='Age transect',fixed=False)
            self.parameters = plot_thermo.transect('age transect',self.thermochronometers,CompareAge,inputc,agecol,errname,agename,colores,predname,ageMarker,profdict,self)
            try:
                self.windowParam.layout.addLayout(self.parameters.MainLayout, 0, 0) 
            except:
                return
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_transect()) 
            
        # 2D maps
        elif d.currentText() == "Map cooling rates":
            self.filePath = [os.path.join(self.PecubePath,self.folderName,"output","CoolingRates.csv")]
            self.loadData(0)
            
        elif d.currentText() == "Tt paths":
            self.filePath = [os.path.join(self.PecubePath,self.folderName,"output","TimeTemperaturePaths.csv")]
            self.loadData(0)
        
        elif d.currentText() == "nN vs time":
            # Create extra window to ask information from the user
            self.windowParam = U_interface.WinExtraParameters(300, 400, text='Age evolution',fixed=False)
            outputPath = os.path.join(self.PecubePath,self.folderName,"output")
            self.parameters = plot_thermo.AgeTime(os.path.join(outputPath),self)
            try:
                self.windowParam.layout.addLayout(self.parameters.MainLayout, 0, 0) 
            except:
                return
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_AgeTime()) 
            
        elif d.currentText() == "Map temperatures":
            # The user has to choose a PecubeXXX.vtk file
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            dialog = QFileDialog()
            dialog.setDirectory(os.path.join(self.PecubePath,self.folderName,"VTK"))
            self.filePath, _ = dialog.getOpenFileNames(
                self, "QFileDialog.getOpenFileName()", "", "(Pecube*.vtk)",options=options)
            self.loadData(0)
            
        elif d.currentText() == "Map ages":
            # The user has to choose a AgeXXX.vtk file
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            dialog = QFileDialog()
            dialog.setDirectory(os.path.join(self.PecubePath,self.folderName,"VTK"))
            self.filePath, _ = dialog.getOpenFileNames(
                self, "QFileDialog.getOpenFileName()", "", "(Age*.vtk)",options=options)
            self.loadData(0)
            
        elif d.currentText() == '4He/3He data':
            Compare43 = self.dataset.MainDataset.filter_by_attrs(description="Compare43He")
            inputc = self.dataset.MainDataset.filter_by_attrs(description="input43He")
            # Create extra window to ask information from the user
            self.windowParam = U_interface.WinExtraParameters(500, 400, text='43He plot',fixed=False)
            self.parameters = plot_thermo.He43_plot('4He/3He',Compare43,inputc,self)
            try:
                self.windowParam.layout.addLayout(self.parameters.MainLayout, 0, 0) 
            except:
                return
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_43He()) 
        
        elif d.currentText() == 'Inversion results':
            # Create extra window to ask information from the user
            self.windowParam = U_interface.WinExtraParameters(200, 400, text='Inversion results',fixed=False)
            path = os.path.join(self.PecubePath,self.folderName,"NA","NA_int_res.csv")
            path2 = os.path.join(self.PecubePath,self.folderName,"NA","nab.out")
            FolderNA = os.path.join(self.PecubePath,self.folderName,"NA")
            self.parameters = plot_NA.plotPecubeNA(path,path2,FolderNA,self)
            self.windowParam.layout.addLayout(self.parameters.VBox, 0, 0)   
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_Inversion()) 
        
        elif d.currentText() == 'Batch results':
            # Create extra window to ask information from the user
            self.windowParam = U_interface.WinExtraParameters(300, 300, text='Batch results',fixed=False)
            ProjectFolder = os.path.join(self.PecubePath,self.folderName)
            self.parameters = batch.AgeBatch(ProjectFolder,self)
            self.windowParam.layout.addLayout(self.parameters.Layout, 0, 0)
            Q = QWidget()
            Q.setLayout(self.windowParam.layout)
            self.windowParam.setCentralWidget(Q)
            conf.WindowsOpen.append(self.windowParam)
            self.windowParam.show()
            # 3) Once the user chose, extract the data from the file
            self.windowParam.OkButton.clicked.connect(lambda: self.plot_Batch()) 
            
    #----------------------------------------------------------
    def plot_transect(self):
        """Call out.transect to plot age transect, then close the window """
        try:
            self.parameters.plotTransect()
        except Exception as E:
            QErrorMessage(self).showMessage("Something wrong happened with age transect: <br>" +  str(E))
        self.windowParam.close()
    
    #----------------------------------------------------------
    def plot_AgeTime(self):
        """ plot age vs time, according to user's choice of thermochronometers. """
        try:
            self.parameters.plot_ageTime()
        except Exception as E:
            QErrorMessage(self).showMessage("Something wrong happened with age vs time: <br>" +  str(E))
        self.windowParam.close()
        
    #----------------------------------------------------------
    def plot_43He(self):
        """Call out.He43_plot to plot 43He profiles, then close the window """
        try:
            self.parameters.plotHe43()
        except Exception as E:
            raise
            QErrorMessage(self).showMessage("Something wrong happened with 43He plot: <br>" +  str(E))
        self.windowParam.close()
        
    #----------------------------------------------------------
    def plot_Batch(self):
        """ Plot results from batch models """
        
        self.parameters.plot_Batch_results()
        self.windowParam.close()
        
    #----------------------------------------------------------
    def plot_Inversion(self):
        """Call out.plotPecubeNA to plot inversion results, then close the window """
        path_NAB = os.path.join(self.PecubePath,self.folderName,"NA","nab.out")
        try:
            if self.parameters.ChoosePlotCombo.currentIndex() == 2: #Plot 2D parameter space with 1d pdf 
                if not os.path.exists(path_NAB):
                    QErrorMessage(self).showMessage("File not found: nab.out. Plot with NA_results.py.")
                   
                self.parameters.plot_results(PDF_1D=True,PDF_2D=True)
                
            elif self.parameters.ChoosePlotCombo.currentIndex() == 3: # Plot 1D pdf only
                if not os.path.exists(path_NAB):
                    QErrorMessage(self).showMessage("File not found: nab.out. Plot with NA_results.py.")
                    
                self.parameters.plot_results(PDF_1D=True,PDF1D_only=True)
                
            elif self.parameters.ChoosePlotCombo.currentIndex() == 1: #Plot 2D parameter space only
                self.parameters.plot_results()
            
            else: # Misfit evolution$
                self.parameters.plot_results(Misfit_evol = True)
                
        except Exception as E:
            QErrorMessage(self).showMessage("Something wrong happened with inversion parameters:<br>" + str(E))
            raise
        self.windowParam.close()
        
    #----------------------------------------------------------
    def plot_MFTL(self,inputdata,data2plot):
        """ plot MFTL data or STD MFTL data."""
        
        # Get dic Pecube
        agecol, errname, agename, predname, colores, profdict, He43dict,coloramp,ageMarker =  Thermo_settings.dict_pecube()
        CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
        inputData = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
        if inputdata == "MFTL comparison":
            try:
                if data2plot == 0: #MFTL
                    self.plotComparisons('MFTL comparison',['MFTL'],CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                     ageMarker=ageMarker)
                else: #STD MFTL
                    self.plotComparisons('MFTL comparison',['DMFTL'],CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                     ageMarker=ageMarker)
            except Exception as E:
                QErrorMessage(self).showMessage('Something wrong happen with plot obs vs pred: ' + str(E))
                return
        elif inputdata == "MFTL elevation":
            # Get dic Pecube
            agecol, errname, agename, predname, colores, profdict, He43dict,coloramp, ageMarker  =  Thermo_settings.dict_pecube()
            CompareAge = self.dataset.MainDataset.filter_by_attrs(description="CompareAges")
            inputData = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
            # Plot Age-elevation profile
            try:
                if data2plot == 0: #MFTL
                    self.plot_AgeData('MFTL elevation',['MFTL'],CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                        predname=predname,ageMarker=ageMarker)
                else:
                    self.plot_AgeData('MFTL elevation',['DMFTL'],CompareAge,inputc=inputData,agecol=agecol,errname=errname,colores=colores,
                                        predname=predname,ageMarker=ageMarker)
            except Exception as E:
                QErrorMessage(self).showMessage('Something wrong happen with plot MFTL vs elevation:' + str(E))
                return
        self.windowParam.close()
        
        
    #----------------------------------------------------------
    def show_Properties(self,tab):
        #to set the properties tab according to whether is a 2D or 3D plot
        if tab.tabText(tab.currentIndex()) == "3D Viewer":
            self.Dock2dProperties.hide()
            index = tab.currentIndex()
            self.List3DProp[index-1].show()
            for i in range(tab.count()):
                if tab.tabText(i) == "3D Viewer" and not i == index:
                    self.List3DProp[i-1].hide()
        else:
            self.Dock2dProperties.show()
            for i in range(tab.count()):
                if tab.tabText(i) == "3D Viewer":
                    self.List3DProp[i-1].hide()
    
    #----------------------------------------------------------  
    def close_current_tab(self,i):
        #To newly created tabs
        #if only one tab
        if self.TabPlot.count() <2:
            #do nothing
            return
        #else remove tab
        if self.TabPlot.tabText(i)=="3D Viewer":
            self.Canvas3D[i-1].plotter.remove_actor(self.MainMesh)
            try:
                self.Canvas3D[i-1].plotter.remove_actor('Sample_location')
            except:
                pass
            self.Canvas3D[i-1].plotter.clear()
            self.Canvas3D[i-1].plotter.close()
            self.Canvas3D[i-1].plotter.deep_clean()
            self.List3DProp[i-1].hide()
            self.List3DProp.pop(i-1)
            self.Canvas3D.pop(i-1)
            self.treeModel3D.removeRow(i-1)
        self.TabPlot.removeTab(i)
    
    #----------------------------------------------------------
    def getObsAges(self):
        # To read observed ages from Samples_settings file
        AgeObs_array = []
        ErrorObs_array = []
        self.elevObs = []
        self.MarkerShapeObs = []
        self.ThermoSysObs = []
        self.thermochronometers = []
        if  int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2:
            CAges_array = self.dataset.MainDataset.filter_by_attrs(description="InputSample")
            varNames = CAges_array.keys()
            # Is there AHe ages observed?
            if "AHE" in varNames:
                AgeObs_array.append(CAges_array.AHE.values)
                ErrorObs_array.append(CAges_array.DAHE.values)
                self.MarkerShapeObs.append(Thermo_settings.Marker_dict['AHE'])
                self.ThermoSysObs.append('AHeObs')
                self.thermochronometers.append("AHE")
            # Is there AFT ages observed?
            if "AFT" in varNames:
                AgeObs_array.append(CAges_array.AFT.values)
                ErrorObs_array.append(CAges_array.DAFT.values)
                self.MarkerShapeObs.append(Thermo_settings.Marker_dict['AFT'])
                self.ThermoSysObs.append('AFTObs')
                self.thermochronometers.append("AFT")
            # Is there ZHe ages observed?
            if "ZHE" in varNames:
                AgeObs_array.append(CAges_array.ZHE.values)
                ErrorObs_array.append(CAges_array.DZHE.values)
                self.MarkerShapeObs.append(Thermo_settings.Marker_dict['ZHE'])
                self.ThermoSysObs.append('ZHeObs')
                self.thermochronometers.append("ZHE")
             # Is there ZFT ages observed?
            if "ZFT" in varNames:
                 AgeObs_array.append(CAges_array.ZFT.values)
                 ErrorObs_array.append(CAges_array.DZFT.values)
                 self.MarkerShapeObs.append(Thermo_settings.Marker_dict['ZFT'])
                 self.ThermoSysObs.append('ZFTObs')
                 self.thermochronometers.append("ZFT")
             # Is there ZHe ages observed?
            if 'KAR' in varNames:
                 AgeObs_array.append(CAges_array.KAR.values)
                 ErrorObs_array.append(CAges_array.DKAR.values)
                 self.MarkerShapeObs.append(Thermo_settings.Marker_dict['KAR'])
                 self.ThermoSysObs.append('KArObs')
                 self.thermochronometers.append("KAR")
              # Is there ZHe ages observed?
            if "BAR" in varNames:
                  AgeObs_array.append(CAges_array.BAR.values)
                  ErrorObs_array.append(CAges_array.DBAR.values)
                  self.MarkerShapeObs.append(Thermo_settings.Marker_dict['BAR'])
                  self.ThermoSysObs.append('BARObs')
                  self.thermochronometers.append("BAR")
              # Is there ZHe ages observed?
            if "MAR" in varNames:
                  AgeObs_array.append(CAges_array.MAR.values)
                  ErrorObs_array.append(CAges_array.DMAR.values)
                  self.MarkerShapeObs.append(Thermo_settings.Marker_dict['MAR'])
                  self.ThermoSysObs.append('MARObs')
                  self.thermochronometers.append("MAR")
              # Is there ZHe ages observed?
            if "HAR" in varNames:
                  AgeObs_array.append(CAges_array.HAR.values)
                  ErrorObs_array.append(CAges_array.DHAR.values)
                  self.MarkerShapeObs.append(Thermo_settings.Marker_dict['HAR'])
                  self.ThermoSysObs.append('HARObs')
                  self.thermochronometers.append("HAR")

        else:
            for i in self.list_sys:
                if i =='HeApatite':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeAHE")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDAHE")])
                elif i == 'HeZircon':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeZHE")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDZHE")])
                elif i == 'FTApatite':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeAFT")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDAFT")])
                elif i == 'FTZircon':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeZFT")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDZFT")])
                elif i == 'ArKFeldspar':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeKAR")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDKAR")])
                elif i == 'ArBiotite':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeBAR")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDBAR")])
                elif i == 'ArMuscovite':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeMAR")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDMAR")])
                elif i == 'ArHornblend':
                    AgeObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("AgeHAR")])
                    ErrorObs_array.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("ErrorDHAR")])

        AgeObs_array = np.asarray(AgeObs_array)
        ErrorObs_array = np.asarray(ErrorObs_array)
        # Get elevations
        self.elevObs.append([float(v) for k, v in self.ParametersInput.TableParameters.items() if k.startswith("Elev_grain")])
        self.elevObs = np.asarray(self.elevObs)
        return AgeObs_array, ErrorObs_array
        
    #----------------------------------------------------------
    def get_eU_array(self):
        """
        To read eU values from Samples_settings file.
        
        """
        if os.path.exists(os.path.join(self.PecubePath,self.folderName,"Data","Samples_settings.txt")):
            eU_array = np.zeros((len(self.ThermSys),self.Ages.shape[1]))
            # for j in range(eU_array.shape[0]):
            #     count = 0
            #     if self.ThermSys[j] == 'AHe':
            #         mineral = 'Apatite'
            #     elif self.ThermSys[j] == 'ZHe':
            #         mineral = 'Zircon'
                    
            # Get U and Th ppm apatite
            try:
                U = self.dataset.MainDataset.AUPPM.values
                Th = self.dataset.MainDataset.ATHPPM.values
                eU_array[0][:] = U + 0.235*Th
            except:
                print('No Uppm, Thppm found for AHe')
            # Get U and Th ppm zircon
            try:
                U = self.dataset.MainDataset.ZUPPM.values
                Th = self.dataset.MainDataset.ZTHPPM.values
                eU_array[1][:] = U + 0.235*Th #Put in list, for now only ZHE have Uppm Thppm
            except:
                print('No Uppm, Thppm found for ZHe')
                # for i in range(int(self.ParametersInput.DParameters['Nb_Samples'])): 
                #    nb_grains = self.ParametersInput.TableParameters['nb_grains'][i][0]
                #    for ii in range(int(nb_grains)):
                #        U = float(self.ParametersInput.TableParameters['U_grain'+mineral+str(i+1)+"_"+str(ii+1)])
                #        Th = float(self.ParametersInput.TableParameters['Th_grain'+mineral+str(i+1)+"_"+str(ii+1)])
                #        eU_array[j][count] = U + 0.235*Th
                #        count+=1
            return eU_array
        else:
            return None
        
    #----------------------------------------------------------
    def get_sample_location(self):
        """
        To get the locations of the samples along X or Y axes.
        
        """
        folderName = self.ParametersInput.DParameters['data_folder']
        if os.path.exists(os.path.join(self.PecubePath,self.folderName,"Data",folderName,folderName+".csv")):
            nb_ThermoSys = len(self.ThermSys)
            Distance_array = np.zeros((nb_ThermoSys,int(self.ParametersInput.DParameters['Nb_TotGrains'])))
            lat_grain =np.zeros((int(self.ParametersInput.DParameters['Nb_TotGrains']),1))
            lon_grain = np.zeros((int(self.ParametersInput.DParameters['Nb_TotGrains']),1))
            for i in range(int(self.ParametersInput.DParameters['Nb_Samples'])):
                sampleName =  self.ParametersInput.TableParameters['SampleName'+str(i+1)]
                lat_grain[i] = self.ParametersInput.TableParameters[sampleName+'_lat']
                lon_grain[i] = self.ParametersInput.TableParameters[sampleName+'_lon']
            # Get coordinates
            y0 = float(self.MainWin.ParamTable.ParamTable.lat0Edit5.text())
            dy = float(self.MainWin.ParamTable.ParamTable.dLatEdit7.text())
            dx = float(self.MainWin.ParamTable.ParamTable.dLonEdit6.text())
            ny = int(self.MainWin.ParamTable.ParamTable.nyEdit3.text())
            nx = int(self.MainWin.ParamTable.ParamTable.nxEdit2.text())
            x0 = float(self.MainWin.ParamTable.ParamTable.lon0Edit4.text())
            #Calculate distance between points and left bottom corner
            xl = x0 + (nx-1)*dx
            yl = y0 + (ny-1)*dy
            xlKm = dx*(nx-1)*111.11*cos((y0+dy*ny/2)*3.141592654/180)
            ylKm = dy*(ny-1)*111.11
            lat_grain = (lat_grain-y0)/(yl-y0)*ylKm
            lon_grain = (lon_grain-x0)/(xl-x0)*xlKm
            for i in range(nb_ThermoSys):
                k=0
                for j in range(len(self.nbGrain_per_sample)):
                    ng = int(self.nbGrain_per_sample[j][0])
                    Distance_array[i][k:k+ng] = lat_grain[j][0]
                    k += ng
            return Distance_array
        else:
            return None
        

    #----------------------------------------------------------
    def plot_AgeElevation(self,inputdata, dataplot,
                        agecol = None, errname = None, colores = None, agename = None,
                        predname = None,ageMarker=None, graphtitle = None):
        """
        Author: Maxime bernard, University of Potsdam. 04 May 2023
        based on Xavier's plotComparisons function
        
        Function to plot:
            - age vs elevation for all nodes models

        Args:
            inputdata (string): path file 'AgeXXX.csv'

            dataplot (list): List of data to plot ; By default, the altitude will be plotted
                             ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
                             Rem: For the moment, MTL and TTp not implemented.
                             Defaults to ['AHe','AFT'].

            agecol (dict, optional): respective column number of the data in the file comparison.txt.
                                     Defaults to None.
            
            errname (dict, optional): respective column number of the error on data in the data input file. 
                                      Defaults to None.
            
            colores (dict, optional): Colors used for the different age system. 
                                      Defaults to None.
            
            agename (dict, optional): legend of each data system. 
                                      Defaults to None.

            graphpath (str, optional): name of the folder where the plot will be written
                                       Usually you do not have to change it. 
                                       Defaults to 'Graphs'.

            graphtitle (str, optional): title to write on the graph. 
                                        Defaults to None.

        """

        # Create canva and toolbar
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        self.toolbar = NavigationToolbar(self.plotSpace,self)
        self.plotSpace.axes.set_box_aspect(1)
        # Read file 
        
        for file in inputdata:
            dataFrame = pd.read_csv(file)
            # Transform to xarray dataset
            Dataset = dataFrame.to_xarray()
            # Rename variables
            datac = Dataset.rename({"HeApatite": 'AHE',
                                                "HeZircon":"ZHE",
                                                "FTApatite":"AFT",
                                                "FTZircon":"ZFT",
                                                "ArKFeldspar":"KAR",
                                                "ArBiotite":"BAR",
                                                "ArMuscovite":"MAR",
                                                "ArHornblend":"HAR",
                                                "nNtl":agecol['THL'],
                                                "nNosl":agecol['OSL'],
                                                'nNesr':agecol['ESR'],
                                                "Height":"HEIGHT"})

            for item in dataplot:
                print(u'\tPlotting %s age elevation' %(str(item)))
                
                try:
                    pred = datac[agecol[item]][:]
                    ElevPred = datac[agecol['alt']][:]
                    # Plot prediction
                    self.plotSpace.axes.plot(pred, ElevPred,
                              marker = ageMarker[item], linestyle = 'None', 
                              label = predname[item], color = colores[item])
      
                except ValueError:
                    print("\n Issue with age-elevation plot, do not plot")

        self.plotSpace.axes.set_xlabel(u'Ages (Ma)')
        self.plotSpace.axes.set_ylabel(u'Elevation (m)') 
                    
        self.plotSpace.axes.legend(loc = 'best')

        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        self.add_plot2win(Q,self.plotSpace,self.toolbar,'Age-elevation')
        
        
    #----------------------------------------------------------
    def plot_AgeData(self,inputdata, dataplot,
                        datac, inputc = None,
                        agecol = None, errname = None, colores = None, agename = None,
                        predname = None, ageMarker=None, graphtitle = None):
        """
        Author: Maxime bernard, University of Potsdam. 04 May 2023
        based on Xavier's plotComparisons function
        
        Function to plot:
            - age vs elevation 
            - Date vs eU
            - age vs grain size
        
        with data and predictions

        Args:
            inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars

            dataplot (list): List of data to plot ; By default, the altitude will be plotted
                             ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
                             Rem: For the moment, MTL and TTp not implemented.
                             Defaults to ['AHe','AFT'].

            datac (xarray): results of Pecube forward modeling 

            inputc (array, optional): Array of data to plot.
                                      Defaults to None.

            agecol (dict, optional): respective column number of the data in the file comparison.txt.
                                     Defaults to None.
            
            errname (dict, optional): respective column number of the error on data in the data input file. 
                                      Defaults to None.
            
            colores (dict, optional): Colors used for the different age system. 
                                      Defaults to None.
            
            agename (dict, optional): legend of each data system. 
                                      Defaults to None.

            graphpath (str, optional): name of the folder where the plot will be written
                                       Usually you do not have to change it. 
                                       Defaults to 'Graphs'.

            graphtitle (str, optional): title to write on the graph. 
                                        Defaults to None.

        """

        # Create canva and toolbar
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        self.toolbar = NavigationToolbar(self.plotSpace,self)
        self.plotSpace.axes.set_box_aspect(1)
        
        for item in dataplot:
            if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
                if item == 'Altitude':
                    # inputdata = None
                    print(u'\tPlotting %s elevation' %(str(item)))
                    item = 'alt'
                else:
                    item = item.upper()
                    print(u'\tPlotting %s age elevation' %(str(item)))
                # fig2 = plt.figure()
                # Plot the age comparison
                # And add error bars on data
                if inputdata == 'Age elevation':
                    # Because age uncertainties are in an input file, ensure we target the right uncertainties to their age
                    # Get sample name of data > -9999
                    sampleName = inputc["SAMPLE"][np.logical_not(inputc[agecol[item]] != inputc[agecol[item]])]
                    sampleNamePred = datac['SID'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    
                    # sort data by sample (alphabetically)
                    sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                    sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                    sampleNameIndexes = np.argsort(sampleNametemp)
                    sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                    index_samplesObs = sampleName.coords
                    index_samplesPred = sampleNamePred.coords
                    errortemp = inputc[errname[item]][np.logical_not(inputc[errname[item]] != inputc[errname[item]])]
                    errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                    errortemp = errortemp2[sampleNamePredIndexes]
                    try:
                        obs = datac[agecol[item]+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        obs = obs[sampleNameIndexes]
                        pred = datac[agecol[item]+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        pred = pred[sampleNameIndexes]
                        
                        ElevPred = datac[agecol['alt']+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        ElevPred = ElevPred[sampleNameIndexes]
                        ElevObs = datac[agecol['alt']+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        ElevObs = ElevObs[sampleNameIndexes]
                        if self.UsePredictedElevation.currentIndex() == 1:
                            ElevPred = ElevObs
                        # If we have negative elevation (i.e. tunnel samples) put negative values for observations
                        # ElevPred2 = np.asarray([-s if pgu.isNeg(s) == True else s for s in ElevObs])
                        # ElevPred = ElevPred2
                        # ElevObs2 = np.asarray([-s if pgu.isNeg(s) == True else s for s in ElevObs])
                        # ElevObs = ElevObs2
                        
                        # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                        # of the array
                        print(obs.shape,pred.shape,errortemp.shape)
                        error = np.empty(pred.shape)
                        error[:] = np.nan
                        count = 0
                        for n in range(len(pred)):
                            if not pgu.isNaN(errortemp[n]):
                                error[count] = errortemp[n]
                                count += 1
                        # Plot observation with error bars
                        segments = [
                            [(xi - err, yi), (xi + err, yi)]
                            for xi, yi, err in zip(obs, ElevObs, error)
                        ]

                        # Create LineCollection for all error bars
                        error_lines = LineCollection(segments, colors=colores[item], linewidths=1)
                        self.plotSpace.axes.add_collection(error_lines)
                        self.plotSpace.axes.plot(obs, ElevObs, marker = ageMarker[item],
                                    linestyle = 'None',  label = 'observed_'+item, color = colores[item],zorder=9)

                        # Plot prediction
                        self.plotSpace.axes.plot(pred, ElevPred,
                                  marker = ageMarker[item], linestyle = 'None', 
                                  label = 'predicted_'+predname[item], color = colores[item], alpha = 0.3,zorder = 10)
  
                    except ValueError:
                        print("\n Issue with errors, do not plot")
                        self.plotSpace.axes.plot(datac[agecol[item]+'OBS'], datac[agecol['alt']+'OBS'], 
                                      marker = ageMarker[item], linestyle = 'None',
                                     label = 'observed_'+item, 
                                     color = colores[item],alpha = 0.3)
                        self.plotSpace.axes.plot(datac[agecol[item]+'PRED'], datac[agecol['alt']+'PRED'], 
                                      marker = ageMarker[item], linestyle = 'None',
                                     label = 'predicted_'+predname[item], 
                                     color = colores[item],alpha = 0.3)
 
                    self.plotSpace.axes.set_xlabel(u'Ages (Ma)')
                    self.plotSpace.axes.set_ylabel(u'Elevation (m)') 
                    self.plotSpace.axes.legend(loc = 'best')
                    
                # -----------------------------------------------
                if inputdata == 'MFTL elevation':
                    item = dataplot[0]
                    # Because age uncertainties are in an input file, ensure we target the right uncertainties to their age
                    # Get sample name of data > -9999
                    sampleName = datac["SID"][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    sampleNamePred = inputc['SAMPLE'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    # sort data by sample (alphabetically)
                    sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                    sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                    sampleNameIndexes = np.argsort(sampleNametemp)
                    sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                    index = sampleName.coords
                    # errortemp = inputc[errname[item]][np.logical_not(inputc[errname[item]] == -9999)]
                    # errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                    # errortemp = errortemp2[sampleNamePredIndexes]
                    
                    # Error on MFTL and std are assumed constant
                    if item == conf.Variable_names["Mean Fission Track Length"]:
                        error_temp = float(self.ParametersInput.DParameters[conf.Variable_names['Mean_fission_track_length_error_inversion']])
                        xlabel = "Mean fission track lengths (µm)"
                    elif item == conf.Variable_names["Mean Fission Track Length error"]:
                        error_temp = float(self.ParametersInput.DParameters[conf.Variable_names['Mean_fission_track_length_std_error_inversion']])
                        xlabel = "std Mean fission track lengths (µm)"
                    try:
                        obs = datac[agecol[item]+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        obs = obs[sampleNameIndexes]
                        pred = datac[agecol[item]+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        pred = pred[sampleNameIndexes]
                        
                        ElevPred = datac[agecol['alt']+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        ElevPred = ElevPred[sampleNameIndexes]
                        ElevObs = datac[agecol['alt']+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        ElevObs = ElevObs[sampleNameIndexes]
                        if self.UsePredictedElevation.currentIndex() == 1:
                            ElevPred = ElevObs
                            
                        error = [error_temp]*len(pred)
                        # If we have negative elevation (i.e. tunnel samples) put negative values for observations
                        # ElevPred2 = np.asarray([-s if pgu.isNeg(s) == True else s for s in ElevObs])
                        # ElevPred = ElevPred2
                        # ElevObs2 = np.asarray([-s if pgu.isNeg(s) == True else s for s in ElevObs])
                        # ElevObs = ElevObs2
                        
                        # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                        # of the array
                        # error = np.empty(pred.shape)
                        # error[:] = np.nan
                        # count = 0
                        # for n in range(len(pred)):
                        #     if not pgu.isNaN(errortemp[n]):
                        #         error[count] = errortemp[n]
                        #         count += 1
                        # Plot observation
                        segments = [
                            [(xi - err, yi), (xi + err, yi)]
                            for xi, yi, err in zip(obs, ElevObs, error)
                        ]

                        # Create LineCollection for all error bars
                        error_lines = LineCollection(segments, colors=colores[item], linewidths=1)
                        self.plotSpace.axes.add_collection(error_lines)
                        self.plotSpace.axes.plot(obs, ElevObs, marker = ageMarker[item],
                                    linestyle = 'None',  label = 'observed_'+item, color = colores[item],zorder=9)
                        
                        # Plot prediction
                        self.plotSpace.axes.plot(pred, ElevPred,
                                  marker = ageMarker[item], linestyle = 'None', 
                                  label = 'predicted_'+predname[item], color = colores[item], alpha = 0.3,zorder = 10)
  
                    except ValueError:
                        print("\n Issue with errors, do not plot")
                        self.plotSpace.axes.plot(datac[agecol[item]+'OBS'], datac[agecol['alt']+'OBS'], 
                                      marker = ageMarker[item], linestyle = 'None',
                                     label = 'observed_'+item, 
                                     color = colores[item],alpha = 0.3)
                        self.plotSpace.axes.plot(datac[agecol[item]+'PRED'], datac[agecol['alt']+'PRED'], 
                                      marker = ageMarker[item], linestyle = 'None',
                                     label = 'predicted_'+predname[item], 
                                     color = colores[item],alpha = 0.3)
 
                    self.plotSpace.axes.set_xlabel(xlabel)
                    self.plotSpace.axes.set_ylabel(u'Elevation (m)') 
                    self.plotSpace.axes.legend(loc = 'best')
                    
                    
                if inputdata == 'Trapped charge':
                    for item in dataplot:
                        if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
                            # Get sample name of data > -9999
                            sampleName = inputc["SAMPLE"+item]
                            sampleNamePred = inputc['SAMPLE'+item]
                            # sort data by sample (alphabetically)
                            sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                            sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                            sampleNameIndexes = np.argsort(sampleNametemp)
                            sampleName_sorted = np.sort(sampleNametemp)
                            sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                            index = sampleName.coords
                            errortemp = inputc[errname[item]]
                            errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                            errortemp = errortemp2[sampleNamePredIndexes]
                            try:
                                obs = inputc[agecol[item]+'OBS_input']
                                pred = datac[agecol[item]+'PRED']
                                ElevPred = datac[agecol['alt']+'PRED_'+item]
                                ElevPred = ElevPred[sampleNameIndexes]
                                ElevObs = datac[agecol['alt']+'OBS_'+item]
                                ElevObs = ElevObs[sampleNameIndexes]

                                #Remove NaNs
                                pred2 = np.asarray([s for s in pred.values if pgu.isNaN(s) == False])
                                pred = pred2
                                pred = pred[sampleNameIndexes]
                                obs2 = np.asarray([s for s in obs.values if pgu.isNaN(s) == False])
                                obs = obs2[sampleNameIndexes]
                                ElevPred2 = np.asarray([s for s in ElevPred.values if pgu.isNaN(s) == False])
                                ElevPred = ElevPred2
                                ElevObs2 = np.asarray([s for s in ElevObs.values if pgu.isNaN(s) == False])
                                ElevObs = ElevObs2
                                if self.UsePredictedElevation.currentIndex() == 1:
                                    ElevPred = ElevObs2

                                # If we have negative elevation (i.e. tunnel samples) put negative values for observations
                                # ElevPred2 = np.asarray([-s if pgu.isNeg(s) == True else s for s in ElevObs ])
                                # ElevPred = ElevPred2
                                # ElevObs2 = np.asarray([-s if pgu.isNeg(s) == True else s for s in ElevObs])
                                # ElevObs = ElevObs2
                                # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                                # of the array
                                error = np.empty(pred.shape)
                                error[:] = np.nan
                                count = 0
    
                                for n in range(len(pred)):
                                    if not pgu.isNaN(errortemp[n]):
                                        error[count] = errortemp[n]
                                        count += 1
                                # print('error = ',error) 
                                # Plot observation
                                # plot sample individually (to change color)
                                segments = [
                                        [(xi - err, yi), (xi + err, yi)]
                                        for xi, yi, err in zip(obs, ElevObs, error)
                                    ]

                                # Create LineCollection for all error bars
                                error_lines = LineCollection(segments, colors=colores[item], linewidths=1)
                                self.plotSpace.axes.add_collection(error_lines)
                                p1 = self.plotSpace.axes.plot(obs, ElevObs, marker = ageMarker[item],
                                            linestyle = 'None',  label = 'observed_'+item, color = colores[item],zorder=9)
                                # for p in range(len(obs)-1):
                                #     self.plotSpace.axes.errorbar(obs[p+1], ElevObs[p+1],
                                #                   xerr = error[p+1],
                                #                   fmt = ageMarker[item], label = 'observed_'+item+sampleName_sorted[p], color = colores[item],zorder=9)
                                # Plot prediction
                                p2 = self.plotSpace.axes.plot(pred, ElevPred,
                                          marker = ageMarker[item], linestyle = 'None',
                                          label = 'predicted_'+predname[item], color = colores[item], alpha = 0.3,zorder=10)
                
          
                            except Exception as E:
                                print("\n Issue with errors, do not plot")
                                print(str(E))
                                self.plotSpace.axes.plot(obs, ElevObs, 
                                              marker = ageMarker[item], linestyle = 'None',
                                             label = 'observed_'+item, 
                                             color = colores[item],zorder=9)
                                self.plotSpace.axes.plot(pred,ElevPred, 
                                              marker = ageMarker[item], linestyle = 'None',
                                             label = 'predicted_'+predname[item], 
                                             color = colores[item],alpha = 0.3,zorder=10)
                                
                        if int(self.ParametersInput.DParameters[conf.Variable_names['ESR_misfit_target']]):
                              self.plotSpace.axes.set_xlabel(u'n/N')      
                        else:
                            self.plotSpace.axes.set_xlabel(u'Age (Ma)')
                        self.plotSpace.axes.set_ylabel(u'Elevation (m)') 
                        self.plotSpace.axes.legend(loc = 'best')
                        self.plotSpace.draw()       
                
                # Date vs eU plot
                if inputdata == 'Date vs eU' and (item == 'AHE' or item == 'ZHE'):
                    # Because age uncertainties are in an input file, ensure we target the right uncertainties to their age
                    # Get sample name of data > -9999
                    sampleName = datac["SID"][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    sampleNamePred = inputc['SAMPLE'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    index_samplesObs = sampleName.coords
                    index_samplesPred = sampleNamePred.coords
                    # sort data by sample (alphabetically)
                    sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                    sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                    sampleNameIndexes = np.argsort(sampleNametemp)
                    sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                    errortemp = inputc[errname[item]][np.logical_not(inputc[errname[item]] == -9999)]
                    errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                    errortemp = errortemp2[sampleNamePredIndexes]
                    try:
                        obs = datac[agecol[item]+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        obs = obs[sampleNameIndexes]
                        pred = datac[agecol[item]+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        pred = pred[sampleNameIndexes]
                        if item == 'AHE':
                            Utemp = inputc['AUPPM'][index_samplesObs]
                            Thtemp = inputc['ATHPPM'][index_samplesObs]
                            Utemp = Utemp[sampleNameIndexes]
                            Thtemp = Thtemp[sampleNameIndexes]
                        elif item == 'ZHE':
                            Utemp = inputc['ZUPPM'][index_samplesObs]
                            Thtemp = inputc['ZTHPPM'][index_samplesObs]
                            Utemp = Utemp[sampleNameIndexes]
                            Thtemp = Thtemp[sampleNameIndexes]
            
                        eUtemp = Utemp + 0.235*Thtemp
                        # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                        # of the array
                        error = np.empty(pred.shape)
                        error[:] = np.nan
                        eU = np.empty(pred.shape)
                        eU[:] = np.nan
                        count = 0
                        for n in range(len(pred)):
                            if not pgu.isNaN(errortemp[n]):
                                error[count] = errortemp[n]
                                eU[count] = eUtemp[n]
                                count += 1
                        # Plot observation
                        segments = [
                                [(xi, yi - err), (xi, yi + err)]
                                for xi, yi, err in zip(eU, obs, error)
                            ]

                        # Create LineCollection for all error bars
                        error_lines = LineCollection(segments, colors=colores[item], linewidths=1)
                        self.plotSpace.axes.add_collection(error_lines)
    
                        self.plotSpace.axes.plot(eU, obs, marker = ageMarker[item],
                                    linestyle = 'None',  label = 'observed_'+item, color = colores[item],zorder=9)
                        
                        # Plot prediction
                        self.plotSpace.axes.plot(eU,pred,
                                  marker = ageMarker[item], linestyle = 'None', 
                                  label = 'predicted_'+predname[item], color = colores[item], alpha = 0.3)
                        
                    except ValueError:
                        print("\n Issue with errors, do not plot")
                        self.plotSpace.axes.plot(eU,obs,
                                  marker = ageMarker[item], linestyle = 'None', 
                                  label = 'observed_'+item, color = colores[item], alpha = 0.3)
                        self.plotSpace.axes.plot(eU,pred,
                                  marker = ageMarker[item], linestyle = 'None', 
                                  label = 'predicted_'+predname[item], color = colores[item], alpha = 0.3)
                    
                    self.plotSpace.axes.set_ylabel(u'Ages (Ma)')
                    self.plotSpace.axes.set_xlabel(u'eU (ppm)') 
                    self.plotSpace.axes.legend(loc = 'best')
                    
        
        #Add toolbar with 2D plot
        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        
        self.add_plot2win(Q,self.plotSpace,self.toolbar,inputdata)
                    
                        
    #----------------------------------------------------------
    def plotComparisons(self,inputdata, dataplot,
                        datac, inputc = None,
                        agecol = None, errname = None, colores = None, agename = None,
                        ageMarker=None,
                        graphpath = 'Graphs', graphtitle = None):
        """
        Author: Xavier Robert, University of Grenoble
        adapted by: Maxime bernard, University of Potsdam. 04 May 2023
        
        Function to plot the comparison between data and predictions

        Args:
            inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars

            dataplot (list): List of data to plot ; By default, the altitude will be plotted
                             ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
                             Rem: For the moment, MTL and TTp not implemented.
                             Defaults to ['AHe','AFT'].

            datac (string): results of Pecube forward modeling

            inputc (array, optional): Array of data to plot.
                                      Defaults to None.

            agecol (dict, optional): respective column number of the data in the file comparison.txt.
                                     Defaults to None.
            
            errname (dict, optional): respective column number of the error on data in the data input file. 
                                      Defaults to None.
            
            colores (dict, optional): Colors used for the different age system. 
                                      Defaults to None.
            
            agename (dict, optional): legend of each data system. 
                                      Defaults to None.

            graphpath (str, optional): name of the folder where the plot will be written
                                       Usually you do not have to change it. 
                                       Defaults to 'Graphs'.

            graphtitle (str, optional): title to write on the graph. 
                                        Defaults to None.

        """

        # Create canva and toolbar
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        self.plotSpace.axes.set_aspect('equal')
        # plt.gca().set_aspect('equal')
        # Plot 1:1 line
        self.plotSpace.axes.plot([0,8000], [0,8000], marker = 'None', linestyle = '-', 
                 label = '1:1 line', color = 'lightgrey', alpha = 1)
        
        # Initiate LL
        LL_array = [] #LL for each thermochronometers
        sigma_array = [] # misfit for each thermochronometers
        N_LL = [] #LL for each thermochronometers
        Item_list = [] # Where to store thermochronometer
        minValues = []
        maxValues = []
        shift = 2 # shift is applied to x and y lim
        ErrorSuccess = 0
        if inputdata == 'Age comparison':
            errorpred = float(self.ParametersInput.DParameters['error_predictions']) / 100
            for item in dataplot:
                if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
                    if item == 'Altitude':
                        # inputdata = None
                        print(u'\tPlotting %s comparison' %(str(item)))
                        item = 'alt'
                    else:
                        item = item.upper()
                        print(u'\tPlotting %s age comparison' %(str(item)))
                    # fig2 = plt.figure()
                    # Plot the age comparison
                    # And add error bars on data
                   
                    # Because age uncertainties are in an input file, ensure we target the right uncertainties to their age
                    # Get sample name of data > -9999
                    sampleName = datac["SID"][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    sampleNamePred = inputc['SAMPLE'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                    # sort data by sample (alphabetically)
                    sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                    sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                    sampleNameIndexes = np.argsort(sampleNametemp)
                    sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                    index = sampleName.coords
                    errortemp = inputc[errname[item]][np.logical_not(inputc[errname[item]] == -9999)]
                    errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                    errortemp = errortemp2[sampleNamePredIndexes]
                    
                    try:
                        obs = datac[agecol[item]+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        obs = obs[sampleNameIndexes]
                        pred = datac[agecol[item]+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        pred = pred[sampleNameIndexes]
                        errorpred_temp = errorpred*pred
                        errorpred_temp = np.asarray([s for s in errorpred_temp.values if pgu.isNaN(s) == False])

                        # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                        # of the array
                        error = np.empty(pred.shape)
                        error[:] = np.nan
                        count = 0
                        for n in range(len(pred)):
                            if not pgu.isNaN(errortemp[n]):
                                error[count] = errortemp[n]
                                count += 1

                        # Plot
                        self.plotSpace.axes.errorbar(obs, pred,
                                  xerr = error, yerr = errorpred_temp,
                                  fmt = ageMarker[item], label = item, color = colores[item])
                        
                        # Compute log likelyhood
                        LL, N, sigma = misfits.compute_misfit(int(self.ParametersInput.DParameters['misfit_function']),obs.values,pred.values,error,errorpred_temp)
                        LL_array.append(LL)
                        sigma_array.append(sigma)
                        N_LL.append(N)
                        Item_list.append(item)
                        
                        ErrorSuccess = 1
                    except ValueError as VE:
                        print("\n Issue with errors, do not plot.", VE)
                        self.plotSpace.axes.plot(datac[agecol[item]+'OBS'], datac[agecol[item]+'PRED'], 
                                      marker = ageMarker[item], linestyle = 'None',
                                     label = item, 
                                     color = colores[item])
                    minValues.append(min(abs(obs.values)))
                    minValues.append(min(abs(pred.values)))
                    maxValues.append(max(abs(obs.values))+shift)
                    maxValues.append(max(abs(pred.values))+shift)
               
            # Write the 1:1 text over the line
            # We may need to revise the way to find the x/y
            self.plotSpace.axes.text((max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2,
                     (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2 - 
                          (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/20,
                     '1:1', rotation = 45, color = 'lightgrey', alpha = 1, ha = 'center', va = 'center')

            self.plotSpace.axes.set_ylim(bottom = min(minValues[:]), top = max(maxValues[:]))
            self.plotSpace.axes.set_xlim(left = min(minValues[:]), right = max(maxValues[:]))

            #plt.xlabel(u'%s Observed (Ma)' %(str(item)))
            #plt.ylabel(u'%s Predicted (Ma)' %(str(item)))
            self.plotSpace.axes.set_xlabel(u'Observed ages (Ma)')
            self.plotSpace.axes.set_ylabel(u'Predicted ages (Ma)') 
                            
        elif inputdata == 'MFTL comparison':
            errorpred = float(self.ParametersInput.DParameters['error_predictions']) / 100
            item = dataplot[0]
            # Because age uncertainties are in an input file, ensure we target the right uncertainties to their age
            # Get sample name of data > -9999
            sampleName = datac["SID"][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
            sampleNamePred = inputc['SAMPLE'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
            # sort data by sample (alphabetically)
            sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
            sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
            sampleNameIndexes = np.argsort(sampleNametemp)
            sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
            index = sampleName.coords
            # errortemp = inputc[errname[item]][np.logical_not(inputc[errname[item]] == -9999)]
            # errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
            # errortemp = errortemp2[sampleNamePredIndexes]
            
            # Error on MFTL and std are assumed constant
            if item == conf.Variable_names["Mean Fission Track Length"]:
                error_temp = float(self.ParametersInput.DParameters[conf.Variable_names['Mean_fission_track_length_error_inversion']])
            elif item == conf.Variable_names["Mean Fission Track Length error"]:
                error_temp = float(self.ParametersInput.DParameters[conf.Variable_names['Mean_fission_track_length_std_error_inversion']])
            
            try:
                obs = datac[agecol[item]+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                obs = obs[sampleNameIndexes]
                pred = datac[agecol[item]+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                pred = pred[sampleNameIndexes]
                errorpred_temp = errorpred*pred
                errorpred_temp = np.asarray([s for s in errorpred_temp.values if pgu.isNaN(s) == False])
                error = [error_temp]*len(pred)
                # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                # of the array
                # error = np.empty(pred.shape)
                # error[:] = np.nan
                # count = 0
                # for n in range(len(pred)):
                #     if not pgu.isNaN(errortemp[n]):
                #         error[count] = errortemp[n]
                #         count += 1
                # Plot
                self.plotSpace.axes.errorbar(obs, pred,
                          xerr = error, yerr = errorpred*pred,
                          fmt = ageMarker[item], label = item, color = colores[item])
                
                # Compute log likelyhood
                LL, N, sigma = misfits.compute_misfit(int(self.ParametersInput.DParameters['misfit_function']),obs.values,pred.values,error,errorpred_temp)
                LL_array.append(LL)
                sigma_array.append(sigma)
                N_LL.append(N)
                Item_list.append(item)
                
                ErrorSuccess = 1
            except ValueError as VE:
                print("\n Issue with errors, do not plot.", VE)
                self.plotSpace.axes.plot(datac[agecol[item]+'OBS'], datac[agecol[item]+'PRED'], 
                              marker = ageMarker[item], linestyle = 'None',
                             label = item, 
                             color = colores[item])
            minValues.append(min(abs(obs.values)))
            minValues.append(min(abs(pred.values)))
            maxValues.append(max(abs(obs.values))+0.3)
            maxValues.append(max(abs(pred.values))+0.3)
           
            # Write the 1:1 text over the line
            # We may need to revise the way to find the x/y
            self.plotSpace.axes.text((max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2,
                     (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2 - 
                          (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/20,
                     '1:1', rotation = 45, color = 'lightgrey', alpha = 1, ha = 'center', va = 'center')
    
            self.plotSpace.axes.set_ylim(bottom = min(minValues[:]), top = max(maxValues[:]))
            self.plotSpace.axes.set_xlim(left = min(minValues[:]), right = max(maxValues[:]))
    
            #plt.xlabel(u'%s Observed (Ma)' %(str(item)))
            #plt.ylabel(u'%s Predicted (Ma)' %(str(item)))
            self.plotSpace.axes.set_xlabel(u'Observed '+ item +' (µm)')
            self.plotSpace.axes.set_ylabel(u'Predicted '+ item +' (µm)') 
                
        elif inputdata == 'Elevation comparison':
                
                for item in dataplot:
                    if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
                        height_obs = datac[agecol['alt']+'OBS'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        height_pred =  datac[agecol['alt']+'PRED'][np.logical_not(datac[agecol[item]+'OBS'] == -9999)]
                        
                        self.plotSpace.axes.plot(height_obs,height_pred, 
                              marker = ageMarker[item], linestyle = 'None',
                              label = item, color = colores[item])
                
                        # Write the 1:1 text over the line
                        # We may need to revise the way to find the x/y
                        item = 'alt'
                        self.plotSpace.axes.text((max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2,
                                (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2 - 
                                    (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/20,
                                '1:1', rotation = 45, color = 'lightgrey', alpha = 1, ha = 'center', va = 'center')
                        self.plotSpace.axes.set_xlabel(u'Observed elevation (m)')
                        self.plotSpace.axes.set_ylabel(u'Predicted elevation (m)') 
                        min_y = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']))
                        min_y = min_y - 0.20*min_y
                        max_y = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED']))
                        max_y = max_y + 0.20*max_y
                        self.plotSpace.axes.set_ylim(bottom = min_y, top = max_y)
                        min_x = min(min(datac[agecol['alt']+'OBS']), min(datac[agecol['alt']+'PRED']))
                        min_x = min_x - 0.20 * min_x
                        max_x = max(max(datac[agecol['alt']+'OBS']), max(datac[agecol['alt']+'PRED']))
                        max_x = max_x + 0.20*max_x
                        self.plotSpace.axes.set_xlim(left = min_x, right = max_x)
                
        elif inputdata == 'Trapped-charge comparison':
            errorpred = float(self.ParametersInput.DParameters['error_predictions']) / 100
            for item in dataplot:
                if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
                    # Get sample name of data > -9999
                    sampleName = inputc["SAMPLE"+item]
                    sampleNamePred = inputc['SAMPLE'+item]
                    # sort data by sample (alphabetically)
                    sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                    sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                    sampleNameIndexes = np.argsort(sampleNametemp)
                    sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                    index = sampleName.coords
                    errortemp = inputc[errname[item]]
                    errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                    errortemp = errortemp2[sampleNamePredIndexes]
                    try:
                        obs = inputc[agecol[item]+'OBS_input']
                        pred = datac[agecol[item]+'PRED']
                        
                        #Remove NaNs
                        pred2 = np.asarray([s for s in pred.values if pgu.isNaN(s) == False])
                        pred = pred2[sampleNamePredIndexes]
                        errorpred_temp = errorpred*pred
                        errorpred_temp = np.asarray([s for s in errorpred_temp if pgu.isNaN(s) == False])
                        obs2 = np.asarray([s for s in obs.values if pgu.isNaN(s) == False])
                        obs = obs2[sampleNamePredIndexes]
                        # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                        # of the array
                        error = np.empty(pred.shape)
                        error[:] = np.nan
                        count = 0

                        for n in range(len(pred)):
                            if not pgu.isNaN(errortemp[n]):
                                error[count] = errortemp[n]
                                count += 1
                        # Plot
                        self.plotSpace.axes.errorbar(obs, pred,
                                  xerr = error, yerr = errorpred*pred,
                                  fmt = ageMarker[item], label = item, color = colores[item])
                        
                        # Compute log likelyhood
                        LL, N, sigma = misfits.compute_misfit(int(self.ParametersInput.DParameters['misfit_function']),obs[:],pred[:],error,errorpred_temp)
                        LL_array.append(LL)
                        sigma_array.append(sigma)
                        N_LL.append(N)
                        Item_list.append(item)
                        
                        ErrorSuccess = 1
                    except Exception as E:
                        print("\n Issue with errors, do not plot")
                        print(str(E))
                        self.plotSpace.axes.plot(datac[agecol[item]+'OBS'], datac[agecol[item]+'PRED'], 
                                      marker = ageMarker[item], linestyle = 'None',
                                     label = item, 
                                     color = colores[item])
                    minValues.append(min(abs(obs)))
                    minValues.append(min(abs(pred)))
                    maxValues.append(max(abs(obs))+shift)
                    maxValues.append(max(abs(pred))+shift)
               
            # Write the 1:1 text over the line
            # We may need to revise the way to find the x/y
            self.plotSpace.axes.text((max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2,
                     (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/2 - 
                          (max(max(datac[agecol[item]+'OBS']), max(datac[agecol[item]+'PRED']))-0)/20,
                     '1:1', rotation = 45, color = 'lightgrey', alpha = 1, ha = 'center', va = 'center')

            self.plotSpace.axes.set_ylim(bottom = min(minValues[:]), top = max(maxValues[:]))
            self.plotSpace.axes.set_xlim(left = min(minValues[:]), right = max(maxValues[:]))

            #plt.xlabel(u'%s Observed (Ma)' %(str(item)))
            #plt.ylabel(u'%s Predicted (Ma)' %(str(item)))
            self.plotSpace.axes.set_xlabel(u'Observed n/N')
            self.plotSpace.axes.set_ylabel(u'Predicted n/N') 

        
        self.plotSpace.axes.legend(loc = 'best')
        # self.plotSpace.axes.title(graphtitle)

        # plt.savefig(graphpath + '/Compare/' + graphtitle +'_Compare_' +str(item) + '.pdf')
        # fig2.clear()
                
        #Compute Log-likelyhood
        if ErrorSuccess:
            try:
                self.plotSpace.fig.suptitle("misfit = " +str([Item_list[v] + ' : '+ str(round(sigma_array[v]/N_LL[v],3)) for v in range(len(sigma_array))]) + "\n" +
                                        "LL = "+str([Item_list[v] + ' : '+ str(round(LL_array[v],3)) for v in range(len(LL_array)) ]),fontsize=8)
            except:
                print("ZeroDivisionError: one or more of the data has no error.")
        #Add toolbar with 2D plot
        self.toolbar = NavigationToolbar(self.plotSpace,self)
        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        self.add_plot2win(Q,self.plotSpace,self.toolbar,inputdata)
    
    #----------------------------------------------------------
    def add_plot2win(self, Q,plotSpace,toolbar,title):
        """ To add plot in the interface """
        
        
        #Check number of plots
        if self.CountPlotsColumn == 2 and self.CountPlotsRow < 1:
            self.CountPlotsRow += 1
            self.CountPlotsColumn = 0
        elif self.CountPlotsRow == 1 and self.CountPlotsColumn == 2:
            self.CountPlotsRow = 0
            self.CountPlotsColumn = 0
        
        if self.MainLay.count() > 4:
            U_interface.MessageBox(text1='The maximum number of plots is reached. Please, remove one of the plot to be able to load another one.')
            return
        if self.CountPlotsColumn == 0 and self.CountPlotsRow == 0:
            self.Splitter1.replaceWidget(0,Q)
            self.Splitter1.widget(0).setVisible(True)
            self.MainLay.refresh()
        elif self.CountPlotsRow == 0 and self.CountPlotsColumn == 1:
            self.Splitter1.replaceWidget(1,Q)
            self.Splitter1.widget(1).setVisible(True)
            self.MainLay.refresh()
        elif self.CountPlotsRow == 1 and self.CountPlotsColumn == 0:
            self.Splitter2.replaceWidget(0,Q) #Add temporary empty widget
            self.Splitter2.widget(0).setVisible(True)
            self.Splitter2.widget(1).setVisible(True)
            self.Splitter2.setVisible(True)
            self.MainLay.refresh()
        elif self.CountPlotsRow == 1 and self.CountPlotsColumn == 1:
            self.Splitter2.replaceWidget(1,Q)
            self.Splitter2.widget(1).setVisible(True)
            self.MainLay.refresh()
        else:
            print('issue for plotting 2D')
            
        # Create items for the tree view
        itemTree = pgu.StandardItem(title+"_"+conf.ProjectName+"("+str(self.CountPlotsColumn)+")", 10)
        self.treeModel.appendRow(itemTree)
        self.treeView.setModel(self.treeModel)
            
        self.SavePlots.append(plotSpace)
        self.SaveToolbar.append(toolbar)
        self.CountPlotsColumn += 1
        
    #----------------------------------------------------------
    def get_ThermoSys(self,filepath):
        #To read the thermochronological systems chosen by the user 
        #to plot an age-elevation profile
        agecol, errname, agename, predname, colores, profdict,He43dict,coloramp,ageMarker  =  Thermo_settings.dict_pecube()
        self.list_sys = []
        for k in sorted(self.System_Thermo.store_systems):
            if self.System_Thermo.store_systems[k] == "HeApatite": 
                self.list_sys.append(agecol['AHE'])
            elif self.System_Thermo.store_systems[k] == "HeZircon": 
                self.list_sys.append(agecol['ZHE'])
            elif self.System_Thermo.store_systems[k] == "FTZircon": 
                self.list_sys.append(agecol['ZFT'])
            elif self.System_Thermo.store_systems[k] == "FTApatite": 
                self.list_sys.append(agecol['AFT'])
            elif self.System_Thermo.store_systems[k] == "ArKFelsdpar": 
                self.list_sys.append(agecol['KAR'])
            elif self.System_Thermo.store_systems[k] == "ArBiotite": 
                self.list_sys.append(agecol['BAR'])
            elif self.System_Thermo.store_systems[k] == "ArMuscovite": 
                self.list_sys.append(agecol['MAR'])
            elif self.System_Thermo.store_systems[k] == "ArHornblend": 
                self.list_sys.append(agecol['HAR'])
            elif self.System_Thermo.store_systems[k] == "ThL" : 
                self.list_sys.append('THL')
            elif self.System_Thermo.store_systems[k] == "OSL": 
                self.list_sys.append('OSL')
            elif self.System_Thermo.store_systems[k] == "ESR": 
                self.list_sys.append('ESR')
              
        # Get dic Pecube
        print(self.list_sys,filepath)
        inputdata = filepath
        self.plot_AgeElevation(inputdata, self.list_sys,
                            agecol = agecol, errname = errname, colores = colores, agename = agename,
                            predname = predname,ageMarker=ageMarker)
        self.wThermo.close()
    
    #----------------------------------------------------------
    def get_AgeElevation(self,filename,headers):
        #To get Age and elevation for all node from Pecube model
        try:
            nfiles = len(filename)
            ages = []
            elevations = []
            for i in range(nfiles):
                df = pd.read_csv(filename[i],usecols=headers)
                AgeElevations = np.asarray(df)
                for j in range(len(headers)-1): #for nb Thermochronometers
                    ages.append(AgeElevations[:,j+1])
                    elevations.append(AgeElevations[:,0])
        except FileNotFoundError:
            QErrorMessage().showMessage("File "+filename+" not found.")

        ages_temp = np.asarray(ages)
        elevations_temp = np.asarray(elevations)
        
        #Remove any Nan column
        ages = ages_temp[:,~np.isnan(ages_temp).any(axis=0)]
        elevations = elevations_temp[:,~np.isnan(elevations_temp).any(axis=0)]
        return ages,elevations
    
    #----------------------------------------------------------
    def get_elevation(self,filename):
        #Get elevations from the CompareAGE.csv file
        try:
            df = pd.read_csv(filename,usecols=['HEIGHTPRED'])
        except FileNotFoundError:
            QErrorMessage().showMessage("File 'CompareAGE.csv' not found in ./output")
            return
        elevations = np.asarray(df)
        #How many grains?
        try:
            nbGrain = int(self.ParametersInput.DParameters['Nb_TotGrains'])
        except KeyError:
            reply = U_interface.MessageBox(text1="No sample-specific location is found. Do you want to load some?")
            if reply.QMessage == QMessageBox.Yes:
                #load file
                options = QFileDialog.Options()
                options |= QFileDialog.DontUseNativeDialog
                self.filePath, _ = QFileDialog.getOpenFileNames(
                    self, "QFileDialog.getOpenFileName()", "", "(*.txt)", options=options)
                if self.filePath:
                    self.FolderPath = os.path.dirname(os.path.dirname(self.filePath[0]))
                    self.folderName = os.path.basename(self.FolderPath)
                    self.readOldInput()
                    filename = os.path.join(self.PecubePath,self.FolderPath,"output","CompareAGE.csv")
                    try:
                        df = pd.read_csv(filename,usecols=['HEIGHTPRED'])
                    except FileNotFoundError:
                        QErrorMessage().showMessage("File 'CompareAGE.csv' not found in ./output")
                        return
                    elevations = np.asarray(df)
                    nbGrain = int(self.ParametersInput.DParameters['Nb_TotGrains'])
                else:
                    return
            else:
                return
        #How many samples?
        nbSamples = len(elevations)
        #How many grain per sample?
        GrainList = self.ParametersInput.TableParameters['nb_grains']
        #Assign elevation to each grain
        elevationArray = np.zeros((nbGrain,1))
        count = 0
        for i in range(nbSamples):
            elevationArray[count] = elevations[i]
            count += 1
        return elevationArray
    
    #----------------------------------------------------------
    def get_ThermoSysSpecific(self):
        #get the thermo system(s) chosen by the user
        #to plot data for specific-sample locations
        outputDir = os.path.join(self.PecubePath,self.folderName,"output")
        #load data
        self.list_sys = [v for k,v in self.System_Thermo.store_systems.items()]
        self.listTimeAgeFiles = []
        for k in self.list_sys:
            if k == 'HeApatite':
                self.listTimeAgeFiles.append('TimeAge.csv')
            elif k == 'FTApatite':
                self.listTimeAgeFiles.append('TimeAgeAFT.csv')
            
    #----------------------------------------------------------
    def get_specific_ages(self):
        #Find all the grain*** files
        outputDir = os.path.join(self.PecubePath,self.folderName,"output")
        fileNames = []
        #From Monte Carlo model or FD ?
        # if int(self.ParametersInput.DParameters['PDModel_signal']): #FD model
        #Is there several thermo systems predicted?
        Ages = []
        self.ThermSys = []
        # Get AHe ages (last column)
        try:
            temp = self.dataset.MainDataset.AHEPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.AHEOBS.values
            #index = np.where(obs < 1e-6)
            #temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("AHE")
        except:
             print("No AHe ages found.")
        # Get AFT ages (last column)
        try:
            temp = self.dataset.MainDataset.AFTPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.AFTOBS.values
            #index = np.where(obs < 1e-6)
            #temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("AFT")
        except:
            print("No AFT ages found.")
        # Get ZHe ages (last column)
        try:
            temp = self.dataset.MainDataset.ZHEPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.ZHEOBS.values
            # index = np.where(obs < 1e-6)
            # temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("ZHE")
        except:
            print("No ZHe ages found.")
        # Get ZFT ages 
        try:
            temp = self.dataset.MainDataset.ZFTPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.ZFTOBS.values
            # index = np.where(obs < 1e-6)
            # temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("ZFT")
        except:
            print("No ZFT ages found.")
        # Get KAr ages 
        try:
            temp = self.dataset.MainDataset.KARPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.KAROBS.values
            # index = np.where(obs < 1e-6)
            # temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("KAR")
        except:
            print("No KAR ages found.")
        # Get BAr ages 
        try:
            temp = self.dataset.MainDataset.BARPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.BAROBS.values
            # index = np.where(obs < 1e-6)
            # temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("BAR")
        except:
            print("No BAR ages found.")
        # Get MAr ages 
        try:
            temp = self.dataset.MainDataset.MARPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.MAROBS.values
            # index = np.where(obs < 1e-6)
            # temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("MAR")
        except:
            print("No MAR ages found.")
        # Get HAr ages 
        try:
            temp = self.dataset.MainDataset.HARPRED.values
            # Get prediction where there is observations
            obs = self.dataset.MainDataset.HAROBS.values
            # index = np.where(obs < 1e-6)
            # temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Ages.append(temp)
                self.ThermSys.append("HAR")
        except:
            print("No HAR ages found.")
        self.Ages = np.asarray(Ages)
        self.AgeError_array = self.Ages*0 # 10% error on predited ages
        return self.Ages, self.AgeError_array 
        # else:   #Monte Carlo model
        #     for file in os.listdir(outputDir):
        #         if file.startswith("Grain"):
        #             fileNames.append(os.path.join(outputDir,file))
        #     #Then read through files to find the age and its uncertainty
        #     self.Ages = np.zeros(len(fileNames))
        #     self.AgeError_array = np.zeros(len(fileNames))
        #     for i in range(len(fileNames)):
        #         try:
        #             file = pd.read_csv(fileNames[i],header=None,
        #                                 squeeze = True,encoding=conf.encoding_label)
        #         except FileNotFoundError:
        #             win = QMessageBox()
        #             win.setText("File not found when trying to read ages.\n"+
        #                           "The process has been stopped.")
        #             win.setIcon(QMessageBox.Information)
        #             win.setStandardButtons(QMessageBox.Ok)
        #             retval = win.exec_()
        #             if retval == QMessageBox.Ok:
        #                 return
        #             else:
        #                 return
                
        #         #Extract some data
        #         index = file.loc[file.str.startswith('age=')].index[0]
        #         data = pd.DataFrame(pd.read_csv(fileNames[i],sep='\  ', nrows=1,
        #                                    skiprows=index,skipinitialspace=True,header=None,encoding=conf.encoding_label)).values
                
        #         #Get age and uncertainty
        #         AgeData = data[0][1].split("+-")
        #         AgeValue = AgeData[0].split("=")
        #         self.Ages[i] = float(AgeData[0])
        #         self.AgeError_array[i] = float(AgeData[1])
                
        #         #Ft Factor 
        #         index = file.loc[file.str.startswith('FTg=')].index[0]
        #         data = pd.DataFrame(pd.read_csv(fileNames[i],sep='\ ', nrows=1,
        #                                    skiprows=index,skipinitialspace=True,header=None,encoding=conf.encoding_label)).values
        #         FTData = data[0][1].split("+-")
        #         FTValue = FTData[0].split("=")
        #         self.Ages[i] = self.Ages[i]/float(FTData[0])
                
        #         return self.Ages, self.AgeError_array
    
    #----------------------------------------------------------
    def get_Trapped_charge(self):
        #Find all the files
        outputDir = os.path.join(self.PecubePath,self.folderName,"output")
        fileNames = []
        Obs = []
        Pred = []
        self.TrapChargSys = []
        # Get ThL
        try:
            TL_output = self.dataset.MainDataset.filter_by_attrs(description="ThL_data")
            pred = TL_output.NNPRED.values
            # Get prediction where there is observations
            obs = TL_output.NNOBS.values
            #index = np.where(obs < 1e-6)
            #temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Obs.append(obs)
                Pred.append(pred)
                self.TrapChargSys.append("ThL")
        except:
             print("No ThL data found.")
        # Get OSL
        try:
            output = self.dataset.MainDataset.filter_by_attrs(description="OSL_data")
            pred = output.NNPRED.values
            # Get prediction where there is observations
            obs = output.NNOBS.values
            #index = np.where(obs < 1e-6)
            #temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Obs.append(obs)
                Pred.append(pred)
                self.TrapChargSys.append("OSL")
        except:
             print("No OSL data found.")
        # Get ESR
        try:
            output = self.dataset.MainDataset.filter_by_attrs(description="ESR_data")
            pred = output.NNPRED.values
            # Get prediction where there is observations
            obs = output.NNOBS.values
            #index = np.where(obs < 1e-6)
            #temp[index[0]] = nan
            if np.any(temp > 1e-6):
                Obs.append(obs)
                Pred.append(pred)
                self.TrapChargSys.append("ESR")
        except:
             print("No ESR data found.")
        
        self.Obs = np.asarray(Obs)
        self.Pred = Pred # 10% error on predited ages
        
        return self.Obs, self.Pred
    
    #----------------------------------------------------------
    def getHeProfile(self,label):
        
        if label.currentText() == "4He/3He spectrum":
            self.He4He3 = 1
        elif label.currentText() == "set profile":
            self.He4He3 = 0
            
    #----------------------------------------------------------
    def getFitFunction(self,label):
        #To set which function to use to fit the 4He/3He data
        if label.currentText() == "quadratic":
                self.FitFunct = 1
        elif label.currentText() == "power law":
                self.FitFunct = 2
        elif label.currentText() == "poly3":
                self.FitFunct = 3
        elif label.currentText() == "poly4":
                self.FitFunct = 4
        elif label.currentText() == "poly5":
                self.FitFunct = 5
            
    #----------------------------------------------------------
    def loadData(self,click):
        """ 
         To load data
         click=1, The user clicked on "Load results", Select file to open
         else load from known file directly.
         
         """
        if click:
            dialog = QFileDialog()
            dialog.setDirectory(os.path.join(self.PecubePath))
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            self.filePath, _ = dialog.getOpenFileNames(
                self, "QFileDialog.getOpenFileName()", "", "(*.in)", options=options)
        if self.filePath:
            # How many files?
            nfiles = len(self.filePath)
            
            # extract file name
            extracted_Fnames = self.filePath[0].split("/")
            str1 = ''
            strlist = extracted_Fnames[-2:]
            str2 = str1.join(strlist)
            Nb_char = len(str2)+2  # count number of character after folder name
            self.fileName = os.path.basename(self.filePath[0])
            self.FolderPath = os.path.dirname(os.path.dirname(self.filePath[0]))
            self.folderName = os.path.basename(self.FolderPath)
            conf.ProjectName = self.folderName
            
            # Extract file type
            file_split = os.path.splitext(self.filePath[0])
            self.file_extension = file_split[1]

            #Update data Combo box to plot Age-elevation profile
            print(self.CurrentFolder,self.folderName)
            self.CurrentFolder = self.folderName #store the current project
            if click == 1:
                self.get_dataList() 
            
                
            #----- Extract data from file(s) -------
            if nfiles > 1: # Several files have been chosen, to summary age predictions
                if self.file_extension == '.out':
                    print("wait")
                else:
                    win = U_interface.MessageBox(text1='PecubeGUI cannot handle this serie of files. The process is cancelled.')
                    if win.QMessage == QMessageBox.Ok or win.Qmessage == QMessageBox.Cancel or win.Qmessage == QMessageBox.No:
                        return
            else:
                #Extract according to file extension
                if self.file_extension == '.vtk': # --------------------------------------------------
                    #Read data from vtk file
                    # Load the vtk file and read the data array names
                    self.file_vtk = pv.read(self.filePath[0])
                    if self.DataCombo.currentText() == "Map ages":
                        self.data_names = self.file_vtk.array_names
                        self.data_names.append("4He/3He res")
                    else:
                        self.data_names = self.file_vtk.array_names
                    # 2) ask the user the data he/she wishes to plot 
                    self.wVTK = U_interface.WinExtraParameters(500, 500, text='Pecube outputs',fixed=True)
                    self.dataVTK_name = self.data_names[0] #by default the first set of data
                    self.parametersVTK = plot_2D.plot2DVTK(self,self.data_names)
                    self.wVTK.layout.addLayout(self.parametersVTK.layout, 0, 0)   
                    Q = QWidget()
                    Q.setLayout(self.wVTK.layout)
                    self.wVTK.setCentralWidget(Q)
                    conf.WindowsOpen.append(self.wVTK)
                    self.wVTK.show()
                    # 3) Once the user chose, extract the data from the file
                    self.wVTK.OkButton.clicked.connect(lambda: self.plotData()) 
                    
                elif self.file_extension == '.csv': # --------------------------------------------------
                    if self.fileName == 'CoolingRates.csv':
                        dataset = pd.read_csv(self.filePath[0],header=None)
                        nb_headers = len(dataset.columns)
                        self.HeadersLabel = list(str(dataset.columns.values))
                        self.dataFrame = pd.DataFrame(dataset)
                        # Load parameters within a window
                        self.CoolingRateObj = plot_2D.CoolingRates_Parameters(self)
                        self.TableWindow = U_interface.WinExtraParameters(500, 300,'Create 2D map of cooling rates...')
                        self.TableWindow.layout.addLayout(self.CoolingRateObj.vBox,0,0)
                        Q = QWidget()
                        Q.setLayout(self.TableWindow.layout)
                        self.TableWindow.setCentralWidget(Q)
                        self.TableWindow.show()
                        #Get the parameter values from the user and load the map
                        self.TableWindow.OkButton.clicked.connect(self.plotData)
                        #load data to table
                        self.tableData = pgu.TableRead(self.dataFrame,self.HeadersLabel)
                        
                    else:
                        try:
                            dataset = pd.read_csv(self.filePath[0])
                        except FileNotFoundError as E:
                            QErrorMessage(self).showMessage("CSV file not found. Please check if it exists: <br>" +
                                                            str(E))
                        nb_headers = len(dataset.columns)
                        self.dataFrame = pd.DataFrame(dataset)
                        #Drop columns with NaN values
                        self.dataFrame = self.dataFrame.dropna(axis=1)
                        self.HeadersLabel = list(self.dataFrame.columns.values)
                        #load data to table to show up
                        self.tableData = pgu.TableRead(self.dataFrame,self.HeadersLabel)
                        #open window to show the data
                        self.TableWindow = U_interface.WinExtraParameters(1000, 800,'Load data...')
                        #Help sentence 
                        Lay = QVBoxLayout()
                        if self.fileName == 'Compare43He.csv': 
                            self.He43_object = out.Data_43He(self.tableData.Table)
                            Lay.addWidget(self.He43_object.Splitter)
                        else:
                            HelpLabel = QLabel("Select a range to plot, or press 'Ok' to plot all data.")
                            HelpLabel.setFont(conf.font12)
                            Lay.addWidget(HelpLabel)
                            Lay.addWidget(self.tableData.Table)
                        self.TableWindow.layout.addLayout(Lay, 0, 0, 1, 3)
                        Q = QWidget()
                        Q.setLayout(self.TableWindow.layout)
                        self.TableWindow.setCentralWidget(Q)
                        conf.WindowsOpen.append(self.TableWindow)
                        self.TableWindow.show()
                        #Get the selected cells to plot
                        self.TableWindow.OkButton.clicked.connect(self.plotData)
                
            #     elif self.file_extension == '.out': # --------------------------------------------------
            #         #### Files from production-diffusion model #####
            #         self.He4He3 = 0
            #         # Read files
            #         file = pd.read_csv(self.filePath[0],header=None,
            #                            squeeze = True,encoding=conf.encoding_label)
        
            #         # Locate the time-Age table evolution
            #         index = file.loc[file.str.startswith('time (My)')].index[0]
            #         #Get Table, headers are : time (My), age (My), error (My), T (°C), Fract
            #         TimeAgeTable = pd.read_csv(self.filePath[0],sep='\t',
            #                                    skiprows=index+1,header=None,encoding=conf.encoding_label)
            #         self.HeadersLabel = ['time (My)', 'age (My)', 'error (My)', 'T (°C)', 'Fract']
            #         self.dataFrame = pd.DataFrame(TimeAgeTable)
            #         self.dataFrame[0] = abs(self.dataFrame[0]) #Time is positive
            #         #load data to table to show up
            #         self.tableData = TableRead(self.dataFrame,self.HeadersLabel)
            #         #open window to show the data
            #         self.TableWindow = WinExtraParameters(800, 1000,'Load data...')
            #         #Extract some data
            #         index = file.loc[file.str.startswith('Tc=')].index[0]
            #         data = pd.DataFrame(pd.read_csv(self.filePath[0], nrows=1,
            #                                    skiprows=index,header=None,encoding=conf.encoding_label)).values
            #         TcRead = data[0][0].split('=')
            #         Tc = TcRead[1]
            #         index = file.loc[file.str.startswith('tc=')].index[0]
            #         data = pd.DataFrame(pd.read_csv(self.filePath[0], nrows=1,
            #                                    skiprows=index,header=None,encoding=conf.encoding_label)).values
            #         tcRead = data[0][0].split('=')
            #         tc = tcRead[1]
            #         index = file.loc[file.str.startswith('FTg=')].index[0]
            #         data = pd.DataFrame(pd.read_csv(self.filePath[0], nrows=1,
            #                                    skiprows=index,header=None,encoding=encoding_label)).values
            #         FTgRead = data[0][0].split('=')
            #         FTg = FTgRead[1]
                    
            #         #Is there any 4He/3He data?
            #         try:
            #             index = file.loc[file.str.startswith('He4_fract')].index[0]
            #             He_fract = pd.DataFrame(pd.read_csv(self.filePath[0],sep='\  ', nrows=2,
            #                                        skiprows=index,header=None,encoding=conf.encoding_label)).values
            #             He4_step = He_fract[0,1:]
            #             He3_step = He_fract[1,1:]
            #             #Handle division by zero
            #             index3 = np.where(He3_step==0)
            #             index4 = np.where(He4_step==0)
            #             if len(index3[0]) < len(index4[0]):
            #                 index3 = index4
            #             elif len(index3[0]) > len(index4[0]):
            #                 index4 = index3
            #             He3_step = np.delete(He3_step,index3)
            #             He4_step = np.delete(He4_step,index4)
            #             #4He/3He degassing spectrum
            #             try:
            #                 He_fract = He4_step/He3_step 
            #             except ZeroDivisionError:
            #                 index3 = np.where(He3_step==0)
            #                 index4 = np.where(He4_step==0)
            #                 He3_step[index3] = 1
            #                 He4_step[index4] = 1
            #                 He_fract = He4_step/He3_step 
            #             bulk = sum(He_fract*He3_step)/sum(He3_step)
            #             self.He_ratio = He_fract/bulk
            #             self.cum3He = np.asarray(np.cumsum(He3_step)/sum(He3_step))
            #             HeProfile = 1
            #         except IndexError: #No 4He/3He data
            #             HeProfile = 0
                    
            #         #Set the window 
            #         Label = QLabel("Model results:")
            #         Label.setFont(fontBold12)
            #         CorrectedAgeLabel = QLabel("Corr. Age (Ma)")
            #         CorrectedAgeLabel.setFont(font10)
            #         CorrectedAgeLabel.setAlignment(Qt.AlignCenter)
            #         FTLabel = QLabel('FT factor')
            #         FTLabel.setFont(font10)
            #         FTLabel.setAlignment(Qt.AlignCenter)
            #         ClosureTempLabel = QLabel('Closure temp. (°C)')
            #         ClosureTempLabel.setFont(font10)
            #         ClosureTempLabel.setAlignment(Qt.AlignCenter)
            #         CorrectedAge = QLineEdit(tc)
            #         CorrectedAge.setEnabled(False)
            #         CorrectedAge.setAlignment(Qt.AlignCenter)
            #         FT = QLineEdit(FTg)
            #         FT.setEnabled(False)
            #         FT.setAlignment(Qt.AlignCenter)
            #         ClosureTemp = QLineEdit(Tc)
            #         ClosureTemp.setEnabled(False)
            #         ClosureTemp.setAlignment(Qt.AlignCenter)
            #         TitleLabel = QLabel("What data to plot?")
            #         TitleLabel.setFont(fontBold12)
            #         HeliumLabel = QLabel("4He/3He data:")
            #         HeliumLabel.setAlignment(Qt.AlignRight)
            #         HeliumLabel.setFont(font10)
            #         HeliumData = QComboBox()
            #         HeliumData.setMinimumContentsLength(15)
            #         items = ['set profile','4He/3He spectrum']
            #         HeliumData.addItems(items)
            #         HeliumData.currentIndexChanged.connect(lambda: self.getHeProfile(HeliumData))
            #         FitCurveLabel = QLabel('Fitting function:')
            #         FitCurveLabel.setAlignment(Qt.AlignRight)
            #         FitCurveLabel.setFont(font10)
            #         FitCurve = QComboBox()
            #         FitCurve.setMinimumContentsLength(15)
            #         items = ['quadratic','power law','poly3','poly4','poly5']
            #         self.FitFunct = 1
            #         FitCurve.addItems(items)
            #         FitCurve.currentIndexChanged.connect(lambda: self.getFitFunction(FitCurve))
            #         if HeProfile:
            #             HeliumData.setEnabled(True)
            #         else:
            #             HeliumData.setEnabled(False)
            #         HelpLabel = QLabel("Select from the table:")
            #         HelpLabel.setFont(font11)
            #         Box1 = QGroupBox()
            #         Lay1 = QGridLayout()
            #         Lay1.setSpacing(10)
            #         Lay1.addWidget(ClosureTempLabel,1,1)
            #         Lay1.addWidget(CorrectedAgeLabel,1,2)
            #         Lay1.addWidget(FTLabel,1,3)
            #         Lay1.addWidget(ClosureTemp,2,1)
            #         Lay1.addWidget(CorrectedAge,2,2)
            #         Lay1.addWidget(FT,2,3)
            #         Lay1.setColumnStretch(0,1)
            #         Lay1.setColumnStretch(4,1)
            #         Box1.setLayout(Lay1)
            #         Lay = QVBoxLayout()
            #         Lay.setSpacing(10)
            #         Lay.addWidget(Label)
            #         Lay.addWidget(Box1)
            #         Lay.addSpacing(20)
            #         Separator = QFrame()
            #         Separator.setFrameShape(QFrame.HLine)
            #         Separator.setFrameShadow(QFrame.Raised)
            #         Separator.setLineWidth(1)
            #         Separator.setMidLineWidth(1)
            #         Lay.addWidget(Separator)
            #         Lay.addWidget(TitleLabel)
            #         Lay.addSpacing(10)
            #         HLay = QHBoxLayout()
            #         HLay.addWidget(HeliumLabel)
            #         HLay.addWidget(HeliumData)
            #         HLay.addWidget(FitCurveLabel)
            #         HLay.addWidget(FitCurve)
            #         HLay.addStretch(1)
            #         HLay.setSpacing(20)
            #         Lay.addLayout(HLay)
            #         Lay.addSpacing(20)
            #         Lay.addWidget(HelpLabel)
            #         Lay.addWidget(self.tableData.Table)
            #         Lay.setSpacing(5)
            #         self.TableWindow.layout.addLayout(Lay, 0, 0)
            #         Q = QWidget()
            #         Q.setLayout(self.TableWindow.layout)
            #         self.TableWindow.setCentralWidget(Q)
            #         conf.WindowsOpen.append(self.TableWindow)
            #         self.TableWindow.show()
            #         #Get the selected cells to plot
            #         self.TableWindow.OkButton.clicked.connect(self.plotData)
                else:
                    return
          
    #----------------------------------------------------------
    def plotData(self):
        #Read selected cells
        try:
            nbSamples = int(self.ParametersInput.DParameters['Nb_Samples'])
            nbGrain_per_sample = self.ParametersInput.TableParameters['nb_grains']
            nbTotGrain = int(self.ParametersInput.DParameters['Nb_TotGrains'])
            #Get name of grains
            SampleName = [v  for k,v in self.ParametersInput.TableParameters.items() if 'SampleName' in k]
            if len(SampleName) > 0:
                textlegend = [SampleName[i]+'_'+str(j+1) for i in range(len(nbGrain_per_sample)) for j in range(int(nbGrain_per_sample[i][0]))]
            else:
                textlegend = ['Sample'+str(i+1)+'_1' for i in range(nbSamples)]
            #Get Label for each object on plots (to handle update in plots due to change in axes properties))
            numLabel = []
            nbdata = 1000
            for i in range(nbdata):
                if i < 10:
                    numLabel.append('000'+str(i))
                elif i >= 10 and i < 100:
                    numLabel.append('00'+str(i))
                elif i >= 100 and i < 1000:
                    numLabel.append('0'+str(i))
        except Exception as E:
            print("No samples data: \n" +str(E))
        countLabel = 0
        if self.file_extension == '.csv' or self.file_extension == '.out': # ----------------------------------------------
            dataSelected_temp = [float(item.text()) for item in self.tableData.Table.selectedItems()]
            dataIndex = np.unique([int(item.column()) for item in self.tableData.Table.selectedItems()])
            dataIndex = dataIndex[1:]-1 #The first index is time, remove it
            if np.array(dataSelected_temp).any():
                #Get number of rows and Column
                nb_Column = max(set(index.column() for index in self.tableData.Table.selectedIndexes()))+1
                nb_Rows = len(set(index.row() for index in self.tableData.Table.selectedIndexes()))
                #Create array to store selected data
                dataSelected = np.zeros((nb_Rows,nb_Column))
                indexColumn = [int(item.column()) for item in self.tableData.Table.selectedIndexes()]
                indexRow = [int(item.row()) for item in self.tableData.Table.selectedIndexes()]
                #Fill array according to indexes of items
                for k in range(len(indexColumn)):
                    indexR = indexRow[k]
                    indexC = indexColumn[k]
                    dataSelected[indexR][indexC] = dataSelected_temp[k]
                #Remove empty Column
                FilledCol_temp = np.all((abs(dataSelected) == 0.0), axis=0)
                FilledCol = []
                for i in range(len(FilledCol_temp)):
                    if FilledCol_temp[i] == False:
                        FilledCol.append(i)
                dataSelected = dataSelected[:,FilledCol]
                dataSelected = pd.DataFrame(dataSelected)
                # Get headers of selected Columns
                index_headers = np.asarray(np.unique(indexColumn), dtype = int)
                index_headers = np.delete(index_headers, 0) # Handle time column for PTT
            else:
                #assume want to plot all data
                dataSelected = self.dataFrame
                # Get headers of all Columns
                index_headers = range(np.size(self.HeadersLabel))
                dataIndex = np.transpose([i for i in range(len(index_headers)-1)])
                
            #### Then plot the data ####
            data_cleaned = np.array(dataSelected.values)
        
            if self.fileName == 'Compare43He.csv':
                cmap = matplotlib.cm.get_cmap('hsv')
                colors = [cmap(i*1/nbSamples) for i in range(len(nbGrain_per_sample)) for j in range(int(nbGrain_per_sample[i][0]))]
                Headers = ['\u03A3F3He', 'Rstep/Rbulk']
                nb_grains = int(len(data_cleaned[1, :])/2)
                textlegend = ['Grain0'+str(i+1) if i+1 <10 else 'Grain' +str(i+1)
                              for i in range(nb_grains)]
                self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
                self.toolbar = NavigationToolbar(self.plotSpace,self)
                #Get ages
                Ages, AgeError = self.get_specific_ages()
                indexAges = [Even for Even in dataIndex if np.mod(Even,2) == 0]
                for i in range(nb_grains):
                    #check for zeros at the end of the array
                    xdata = data_cleaned[:, i*2+1]
                    xdata = pgu.removeZeros(xdata,len(xdata))
                    ydata = data_cleaned[:, i*2]
                    ydata = pgu.removeZeros(ydata,len(ydata))
                    try:
                        if len(xdata) < len(ydata):
                            xdata =data_cleaned[:,i*2+1]
                            xdata = xdata[0:len(ydata)]
                        elif len(ydata) < len(xdata):
                            ydata =data_cleaned[:,i*2]
                            ydata = ydata[0:len(xdata)]
                    except TypeError:
                        QErrorMessage().showMessage("An error ocurred while reading the data from file. Please, check the data.")
                    #4He/3He spectra or step age ?
                    if self.He43_object.plotStepAge:
                        #Did we already produce the required file ?
                        StepAgeFile = os.path.join(self.PecubePath,self.folderName,"output","43He_stepAge.csv")
                        if os.path.exists(StepAgeFile):
                            #Read into file 
                            dataset = pd.read_csv(StepAgeFile)
                            nb_headers = len(dataset.columns)
                            HeadersLabel = list(dataset.columns.values)
                            dataFrame = pd.DataFrame(dataset)
                            #Get corresponding data from selected range
                            dataStepAge = dataFrame.iloc[:,dataIndex]
                            #Extract values as 4He/3He data
                            index = len(ydata)
                            Ratio43_StepAge = dataStepAge.iloc[:,i*2].values
                            Ratio43_StepAge = Ratio43_StepAge[0:index]
                            #Get selected grain(s) for ages
                            headers = dataset.columns.values[indexAges[i]]
                            grainID = headers[-4:]
                            #remove zeros
                            index = [int(value) for value in range(len(grainID)) if int(grainID[value]) > 0]
                            grainID = [grainID[int(index[0]):]]
                            #Compute Step age
                            StepAge = ydata / Ratio43_StepAge * Ages[0][int(grainID[0])-1]
                            #plot spectra
                            self.plotSpace.axes.plot(xdata,StepAge,'-',color=colors[i],label=numLabel[countLabel]+'43HeSA'+textlegend[i])
                            countLabel += 1
                            self.plotSpace.axes.set_ylabel('Step age (Ma)')
                            self.plotSpace.axes.set_title('Age spectrum')
                            # Create item for the tree view
                            itemTree = pgu.StandardItem('Age spectrum'+"_"+conf.ProjectName, 10)
                    else:
                        #Any observed data?
                        #plot spectra
                        self.plotSpace.axes.plot(xdata,ydata,'-',color=colors[i],label=numLabel[countLabel]+'43He'+textlegend[i])
                        countLabel += 1
                        if self.He43_object.data2Plot:
                            xobs = self.He43_object.data2Plot['SF3He']
                            yobs = self.He43_object.data2Plot['Rstep/Rbulk']
                            obsError = self.He43_object.data2Plot['Error']
                            ###### Goodness of fit #######
                            #Interpolate pred data at obs data location
                            ypred = []
                            for k in range(len(xobs)):
                                ilower = [v for v in range(len(xdata)) if xdata[v] <= xobs[k]]
                                ilower = ilower[-1]
                                ihigher = [v for v in range(len(xdata)) if xdata[v] >= xobs[k]]
                                ihigher = ihigher[0]
                                #linear interpolation
                                ypred.append(ydata[ilower]+(xdata[ihigher]-xdata[ilower])*((ydata[ihigher]-ydata[ilower])/(xdata[ihigher]-xdata[ilower])))
                            ypred[-1] = ydata[ihigher]
                            #Calculate log-likelihood
                            N = len(xobs)
                            ypred = np.asarray(ypred)
                            yobs = np.asarray(yobs)
                            obsError = np.asarray(obsError)
                            LL = -(1/(2*N)) * sum(((yobs-ypred)/obsError)**2)
                            self.plotSpace.fig.suptitle("LL = "+str(round(LL,3)))
                            for p in range(len(obsError)):
                                self.plotSpace.axes.plot([xobs[p], xobs[p]],[yobs[p]-obsError[p], yobs[p]+obsError[p]],'k-',label=numLabel[countLabel]+'43He'+textlegend[i])
                                countLabel += 1
                            self.plotSpace.axes.plot(xobs,yobs,'r.')
                        self.plotSpace.axes.set_ylabel(Headers[1])
                        self.plotSpace.axes.set_ylim([0.0,1.5])
                        self.plotSpace.axes.set_title(self.fileName)
                        # Create item for the tree view
                        itemTree = pgu.StandardItem(str(self.fileName[0:-4])+"_"+conf.ProjectName, 10)
                self.plotSpace.axes.set_xlabel(Headers[0])
                self.plotSpace.axes.set_xlim([0.0,1.0])
                self.plotSpace.axes.legend(textlegend)
                # self.plotSpace.axes.set_aspect(1.0/self.plotSpace.axes.get_data_ratio(), adjustable='box')
                try:
                    self.treeModel.appendRow(itemTree)
                except:
                    QErrorMessage().showMessage("Something wrong happen with the plot. Please check the data selected.")
                    return
                self.treeView.setModel(self.treeModel)
                self.TableWindow.close()
                
            elif self.fileName == 'CoolingRates.csv':    
                Headers = ['Y distance', 'X distance']
                # Get time array, is the first column of dataset
                timeArray = data_cleaned[:,0]
                # Get min and max temp defined bu the user
                minTemp = float(self.CoolingRateObj.minTemp.text())
                maxTemp = float(self.CoolingRateObj.maxTemp.text())
                
                #Get some parameters
                nx = int(self.MainWin.ParamTable.ParamTable.nxEdit2.text())
                ny = int(self.MainWin.ParamTable.ParamTable.nyEdit3.text())
                # Get nskip: topographic points skipped by the user
                nskip = int(self.MainWin.ParamTable.ParamTable.nskipEdit8.text())
                nx = int((nx - 1) /nskip + 1)
                ny = int((ny - 1) /nskip + 1)

                # Now search for each Tt path, the location of temperatures
                # and the associated time.
                # interpolate according to temperature or time limits
                if self.CoolingRateObj.data == 'temperature':
                    # Get min and max temperature in the dataset
                    maxTempObs = min(data_cleaned[0,1:])
                    minTempObs = max(data_cleaned[-1,1:])
                    reply = None
                    if maxTempObs < maxTemp:
                        reply = U_interface.MessageBox(text1="The maximum temperature provided is out of Tt paths observed.\n"+
                                   "Please use min < "+str(maxTempObs)+ " °C to continue.")
                    elif  minTempObs > minTemp:
                        reply = U_interface.MessageBox(text1="The minimum temperature provided is out of Tt paths observed.\n"+
                                   "Please use max > "+str(minTempObs)+ " °C to continue.")
                    if reply:
                        if reply.QMessage == QMessageBox.Yes:
                            return
                        else:
                            self.TableWindow.close()
                            return
                    # We start by the lower temperature, we search for the first temperature that is below minTemp
                    indexMinTemp_temp = [[i for i in range(len(data_cleaned[:,0])) if data_cleaned[i,j+1]<=minTemp] for j in range(len(data_cleaned[0,1:nx*ny+1]))]
                    indexMinTem = [i[0] if len(i) > 1 else i for i in indexMinTemp_temp ]
                    #Double check
                    indexMinTemp = [i[0] if isinstance(i,list) else i for i in indexMinTem]
                    
                    # Now the higher temperature, we search for the first temperature that is above minTemp
                    indexMaxTemp_temp = [[i for i in range(len(data_cleaned[:,0])) if data_cleaned[i,j+1]>=maxTemp] for j in range(len(data_cleaned[0,1:nx*ny+1]))]
                    indexMaxTem = [i[-1] if len(i) > 1 else i for i in indexMaxTemp_temp ]
                    #Double check
                    indexMaxTemp = [i[-1] if isinstance(i,list) else i for i in indexMaxTem]
                
                    # Find the age limits
                    coef = [(data_cleaned[indexMinTemp[i]-1,0]-data_cleaned[indexMinTemp[i],0])/
                            (data_cleaned[indexMinTemp[i]-1,i+1]-data_cleaned[indexMinTemp[i],i+1]) for i in range(len(indexMinTemp))]
                    # Calculate age min
                    minAge = np.asarray([data_cleaned[indexMinTemp[i]-1,0]-coef[i]*(data_cleaned[indexMinTemp[i]-1,i+1]-minTemp)
                              for i in range(len(indexMinTemp))])
                    # for upper age
                    coef = [(data_cleaned[indexMaxTemp[i],0]-data_cleaned[indexMaxTemp[i]+1,0])/
                            (data_cleaned[indexMaxTemp[i],i+1]-data_cleaned[indexMaxTemp[i]+1,i+1]) for i in range(len(indexMaxTemp))]
                    # Calculate age min
                    maxAge = np.asarray([data_cleaned[indexMaxTemp[i],0]-coef[i]*(data_cleaned[indexMaxTemp[i],i+1]-maxTemp)
                              for i in range(len(indexMaxTemp))])
                    # Calculate cooling rates
                    CoolingRates = np.asarray([maxTemp-minTemp]*len(indexMaxTemp))/(maxAge-minAge)
                    
                else: # Mean the user provided time range
                    # Get min and max time in the dataset
                    maxTempObs = max(data_cleaned[:,0])
                    minTempObs = min(data_cleaned[:,0])
                    reply = None
                    if maxTempObs < maxTemp:
                        reply = U_interface.MessageBox(text1="The maximum time provided is out of Tt paths observed.\n"+
                                   "Please use max < "+str(maxTempObs)+ " Myr to continue.")
                    elif  minTempObs > minTemp:
                        reply = U_interface.MessageBox(text1="The minimum time provided is out of Tt paths observed.\n"+
                                   "Please use min > "+str(minTempObs)+ " Myr to continue.")
                    if reply:
                        if reply.Qmessage == QMessageBox.Yes:
                            return
                        else:
                            self.TableWindow.close()
                            return
                    # We start by the youngest time, we search for the first time that is below minTemp
                    indexMinTemp = [i for i in range(len(data_cleaned[:,0])) if data_cleaned[i,0]<=minTemp]
                    indexMinTemp = indexMinTemp[0]
                    # Now the oldest time, we search for the first time that is above minTemp
                    indexMaxTemp = [i for i in range(len(data_cleaned[:,0])) if data_cleaned[i,0]>=maxTemp]
                    indexMaxTemp = indexMaxTemp[-1]
                    
                    # Find the age limits
                    coef = [(data_cleaned[indexMinTemp-1,j+1]-data_cleaned[indexMinTemp,j+1])/
                            (data_cleaned[indexMinTemp-1,0]-data_cleaned[indexMinTemp,0]) for j in range(len(data_cleaned[0,1:nx*ny+1]))]

                    # Calculate temp min
                    minAge = np.asarray([data_cleaned[indexMinTemp-1,i+1]-coef[i]*(data_cleaned[indexMinTemp-1,0]-minTemp)
                              for i in range(len(data_cleaned[0,1:nx*ny+1]))])
                    # for upper temp
                    coef = [(data_cleaned[indexMaxTemp,j+1]-data_cleaned[indexMaxTemp+1,j+1])/
                            (data_cleaned[indexMaxTemp,0]-data_cleaned[indexMaxTemp+1,0]) for j in range(len(data_cleaned[0,1:nx*ny+1]))]
                    # Calculate age min
                    maxAge = np.asarray([data_cleaned[indexMaxTemp,i+1]-coef[i]*(data_cleaned[indexMaxTemp,0]-maxTemp)
                              for i in range(len(data_cleaned[0,1:nx*ny+1]))])
                    # Calculate cooling rates
                    CoolingRates = (maxAge-minAge)/np.asarray([maxTemp-minTemp]*len(coef)) if maxTemp - minTemp != 0 else 0
                        
                # Create 2D map
                dx = float(self.MainWin.ParamTable.ParamTable.dLonEdit6.text())
                dy = float(self.MainWin.ParamTable.ParamTable.dLatEdit7.text())
                x0 = float(self.MainWin.ParamTable.ParamTable.lon0Edit4.text())
                y0 = float(self.MainWin.ParamTable.ParamTable.lat0Edit5.text())
                xl = (nx-1)*dx*111.11*nskip*cos((y0+dy*ny/2)*3.141592654/180)
                yl = (ny-1)*dy*111.11*nskip
                x = np.linspace(x0,x0+xl,nx)
                y = np.linspace(y0,y0+yl,ny)
                X,Y = np.meshgrid(y,x)
                # Reshape data accordingly
                Z = np.transpose(np.reshape(CoolingRates,(ny,nx)))
                reversed_map = plt.cm.get_cmap('coolwarm')
                self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
                self.toolbar = NavigationToolbar(self.plotSpace,self)
                self.plotSpace.axes.set_title('Cooling rates (°C/Myr) Range: '+str(minTemp)+'-'+str(maxTemp))
                self.plotSpace.axes.set_xlabel(Headers[0])
                self.plotSpace.axes.set_ylabel(Headers[1])

                plot2D = self.plotSpace.axes.pcolormesh(Y,X,Z,cmap=reversed_map)
                cb = self.plotSpace.fig.colorbar(plot2D, ax=self.plotSpace.axes, shrink=0.4, aspect=7)
                # self.plotSpace.axes.set_aspect(1.0/self.plotSpace.axes.get_data_ratio(), adjustable='box')
                # Create item for the tree view
                itemTree = pgu.StandardItem(str(self.fileName[0:-4])+"_"+conf.ProjectName, 10)
                self.treeModel.appendRow(itemTree)
                self.treeView.setModel(self.treeModel)
                #Close the Window
                self.TableWindow.close()
                
            else:
                # Define Headers
                if self.fileName == 'TimeTemperaturePaths.csv':
                    Headers = ['Time (Ma)', 'Temperature (°C)']
                    textTemp = self.HeadersLabel
                    cmap = matplotlib.cm.get_cmap('hsv')
                    colors = [cmap(i*1/nbSamples) for i in range(len(nbGrain_per_sample)) for j in range(int(nbGrain_per_sample[i][0]))]
                    textlegend = [textTemp[i+1] for i in dataIndex[:]]
                    itemText = 'tT_paths'
                elif self.fileName == 'TimeAge.csv':
                    Headers = ['Time (Ma)', 'Age (Ma)']
                    textlegend = ['Grain'+str(i+1)
                                  for i in range(len(data_cleaned[1, :])-2)]
                    itemText = str(self.fileName[0:-4])
                elif self.fileName.startswith("Grain"):
                    if self.He4He3:
                        Headers = ['\u03A3F$^3$He', 'Rstep/Rbulk']
                        cmap = matplotlib.cm.get_cmap('hsv')
                        colors = [cmap(i*1/nbSamples) for i in range(len(nbGrain_per_sample)) for j in range(int(nbGrain_per_sample[i][0]))]
                        data_cleaned = pd.DataFrame([self.cum3He,self.He_ratio]).T
                        data_cleaned = np.array(data_cleaned.values)
                        itemText = "4He/3He spectrum "+str(self.fileName[-8:-4])
                    else:
                        Headers = np.asarray(self.HeadersLabel)[index_headers]
                        itemText = str(self.fileName[0:-4])
                        textlegend = ['']
                    
                self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
                self.toolbar = NavigationToolbar(self.plotSpace,self)
                if 'error (My)' in Headers:
                    indexError = np.where(Headers == 'error (My)')[0][0]
                    indexAge = np.where(Headers == 'age (My)')[0][0]
                    xdata = data_cleaned[:, 0]
                    ydata1 = data_cleaned[:, indexAge] + data_cleaned[:, indexError] 
                    ydata2 = data_cleaned[:, indexAge] - data_cleaned[:, indexError]
                    ydata = data_cleaned[:, indexAge]
                    self.plotSpace.axes.plot(xdata,ydata1, 'k--')
                    self.plotSpace.axes.plot(xdata,ydata2, 'k--')
                    self.plotSpace.axes.plot(xdata,ydata, 'k-')
                    self.plotSpace.axes.set_xlim([max(data_cleaned[:,0]),0])
                    # self.plotSpace.axes.set_aspect(1.0/self.plotSpace.axes.get_data_ratio(), adjustable='box')
                    textlegend = ['+1$\sigma$','-1$\sigma$','age']
                else:
                    for i in range(len(data_cleaned[1, :])-1):
                        #check for zeros at the en of the array
                        xdata = data_cleaned[:, 0]
                        # xdata = pgu.removeZeros(xdata,len(xdata))
                        ydata = data_cleaned[:, i+1]
                        # ydata = pgu.removeZeros(ydata,len(ydata))
                        if len(xdata) < len(ydata):
                            xdata =data_cleaned[:,0]
                            xdata = xdata[0:len(ydata)]
                        elif len(ydata) < len(xdata):
                            ydata =data_cleaned[:,i+1]
                            ydata = ydata[0:len(xdata)]
                        if self.He4He3 and self.fileName.startswith("Grain"):
                            #Compute fitting polynomial curve
                            if self.FitFunct == 1: #quadratic function
                                He43Fit, pcov = curve_fit(pgu.objective_poly, xdata, ydata,p0=[ydata[0],1,1])
                                #Summarize fitting parameters
                                a,b,c = He43Fit
                                std_error = np.sqrt(np.diag(pcov))
                                annotationText = 'y =  %.3f + %.3fx + %.3fx$^2$ \n 1$\sigma$: a = %.2f  b = %.2f c = %.2f' % (a, b, c, std_error[0], std_error[1], std_error[2])
                                yFit = pgu.objective_poly(xdata, a, b, c)
                            elif self.FitFunct == 2: #Power law function
                                He43Fit, pcov = curve_fit(pgu.objective_power, xdata, ydata,p0=[1,0.5])
                                #Summarize fitting parameters
                                a,b = He43Fit
                                std_error = np.sqrt(np.diag(pcov))
                                annotationText = 'y =  %.3fx$^{%.3f}$ \n 1$\sigma$: a = %.2f  b = %.2f' % (a, b, std_error[0], std_error[1])
                                yFit = pgu.objective_power(xdata, a, b)
                            elif self.FitFunct == 3: #3th order polynomial function
                                He43Fit, pcov = curve_fit(pgu.objective_poly3, xdata, ydata,p0=[ydata[0],1,1,1])
                                #Summarize fitting parameters
                                a,b,c,d = He43Fit
                                std_error = np.sqrt(np.diag(pcov))
                                annotationText = ('y =  %.3f + %.3fx + %.3fx$^2$ + %.3fx$^3$ \n 1$\sigma$: a = %.2f  b = %.2f c = %.2f d = %.2f' %
                                                  (a, b, c, d, std_error[0], std_error[1], std_error[2], std_error[3]))
                                yFit = pgu.objective_poly3(xdata, a, b, c, d)
                            elif self.FitFunct == 4: #4th order polynomial function
                                He43Fit, pcov = curve_fit(pgu.objective_poly4, xdata, ydata,p0=[ydata[0],1,1,1,1])
                                #Summarize fitting parameters
                                a,b,c,d,e = He43Fit
                                std_error = np.sqrt(np.diag(pcov))
                                annotationText = ('y =  %.3f + %.3fx + %.3fx$^2$ + %.3fx$^3$ + %.3fx$^4$ \n 1$\sigma$: a = %.2f  b = %.2f c = %.2f d = %.2f e = %.2f' %
                                                  (a, b, c, d, e, std_error[0], std_error[1], std_error[2], std_error[3],
                                                     std_error[4]))
                                yFit = pgu.objective_poly4(xdata, a, b, c, d, e)
                            elif self.FitFunct == 5: #5th order polynomial function
                                He43Fit, pcov = curve_fit(pgu.objective_poly5, xdata, ydata,p0=[ydata[0],1,1,1,1,1])
                                #Summarize fitting parameters
                                a,b,c,d,e,f = He43Fit
                                std_error = np.sqrt(np.diag(pcov))
                                annotationText = ('y =  %.3f + %.3fx + %.3fx$^2$ + %.3fx$^3$ + %.3fx$^4$ + %.3fx$^5$ \n 1$\sigma$: a = %.2f  b = %.2f c = %.2f d = %.2f e = %.2f f = %.2f' %
                                                  (a, b, c, d, e, f, std_error[0], std_error[1], std_error[2], std_error[3],
                                                     std_error[4], std_error[5]))
                                yFit = pgu.objective_poly5(xdata, a, b, c, d, e, f)

                            textlegend = ['predictions',annotationText]
                            self.plotSpace.axes.plot(xdata,ydata, 'ro',label=numLabel[countLabel]+textlegend[i])
                            countLabel += 1
                            self.plotSpace.axes.plot(xdata,yFit, 'k-',label=numLabel[countLabel]+textlegend[i])
                            countLabel += 1
                            
                            self.plotSpace.axes.set_xlim([0,1])
                        else:
                            self.plotSpace.axes.plot(xdata,ydata, '-',color=colors[i],label=numLabel[countLabel]+textlegend[i])
                            countLabel += 1
                    self.plotSpace.axes.set_title(itemText)
                    self.plotSpace.axes.set_xlabel(Headers[0])
                    self.plotSpace.axes.set_ylabel(Headers[1])
                    self.plotSpace.fig.legend(textlegend)
                    # self.plotSpace.axes.set_aspect(1.0/self.plotSpace.axes.get_data_ratio(), adjustable='box')
                #Close the Window
                self.TableWindow.close()
                # Create item for the tree view
                itemTree = pgu.StandardItem(itemText+"_"+conf.ProjectName, 10)
                self.treeModel.appendRow(itemTree)
                self.treeView.setModel(self.treeModel)
                
        elif self.file_extension == '.vtk': # ---------------------------------------------
            # First check that all the minimum required parameters have bee provided
            if self.parametersVTK.contours and not self.parametersVTK.ncon.text():
                QErrorMessage().showMessage("One or several parameters are missing ! Please, provide them before to move on.")
                return
            
            # Get the name of the data chose by the user
            self.dataVTK_name = self.parametersVTK.DataSelection.currentText()
            
            # close the window
            self.wVTK.close()
            
            # Extract the data: is 1D array
            # 2 kinds, AgesXXX or PecubeXXX
            if self.dataVTK_name == '4He/3He res':
                data_EdgeAge = self.file_vtk.get_array(conf.Variable_names['Edge Age for 4He/3He']) # Get Edge Age
                dataAHe = self.file_vtk.get_array('ApatiteHeAge') 
                data_vtk = (1 - (data_EdgeAge/dataAHe)) *100 # IRP
            else:
                data_vtk = self.file_vtk.get_array(self.dataVTK_name)
            
            # Get nskip: topographic points skipped by the user
            nskip = int(self.MainWin.ParamTable.ParamTable.nskipEdit8.text())
            
            # Get nx and ny
            nx = int((int(self.MainWin.ParamTable.ParamTable.nxEdit2.text())-1)/nskip+1)
            ny = int((int(self.MainWin.ParamTable.ParamTable.nyEdit3.text())-1)/nskip+1)
            
            # The file PecubeXXX contains data in 3D, but we only want surface nodes
            # The surface nodes correspond to every nz node in the data
            if self.fileName.startswith('Pecube'):
                # Get nz
                nz = int(self.MainWin.ParamTable.ParamTable.nzEdit3.text())
                
                # Get coordinates
                coord = self.file_vtk.points
                coordx_temp = coord[:,0]
                coordy_temp = coord[:,1]
                coordz_temp = coord[:,2]
                
                # Extract every nz point from coordinates and data
                coordx = coordx_temp[nz-1::nz]
                coordy = coordy_temp[nz-1::nz]
                coordz = coordz_temp[nz-1::nz]
                
                # Get number of points (handle nskip provided by the user)
                npoints = coordx.size
                
                # Do we want to interpolate any depth or isotherm?
                if self.parametersVTK.interpolate:
                    # Which data ?
                    data_name = self.parametersVTK.data
                    
                    # Get value, can be depth or isotherm
                    value = float(self.parametersVTK.isotherm.text())
                    
                    if data_name == 'Depth':
                        # Depth from surface
                        depth_from_surface = value
                        # Extract each z column of the mesh and interpolate value
                        temperature_at_depth = np.zeros((npoints))
                        for i in range(npoints):
                            z_profile = coordz_temp[i*nz:(i+1)*nz]
                            T_profile = np.flipud(data_vtk[i*nz:(i+1)*nz])
                            xi = np.array(z_profile)
                            yi = np.array(T_profile)
                            f = interp1d(z_profile, T_profile, kind = self.parametersVTK.interpolation)
                            new_temp = f(depth_from_surface)
                            temperature_at_depth[i] = new_temp
                        data_vtk = temperature_at_depth
                        self.dataVTK_name = "Temperature at depth = " + str(value) + " km"
                        
                    elif data_name == 'Isotherm':
                        ### Interpolate depth of isotherm ###
                        # Extract columns of temperatures
                        depth_at_temperature = np.zeros((npoints))
                        for i in range(npoints):
                            T_profile = np.flipud(data_vtk[i*nz:(i+1)*nz]) #increasing temperatures
                            # Extract depth profile
                            z_profile = coordz_temp[i*nz:(i+1)*nz]
                            depth_from_surface = np.flipud(coordz[i] - z_profile)
                            f = interp1d(T_profile, z_profile, kind = self.parametersVTK.interpolation)
                            new_depth = f(value)
                            depth_at_temperature[i] = new_depth
                        data_vtk = depth_at_temperature
                        self.dataVTK_name = "Depth of isotherm " + str(value) + " °C from the surface"
                        
                else: # if no interpolation required
                    # Extract data at the surface of the mesh
                    data_vtk = data_vtk[nz-1::nz]
                # Set coloramp
                if self.dataVTK_name == 'Temperature' or self.dataVTK_name == "Temperature at depth = " + str(value) + " km":
                    reversed_map = plt.cm.get_cmap('coolwarm')
                elif self.dataVTK_name == "Depth of isotherm " + str(value) + " °C":
                    reversed_map = plt.cm.get_cmap('copper_r')
                else:
                    colmap = plt.cm.get_cmap('Oranges')
                    reversed_map = colmap.reversed()
                            
            else: # if AgesXXX.vtk
                # Get coordinates
                coord = self.file_vtk.points
                coordx = coord[:,0]
                coordy = coord[:,1]
                coordy = np.flipud(coordy)
                coordz = coord[:,2]
                
                colmap = plt.cm.get_cmap('Oranges')
                reversed_map = colmap.reversed()
                
            # Want normalized value
            if self.parametersVTK.normalized:
                mindata = np.min(data_vtk)
                maxdata = np.max(data_vtk)
                data_vtk = (data_vtk-mindata) / (maxdata-mindata)
                
            # Get coordinates
            try:
                X = np.transpose(np.reshape(coordx,(ny,nx)))
                Y = np.transpose(np.reshape(coordy,(ny,nx)))
            except ValueError:
                QErrorMessage().showMessage("These data cannot be plot. The shape of the array do not correspond to the shape of the mesh grid."+
                             "\nSee more details in the console...")
                print('Shape X: ', np.shape(coordx))
                print('nx: ', nx)
                print('ny: ', ny)
                return
            # Reshape data accordingly
            Z = np.reshape(data_vtk,(ny,nx))
            if Z.shape[0] == X.shape[1]:
                Z = np.transpose(Z)
            # Create item for the tree view
            itemTree = pgu.StandardItem(self.dataVTK_name+"_"+conf.ProjectName, 10)
            self.treeModel.appendRow(itemTree)
            self.treeView.setModel(self.treeModel)
            
            # Plot all the data
            self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
            self.toolbar = NavigationToolbar(self.plotSpace,self)
            plot2D = self.plotSpace.axes.pcolormesh(X,Y,Z,cmap=reversed_map)
            cb = self.plotSpace.fig.colorbar(plot2D, ax=self.plotSpace.axes, shrink=0.4, aspect=7)
            self.plotSpace.axes.set_title(self.dataVTK_name)
            self.plotSpace.axes.set_xlabel('Distance (km)')
            self.plotSpace.axes.set_ylabel('Distance (km)')  
            self.plotSpace.axes.set_aspect(1.0/self.plotSpace.axes.get_data_ratio(), adjustable='box')
            self.plotSpace.axes.invert_yaxis()
            # Does the user want contour plot?
            if self.parametersVTK.contours:
                # First load elevations and correct for the thickness of the crust
                elevations =  np.reshape(coordz,(ny,nx))
                if elevations.shape[0] == X.shape[1]:
                    elevations = np.transpose(elevations)
                self.plotSpace.axes.contour(X,Y,elevations,int(self.parametersVTK.ncon.text()),
                                            colors='k',linewidths=1)
            self.plotSpace.axes.set_aspect('equal',anchor='C')
        #Check number of plots
        if self.CountPlotsColumn == 2 and self.CountPlotsRow < 1:
            self.CountPlotsRow += 1
            self.CountPlotsColumn = 0
        elif self.CountPlotsRow == 1 and self.CountPlotsColumn == 2:
            self.CountPlotsRow = 0
            self.CountPlotsColumn = 0
            
        #Add toolbar with 2D plot
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        if self.MainLay.count() > 4:
            U_interface.MessageBox(text1='The maximum number of plots is reached. Please, remove one of the plot to be able to load another one.')
            return
        if self.CountPlotsColumn == 0 and self.CountPlotsRow == 0:
            self.Splitter1.replaceWidget(0,Q)
            self.Splitter1.widget(0).setVisible(True)
            self.MainLay.refresh()
        elif self.CountPlotsRow == 0 and self.CountPlotsColumn == 1:
            self.Splitter1.replaceWidget(1,Q)
            self.Splitter1.widget(1).setVisible(True)
            self.MainLay.refresh()
        elif self.CountPlotsRow == 1 and self.CountPlotsColumn == 0:
            self.Splitter2.replaceWidget(0,Q) #Add temporary empty widget
            self.Splitter2.widget(0).setVisible(True)
            self.Splitter2.widget(1).setVisible(True)
            self.Splitter2.setVisible(True)
            self.MainLay.refresh()
        elif self.CountPlotsRow == 1 and self.CountPlotsColumn == 1:
            self.Splitter2.replaceWidget(1,Q)
            self.Splitter2.widget(1).setVisible(True)
            self.MainLay.refresh()
        else:
            print('issue for plotting 2D')

        # #Add plot space
        self.SavePlots.append(self.plotSpace)
        self.SaveToolbar.append(self.toolbar)
        self.CountPlotsColumn += 1
        self.indexSplitter += 1
          
    #----------------------------------------------------------
    def load3DData(self):
        # Import 3D file
        dialog = QFileDialog()
        dialog.setDirectory(os.path.join(self.PecubePath))
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.filePath, _ = dialog.getOpenFileNames(
            self, "QFileDialog.getOpenFileName()", "", "VTK Files (*.vtk)", options=options)
        if self.filePath:
            
            # extract file name
            extracted_Fnames = self.filePath[0].split("/")
            str1 = ''
            strlist = extracted_Fnames[-2:]
            str2 = str1.join(strlist)
            Nb_char = len(str2)+2  # count number of character after folder name
            self.fileName3D = extracted_Fnames[-1]
            self.folderName = self.filePath[0][0:-Nb_char]
            
            #Update data Combo box
            if self.DataCombo.count() <= 1 or self.CurrentFolder != self.folderName:
                self.CurrentFolder = self.folderName #store the current project
                self.get_dataList()
                   
            #Read VTK file with Pyvista
            self.grid = pv.read(os.path.join(self.PecubePath,self.folderName,"VTK",self.fileName3D))
            #Get min and max values
            DataRange = self.grid.get_data_range()
    
            # set 3D interface
            self.Object3D = pgu.Canvas3D()
            self.MainMesh= self.Object3D.plotter.add_mesh(
            self.grid,show_edges=False, cmap='RdBu', flip_scalars=True, scalar_bar_args=self.sargs)
            #Add axes orientation (x,y,z)
            self.Object3D.plotter.add_axes(line_width=2,interactive=True)
            
            #Get Data arrays
            self.VTKData = self.grid.array_names       
      
            ####### Check for sample location #######
            for t in range(2):
                try:
                    file = os.path.join(self.PecubePath,self.folderName,"data",self.ParametersInput.DParameters['data_folder'] ,self.ParametersInput.DParameters['data_folder']+".csv")
                    graindata = pd.read_csv(file, sep=',', skiprows=0)
                    datagrain = np.asarray(graindata)
                    # Get nskip: topographic points skipped by the user
                    nskip = int(self.MainWin.ParamTable.ParamTable.nskipEdit8.text())
                    y0 = float(self.MainWin.ParamTable.ParamTable.lat0Edit5.text())
                    dy = float(self.MainWin.ParamTable.ParamTable.dLatEdit7.text())
                    dx = float(self.MainWin.ParamTable.ParamTable.dLonEdit6.text())
                    ny = int(self.MainWin.ParamTable.ParamTable.nyEdit3.text())
                    nx = int(self.MainWin.ParamTable.ParamTable.nxEdit2.text())
                    x0 = float(self.MainWin.ParamTable.ParamTable.lon0Edit4.text())
                    #Calculate distance between points and left bottom corner
                    grainlon = (datagrain[:, 1])
                    grainlat = (datagrain[:, 2])
                    xl = x0 + (nx-1)*dx
                    yl = y0 + (ny-1)*dy
                    xlKm = dx*(nx-1)*111.11*cos((y0+dy*ny/2)*3.141592654/180)
                    ylKm = dy*(ny-1)*111.11
                    grainlat = (grainlat-y0)/(yl-y0)*ylKm
                    grainlon = (grainlon-x0)/(xl-x0)*xlKm
                    maxaltitude = self.read_elevation()
                    filename = os.path.join(self.PecubePath,self.folderName,"output","CompareAGE.csv")
                    if (os.path.exists(os.path.join(self.PecubePath,self.folderName,"data","Samples_settings.txt")) and
                        int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']])==2):
                        #read final elevation of samples (in m)
                        elevations = []
                        nbSamples = int(self.ParametersInput.DParameters['Nb_Samples'])
                        nbGrain_per_sample = self.ParametersInput.TableParameters['nb_grains']
                        #Get name of grains
                        SampleName = [v  for k,v in self.ParametersInput.TableParameters.items() if 'SampleName' in k]
                        textlegend = [SampleName[i]+'_'+str(j+1) for i in range(len(nbGrain_per_sample)) for j in range(int(nbGrain_per_sample[i][0]))]
                        elevations = self.get_elevation(filename)
                        if elevations.any():
                            elevations = elevations/1000
                            array_z =elevations+maxaltitude #Elevation of samples
                            array = np.ones((len(grainlat),3))
                            for i in range(len(grainlat)):
                                array[i,0] = grainlon[i]
                                array[i,1] = grainlat[i]
                                array[i,2] = array_z[i]
                            #Use Pyvista to plot samples as spheres  
                            self.poly = pv.PolyData(array)
                            self.poly["My Labels"] = textlegend
                            self.Samples_loc = self.Object3D.plotter.add_point_labels(self.poly,"My Labels",point_size=10,
                                                                font_size=10,name='Sample_location')
                    else:
                        pass
                except FileNotFoundError:
                    print(self.PecubePath,self.folderName,"data","Sample","Sample.csv")
                    print('No sample specific location')
                    pass
                except AttributeError:
                    print("Attribute Error, load old input and retry.")
      
            #Create properties tab
            #store 3D canvas
            self.Canvas3D.append(self.Object3D)
            self.activePlot = 1
            self.Prop3D = pgu.Properties3D(self,self.MainWin,self.Object3D,self.MainMesh,self.Samples_loc,
                                        self.grid,self.activePlot,self.sargs,self.poly)
            #set min and max value in properties tab
            self.Prop3D.MinVal3D.setText(str(DataRange[0]))
            self.Prop3D.MaxVal3D.setText(str(DataRange[1]))
            #Update the combobox to set properties in the pecube interface
            self.Prop3D.Select3dData.clear()
            self.Prop3D.Select3dData.addItems(self.VTKData)
            self.Prop3D.Select3dData.activated.connect(self.Prop3D.DataChanged) 
            if self.Samples_loc: #enable show sample locations
                self.Prop3D.SLocCheck.setEnabled(True)
                self.Prop3D.SLocCheck.setChecked(True)
            
            # add a new tab
            VBox = QVBoxLayout()
            # create  toolbar object
            toolbar = pgu.toolbar_3Dplot(self,self.MainWin,self.Object3D)
            VBox.addWidget(toolbar.toolbar)
            VBox.addWidget(self.Object3D.Frame)
            Q = QWidget()
            Q.setLayout(VBox)
            self.TabPlot.addTab(Q, "3D Viewer")
            # Create item for the tree view
            itemTree = pgu.StandardItem(str(self.fileName3D[0:-4])+"_"+conf.ProjectName, 10)
            self.treeModel3D.appendRow(itemTree)
            self.treeView3D.setModel(self.treeModel3D)
            # self.itemID.append(1)
            #Store newly created 3D plot
            self.List3DProp.append(self.Prop3D.Dock3dProperties)
            self.MainWin.addDockWidget(Qt.LeftDockWidgetArea,self.Prop3D.Dock3dProperties)
            self.MainWin.tabifyDockWidget(self.DockData,self.Prop3D.Dock3dProperties)
            self.Dock2dProperties.hide()
            if self.List3DProp:
                for i in range(self.TabPlot.count()-1):
                    if self.TabPlot.tabText(i) == "3D Viewer":
                        self.List3DProp[i-1].hide()
                        
    #----------------------------------------------------------
    def read_elevation(self):
        #to read the elevation of the sample location
        #Go search the crust thickness in the input file
        input_file = os.path.join(self.PecubePath,self.folderName,"input","Pecube.in")
        Read_parameters = {}
        with open(input_file) as file:
            for line in file:
                if "=" in line:
                    (key, _, par) = line.split()
                    Read_parameters[str(key)] = par
                    
        DefaultParameters = conf.DefaultParameterValues()
        for k, v in Read_parameters.items():
            if k == 'thickness':
                thickness = float(v)
                break
            else:
                thickness = float(DefaultParameters.DParameters['thickness'])
        return thickness
    
        
    #----------------------------------------------------------
    def ShowData(self):
        # to show/hide data in plot
        # #first delete all widget from MainLay
        # for i in reversed(range(self.MainLay.count())):
        #     self.MainLay.widget(i).setVisible(False)
        countCol = 0
        countRow = 0
        i = 0
        while self.treeModel.item(i):
            if self.treeModel.item(i).isCheckable():
                if countCol == 2 and countRow == 0:
                    countCol = 0
                    countRow += 1
                elif countRow == 1 and countCol==2:
                    countRow = 0
                    countCol = 0
                if self.treeModel.item(i).checkState() == Qt.Checked:
                    self.activePlot = 0
                    if i > 1:
                        self.MainLay.widget(1).widget(i-2).setVisible(True)#addWidget(self.SaveToolbar[i],countRow-1,countCol)
                    else:
                        self.MainLay.widget(0).widget(i).setVisible(True)
                    #self.MainLay.addWidget(self.SavePlots[i],countRow,countCol)
                    #self.SaveToolbar[i].show()
                    #self.SavePlots[i].show()
                else:
                    self.activePlot = 0
                    if i > 1:
                        self.MainLay.widget(1).widget(i-2).setVisible(False)
                    else:
                        self.MainLay.widget(0).widget(i).setVisible(False)
                countCol +=1
            i += 1
    
    #----------------------------------------------------------
    def RemoveData(self):
        ### to remove data from tree list ###
        countCol = 0
        countRow = 0
        i = 0
        index = []
        #Delete from tree view and storeed plots and toolbars
        try:
            while self.treeModel.item(i):
                if self.treeModel.item(i).isCheckable():
                    if self.treeModel.item(i).checkState() == Qt.Checked:
                        self.treeModel.removeRow(i)
                        self.SavePlots.pop(i)
                        self.SaveToolbar.pop(i)
                        i -= 1
                i+=1
        except IndexError:
            pass
                
        # Replot the left over
        self.CountPlotsColumn = 0
        self.CountPlotsRow = 0
        for i in range(4):
            if i < len(self.SavePlots):
                VBox = QVBoxLayout()
                VBox.addWidget(self.SaveToolbar[i])
                VBox.addWidget(self.SavePlots[i])
                #add Widget
                Q = QWidget()
                Q.setLayout(VBox)
                if i < 2:
                    self.CountPlotsColumn += 1
                elif i == 2:
                    self.CountPlotsColumn = 0
                    self.CountPlotsRow = 1
                elif i > 2:
                    self.CountPlotsColumn += 1
                    
            else:  # Feed empty widget to the remaining slots
                Q = QWidget()
                
            if i < 2:
                self.Splitter1.replaceWidget(i,Q)
            else:
                self.Splitter2.replaceWidget(i-2,Q)
            
        if self.MainLay.count() < 3:
            self.Splitter2.setVisible(False)
        
        self.MainLay.refresh()
        


# -------------------- End classes for Graphic Window ---------------------------