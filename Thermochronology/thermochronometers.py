#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This module contains parameters related to the thermochronometers settings.

@author: maxime bernard

"""

import matplotlib.backends.backend_svg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QAction,QPushButton,QWidget,QLabel,QCheckBox,QVBoxLayout,
                             QGridLayout,QComboBox,QGroupBox,QStackedWidget,QHBoxLayout,
                             QErrorMessage,QLineEdit,QTableWidgetItem,QRadioButton,QFrame,
                             QTreeView,QListView)
from PyQt5.QtCore import Qt
from PyQt5.Qt import QStandardItemModel, QStandardItem
from PyQt5.QtGui import (QIntValidator)
import configs as conf
import os
import stat
import sys
import numpy as np
from scipy.interpolate import interp1d
import main as pgui
import Utils.PGUI_utils as pgu
import xarray as xr
import Thermochronology.helium as Helium
import Thermochronology.fission_track as FT
import Thermochronology.argon as Argon
import Thermochronology.trapped_charge as TC
import Thermochronology.settings as Thermo_settings


##############################################################################
############################## Classes #######################################
##############################################################################

#----------------------------------------------------------------
class Toolbar_ages:
    """
    This class creates a customized navigation toolbar for 1D age plots.
    It allows to show/hide data according to the thermochronometer.
    
    """
    def __init__(self,plotCanvas,parent):
        super().__init__()
        self.toolbar = NavigationToolbar(plotCanvas, None)
        self.plotCanvas = plotCanvas
        self.Separator = QAction()
        self.Separator.setSeparator(True)
        self.countAHe = 1
        self.countZHe = 1
        self.countAFT = 1
        #AHe
        self.ShowAHeButton = QPushButton("AHe")
        self.ShowAHeButton.setFont(conf.font8)
        self.ShowAHeButton.setFlat(True)
        self.ShowAHeButton.clicked.connect(lambda: self.show_AHe_data())
        self.ShowAHeButton.setEnabled(False)
        self.ShowAHeButton.setToolTip("Show/Hide AHe data")
        #ZHe
        self.ShowZHeButton = QPushButton("ZHe")
        self.ShowZHeButton.setFont(conf.font8)
        self.ShowZHeButton.setFlat(True)
        self.ShowZHeButton.clicked.connect(lambda: self.show_ZHe_data())
        self.ShowZHeButton.setEnabled(False)
        self.ShowZHeButton.setToolTip("Show/Hide ZHe data")
        #AFT
        self.ShowAFTButton = QPushButton("AFT")
        self.ShowAFTButton.setFont(conf.font8)
        self.ShowAFTButton.setFlat(True)
        self.ShowAFTButton.clicked.connect(lambda: self.show_AFT_data())
        self.ShowAFTButton.setEnabled(False)
        self.ShowAFTButton.setToolTip("Show/Hide AFT data")
        
        #fig = plt.figure()
        
        #Update the toolbar with new buttons
        self.toolbar.addAction(self.Separator)
        self.toolbar.addWidget(self.ShowAHeButton)
        self.toolbar.addWidget(self.ShowZHeButton)
        self.toolbar.addWidget(self.ShowAFTButton)
    
    #-----------------------------------------------------------
    def show_AHe_data(self):
        """
        Show/hide AHe data in plots.
        
        """
        
        if self.countAHe == 0:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "AHe" or line.get_gid() == "AHeObs":
                    line.set_visible(True)
            self.countAHe = 1
        elif self.countAHe == 1:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "AHe":
                    line.set_visible(True)
                elif line.get_gid() == "AHeObs":
                    line.set_visible(False)
            self.countAHe = 2
        else:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "AHe":
                    line.set_visible(False)
                elif line.get_gid() == "AHeObs":
                    line.set_visible(False)
            self.countAHe = 0
        self.plotCanvas.fig.canvas.draw()
    
    #-----------------------------------------------------------
    def show_ZHe_data(self):
        """
        Show/hide ZHe data in plots.
        
        """
        if self.countZHe == 0:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "ZHe" or line.get_gid() == "ZHeObs":
                    line.set_visible(True)
            self.countZHe = 1
        elif self.countZHe == 1:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "ZHe":
                    line.set_visible(True)
                elif line.get_gid() == "ZHeObs":
                    line.set_visible(False)
            self.countZHe = 2
        else:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "ZHe":
                    line.set_visible(False)
                elif line.get_gid() == "ZHeObs":
                    line.set_visible(False)
            self.countZHe = 0
        self.plotCanvas.fig.canvas.draw()
        
    #-----------------------------------------------------------
    def show_AFT_data(self):
        """
        Show/hide AFT data in plots.
        
        """
        if self.countAFT== 0:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "AFT" or line.get_gid() == "AFTObs":
                    line.set_visible(True)
            self.countAFT = 1
        elif self.countAFT == 1:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "AFT":
                    line.set_visible(True)
                elif line.get_gid() == "AFTObs":
                    line.set_visible(False)
            self.countAFT = 2
        else:
            for line in self.plotCanvas.axes.findobj():
                if line.get_gid() == "AFT":
                    line.set_visible(False)
                elif line.get_gid() == "AFTObs":
                    line.set_visible(False)
            self.countAFT = 0
        self.plotCanvas.fig.canvas.draw()
        
        
##################################################
class ThermoAges(QWidget):
    """
    This class defines the thermochronometers to plot in an age-elevation profile
    when in "for all nodes" mode.
    
    """
    #see class 'GraphWin'
    def __init__(self):
        super().__init__()
        
        self.store_systems = {}
        self.Luminescence = QLabel("Trapped_charge:")
        self.Luminescence.setFont(conf.fontBold12)
        self.Helium = QLabel("(U-Th)/He:")
        self.Helium.setFont(conf.fontBold12)
        self.FissionTrack = QLabel("Fission tracks:")
        self.FissionTrack.setFont(conf.fontBold12)
        self.Argon = QLabel("Argon:")
        self.Argon.setFont(conf.fontBold12)
        
        self.ThL = QCheckBox('Thermo-Luminescence')
        self.OSL = QCheckBox('OSL')
        self.ESR = QCheckBox('ESR')
        self.AHe = QCheckBox('Apatite')
        self.ZHe = QCheckBox("Zircon")
        self.AFT = QCheckBox("Apatite")
        self.ZFT = QCheckBox("Zircon")
        self.AFeld = QCheckBox("Feldspar")
        self.ABiot = QCheckBox("Biotite")
        self.AMusc = QCheckBox("Muscovite")
        self.AHorn = QCheckBox("Hornblende")
        
        self.ThL.stateChanged.connect(lambda: self.updateSystem(self.ThL,9))
        self.OSL.stateChanged.connect(lambda: self.updateSystem(self.OSL,10))
        self.ESR.stateChanged.connect(lambda: self.updateSystem(self.ESR,11))
        self.AHe.stateChanged.connect(lambda: self.updateSystem(self.AHe,1))
        self.ZHe.stateChanged.connect(lambda: self.updateSystem(self.ZHe,2))
        self.AFT.stateChanged.connect(lambda: self.updateSystem(self.AFT,3))
        self.ZFT.stateChanged.connect(lambda: self.updateSystem(self.ZFT,4))
        self.AFeld.stateChanged.connect(lambda: self.updateSystem(self.AFeld,5))
        self.ABiot.stateChanged.connect(lambda: self.updateSystem(self.ABiot,6))
        self.AMusc.stateChanged.connect(lambda: self.updateSystem(self.AMusc,7))
        self.AHorn.stateChanged.connect(lambda: self.updateSystem(self.AHorn,8))
        
        self.VMainBox = QVBoxLayout()
        Grid1 = QGridLayout()
        Grid1.setSpacing(10)
        Grid1.addWidget(self.Luminescence,0,0,1,2)
        Grid1.addWidget(self.ThL,0,2,1,2)
        Grid1.addWidget(self.OSL,1,2,1,2)
        Grid1.addWidget(self.ESR,2,2,1,2)
        Grid1.addWidget(self.Helium,3,0,1,2)
        Grid1.addWidget(self.AHe,3,2,1,2)
        Grid1.addWidget(self.ZHe,4,2,1,2)
        Grid1.addWidget(self.FissionTrack,5,0,1,2)
        Grid1.addWidget(self.AFT,5,2,1,2)
        Grid1.addWidget(self.ZFT,6,2,1,2)
        Grid1.addWidget(self.Argon,7,0,1,2)
        Grid1.addWidget(self.AFeld,7,2,1,2)
        Grid1.addWidget(self.ABiot,8,2,1,2)
        Grid1.addWidget(self.AMusc,9,2,1,2)
        Grid1.addWidget(self.AHorn,10,2,1,2)
        
        Grid1.setColumnStretch(4,1)
        Grid1.setRowStretch(11,1)
        self.VMainBox.addLayout(Grid1)
    
    #----------------------------------------------------------
    def updateSystem(self,sys,flag):
        #flag : 1 = AHe
        #2 = ZHe
        #3 = AFT
        #4 = ZFT
        #5 = AFeld
        #6 = ABiot
        #7 = AMusc
        #8 = AHorn
        if sys.isChecked(): #alphabetic order
            if flag == 1:
                self.store_systems[flag] = 'HeApatite'
            elif flag == 2:
                self.store_systems[flag] = 'HeZircon'
            elif flag == 3:
                self.store_systems[flag] = 'FTApatite'
            elif flag == 4:
                self.store_systems[flag] = 'FTZircon'
            elif flag == 5:
                self.store_systems[flag] = 'ArKFeldspar'
            elif flag == 6:
                self.store_systems[flag] = 'ArBiotite'
            elif flag == 7:
                self.store_systems[flag] = 'ArMuscovite'
            elif flag == 8:
                self.store_systems[flag] = 'ArHornblend'
            elif flag == 9:
                self.store_systems[flag] = 'ThL'
            elif flag == 10:
                self.store_systems[flag] = 'OSL'
            elif flag == 11:
                self.store_systems[flag] = 'ESR'
        else:
            del self.store_systems[flag]
            

#############################################################
class set_Thermochron_Parameters(QWidget):
    """ 
      This class set the parameters for each thermochronometers when using
      sample-specific prediction.
      This is what is shown in tab ('Ages')
      Parent is WindowsParameters.
      
    """
    def __init__(self, d, Param, PFolder):
        super().__init__()
        
        # Project folder name
        self.PFolder = PFolder
        self.Param = Param
        self.d = d
        self.parent = d
        self.Apatites = {}
        self.Zircons = {}
        self.AFTParam = 0
        
        # List to store zontaion profiles (see class GrainParameters)
        self.zonationUxList = []
        self.zonationUyList = []
        self.zonationThxList = []
        self.zonationThyList = []
        # Flag for active zonation
        self.zonation = 0
        
        # Parameters
        self.ParamTitle = QLabel('Thermochronological system:')
        self.ParamTitle.setFont(conf.font12)
        self.AgesCombo = QComboBox()
        self.AgesCombo.setMinimumContentsLength(5)
        self.groupBox = QGroupBox("Parameters")
        # Grain characterisics
        # self.GrainParam = GrainsParameters(self, self.Param,self.parent.nbSampleValue,self.parent.Samples,"AHE")
        # self.AFTGrainParam = GrainsParameters(self, self.Param,self.parent.nbSampleValue,self.parent.Samples,"AFT")
        # self.ZirconParam = ZirconGrainsParameters(self, self.Param,self.parent.nbSampleValue,self.parent.Samples,"ZHE")
        self.TLParam = TC.ThermoLuminescence(self,self.parent,self.Param,self.PFolder)
        self.OSLParam = TC.OSL(self,self.parent,self.Param,self.PFolder)
        self.ESRParam = TC.ESR(self,self.parent,self.Param,self.PFolder)
        self.apatiteHelium = Helium.apatite(self,self.parent,self.Param,self.PFolder)
        self.ZHeParam = Helium.zircon(self,self.parent,self.Param,self.PFolder)
        self.AFTParam = FT.apatite(self,self.parent,self.Param,self.PFolder)
        self.ZFTParam = FT.zircon(self,self.parent,self.Param,self.PFolder)
        self.KArParam = Argon.K_Feldspar(self,self.parent,self.Param,self.PFolder)
        self.BArParam = Argon.Biotite(self,self.parent,self.Param,self.PFolder)
        self.MArParam = Argon.Muscovite(self,self.parent,self.Param,self.PFolder)
        self.HArParam = Argon.Hornblende(self,self.parent,self.Param,self.PFolder)
        self.ParamLayout = QVBoxLayout()
        
        # Add in layout
        self.Stack = QStackedWidget()
        self.Stack.addWidget(self.TLParam.Q)
        self.Stack.addWidget(self.OSLParam.Q)
        self.Stack.addWidget(self.ESRParam.Q)
        self.Stack.addWidget(self.apatiteHelium.Q)
        self.Stack.addWidget(self.ZHeParam.Q)
        self.Stack.addWidget(self.AFTParam.Q)
        self.Stack.addWidget(self.ZFTParam.Q)
        self.Stack.addWidget(self.KArParam.Q)
        self.Stack.addWidget(self.BArParam.Q)
        self.Stack.addWidget(self.MArParam.Q)
        self.Stack.addWidget(self.HArParam.Q)
        self.VBox1 = QVBoxLayout()
        HBox = QHBoxLayout()
        HBox.addWidget(self.ParamTitle)
        HBox.addWidget(self.AgesCombo)
        HBox.addStretch(2)
        self.VBox1.addWidget(self.Stack)
        self.groupBox.setLayout(self.VBox1)
        self.Layout = QVBoxLayout()
        self.Layout.addLayout(HBox)
        self.Layout.addWidget(self.groupBox)
        
        # Signals
        self.AgesCombo.currentIndexChanged.connect(lambda: self.updateThermoSysParam())
    
    #-----------------------------------------------------------
    def updateTables(self):
        """To updates Tables in all thermo tabs """
        
        # Update number of rows
        self.TLParam.updateTableGrain(self.TLParam.GrainTable, 2, self.Param)
        self.OSLParam.updateTableGrain(self.OSLParam.GrainTable, 2, self.Param)
        self.ESRParam.updateTableGrain(self.ESRParam.GrainTable, 2, self.Param)
        self.apatiteHelium.updateTableGrain(self.apatiteHelium.GrainTable, 2, self.Param)
        self.ZHeParam.updateTableGrain(self.ZHeParam.GrainTable, 2, self.Param)
        self.AFTParam.updateTableGrain(self.AFTParam.GrainTable, 2, self.Param)
        self.ZFTParam.updateTableGrain(self.ZFTParam.GrainTable, 2, self.Param)
        self.KArParam.updateTableGrain(self.KArParam.GrainTable, 2, self.Param)
        self.BArParam.updateTableGrain(self.BArParam.GrainTable, 2, self.Param)
        self.MArParam.updateTableGrain(self.MArParam.GrainTable, 2, self.Param)
        self.HArParam.updateTableGrain(self.HArParam.GrainTable, 2, self.Param)
        
        # self.apatiteHelium.GrainParam.updateTableGrain(self.apatiteHelium.GrainParam.GrainTable, 2, self.Param, self)
        # self.ZHeParam.GrainParam.updateTableGrain(self.ZHeParam.GrainParam.GrainTable, 2, self.Param, self)
        
    #-----------------------------------------------------------
    def updateWidgets(self):
        """To updates widgets in all thermo tabs """
        self.apatiteHelium.updateWidgets()

    #-----------------------------------------------------------
    def updateComboBox(self):
        """ 
          To update the combo box according to the selected thermochronometers.
          
        """
        self.AgesCombo.clear()
        if self.parent.ageTLb.isChecked() or int(self.parent.Samples['TL_TotGrain']) > 0:
            self.AgesCombo.addItem('TL')
        if self.parent.ageOSLb.isChecked() or int(self.parent.Samples['OSL_TotGrain']) > 0:
            self.AgesCombo.addItem('OSL')
        if self.parent.ageESRb.isChecked() or int(self.parent.Samples['ESR_TotGrain']) > 0:
            self.AgesCombo.addItem('ESR')
        if self.parent.ageAHeb.isChecked() or int(self.parent.Samples['AHe_TotGrain']) > 0:
            self.AgesCombo.addItem('AHe')
        if self.parent.ageZHeb.isChecked() or int(self.parent.Samples['ZHe_TotGrain']) > 0:
            self.AgesCombo.addItem('ZHe')
        if self.parent.ageAFTb.isChecked() or int(self.parent.Samples['AFT_TotGrain']) > 0:
            self.AgesCombo.addItem('AFT')
        if self.parent.ageZFTb.isChecked() or int(self.parent.Samples['ZFT_TotGrain']) > 0:
            self.AgesCombo.addItem('ZFT')
        if self.parent.ageKArb.isChecked() or int(self.parent.Samples['KAr_TotGrain']) > 0:
            self.AgesCombo.addItem('KAr')
        if self.parent.ageBArb.isChecked() or int(self.parent.Samples['BAr_TotGrain']) > 0:
            self.AgesCombo.addItem('BAr')
        if self.parent.ageMArb.isChecked() or int(self.parent.Samples['MAr_TotGrain']) > 0:
            self.AgesCombo.addItem('MAr')
        if self.parent.ageHArb.isChecked() or int(self.parent.Samples['HAr_TotGrain']) > 0:
            self.AgesCombo.addItem('HAr')
        
    #-----------------------------------------------------------
    def updateThermoSysParam(self):
        """
          To update the displayed parameters according to the selected thermo system.
        
        """
        if self.AgesCombo.currentText() == 'TL': #Thermoluminescence
            self.Stack.setCurrentWidget(self.TLParam.Q)
        elif self.AgesCombo.currentText() == 'OSL': #OSL
            self.Stack.setCurrentWidget(self.OSLParam.Q)
        if self.AgesCombo.currentText() == 'ESR': #ESR
            self.Stack.setCurrentWidget(self.ESRParam.Q)
        elif self.AgesCombo.currentText() == 'AHe': #AHe
            self.Stack.setCurrentWidget(self.apatiteHelium.Q)
        elif self.AgesCombo.currentText() == 'ZHe': #ZHe
            self.Stack.setCurrentWidget(self.ZHeParam.Q)
        elif self.AgesCombo.currentText() == 'AFT': #AFT
            self.Stack.setCurrentWidget(self.AFTParam.Q)
        elif self.AgesCombo.currentText() == 'ZFT': 
            self.Stack.setCurrentWidget(self.ZFTParam.Q)
        elif self.AgesCombo.currentText() == 'KAr': 
            self.Stack.setCurrentWidget(self.KArParam.Q)
        elif self.AgesCombo.currentText() == 'BAr': 
            self.Stack.setCurrentWidget(self.BArParam.Q)
        elif self.AgesCombo.currentText() == 'MAr': 
            self.Stack.setCurrentWidget(self.MArParam.Q)
        elif self.AgesCombo.currentText() == 'HAr': 
            self.Stack.setCurrentWidget(self.HArParam.Q)    
    
    #----------------------------------------------------------
    def Merge_dico(self,dict1,dict2):
        """
          To merge two dictionaries with key and values one by one
        
        """
        res = {}
        list1 = []
        list2 = []
        for i in zip(dict1.keys(),dict2.keys()):
            list1.append(dict1[str(i[0])])
            list2.append(dict2[str(i[1])])
        return list1, list2
    
    #----------------------------------------------------------
    def savefile(self):
        """
         To save sample-specific input parameters in "Samples_settings.txt".
        
        """
        # if self.parent.AgeCombo.currentIndex() == 2:#Means we want grain-specific
        name = os.path.join(self.PFolder,"data","Samples_settings.txt")
        try: 
            file = open(name, 'w+')
        except PermissionError:
            if os.path.exists(name):
                os.chmod(name,stat.S_IWRITE)
                if sys.platform == 'win32' or sys.platform =='cygwin':
                    conf.win32api.SetFileAttributes(name,conf.win32con.FILE_ATTRIBUTE_NORMAL)
                os.remove(name)
            file = open(name, 'w+')
                
        ### Start of writing sample_setting file ####
        # First Line, Name of the data folder
        file.write(self.parent.dataFolderEdit1.text()+'\n')
        # Second line - Number of samples
        nSamples =str(self.parent.nbSampleValue.text())
        file.write(nSamples+'\n')
        # Third line - number of total grains = sum(max(nbgrain per sample))
        file.write(str(self.parent.Samples['Total_maxGrain'])+"\n")
        
        # Fourth line - thermochronometers
        try:
            file.write(str(self.parent.input_parameters['age_AHe_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_AHe_flag'])+'\t')
        try:
            file.write(str(self.parent.input_parameters['age_ZHe_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_ZHe_flag'])+'\t') 
        try:
            file.write(str(self.parent.input_parameters['age_AFT_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_AFT_flag'])+'\t') 
        try:
            file.write(str(self.parent.input_parameters['age_KAr_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_KAr_flag'])+'\t') 
        try:
            file.write(str(self.parent.input_parameters['age_BAr_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_BAr_flag'])+'\t') 
        try:
            file.write(str(self.parent.input_parameters['age_MAr_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_MAr_flag'])+'\t') 
        try:
            file.write(str(self.parent.input_parameters['age_HAr_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_HAr_flag'])+'\t') 
        try:
            file.write(str(self.parent.input_parameters['age_TL_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_TL_flag'])+'\t')
        try:
            file.write(str(self.parent.input_parameters['age_OSL_flag'])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters['age_OSL_flag'])+'\t')
        try:
            file.write(str(self.parent.input_parameters['age_ESR_flag'])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters['age_ESR_flag'])+'\n')
            
        # Fifth Line - 4He/3He?
        try:
            file.write(str(self.parent.Samples[conf.Variable_names['4He/3He']]))
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['4He/3He']]))         
        
        # Seventh line - zonation flag
        file.write('\n')
        file.write(str(self.zonation))
        file.write('\n')
        # Write sample names + Lat, lon, Elevation
        for i in range(int(self.parent.nbSampleValue.text())):
            sampleName_temp = self.parent.Samples['SampleName'+str(i+1)]
            file.write(self.parent.Samples['SampleName'+str(i+1)]+"\t")
            file.write(self.parent.Samples[sampleName_temp+'_lon']+"\t")
            file.write(self.parent.Samples[sampleName_temp+'_lat']+"\t")
            file.write(self.parent.Samples[sampleName_temp+'_Elev']+"\n")
        # from eigth line - write number of grain for each sample
        # TotGrain = 0
        for i in range(int(self.parent.nbSampleValue.text())):
            nb_grains = self.parent.Samples['MaxNbGrain_'+str(i+1)]
            # TotGrain += int(nb_grains)
            file.write(str(nb_grains)+"\n")
        
        # # Write grain characteristics (radius, U, Th) Apatite
        # try:
        #     for i in range(int(self.parent.nbSampleValue.text())):
        #         nb_grains = self.parent.Samples['ngrains'+str(i+1)]
        #         for j in range(int(nb_grains)):
        #             file.write(self.Apatites['grainsizeApatite'+str(i+1)+"_"+str(j+1)]+
        #                        '\t'+self.Apatites['U_grainApatite'+str(i+1)+"_"+str(j+1)]+
        #                        '\t'+self.Apatites['Th_grainApatite'+str(i+1)+"_"+str(j+1)]+
        #                        '\t'),#Radius, U, Th
        #             file.write(self.Zircons['grainsizeZircon'+str(i+1)+"_"+str(j+1)]+
        #                                     '\t'+self.Zircons['U_grainZircon'+str(i+1)+"_"+str(j+1)]+
        #                                     '\t'+self.Zircons['Th_grainZircon'+str(i+1)+"_"+str(j+1)]+
        #                                     '\n')
        # except IndexError:
        #     pass
        # except KeyError as KE:
        #     QErrorMessage().showMessage("Cannot find " + str(KE) + '.')
        #     return
        # Write parameters for AHe
        file.write(str(self.apatiteHelium.alphaCombo.currentIndex())+'\t')
        try:
            file.write(str(self.apatiteHelium.alphaDistCombo.currentIndex())+'\t')
        except RuntimeError: #Means use Monte Carlo model
            file.write('1\t')
        file.write(str(self.apatiteHelium.MDCombo.currentIndex())+'\t')
        file.write(str(self.apatiteHelium.D0Value.text())+'\t')
        file.write(str(self.apatiteHelium.ActEnerValue.text())+'\n')
        
        # Write parameters for ZHe
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['Diffusion_model_ZHe']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['Diffusion_model_ZHe']])+'\t')
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['D0_Zircon']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['D0_Zircon']])+'\t')
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['Ea_Zircon']])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['Ea_Zircon']])+'\n')
            
        # Write parameters for AFT
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['FissionTrackModel']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['FissionTrackModel']])+'\t') # By default Ketcham et al. (2007)
        try:
            file.write(str(self.parent.input_parameters['rhoST'])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters['rhoST'])+'\n')
            
       
        # Write parameters for KAr
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['D0_Feldspar']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['D0_Feldspar']])+'\t') # 
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['Ea_Feldspar']])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['Ea_Feldspar']])+'\n') # 
        
        # Write parameters for BAr
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['D0_Biotite']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['D0_Biotite']])+'\t') # 
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['Ea_Biotite']])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['Ea_Biotite']])+'\n') # 
            
        # Write parameters for MAr
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['D0_Muscovite']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['D0_Muscovite']])+'\t') # 
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['Ea_Muscovite']])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['Ea_Muscovite']])+'\n') # 
            
        # Write parameters for HAr
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['D0_Hornblende']])+'\t')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['D0_Hornblende']])+'\t') # 
        try:
            file.write(str(self.parent.input_parameters[conf.Variable_names['Ea_Hornblende']])+'\n')
        except KeyError:
            file.write(str(self.Param.DParameters[conf.Variable_names['Ea_Hornblende']])+'\n') # 
            
        # # Write Number of heating steps for 43He
        # file.write('\n')
        # try:
        #     file.write(str(self.apatiteHelium.He43Param.NbstepsEdit.text())+'\t')
        # except:
        #     file.write(str(self.Param.DParameters['NstepsHe43'])+'\t')
        # # Write Number of column for heating schedule
        # nb_col = self.apatiteHelium.He43Param.HeatingTable.columnCount()
        # nb_row = self.apatiteHelium.He43Param.HeatingTable.rowCount()
        # file.write(str(nb_col)+'\t')
        # file.write(str(nb_row)+'\n')
        # # Write the heating schedule(s)
        # count = 0
        # for k, v in self.apatiteHelium.He43Param.He43Storage.items():
        #     if count == nb_col:
        #         file.write('\n')
        #         count = 0
        #     file.write(str(v) + '\t')
        #     count += 1
                
        # Write Zonation profiles for each grain
        # if self.zonation:
        #     # First write number of data for grain, (evenly spaced = 49)
        #     file.write('\n')
        #     file.write(self.GrainParam.setNbLayers.text())
        #     # Then write Tables
        #     count = 0
        #     for i in range(int(self.parent.nbSampleValue.text())):
        #         nb_grains = self.parent.Samples['ngrains'+str(i+1)]
        #         for j in range(int(nb_grains)):
        #             lenData = max(len(self.zonationUxList[i]),len(self.zonationThxList[i]))
        #             # then write 2 columns with depth of the zonation layer, and Thorium/Uranium ratio
        #             # get the size of the grain
        #             gs = int(self.grains['grainsize'+str(i+1)+"_"+str(j)])
        #             ### Calculate the concentration ratio of each zonation layer ###
        #             # we interpolate the array that has less elements according to the other 
        #             if lenData == len(self.zonationUxList[count]):
        #                 xdata = np.array(self.zonationUxList[count]) 
        #                 Udata = np.array(self.zonationUyList[count]) #Uranium
        #                 # Interpolate data for thorium in regularly spaced array
        #                 f = interp1d(self.zonationThxList[count],self.zonationThyList[count])
        #                 ydata = Thdata = f(xdata)
        #                 # Compute ratio
        #                 if not  len(xdata) == len(ydata): #check whether Th and U vector have the same size
        #                    QErrorMessage().showMessage("Error in zonation profiles: Uranium and Thorium arrays don't have the same size...")
        #                    file.close()
        #                    return
        #                 # ratioPPM = np.array(Thydata)/np.array(self.zonationUyList[i])
        #             else: # Thorium array is smaller than the Uranium array
        #                 xdata = np.array(self.zonationThxList[count]) 
        #                 Thdata = np.array(self.zonationThyList[count]) #Thorium
        #                 # interpolate data for uranium
        #                 f = interp1d(self.zonationUxList[count],self.zonationUyList[count])
        #                 ydata = Udata = f(xdata)
        #                 # Compute ratio
        #                 if not  len(xdata) == len(ydata): #check whether Th and U vector have the same size
        #                    QErrorMessage().showMessage("Error in zonation profiles: Uranium and Thorium arrays don't have the same size...")
        #                    file.close()
        #                    return
                    
        #             ### Calculate the thickness of each zonation layer ###
        #             xdata_reg = np.linspace(0,1,50)
        #             # interpolate values along xdata_reg
        #             interpolant = interp1d(xdata, Udata)
        #             Udata_interp = interpolant(xdata_reg)
        #             interpolant = interp1d(xdata, Thdata)
        #             Thdata_interp = interpolant(xdata_reg)
        #             depthLayer = xdata_reg[1:] - xdata_reg[0:-1]
        #             depthLayer = np.array(depthLayer)*gs, #distance in µm
        #             depthLayer = depthLayer[0]
     
        #             # Get the mean concentration for each layer of the grain
        #             eU_array = np.zeros((len(depthLayer),))
        #             for k in range(len(xdata_reg)-1): # Loop over interpolation array
        #                 left = xdata_reg[k]
        #                 right = xdata_reg[k+1]
        #                 extracted_U = []
        #                 extracted_Th = []
        #                 for l in range(len(xdata)):# Loop over U and Th profiles
        #                     if xdata[l] > right:
        #                         break
        #                     if left < xdata[l] < right: # value in profiles between interpolation interval
        #                         extracted_U.append(Udata[l])
        #                         extracted_Th.append(Thdata[l])
        #                 # Add the y value of the limits of the interval
        #                 extracted_U.append(Udata_interp[k])
        #                 extracted_U.append(Udata_interp[k+1])
        #                 extracted_Th.append(Thdata_interp[k])
        #                 extracted_Th.append(Thdata_interp[k+1])
        #                 # Compute the means for each interval
        #                 mean_U = sum(extracted_U) / len(extracted_U)
        #                 mean_Th = sum(extracted_Th) / len(extracted_Th)
        #                 eU_array[k] = mean_U + 0.235 * mean_Th
        
        #             # Normalize to its maximal value
        #             eU_array = eU_array / np.max(eU_array)
       
        #             # Write the profiles in the file
        #             for j in range(len(eU_array)):
        #                 file.write('\n')
        #                 file.write(str(round(depthLayer[j],3)) + '\t')
        #                 file.write(str(round(eU_array[j],3)))
                        
        #             count += 1
                    
        # Write observed ages
        file.close()
        ### End of writing sample_setting file ####
        
    #----------------------------------------------------------
    def save_obs_file(self):
        
        """
            Create file for observations.
            Observations always have at minimum:
                sample name,LON,LAT,HEIGHT
                
            other observations include:
                - ages and error
                - grain size
                - U and Th
                - rmr0
                - 43 He data
            
            All the data are stored in xarray dataset.
        """
        
        # Get the name of the directory where observations will be stored
        try:
            Fname = self.parent.input_parameters['data_folder']
        except KeyError:
            Fname = 'Samples'
            pgui.store_input(self.parent, {'data_folder': Fname})
            
        # Set the path
        name = os.path.join(self.PFolder,"data",Fname)
        if os.path.exists(name):
            pass
        else:
            os.mkdir(name)
        name = os.path.join(self.PFolder, "data",Fname,Fname+".csv")
        
        # Get data IDs and ages/error of observations
        # Set dictionary for the dataset
        # All these data hae the same dimension
        dictionary = {}
        Headers = ['SAMPLE','LON','LAT','HEIGHT']
        # Get sample name and coordinates
        count = 0
        # Each thermochronometer can have different amount of grains per sample
        # We need to iterate through thermochronometers and check maximum number of grain per sample
        # The nb grains per sample is recorded in the self.parent.Samples dictionary that reads the samples table
        SAMPLE = ([self.parent.Samples['SampleName'+str(i+1)]+'_'+str(j+1) for i in range(int(self.parent.nbSampleValue.text()))
                   for j in range(int(self.parent.Samples['MaxNbGrain_'+str(i+1)]))])
                   
        LON = ([float(self.parent.SampleTable.item(i, 1).text()) for i in range(int(self.parent.nbSampleValue.text()))
                for j in range(int(self.parent.Samples['MaxNbGrain_'+str(i+1)]))])
                   
        LAT = ([float(self.parent.SampleTable.item(i, 2).text()) for i in range(int(self.parent.nbSampleValue.text()))
                for j in range(int(self.parent.Samples['MaxNbGrain_'+str(i+1)]))])
        HEIGHT = ([float(self.parent.SampleTable.item(i, 3).text()) for i in range(int(self.parent.nbSampleValue.text()))
                for j in range(int(self.parent.Samples['MaxNbGrain_'+str(i+1)]))])
        
        dictionary['SAMPLE'] = (['x'], SAMPLE)
        dictionary['LON'] = (['x'], LON)
        dictionary['LAT'] = (['x'], LAT)
        dictionary['HEIGHT'] = (['x'], HEIGHT)
        
        # for i in range(int(self.parent.nbSampleValue.text())):
        #     nb_grains = self.parent.Samples['ngrains'+str(i+1)]
        #     for j in range(int(nb_grains)):
        #         # Write sample names
        #         file.write(self.parent.Samples['SampleName'+str(i+1)]+'_'+str(j+1)+",")
        #         file.write(str(self.parent.SampleTable.item(i, 1).text()) + ',')
        #         file.write(str(self.parent.SampleTable.item(i, 2).text()) + ',')
        #         file.write(str(self.parent.SampleTable.item(i, 3).text()) + ',')
        #         for t in range(count_thermo):
        #             [file.write(str(ObservedData[t,count])+',') if float(ObservedData[t,count])>0.0 else file.write(' ,')]
        #             [file.write(str(ObservedError[t,count])+',') if float(ObservedData[t,count])>0.0 else count if float(ObservedData[t,count])==9999 else file.write(" ,")]
        #         count += 1
        #         file.write('\n')
        
        # Set thermochronometers flag
        if int(self.parent.Samples['TL_TotGrain']):
            self.parent.input_parameters['age_TL_flag'] = '1'
        else:
            self.parent.input_parameters['age_TL_flag'] = '0'
        if int(self.parent.Samples['OSL_TotGrain']):
            self.parent.input_parameters['age_OSL_flag'] = '1'
        else:
            self.parent.input_parameters['age_OSL_flag'] = '0'
        if int(self.parent.Samples['AHe_TotGrain']):
            self.parent.input_parameters['age_AHe_flag'] = '1'
        else:
            self.parent.input_parameters['age_AHe_flag'] = '0'
        if int(self.parent.Samples['AFT_TotGrain']):
            self.parent.input_parameters['age_AFT_flag'] = '1'
        else:
            self.parent.input_parameters['age_AFT_flag'] = '0'
        if int(self.parent.Samples['ZHe_TotGrain']):
            self.parent.input_parameters['age_ZHe_flag'] = '1'
        else:
            self.parent.input_parameters['age_ZHe_flag'] = '0'
        if int(self.parent.Samples['ZFT_TotGrain']):
            self.parent.input_parameters['age_ZFT_flag'] = '1'
        else:
            self.parent.input_parameters['age_ZFT_flag'] = '0'
        if int(self.parent.Samples['KAr_TotGrain']):
            self.parent.input_parameters['age_KAr_flag'] = '1'
        else:
            self.parent.input_parameters['age_KAr_flag'] = '0'
        if int(self.parent.Samples['BAr_TotGrain']):
            self.parent.input_parameters['age_BAr_flag'] = '1'
        else:
            self.parent.input_parameters['age_BAr_flag'] = '0'
        if int(self.parent.Samples['MAr_TotGrain']):
            self.parent.input_parameters['age_MAr_flag'] = '1'
        else:
            self.parent.input_parameters['age_MAr_flag'] = '0'
        if int(self.parent.Samples['HAr_TotGrain']):
            self.parent.input_parameters['age_HAr_flag'] = '1'
        else:
            self.parent.input_parameters['age_HAr_flag'] = '0'
        if int(self.parent.Samples['ESR_TotGrain']):
            self.parent.input_parameters['age_ESR_flag'] = '1'
        else:
            self.parent.input_parameters['age_ESR_flag'] = '0'
        
        count_thermo = 0
        Ahe_flag = 0
        Zhe_flag = 0
        Aft_flag = 0
        RMR0AHE = 0
        RMR0AFT = 0
        for km, vm in self.parent.input_parameters.items():
            # if km == 'age_TL_flag' and vm == '1':
            #     temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL')}
            #     TL = [float(temp['TL_'+s]) if 'TL_'+s in temp.keys() else None for s in SAMPLE]
            #     temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('DTL')}
            #     DTL = [float(temp['DTL_'+s]) if 'DTL_'+s in temp.keys() else None for s in SAMPLE]
            #     dictionary['TL'] = (['x'], TL)
            #     dictionary['DTL'] = (['x'], DTL)
            # elif km == 'age_OSL_flag' and vm == '1':
            #     temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL')}
            #     OSL = [float(temp['OSL_'+s]) if 'OSL_'+s in temp.keys() else None for s in SAMPLE]
            #     temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('DOSL')}
            #     DOSL = [float(temp['DOSL_'+s]) if 'DOSL_'+s in temp.keys() else None for s in SAMPLE]
            #     dictionary['OSL'] = (['x'], OSL)
            #     dictionary['DOSL'] = (['x'], DOSL)
            # elif km == 'age_ESR_flag' and vm == '1':
            #     temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR')}
            #     ESR = [float(temp['ESR_'+s]) if 'ESR_'+s in temp.keys() else None for s in SAMPLE]
            #     temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('DESR')}
            #     DESR = [float(temp['DESR_'+s]) if 'DESR_'+s in temp.keys() else None for s in SAMPLE]
            #     dictionary['ESR'] = (['x'], ESR)
            #     dictionary['DESR'] = (['x'], DESR)
            if km == 'age_AHe_flag' and vm == '1' :
                Ahe_flag = 1
                # take only ages
                temp = {k:v for k,v in self.apatiteHelium.grains.items() if k.startswith('AHE')}
                AHE = [float(temp['AHE_'+s]) if 'AHE_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.apatiteHelium.grains.items() if k.startswith('DAHE')}
                DAHE = [float(temp['DAHE_'+s]) if 'DAHE_'+s in temp.keys() else None for s in SAMPLE]
                # if zero -> none
                # DAHE = np.asarray([DAHE[i] if AHE[i] != 0.0 else None for i in range(len(AHE))])
                # AHE = np.asarray([AHE[i] if AHE[i] != 0.0 else None for i in range(len(AHE))])
                dictionary['AHE'] = (['x'], AHE)
                dictionary['DAHE'] = (['x'], DAHE)
                # Grain size
                temp = {k:v for k,v in self.apatiteHelium.grains.items() if k.startswith('ASIZE')}
                SAHE = [float(temp['ASIZE_'+s]) if 'ASIZE_'+s in temp.keys() else None for s in SAMPLE]
                # Uranium (ppm)
                temp = {k:v for k,v in self.apatiteHelium.grains.items() if k.startswith('AUPPM')}
                AUPPM = [float(temp['AUPPM_'+s]) if 'AUPPM_'+s in temp.keys() else None for s in SAMPLE]
                # Thorium (ppm)
                temp = {k:v for k,v in self.apatiteHelium.grains.items() if k.startswith('ATHPPM')}
                ATHPPM = [float(temp['ATHPPM_'+s]) if 'ATHPPM_'+s in temp.keys() else None for s in SAMPLE]
                # RMR0
                temp = {k:v for k,v in self.apatiteHelium.grains.items() if k.startswith(conf.Variable_names['FTL_kinetic_parameter_value_AHe'])}
                RMR0AHE = [float(temp[conf.Variable_names['FTL_kinetic_parameter_value_AHe']+'_'+s]) if conf.Variable_names['FTL_kinetic_parameter_value_AHe']+'_'+s in temp.keys() else None for s in SAMPLE]
                count_thermo += 1
            elif km == 'age_ZHe_flag' and vm == '1':
                Zhe_flag = 1
                temp = {k:v for k,v in self.ZHeParam.grains.items() if k.startswith('ZHE')}
                ZHE = [float(temp['ZHE_'+s]) if 'ZHE_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.ZHeParam.grains.items() if k.startswith('DZHE')}
                DZHE = [float(temp['DZHE_'+s]) if 'DZHE_'+s in temp.keys() else None for s in SAMPLE]
                # DZHE = np.asarray([DZHE[i] if ZHE[i] != 0.0 else None for i in range(len(ZHE))])
                # ZHE = np.asarray([ZHE[i] if ZHE[i] != 0.0 else None for i in range(len(ZHE))])
                dictionary['ZHE'] = (['x'], ZHE)
                dictionary['DZHE'] = (['x'], DZHE)
                # Grain size
                temp = {k:v for k,v in self.ZHeParam.grains.items() if k.startswith('ZSIZE')}
                SZHE = [float(temp['ZSIZE_'+s]) if 'ZSIZE_'+s in temp.keys() else None for s in SAMPLE]
                # Uranium
                temp = {k:v for k,v in self.ZHeParam.grains.items() if k.startswith('ZUPPM')}
                ZUPPM = [float(temp['ZUPPM_'+s]) if 'ZUPPM_'+s in temp.keys() else None for s in SAMPLE]
                # Thorium
                temp = {k:v for k,v in self.ZHeParam.grains.items() if k.startswith('ZTHPPM')}
                ZTHPPM = [float(temp['ZTHPPM_'+s]) if 'ZTHPPM_'+s in temp.keys() else None for s in SAMPLE]
                count_thermo += 1
            elif km == 'age_AFT_flag' and vm == '1':
                Aft_flag=1
                temp = {k:v for k,v in self.AFTParam.grains.items() if k.startswith('AFT')}
                AFT = [float(temp['AFT_'+s]) if 'AFT_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.AFTParam.grains.items() if k.startswith('DAFT')}
                DAFT = [float(temp['DAFT_'+s]) if 'DAFT_'+s in temp.keys() else None for s in SAMPLE]
                # DAFT = np.asarray([None if AFT[i] == 0.0 else DAFT[i] for i in range(len(AFT))])
                # AFT = np.asarray([None if AFT[i] == 0.0 else AFT[i] for i in range(len(AFT))])
                dictionary['AFT'] = (['x'], AFT)
                dictionary['DAFT'] = (['x'], DAFT)
                # RMR0
                temp = {k:v for k,v in self.AFTParam.grains.items() if k.startswith(conf.Variable_names['FTL_kinetic_parameter_value_AFT'])}
                RMR0AFT = [float(temp[conf.Variable_names['FTL_kinetic_parameter_value_AFT']+'_'+s]) if conf.Variable_names['FTL_kinetic_parameter_value_AFT']+'_'+s in temp.keys() else None for s in SAMPLE]
                # Mean Track length
                temp = {k:v for k,v in self.AFTParam.grains.items() if k.startswith(conf.Variable_names['Mean Fission Track Length'])}
                MTLAFT = [float(temp[conf.Variable_names['Mean Fission Track Length']+'_'+s]) if conf.Variable_names['Mean Fission Track Length']+'_'+s in temp.keys() else None for s in SAMPLE]
                dictionary[conf.Variable_names['Mean Fission Track Length']] = (['x'], MTLAFT)
                temp = {k:v for k,v in self.AFTParam.grains.items() if k.startswith(conf.Variable_names['Mean Fission Track Length error'])}
                DMTLAFT = [float(temp[conf.Variable_names['Mean Fission Track Length error']+'_'+s]) if conf.Variable_names['Mean Fission Track Length error']+'_'+s in temp.keys() else None for s in SAMPLE]
                dictionary[conf.Variable_names['Mean Fission Track Length error']] = (['x'], DMTLAFT)
                count_thermo += 1
            elif km == 'age_ZFT_flag' and vm == '1':
                temp = {k:v for k,v in self.ZFTParam.grains.items() if k.startswith('ZFT')}
                ZFT = [float(temp['ZFT_'+s]) if 'ZFT_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.ZFTParam.grains.items() if k.startswith('DZFT_')}
                DZFT = [float(temp['DZFT_'+s]) if 'DZFT_'+s in temp.keys() else None for s in SAMPLE]
                dictionary['ZFT'] = (['x'], ZFT)
                dictionary['DZFT'] = (['x'], DZFT)
                count_thermo += 1
            elif km == 'age_KAr_flag' and vm == '1':
                temp = {k:v for k,v in self.KArParam.grains.items() if k.startswith('KAR')}
                KAR = [float(temp['KAR_'+s]) if 'KAR_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.KArParam.grains.items() if k.startswith('DKAR')}
                DKAR = [float(temp['DKAR_'+s]) if 'DKAR_'+s in temp.keys() else None for s in SAMPLE]
                dictionary['KAR'] = (['x'], KAR)
                dictionary['DKAR'] = (['x'], DKAR)
                count_thermo += 1
            elif km == 'age_BAr_flag' and vm == '1':
                temp = {k:v for k,v in self.BArParam.grains.items() if k.startswith('BAR')}
                BAR = [float(temp['BAR_'+s]) if 'BAR_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.BArParam.grains.items() if k.startswith('DBAR')}
                DBAR = [float(temp['DBAR_'+s]) if 'DBAR_'+s in temp.keys() else None for s in SAMPLE]
                dictionary['BAR'] = (['x'], BAR)
                dictionary['DBAR'] = (['x'], DBAR)
                count_thermo += 1
            elif km == 'age_MAr_flag' and vm == '1':
                temp = {k:v for k,v in self.MArParam.grains.items() if k.startswith('MAR')}
                MAR = [float(temp['MAR_'+s]) if 'MAR_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.MArParam.grains.items() if k.startswith('DMAR')}
                DMAR = [float(temp['DMAR_'+s]) if 'DMAR_'+s in temp.keys() else None for s in SAMPLE]
                dictionary['MAR'] = (['x'], MAR)
                dictionary['DMAR'] = (['x'], DMAR)
                count_thermo += 1
            elif km == 'age_HAr_flag' and vm == '1':
                temp = {k:v for k,v in self.HArParam.grains.items() if k.startswith('HAR')}
                HAR = [float(temp['HAR_'+s]) if 'HAR_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.HArParam.grains.items() if k.startswith('DHAR')}
                DHAR = [float(temp['DHAR_'+s]) if 'DHAR_'+s in temp.keys() else None for s in SAMPLE]
                dictionary['HAR'] = (['x'], HAR)
                dictionary['DHAR'] = (['x'], DHAR)
                count_thermo += 1

        # Extra parameters
        if Ahe_flag:
            if SAHE:
                dictionary['ASIZE'] = (['x'], SAHE)
            if AUPPM:
                dictionary['AUPPM'] = (['x'], AUPPM)
            if ATHPPM:
                dictionary['ATHPPM'] = (['x'], ATHPPM)
        if Zhe_flag:
            if SZHE:
                dictionary['ZSIZE'] = (['x'], SZHE)
            if ZUPPM:
                dictionary['ZUPPM'] = (['x'], ZUPPM)
            if ZTHPPM:
                dictionary['ZTHPPM'] = (['x'], ZTHPPM)

        # Write annealing kinetic for AFT and AHe
        if RMR0AFT:
            dictionary[conf.Variable_names['FTL_kinetic_parameter_value_AFT']] = (['x'], RMR0AFT)
        if RMR0AHE:
            dictionary[conf.Variable_names['FTL_kinetic_parameter_value_AHe']] = (['x'], RMR0AHE)
        # RMR0 = []
        # if RMR0AHE or RMR0AFT:
        #     if RMR0AHE == 0:
        #         RMR0 = [RMR0AFT[i] for i in range(len(SAMPLE))]
        #     elif RMR0AFT == 0:
        #         RMR0 = [RMR0AHE[i] for i in range(len(SAMPLE))]
        #     else:
        #         RMR0 = [RMR0AFT[i] if RMR0AFT[i] != None else RMR0AHE[i] for i in range(len(SAMPLE))]
        #     if RMR0:
        #         dictionary[conf.Variable_names['FTL_kinetic_parameter_value']] = (['x'], RMR0)
            
            
        # This is the dataset that stores the samples
        # coordinates and ages with error
        try:
            self.datasetAges = xr.Dataset(
                data_vars=dictionary,
                attrs = dict(description="Ages data")
                )
        except ValueError as VE:
            QErrorMessage(self).showMessage("Something wrong happened with the grain parameters:"+ str(VE))
            raise
            return
        
        # Write csv file - Age observation
        name = os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,Fname+".csv")
        try:
            df = self.datasetAges.to_dataframe().to_csv(name,index=False)
        except PermissionError:
            if os.path.exists(name):
                os.chmod(os.path.join(self.PFolder,"data",Fname),0o777)
                os.chmod(name,0o666)
                file = open(name, 'w+')
        except OSError:
            os.mkdir(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname))
            df = self.datasetAges.to_dataframe().to_csv(name,index=False)
           
        
        #### Is there any ThermoLuminescence/OSL/ESR data ? #######
        for km, vm in self.parent.input_parameters.items():
            if km == 'age_TL_flag' and vm == '1' :
                dictionary = {}
                SAMPLE = ([self.parent.Samples['SampleName'+str(i+1)]+'_'+str(j+1) for i in range(int(self.parent.nbSampleValue.text()))
                           for j in range(int(self.parent.Samples['TL_ngrains'+str(i+1)])) if int(self.parent.Samples['TL_ngrains'+str(i+1)]) > 0 ])
                
                LON = ([float(self.parent.SampleTable.item(i, 1).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['TL_ngrains'+str(i+1)])) if int(self.parent.Samples['TL_ngrains'+str(i+1)]) > 0 ])
                           
                LAT = ([float(self.parent.SampleTable.item(i, 2).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['TL_ngrains'+str(i+1)])) if int(self.parent.Samples['TL_ngrains'+str(i+1)]) > 0 ])
               
                HEIGHT = ([float(self.parent.SampleTable.item(i, 3).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['TL_ngrains'+str(i+1)])) if int(self.parent.Samples['TL_ngrains'+str(i+1)]) > 0 ])
                # take only ages
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_NN')}
                TL_NN = [float(temp['TL_NN_'+s]) if 'TL_NN_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('DTL_NN')}
                DTL_NN = [float(temp['DTL_NN_'+s]) if 'DTL_NN_'+s in temp.keys() else None for s in SAMPLE]
        
                # Doser
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_Doser')}
                TL_DOSER = [float(temp['TL_Doser_'+s]) if 'TL_Doser_'+s in temp.keys() else None for s in SAMPLE]
                # D0
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_D0')}
                D0 = [float(temp['TL_D0_'+s]) if 'TL_D0_'+s in temp.keys() else None for s in SAMPLE]
                # a
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_a')}
                TLA = [float(temp['TL_a_'+s]) if 'TL_a_'+s in temp.keys() else None for s in SAMPLE]
                # b
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_b')}
                TLB = [float(temp['TL_b_'+s]) if 'TL_b_'+s in temp.keys() else None for s in SAMPLE]
                # Et
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_Et')}
                TL_Et = [float(temp['TL_Et_'+s]) if 'TL_Et_'+s in temp.keys() else None for s in SAMPLE]
                # logs
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_logs')}
                TL_logs = [float(temp['TL_logs_'+s]) if 'TL_logs_'+s in temp.keys() else None for s in SAMPLE]
                # rhop
                temp = {k:v for k,v in self.TLParam.grains.items() if k.startswith('TL_logrho_')}
                TL_logrho = [float(temp['TL_logrho_'+s]) if 'TL_logrho_'+s in temp.keys() else None for s in SAMPLE]
                
                # Save in dictionary
                dictionary['SAMPLETL'] = (['x'], SAMPLE)
                dictionary['LON'] = (['x'], LON)
                dictionary['LAT'] = (['x'], LAT)
                dictionary['HEIGHT'] = (['x'], HEIGHT)
                dictionary['DOSER'] = (['x'], TL_DOSER)
                dictionary['D0'] = (['x'], D0)
                dictionary['ATL'] = (['x'], TLA)
                dictionary['BTL'] = (['x'], TLB)
                dictionary['ET'] = (['x'], TL_Et)
                dictionary['LOGS'] = (['x'], TL_logs)
                dictionary['LOGRHO'] = (['x'], TL_logrho)
                dictionary['N/N'] = (['x'], TL_NN)
                dictionary['DNN'] = (['x'], DTL_NN)
                
                # This is the dataset that stores the data
                try:
                    self.dataset = xr.Dataset(
                        data_vars=dictionary,
                        attrs = dict(description="ThL data")
                        )
                except ValueError as VE:
                    QErrorMessage(self).showMessage("Something wrong happened with ThL data.")
                    return
                 
                # Write csv file - Th.L observation
                name = os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,"ThL_data.csv")
                try:
                    df = self.dataset.to_dataframe().to_csv(name,index=False)
                except PermissionError:
                    if os.path.exists(name):
                        os.remove(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,"ThL_data.csv"))
                        # os.chmod(os.path.join(self.PFolder,"data",Fname),0o777)
                        # os.chmod(name,0o666)
                        file = open(name, 'w+')
                except OSError:
                    os.mkdir(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname))
                    df = self.dataset.to_dataframe().to_csv(name,index=False)
                    
            ############# OSL data ? ########################
            if km == 'age_OSL_flag' and vm == '1' :
                dictionary = {}
                SAMPLE = ([self.parent.Samples['SampleName'+str(i+1)]+'_'+str(j+1) for i in range(int(self.parent.nbSampleValue.text()))
                           for j in range(int(self.parent.Samples['OSL_ngrains'+str(i+1)])) if int(self.parent.Samples['OSL_ngrains'+str(i+1)]) > 0 ])
                
                LON = ([float(self.parent.SampleTable.item(i, 1).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['OSL_ngrains'+str(i+1)])) if int(self.parent.Samples['OSL_ngrains'+str(i+1)]) > 0 ])
                           
                LAT = ([float(self.parent.SampleTable.item(i, 2).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['OSL_ngrains'+str(i+1)])) if int(self.parent.Samples['OSL_ngrains'+str(i+1)]) > 0 ])
               
                HEIGHT = ([float(self.parent.SampleTable.item(i, 3).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['OSL_ngrains'+str(i+1)])) if int(self.parent.Samples['OSL_ngrains'+str(i+1)]) > 0 ])
                # take only ages
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_NN')}
                OSL_NN = [float(temp['OSL_NN_'+s]) if 'OSL_NN_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_DNN')}
                OSL_DNN = [float(temp['OSL_DNN_'+s]) if 'OSL_DNN_'+s in temp.keys() else None for s in SAMPLE]
        
                # Doser
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_Doser')}
                OSL_DOSER = [float(temp['OSL_Doser_'+s]) if 'OSL_Doser_'+s in temp.keys() else None for s in SAMPLE]
                # D0
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_D0')}
                D0 = [float(temp['OSL_D0_'+s]) if 'OSL_D0_'+s in temp.keys() else None for s in SAMPLE]
                # Eu or GOK_b
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_Eu')}
                Eu = [float(temp['OSL_Eu_'+s]) if 'OSL_Eu_'+s in temp.keys() else None for s in SAMPLE]
                # Et
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_Et')}
                OSL_Et = [float(temp['OSL_Et_'+s]) if 'OSL_Et_'+s in temp.keys() else None for s in SAMPLE]
                # Logs
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_logs')}
                OSL_logs = [float(temp['OSL_logs_'+s]) if 'OSL_logs_'+s in temp.keys() else None for s in SAMPLE]
                # logrho
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_logrho_')}
                OSL_logrho = [float(temp['OSL_logrho_'+s]) if 'OSL_logrho_'+s in temp.keys() else None for s in SAMPLE]
                # a coefficient
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_a_')}
                OSL_a = [float(temp['OSL_a_'+s]) if 'OSL_a_'+s in temp.keys() else None for s in SAMPLE]
                # imax coefficient
                temp = {k:v for k,v in self.OSLParam.grains.items() if k.startswith('OSL_Lmax_')}
                OSL_Lmax = [float(temp['OSL_Lmax_'+s]) if 'OSL_Lmax_'+s in temp.keys() else None for s in SAMPLE]
                
                # Save in dictionary
                dictionary['SAMPLEOSL'] = (['x'], SAMPLE)
                dictionary['LON'] = (['x'], LON)
                dictionary['LAT'] = (['x'], LAT)
                dictionary['HEIGHT'] = (['x'], HEIGHT)
                dictionary['DOSER'] = (['x'], OSL_DOSER)
                dictionary['D0'] = (['x'], D0)
                dictionary['EU'] = (['x'], Eu)
                dictionary['ET'] = (['x'], OSL_Et)
                dictionary['LOGS'] = (['x'], OSL_logs)
                dictionary['LOGRHO'] = (['x'], OSL_logrho)
                dictionary['N/N'] = (['x'], OSL_NN)
                dictionary['DNN'] = (['x'], OSL_DNN)
                dictionary['a'] = (['x'], OSL_a)
                dictionary['Lmax'] = (['x'], OSL_Lmax)
                
                # This is the dataset that stores the data
                try:
                    self.dataset = xr.Dataset(
                        data_vars=dictionary,
                        attrs = dict(description="OSL data")
                        )
                except ValueError as VE:
                    QErrorMessage(self).showMessage("Something wrong happened with ThL data.")
                    return
                 
                # Write csv file - Th.L observation
                name = os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,"OSL_data.csv")
                try:
                    df = self.dataset.to_dataframe().to_csv(name,index=False)
                except PermissionError:
                    if os.path.exists(name):
                        os.remove(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,"OSL_data.csv"))
                        # os.chmod(os.path.join(self.PFolder,"data",Fname),0o777)
                        # os.chmod(name,0o666)
                        file = open(name, 'w+')
                except OSError:
                    os.mkdir(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname))
                    df = self.dataset.to_dataframe().to_csv(name,index=False)
            
            ############# ESR data ? ########################
            if km == 'age_ESR_flag' and vm == '1' :
                dictionary = {}
                SAMPLE = ([self.parent.Samples['SampleName'+str(i+1)]+'_'+str(j+1) for i in range(int(self.parent.nbSampleValue.text()))
                           for j in range(int(self.parent.Samples['ESR_ngrains'+str(i+1)])) if int(self.parent.Samples['ESR_ngrains'+str(i+1)]) > 0 ])
                
                LON = ([float(self.parent.SampleTable.item(i, 1).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['ESR_ngrains'+str(i+1)])) if int(self.parent.Samples['ESR_ngrains'+str(i+1)]) > 0 ])
                           
                LAT = ([float(self.parent.SampleTable.item(i, 2).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['ESR_ngrains'+str(i+1)])) if int(self.parent.Samples['ESR_ngrains'+str(i+1)]) > 0 ])
               
                HEIGHT = ([float(self.parent.SampleTable.item(i, 3).text()) for i in range(int(self.parent.nbSampleValue.text()))
                        for j in range(int(self.parent.Samples['ESR_ngrains'+str(i+1)])) if int(self.parent.Samples['ESR_ngrains'+str(i+1)]) > 0 ])
                # take only ages
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_NN')}
                ESR_NN = [float(temp['ESR_NN_'+s]) if 'ESR_NN_'+s in temp.keys() else None for s in SAMPLE]
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_DNN')}
                ESR_DNN = [float(temp['ESR_DNN_'+s]) if 'ESR_DNN_'+s in temp.keys() else None for s in SAMPLE]
        
                # Doser
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_Doser')}
                ESR_DOSER = [float(temp['ESR_Doser_'+s]) if 'ESR_Doser_'+s in temp.keys() else None for s in SAMPLE]
                # D0
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_D0')}
                D0 = [float(temp['ESR_D0_'+s]) if 'ESR_D0_'+s in temp.keys() else None for s in SAMPLE]
                # sigma Et
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_sigmaEt')}
                sigmaEt = [float(temp['ESR_sigmaEt_'+s]) if 'ESR_sigmaEt_'+s in temp.keys() else None for s in SAMPLE]
                # a
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_a')}
                GOK_a = [float(temp['ESR_a_'+s]) if 'ESR_a_'+s in temp.keys() else None for s in SAMPLE]
                # b
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_b')}
                GOK_b = [float(temp['ESR_b_'+s]) if 'ESR_b_'+s in temp.keys() else None for s in SAMPLE]
                # Et
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_Et')}
                ESR_Et = [float(temp['ESR_Et_'+s]) if 'ESR_Et_'+s in temp.keys() else None for s in SAMPLE]
                # Logs
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_logs')}
                ESR_logs = [float(temp['ESR_logs_'+s]) if 'ESR_logs_'+s in temp.keys() else None for s in SAMPLE]
                # Lmax
                temp = {k:v for k,v in self.ESRParam.grains.items() if k.startswith('ESR_Lmax')}
                ESR_Lmax = [float(temp['ESR_Lmax_'+s]) if 'ESR_Lmax_'+s in temp.keys() else None for s in SAMPLE]
                
                # Save in dictionary
                dictionary['SAMPLEESR'] = (['x'], SAMPLE)
                dictionary['LON'] = (['x'], LON)
                dictionary['LAT'] = (['x'], LAT)
                dictionary['HEIGHT'] = (['x'], HEIGHT)
                dictionary['DOSER'] = (['x'], ESR_DOSER)
                dictionary['D0'] = (['x'], D0)
                dictionary['SIGMAET'] = (['x'], sigmaEt)
                dictionary['ET'] = (['x'], ESR_Et)
                dictionary['LOGS'] = (['x'], ESR_logs)
                dictionary['N/N'] = (['x'], ESR_NN)
                dictionary['DNN'] = (['x'], ESR_DNN)
                dictionary['a'] = (['x'], GOK_a)
                dictionary['b'] = (['x'], GOK_b)
                dictionary['Lmax'] = (['x'], ESR_Lmax)
                
                # This is the dataset that stores the data
                try:
                    self.dataset = xr.Dataset(
                        data_vars=dictionary,
                        attrs = dict(description="ESR data")
                        )
                except ValueError as VE:
                    QErrorMessage(self).showMessage("Something wrong happened with ThL data.")
                    return
                 
                # Write csv file - Th.L observation
                name = os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,"ESR_data.csv")
                try:
                    df = self.dataset.to_dataframe().to_csv(name,index=False)
                except PermissionError:
                    if os.path.exists(name):
                        os.remove(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname,"ESR_data.csv"))
                        # os.chmod(os.path.join(self.PFolder,"data",Fname),0o777)
                        # os.chmod(name,0o666)
                        file = open(name, 'w+')
                except OSError:
                    os.mkdir(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",Fname))
                    df = self.dataset.to_dataframe().to_csv(name,index=False)
            

        
        
#############################################################
# class GrainsParameters(QWidget):
#     """
#       This class contains the grains characteristics for He computation
#       it also show zonation if asked by the user.
    
#     """
#     def __init__(self, parent, Param, nSamples,Samples,system):
#         super().__init__()
        
#         # To store the grain characteristics
#         self.parent = parent
#         self.old_input = 0 #old input file signal
#         self.toolbarList = []
#         self.Param = Param
        
#         # Define Splitters
#         self.GrainSplitter = pgu.QCustomSplitter(Qt.Vertical)
#         self.TableSplit = QWidget()
#         self.ZonationSplit = QWidget()
#         self.ZonationSplit.setEnabled(False)
        
#         #### Parameters #####
#         # Table
#         if not nSamples: #if number of sample == 0, force to 1 at least
#             nSamples = 1
#         self.nbSampleValue = nSamples
#         self.Samples = Samples
#         self.Param = Param
#         self.Labelw = QLabel('Set grains characteristics:')
#         self.Labelw.setFont(conf.fontLine12)
#         self.Labelw.setAlignment(Qt.AlignLeft)
#         data = {" ID": ['1']*int(Param.DParameters['ngrains']),
#                 "Size (µm)": ["0"]*int(Param.DParameters['ngrains']),
#                 "U (ppm)": ["0"]*int(Param.DParameters['ngrains']),
#                 "Th (ppm)": ["0"]*int(Param.DParameters['ngrains']),
#                 "rmr0 ": ["0"]*int(Param.DParameters['ngrains'])}
#         # Zonation
#         self.ZonationCheckBox = QCheckBox("Zonation")
#         self.ZonationCheckBox.setEnabled(False)
#         if Param.DParameters['zonation_flag'] == '1':
#             self.old_input = 1
#             self.ZonationCheckBox.setChecked(True)
#             self.ZonationSplit.setEnabled(True)
#         self.ZonationLabel = QLabel("Set zonations:")
#         self.ZonationLabel.setFont(conf.fontLine12)
#         self.grainListLabel = QLabel('Select the grain:')
#         self.grainListLabel.setFont(conf.font8)
#         self.grainComboBox = QComboBox()
#         GrainList = (['grain'+str(j+1)+'_'+str(i+1)+'  ' for j in range(int(self.nbSampleValue.text()))
#                       for i in range(int(self.Samples['AHe_ngrains'+str(j+1)]))])
#         self.grainComboBox.addItems(GrainList)
#         self.UserDescription = QLabel("Use the keyboard to set zonation for:\n'u': uranium, 't': thorium, 'r': restore both curves.")
#         self.SaveButton = QPushButton("save")
#         self.SaveButton.clicked.connect(lambda: self.saveZonation())
#         self.setNbLayersLabel = QLabel("Number of layers: ")
#         self.setNbLayersLabel.setFont(conf.font8)
#         self.setNbLayersLabel.setToolTip("Concentrations are averaged in each layer")
#         self.setNbLayers = QLineEdit()
#         self.setNbLayers.setValidator(QIntValidator(0, 50, self))
#         self.setNbLayers.setText(Param.DParameters['nLayers'])
        
#         # Check for user update
#         self.grainComboBox.currentIndexChanged.connect(
#             lambda: self.update_ComboBox(self.grainComboBox))
#         self.ZonationCheckBox.stateChanged.connect(
#             lambda: self.zonationState(self.ZonationCheckBox))
        
#         # How many grains in total?
#         self.TotGrain = 0
#         for j in range(int(self.nbSampleValue.text())):
#                 nb_grains = self.Samples['AHe_ngrains'+str(j+1)]
#                 self.TotGrain += int(nb_grains)
        
#         self.GrainTable = U_interface.WinTable(data, self.TotGrain, 5)
#         self.updateTableGrain(self.GrainTable, 1, Param, self.parent)
#         self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
#             self.GrainTable, 0, Param, self.parent))
#         # self.parent.parent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
#         #     self.GrainTable, 2, Param, self.parent))
        
#         # Relative boundaries of the grain
#         self.x = [0.0, 1.0]
        
#         #Get and store zonation plot
#         # self.getZonationPlots()
        
#         # Add toolbar
#         # self.toolbar = NavigationToolbar(self.ZonationPlot,self)
        
#         #### Build the layouts ####
#         # Table split
#         VBox1 = QVBoxLayout()
#         VBox2 = QVBoxLayout()
#         self.Grid = QGridLayout()
#         self.Grid.setSpacing(10)
#         self.Grid.addWidget(self.ZonationCheckBox, 2, 0)
#         self.Grid.addWidget(self.GrainTable, 3, 0, 1, 4)
#         self.Grid.setRowStretch(3,5)
#         VBox1.addWidget(self.Labelw)
#         VBox1.addLayout(self.Grid)
 
#         # Zonation Split
#         VBox2.addWidget(self.ZonationLabel)
#         VBox2.addStretch(1)
#         Grid = QGridLayout()
#         Grid.addWidget(self.grainListLabel, 0, 1)
#         Grid.addWidget(self.grainComboBox, 0, 2, 1, 1)
#         Grid.addWidget(self.SaveButton, 0, 4)
#         Grid.addWidget(self.setNbLayersLabel, 1, 1)
#         Grid.addWidget(self.setNbLayers, 1, 2)
#         Grid.addWidget(self.UserDescription, 2, 1, 1, 4)
#         Grid.setColumnStretch(6,1)
#         VBox2.addLayout(Grid)
#         VBox2.addStretch(1)
#         # VBox2.addWidget(self.toolbar)
#         # VBox2.addWidget(self.ZonationPlot)
#         VBox2.addStretch(1)
        
#         # Add Layouts to widgets
#         self.TableSplit.setLayout(VBox1)
#         # self.ZonationSplit.setLayout(VBox2)
        
#         # Add widgets to splitter
#         self.GrainSplitter.addWidget(self.TableSplit)
#         # self.GrainSplitter.addWidget(self.ZonationSplit)
        
#         # Add splitter to layout
#         self.Layout = QVBoxLayout()
#         self.Layout.addWidget(self.GrainSplitter)
        
#     #----------------------------------------------------------
#     def zonationState(self,checkBox):
#         """ Enable/disable the zonation area """
#         if checkBox.isChecked():
#             self.ZonationSplit.setEnabled(True)
#             self.parent.zonation = 1
#         else:
#             self.ZonationSplit.setEnabled(False)
#             self.parent.zonation = 0
            
#     #----------------------------------------------------------
#     def getZonationPlots(self):
#         """ Gets and stores the zonation plot for each grain"""
#         if len(self.parent.zonationUxList) == 1 or not self.parent.zonationUxList: 
#             # means first time the user open the window, otherwise the zonation profile are
#             # loaded from previous opening
#             # Store the default values. Each grain has two lists, for U and Th values respectively.
#             count = 0
#             for j in range(int(self.nbSampleValue.text())):
#                 nb_grains = self.Samples['AHe_ngrains'+str(j+1)]
#                 for i in range(int(nb_grains)):
#                     item = self.GrainTable.item(count, 0)
#                     if item is None:
#                         break
#                     sampleName = item.text()
#                     Ug = float(self.grains['AUPPM_'+sampleName])
#                     Th = float(self.grains['ATHPPM_'+sampleName])
#                     self.parent.zonationUxList.append(self.x)
#                     self.parent.zonationUyList.append([Ug, Ug])
#                     self.parent.zonationThxList.append(self.x)
#                     self.parent.zonationThyList.append([Th, Th])  
#                     count += 1
#         if len(self.parent.zonationUxList) > 0:
#             # Create zonation profiles for the first grain
#             self.ZonationPlot = pgu.InteractivePlot(Uxdata=self.parent.zonationUxList[0],
#                                                 Uydata=self.parent.zonationUyList[0],
#                                                 Thxdata= self.parent.zonationThxList[0],
#                                                 Thydata=self.parent.zonationThyList[0])
#         else:
#             self.ZonationPlot = pgu.InteractivePlot(Uxdata=[0],
#                                                 Uydata=[10],
#                                                 Thxdata= [0],
#                                                 Thydata=[10])
                
#     #--------------------------------------------------------- 
#     def update_ComboBox(self, cb):
#         # Get index of the current label shown (grain number)
#         index = cb.currentIndex()
        
#         # First clear the axes, ax1 is for uranium, ax2 is for thorium
#         self.ZonationPlot.ax1.clear()
#         self.ZonationPlot.ax2.clear()
#         # then clear the previous list of stored points
#         self.ZonationPlot.listPoints = []
#         # then update zonation plot
#         self.ZonationPlot.Uxdata = np.array(self.parent.zonationUxList[index])
#         self.ZonationPlot.Uydata = np.array(self.parent.zonationUyList[index])
#         self.ZonationPlot.Thxdata = np.array(self.parent.zonationThxList[index])
#         self.ZonationPlot.Thydata = np.array(self.parent.zonationThyList[index])
#         # Create the plots
#         self.ZonationPlot.create_draggable_points()
#         # Restore the grids
#         self.ZonationPlot.ax1.minorticks_on()
#         # self.ZonationPlot.ax1.tick_params(axis='y', labelcolor='red')
#         self.ZonationPlot.ax1.grid(b=True,which='major', color='0.8', linestyle='-')
#         self.ZonationPlot.ax1.grid(b=True,which='minor', color='0.5', linestyle='--')
#         # self.ZonationPlot.ax2.tick_params(axis='y', labelcolor='blue')
#         self.ZonationPlot.ax2.grid(b=True,which='major', color='0.8', linestyle='-')
#         self.ZonationPlot.ax2.grid(b=True,which='minor', color='0.5', linestyle='--')
#         self.ZonationPlot.ax1.set_xlim((0,1))
#         self.ZonationPlot.ax1.set_ylim((0,max(self.ZonationPlot.Uydata)+30))
#         self.ZonationPlot.ax2.set_ylim((0,max(self.ZonationPlot.Thydata)+30))
#         self.ZonationPlot.ax1.set_xlabel('edge <- Grain distance -> core')
#         self.ZonationPlot.ax1.set_ylabel('U (ppm)', color = 'red')
#         self.ZonationPlot.ax2.set_ylabel('Th (ppm)', color = 'blue')
#         self.ZonationPlot.ax2.set_xlabel('edge <- Grain distance -> core')
#         # update the plots in the interface
#         self.ZonationPlot.draw()
    
#     #--------------------------------------------------------
#     def saveZonation(self):
#         # save the new plots, possibly changed manually by the user
#         index = self.grainComboBox.currentIndex()
#         # Uranium
#         self.parent.zonationUxList[index] = self.ZonationPlot.listPoints[0].line_object[0].get_xdata()
#         self.parent.zonationUyList[index] = self.ZonationPlot.listPoints[0].line_object[0].get_ydata()
#         # Thorium
#         self.parent.zonationThxList[index] = self.ZonationPlot.listPoints[1].line_object[0].get_xdata()
#         self.parent.zonationThyList[index] = self.ZonationPlot.listPoints[1].line_object[0].get_ydata()
    
#     #--------------------------------------------------------
#     def updateTableGrain(self, t, h, p, d):
#         """ 
#           To store input value for the Grain table
#           t = table
#           h : parameter
#           p : parameters
#           d : pecube parameters (window)
#         """
#         self.grains = {}
#         h_temp = 1
#         if h == 2: # Update from Table Samples
#             h_temp = 2
#             try:
#                 self.GrainTable.cellChanged.disconnect() # Shut down update from the user
#             except TypeError:
#                 pass
#             if self.parent.parent.SampleTable.column(self.parent.parent.SampleTable.currentItem()) == 4:   
#                 h = 1
#         t.setRowCount(self.TotGrain)
#         if h:
#             # dParam = DefaultParameterValues()
#             count = 0
#             for j in range(int(self.nbSampleValue.text())):
#                 nb_grains = int(self.Samples['AHe_ngrains'+str(j+1)])
#                 if nb_grains > 0:
#                     for i in range(int(nb_grains)):
#                         # Grain ID
#                         if "SampleName"+str(j+1) in p.TableParameters:
#                             grainID = p.TableParameters['SampleName'+str(j+1)]+"_"+str(i+1)
#                             t.setItem(count, 0, QTableWidgetItem(grainID))
#                         else:
#                             # The sample ID is found in class 'WindowParameters' in PecubeGUI.py
#                             # so go back to parent (AgesParameter) then parent again (WindowParameters)
#                             grainID = self.parent.parent.Samples['SampleName'+str(j+1)]+"_"+str(i+1)
#                             t.setItem(count,0, QTableWidgetItem(grainID))
#                         if "ASIZE_"+grainID in p.TableParameters:
#                             t.setItem(count, 1, QTableWidgetItem(
#                             str(p.TableParameters['ASIZE_'+grainID])))
#                         else:
#                             t.setItem(count, 1, QTableWidgetItem(
#                                 p.DParameters['grainradius']))
#                         if "AUPPM_"+grainID in p.TableParameters:
#                             t.setItem(count, 2, QTableWidgetItem(
#                                 str(p.TableParameters['AUPPM_'+grainID])))
#                         else:
#                             t.setItem(count, 2, QTableWidgetItem(
#                                 p.DParameters['Uppm']))
#                         if "ATHPPM_"+grainID in p.TableParameters:
#                             t.setItem(count, 3, QTableWidgetItem(
#                                 str(p.TableParameters['ATHPPM_'+grainID])))
#                         else:
#                             t.setItem(count, 3, QTableWidgetItem(
#                                 p.DParameters['Thppm']))
#                         if "RMR0_"+grainID in p.TableParameters:
#                             t.setItem(count, 4, QTableWidgetItem(
#                                 str(p.TableParameters['RMR0_'+grainID])))
#                         else:
#                             t.setItem(count, 4, QTableWidgetItem(
#                                 p.DParameters['rmr0_spec']))
#                         count += 1
                    
#         num_rows, num_cols = t.rowCount(), t.columnCount()
#         temp_dict = {}
#         temp_GrainID = []
#         try:
#             for row in range(num_rows):
#                 for col in range(num_cols):
#                     item = t.item(row, col)
#                     if item is None:
#                         break
#                     if col == 0:
#                         sampleName = item.text()
#                         temp_GrainID.append(sampleName)
#                     if col == 1:
#                         value = float(item.text())
#                         temp_dict['ASIZE_'+sampleName] = str(item.text())
#                         p.TableParameters['ASIZE_'+sampleName] = str(item.text())
#                     elif col == 2:
#                         value = float(item.text())
#                         temp_dict['AUPPM_'+sampleName] = str(item.text())
#                         p.TableParameters['AUPPM_'+sampleName] = str(item.text())
#                     elif col == 3:
#                         value = float(item.text())
#                         temp_dict['ATHPPM_'+sampleName] = str(item.text())
#                         p.TableParameters['ATHPPM_'+sampleName] = str(item.text())
#                     elif col == 4:
#                         value = float(item.text())
#                         temp_dict['RMR0_'+sampleName] = str(item.text())
#                         p.TableParameters['RMR0_'+sampleName] = str(item.text())
                        
#                 if item is None:
#                     break
#         except ValueError as VE:
#             print("Error", VE)
#             QErrorMessage(self).showMessage("In Grain table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
#             temp_dict['GrainIDs'] = temp_GrainID
#             self.grains.update(temp_dict)
#             return
#         temp_dict['GrainIDs'] = temp_GrainID
#         self.grains.update(temp_dict)
#         # Update AFT rmr0
#         if self.parent.AFTParam:
#             try:
#                 self.parent.AFTParam.GrainTable.cellChanged.disconnect()
#                 self.parent.AFTParam.updateTableGrain(self.parent.AFTParam.GrainTable, 1, p)
#                 self.parent.AFTParam.GrainTable.cellChanged.connect(lambda:self.parent.AFTParam.updateTableGrain(self.parent.AFTParam.GrainTable, 0, p))
#             except:
#                 print('Line 3092 in Thermochronometers.py: disconnect() failed')
                
#         # We write U, Th for all nodes
#         if self.parent.parent.AgeCombo.currentIndex() == 1:
#             for k, v in self.grains.items():
#                 if k.startswith('AUPPM'):
#                     pgui.store_input(self.parent.parent, {'AUPPM': str(v)})
#                 elif k.startswith('ATHPPM'):
#                     pgui.store_input(self.parent.parent, {'AThPPM': str(v)})
#                 elif k.startswith('ASIZE'):
#                     pgui.store_input(self.parent.parent, {'ASIZE': str(v)})
#                 elif k.startswith('RMR0'):
#                     pgui.store_input(self.parent.parent, {'rmr0_spec': str(v)})
                    
#         if h == 2: #update from Table Samples
#             self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
#                 self.GrainTable, 0, self.Param))
        

#         # Store the default values for zonation plot.
#         # Each grain has two lists, for U and Th values respectively.
#         if not h: #means the user changes the cell of the table manually
#             #means first time the user open the window, otherwise the zonation profile are
#             #loaded from previous opening
#             self.parent.zonationUxList = []
#             self.parent.zonationUyList = []
#             self.parent.zonationThxList = []
#             self.parent.zonationThyList = []
#             for i in range(int(self.nbSampleValue.text())):
#                 nb_grains = self.Samples['AHe_ngrains'+str(i+1)]
#                 for j in range(int(nb_grains)):
#                     Ug = float(self.grains['AUPPM_'+temp_GrainID[j]])
#                     Th = float(self.grains['ATHPPM_'+temp_GrainID[j]])
#                     self.parent.zonationUxList.append([0.0, 1.0])
#                     self.parent.zonationUyList.append([Ug, Ug])
#                     self.parent.zonationThxList.append([0.0, 1.0])
#                     self.parent.zonationThyList.append([Th, Th])

        
##########################################################################
# class ZirconGrainsParameters(QWidget):
#     """
#       This class contains the grains characteristics for He computation
#       it also show zonation if asked by the user.
    
#     """
#     def __init__(self, parent, Param, nSamples,Samples,system):
#         super().__init__()
        
#         # To store the grain characteristics
#         self.parent = parent
#         self.old_input = 0 #old input file signal
#         self.toolbarList = []
#         self.Param = Param
        
#         # Define Splitters
#         self.GrainSplitter = pgu.QCustomSplitter(Qt.Vertical)
#         self.TableSplit = QWidget()
#         self.ZonationSplit = QWidget()
#         self.ZonationSplit.setEnabled(False)
        
#         #### Parameters #####
#         # Table
#         if not nSamples: #if number of sample == 0, force to 1 at least
#             nSamples = 1
#         self.nbSampleValue = nSamples
#         self.Samples = Samples
#         self.Param = Param
#         self.Labelw = QLabel('Set grains characteristics:')
#         self.Labelw.setFont(conf.fontLine12)
#         self.Labelw.setAlignment(Qt.AlignLeft)
#         data = {" ID": ['1']*int(Param.DParameters['ngrains']),
#                 "Size (µm)": ["0"]*int(Param.DParameters['ngrains']),
#                 "U (ppm)": ["0"]*int(Param.DParameters['ngrains']),
#                 "Th (ppm)": ["0"]*int(Param.DParameters['ngrains'])}
#         # Zonation
#         self.ZonationCheckBox = QCheckBox("Zonation")
#         self.ZonationCheckBox.setEnabled(False)
#         if Param.DParameters['zonation_flag'] == '1':
#             self.old_input = 1
#             self.ZonationCheckBox.setChecked(True)
#             self.ZonationSplit.setEnabled(True)
#         self.ZonationLabel = QLabel("Set zonations:")
#         self.ZonationLabel.setFont(conf.fontLine12)
#         self.grainListLabel = QLabel('Select the grain:')
#         self.grainListLabel.setFont(conf.font8)
#         self.grainComboBox = QComboBox()
#         GrainList = (['grain'+str(j+1)+'_'+str(i+1)+'  ' for j in range(int(self.nbSampleValue.text()))
#                       for i in range(int(self.Samples['ZHe_ngrains'+str(j+1)]))])
#         self.grainComboBox.addItems(GrainList)
#         self.UserDescription = QLabel("Use the keyboard to set zonation for:\n'u': uranium, 't': thorium, 'r': restore both curves.")
#         self.SaveButton = QPushButton("save")
#         self.SaveButton.clicked.connect(lambda: self.saveZonation())
#         self.setNbLayersLabel = QLabel("Number of layers: ")
#         self.setNbLayersLabel.setFont(conf.font8)
#         self.setNbLayersLabel.setToolTip("Concentrations are averaged in each layer")
#         self.setNbLayers = QLineEdit()
#         self.setNbLayers.setValidator(QIntValidator(0, 50, self))
#         self.setNbLayers.setText(Param.DParameters['nLayers'])
        
#         # Check for user update
#         self.grainComboBox.currentIndexChanged.connect(
#             lambda: self.update_ComboBox(self.grainComboBox))
#         self.ZonationCheckBox.stateChanged.connect(
#             lambda: self.zonationState(self.ZonationCheckBox))
        
#         # How many grains in total?
#         self.TotGrain = 0
#         for j in range(int(self.nbSampleValue.text())):
#                 nb_grains = self.Samples['ZHe_ngrains'+str(j+1)]
#                 self.TotGrain += int(nb_grains)
        
#         self.GrainTable = U_interface.WinTable(data, self.TotGrain, 4)
#         self.updateTableGrain(self.GrainTable, 1, Param, self.parent)
#         self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
#             self.GrainTable, 0, Param, self.parent))
#         # self.parent.parent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
#         #     self.GrainTable, 2, Param, self.parent))
        
#         # Relative boundaries of the grain
#         self.x = [0.0, 1.0]
        
#         #Get and store zonation plot
#         # self.getZonationPlots()
        
#         # Add toolbar
#         # self.toolbar = NavigationToolbar(self.ZonationPlot,self)
        
#         #### Build the layouts ####
#         # Table split
#         VBox1 = QVBoxLayout()
#         VBox2 = QVBoxLayout()
#         self.Grid = QGridLayout()
#         self.Grid.setSpacing(10)
#         self.Grid.addWidget(self.ZonationCheckBox, 2, 0)
#         self.Grid.addWidget(self.GrainTable, 3, 0, 1, 4)
#         self.Grid.setRowStretch(3,5)
#         VBox1.addWidget(self.Labelw)
#         VBox1.addLayout(self.Grid)
 
#         # Zonation Split
#         VBox2.addWidget(self.ZonationLabel)
#         VBox2.addStretch(1)
#         Grid = QGridLayout()
#         Grid.addWidget(self.grainListLabel, 0, 1)
#         Grid.addWidget(self.grainComboBox, 0, 2, 1, 1)
#         Grid.addWidget(self.SaveButton, 0, 4)
#         Grid.addWidget(self.setNbLayersLabel, 1, 1)
#         Grid.addWidget(self.setNbLayers, 1, 2)
#         Grid.addWidget(self.UserDescription, 2, 1, 1, 4)
#         Grid.setColumnStretch(6,1)
#         VBox2.addLayout(Grid)
#         VBox2.addStretch(1)
#         # VBox2.addWidget(self.toolbar)
#         # VBox2.addWidget(self.ZonationPlot)
#         VBox2.addStretch(1)
        
#         # Add Layouts to widgets
#         self.TableSplit.setLayout(VBox1)
#         # self.ZonationSplit.setLayout(VBox2)
        
#         # Add widgets to splitter
#         self.GrainSplitter.addWidget(self.TableSplit)
#         # self.GrainSplitter.addWidget(self.ZonationSplit)
        
#         # Add splitter to layout
#         self.Layout = QVBoxLayout()
#         self.Layout.addWidget(self.GrainSplitter)
        
#     #----------------------------------------------------------
#     def zonationState(self,checkBox):
#         """ Enable/disable the zonation area """
#         if checkBox.isChecked():
#             self.ZonationSplit.setEnabled(True)
#             self.parent.zonation = 1
#         else:
#             self.ZonationSplit.setEnabled(False)
#             self.parent.zonation = 0
            
#     #----------------------------------------------------------
#     def getZonationPlots(self):
#         """ Gets and stores the zonation plot for each grain"""
#         if len(self.parent.zonationUxList) == 1 or not self.parent.zonationUxList: 
#             # means first time the user open the window, otherwise the zonation profile are
#             # loaded from previous opening
#             # Store the default values. Each grain has two lists, for U and Th values respectively.
#             count = 0
#             for j in range(int(self.nbSampleValue.text())):
#                 nb_grains = self.Samples['ZHe_ngrains'+str(j+1)]
#                 for i in range(int(nb_grains)):
#                     item = self.GrainTable.item(count, 0)
#                     if item is None:
#                         break
#                     sampleName = item.text()
#                     Ug = float(self.grains['ZUPPM_'+sampleName])
#                     Th = float(self.grains['ZTHPPM_'+sampleName])
#                     self.parent.zonationUxList.append(self.x)
#                     self.parent.zonationUyList.append([Ug, Ug])
#                     self.parent.zonationThxList.append(self.x)
#                     self.parent.zonationThyList.append([Th, Th])  
#                     count += 1
#         if len(self.parent.zonationUxList) > 0:
#             # Create zonation profiles for the first grain
#             self.ZonationPlot = pgu.InteractivePlot(Uxdata=self.parent.zonationUxList[0],
#                                                 Uydata=self.parent.zonationUyList[0],
#                                                 Thxdata= self.parent.zonationThxList[0],
#                                                 Thydata=self.parent.zonationThyList[0])
#         else:
#             self.ZonationPlot = pgu.InteractivePlot(Uxdata=[0],
#                                                 Uydata=[10],
#                                                 Thxdata= [0],
#                                                 Thydata=[10])
                
#     #--------------------------------------------------------- 
#     def update_ComboBox(self, cb):
#         # Get index of the current label shown (grain number)
#         index = cb.currentIndex()
        
#         # First clear the axes, ax1 is for uranium, ax2 is for thorium
#         self.ZonationPlot.ax1.clear()
#         self.ZonationPlot.ax2.clear()
#         # then clear the previous list of stored points
#         self.ZonationPlot.listPoints = []
#         # then update zonation plot
#         self.ZonationPlot.Uxdata = np.array(self.parent.zonationUxList[index])
#         self.ZonationPlot.Uydata = np.array(self.parent.zonationUyList[index])
#         self.ZonationPlot.Thxdata = np.array(self.parent.zonationThxList[index])
#         self.ZonationPlot.Thydata = np.array(self.parent.zonationThyList[index])
#         # Create the plots
#         self.ZonationPlot.create_draggable_points()
#         # Restore the grids
#         self.ZonationPlot.ax1.minorticks_on()
#         # self.ZonationPlot.ax1.tick_params(axis='y', labelcolor='red')
#         self.ZonationPlot.ax1.grid(b=True,which='major', color='0.8', linestyle='-')
#         self.ZonationPlot.ax1.grid(b=True,which='minor', color='0.5', linestyle='--')
#         # self.ZonationPlot.ax2.tick_params(axis='y', labelcolor='blue')
#         self.ZonationPlot.ax2.grid(b=True,which='major', color='0.8', linestyle='-')
#         self.ZonationPlot.ax2.grid(b=True,which='minor', color='0.5', linestyle='--')
#         self.ZonationPlot.ax1.set_xlim((0,1))
#         self.ZonationPlot.ax1.set_ylim((0,max(self.ZonationPlot.Uydata)+30))
#         self.ZonationPlot.ax2.set_ylim((0,max(self.ZonationPlot.Thydata)+30))
#         self.ZonationPlot.ax1.set_xlabel('edge <- Grain distance -> core')
#         self.ZonationPlot.ax1.set_ylabel('U (ppm)', color = 'red')
#         self.ZonationPlot.ax2.set_ylabel('Th (ppm)', color = 'blue')
#         self.ZonationPlot.ax2.set_xlabel('edge <- Grain distance -> core')
#         # update the plots in the interface
#         self.ZonationPlot.draw()
    
#     #--------------------------------------------------------
#     def saveZonation(self):
#         # save the new plots, possibly changed manually by the user
#         index = self.grainComboBox.currentIndex()
#         # Uranium
#         self.parent.zonationUxList[index] = self.ZonationPlot.listPoints[0].line_object[0].get_xdata()
#         self.parent.zonationUyList[index] = self.ZonationPlot.listPoints[0].line_object[0].get_ydata()
#         # Thorium
#         self.parent.zonationThxList[index] = self.ZonationPlot.listPoints[1].line_object[0].get_xdata()
#         self.parent.zonationThyList[index] = self.ZonationPlot.listPoints[1].line_object[0].get_ydata()
    
#     #--------------------------------------------------------
#     def updateTableGrain(self, t, h, p, d):
#         """ 
#           To store input value for the Grain table
#           t = table
#           h : parameter
#           p : parameters
#           d : pecube parameters (window)
#         """
#         self.grains = {}
#         h_temp = 1
#         if h == 2: # Update from Table Samples
#             h_temp = 2
#             try:
#                 self.GrainTable.cellChanged.disconnect() # Shut down update from the user
#             except TypeError:
#                 pass
#             if self.parent.parent.SampleTable.column(self.parent.parent.SampleTable.currentItem()) == 4:   
#                 h = 1
#         t.setRowCount(self.TotGrain)
#         if h:
#             # dParam = DefaultParameterValues()
#             count = 0
#             for j in range(int(self.nbSampleValue.text())):
#                 nb_grains = int(self.Samples['ZHe_ngrains'+str(j+1)])
#                 if nb_grains > 0:
#                     for i in range(int(nb_grains)):
#                         # Grain ID
                        
#                         if "SampleName"+str(j+1) in p.TableParameters:
#                             grainID = p.TableParameters['SampleName'+str(j+1)]+"_"+str(i+1)
#                             t.setItem(count, 0, QTableWidgetItem(grainID))
#                         else:
#                             # The sample ID is found in class 'WindowParameters' in PecubeGUI.py
#                             # so go back to parent (AgesParameter) then parent again (WindowParameters)
#                             grainID = self.parent.parent.Samples['SampleName'+str(j+1)]+"_"+str(i+1)
#                             t.setItem(count,0, QTableWidgetItem(grainID))
#                         if "ZSIZE_"+grainID in p.TableParameters:
#                             t.setItem(count, 1, QTableWidgetItem(
#                             str(p.TableParameters['ZSIZE_'+grainID])))
#                         else:
#                             t.setItem(count, 1, QTableWidgetItem(
#                                 p.DParameters['grainradius']))
#                         if "ZUPPM_"+grainID in p.TableParameters:
#                             t.setItem(count, 2, QTableWidgetItem(
#                                 str(p.TableParameters['ZUPPM_'+grainID])))
#                         else:
#                             t.setItem(count, 2, QTableWidgetItem(
#                                 p.DParameters['Uppm']))
#                         if "ZTHPPM_"+grainID in p.TableParameters:
#                             t.setItem(count, 3, QTableWidgetItem(
#                                 str(p.TableParameters['ZTHPPM_'+grainID])))
#                         else:
#                             t.setItem(count, 3, QTableWidgetItem(
#                                 p.DParameters['Thppm']))
                    
#                         count += 1
                    
#         num_rows, num_cols = t.rowCount(), t.columnCount()
#         temp_dict = {}
#         temp_GrainID = []
#         try:
#             for row in range(num_rows):
#                 for col in range(num_cols):
#                     item = t.item(row, col)
#                     if item is None:
#                         break
#                     if col == 0:
#                         sampleName = item.text()
#                         temp_GrainID.append(sampleName)
#                     if col == 1:
#                         value = float(item.text())
#                         temp_dict['ZSIZE_'+sampleName] = str(item.text())
#                         p.TableParameters['ZSIZE_'+sampleName] = str(item.text())
#                     elif col == 2:
#                         value = float(item.text())
#                         temp_dict['ZUPPM_'+sampleName] = str(item.text())
#                         p.TableParameters['ZUPPM_'+sampleName] = str(item.text())
#                     elif col == 3:
#                         value = float(item.text())
#                         temp_dict['ZTHPPM_'+sampleName] = str(item.text())
#                         p.TableParameters['ZTHPPM_'+sampleName] = str(item.text())
                        
#                 if item is None:
#                     break
#         except ValueError as VE:
#             print("Error", VE)
#             QErrorMessage(self).showMessage("In Grain table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
#             temp_dict['GrainIDs'] = temp_GrainID
#             self.grains.update(temp_dict)
#             return
#         temp_dict['GrainIDs'] = temp_GrainID
#         self.grains.update(temp_dict)
#         # Update AFT rmr0
#         if self.parent.AFTParam:
#             try:
#                 self.parent.AFTParam.GrainTable.cellChanged.disconnect()
#                 self.parent.AFTParam.updateTableGrain(self.parent.AFTParam.GrainTable, 1, p)
#                 self.parent.AFTParam.GrainTable.cellChanged.connect(lambda:self.parent.AFTParam.updateTableGrain(self.parent.AFTParam.GrainTable, 0, p))
#             except:
#                 print('Line 3092 in Thermochronometers.py: disconnect() failed')
                
#         # We write U, Th for all nodes
#         if self.parent.parent.AgeCombo.currentIndex() == 1:
#             for k, v in self.grains.items():
#                 if k.startswith('ZUPPM'):
#                     pgui.store_input(self.parent.parent, {'ZUPPM': str(v)})
#                 elif k.startswith('ZTHPPM'):
#                     pgui.store_input(self.parent.parent, {'ZThPPM': str(v)})
#                 elif k.startswith('ZSIZE'):
#                     pgui.store_input(self.parent.parent, {'ZSIZE': str(v)})
                    
#         if h == 2: #update from Table Samples
#             self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
#                 self.GrainTable, 0, self.Param))
        

#         # Store the default values for zonation plot.
#         # Each grain has two lists, for U and Th values respectively.
#         if not h: #means the user changes the cell of the table manually
#             #means first time the user open the window, otherwise the zonation profile are
#             #loaded from previous opening
#             self.parent.zonationUxList = []
#             self.parent.zonationUyList = []
#             self.parent.zonationThxList = []
#             self.parent.zonationThyList = []
#             for i in range(int(self.nbSampleValue.text())):
#                 nb_grains = int(self.Samples['ZHe_ngrains'+str(i+1)])
#                 for j in range(int(nb_grains)):
#                     Ug = float(self.grains['ZUPPM_'+temp_GrainID[j]])
#                     Th = float(self.grains['ZTHPPM_'+temp_GrainID[j]])
#                     self.parent.zonationUxList.append([0.0, 1.0])
#                     self.parent.zonationUyList.append([Ug, Ug])
#                     self.parent.zonationThxList.append([0.0, 1.0])
#                     self.parent.zonationThyList.append([Th, Th])



               
               
