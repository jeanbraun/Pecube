#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This module contains functions and classes used in the interface
for the input parameters for Pecube and CalcAge.

@author: maxime Bernard

"""


import sys
import os

import configs as conf
from PyQt5.QtWidgets import (QWidget, QMainWindow,QPushButton,
                             QVBoxLayout, QHBoxLayout, QGridLayout,
                             QTabWidget, QLabel,QLineEdit,QMessageBox, QCheckBox,
                             QSizePolicy,QFrame,QComboBox,QTableWidgetItem, 
                             QFileDialog,QGroupBox,QSlider, QScrollArea, QRadioButton,
                             QErrorMessage)
from PyQt5.QtGui import (QPixmap, QIcon,QIntValidator)
from PyQt5.QtCore import QProcess, Qt
import re
import shutil
from math import *
from scipy.io import loadmat
import numpy as np
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import pyvista as pv

import Utils.PGUI_utils as pgu
import Utils.interface as U_interface
import Utils.GIS as GIS

import Topography.update_topography as update_topo

try:
    from cmcrameri import cmc #import colormaps from Crameri, F. (2018), Scientific colour-maps, Zenodo,doi:10.5281/zenodo.1243862    
except: # cmap is already registered
    pass
import Thermochronology.thermochronometers as th
from bmi_topography import Topography
try:
    from PyQt5 import QtWebEngineWidgets
except Exception as E:
    print(type(E))
    print(E.args)

#Locate Icones directory
IconPath = conf.IconPath

# Path of executable
FolderPath = conf.FolderPath

#Force dot notation for the decimals
DoubleValidator = conf.DoubleValidator

#################################################################################
#################################################################################
#################################################################################
class ParamWin(QMainWindow, object):
    """ Window for the input parameters tabs
     PName is the name of the Project
     DParam is the default parameters stored
     IParam is the input parameters stored
    """
    def __init__(self, parent, PName, DParam, IParam, oldProject):
        super().__init__()

        self.success = 1.
        self.parent = parent
        self.UsePredictedElevation = self.parent.UsePredictedElevation
        self.aborted = 0 #Signal if run are stopped by the user
        self.layout = QGridLayout(self)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.VLine)
        self.PNAME = PName
        self.IParam = IParam
        self.DParam = DParam
        self.pre_hist_signal = int(self.DParam.DParameters['npreStep'])
        self.p = None  # Default empty value for 'Run Pecube' button
        self.RunWindow = 0
        self.oldProject = oldProject
        # Create a new directory with PName
        self.PecubeFolder = self.parent.PecubePath
        self.project_folder = os.path.join(self.PecubeFolder, self.PNAME)
        self.input_folder = os.path.join(self.project_folder, 'input')
        self.VTK_folder = os.path.join(self.project_folder, 'VTK')
        self.data_folder = os.path.join(self.project_folder, 'data')
        self.NA_folder = os.path.join(self.project_folder, 'NA')
        self.LOG_folder = os.path.join(self.project_folder, 'LOG')
        self.Output_folder = os.path.join(self.project_folder, 'output')
        self.inInversionMode = 0
        
        #Check if folder exists
        # if oldProject != 2:#means we create a new file or load from an old input file
        #     if os.path.exists(self.project_folder):
        #         win = pgu.MessageBox(text1='This Project already exists, do you want to replace it?')
        #         if win.QMessage == QMessageBox.Ok:
        #             if IParam.storeInput: #means old input
        #                 if os.path.exists(self.NA_folder):
        #                     shutil.rmtree(self.VTK_folder)
        #                     shutil.rmtree(self.NA_folder)
        #                     shutil.rmtree(self.LOG_folder)
        #                     os.mkdir(self.VTK_folder)
        #             else:
        #                 shutil.rmtree(self.project_folder,onerror=handleReadonly)
        #                 os.mkdir(self.self.project_folder)
        #                 os.mkdir(self.input_folder)
        #                 os.mkdir(self.VTK_folder)
        #                 os.mkdir(self.data_folder)
        #         elif win.QMessage == QMessageBox.Cancel or win.QMessage == QMessageBox.No:
        #             self.success = 0
        #             return 
        #     else:
        #         os.mkdir(self.project_folder)
        #         os.mkdir(self.input_folder)
        #         os.mkdir(self.VTK_folder)
        #         os.mkdir(self.data_folder)
        # First set Title/instructions for the Table
        self.textDes = QLabel('Navigate through tables to set parameters')
        self.textDes.setFont(conf.font10)
        self.textDes.adjustSize()

        self.setWindowTitle('Create new input file')
        self.setGeometry(100, 100, 500, 800)

        self.ParamTable = WindowParameters(
            self, parent, DParam, IParam, self.project_folder)

        # Second set 'Run Pecube Button'
        IconRunPecube = os.path.join(IconPath,"Run_Pecube.png")
        IconRunAges = os.path.join(IconPath,"Run_Ages.png")
        self.RunButton = QPushButton(QIcon(IconRunPecube), 'Run Pecube')
        self.RunButton.clicked.connect(lambda:self.Run_process())
        self.RunAges = QPushButton(QIcon(IconRunAges), 'run CalcAge')
        self.RunAges.clicked.connect(lambda:self.run_Ages())
        if oldProject == 1 and self.ParamTable.input_parameters[conf.Variable_names['Mode_age_computation']] == '2': #Means old project
            self.RunAges.setEnabled(True)
            # self.RunAges.setEnabled(False)
        else:
            self.RunAges.setEnabled(False)
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.RunButton)
        hbox.addWidget(self.RunAges)
        hbox.addStretch(1)

        # build layout
        # self.layout.addWidget(Separator, 0, 0, 4, 1)
        self.layout.addWidget(self.textDes, 1, 0)
        self.layout.addLayout(hbox, 2, 0, 1, 1)
        self.layout.addWidget(self.ParamTable.tabs, 3, 0, 1, 1)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.widget = QWidget()
        self.widget.setLayout(self.layout)
    
    #----------------------------------------------------------
    def Run_process(self):
        """ This function serves to make the 'Run Pecube' button active. """
        
        if self.p is None:  # No process running
            # 1) Check if VTK folder and Write input file
            if not os.path.exists(self.VTK_folder):
                os.mkdir(self.VTK_folder)
            filePath = self.input_folder
            fileName = os.path.join(filePath, 'Pecube.in')
            if os.path.exists(fileName):
                os.remove(fileName)
            file = open(fileName, 'w')
            inputPar = self.ParamTable.input_parameters
            # first clean stored parameter values
            defaultParam = conf.DefaultParameterValues()
            inputPar = pgu.clean_input_file(inputPar, defaultParam)
            for k, v in sorted(inputPar.items()):
                file.write(str(k) + ' = ' + str(v) + '\n')
                # check for parameter to invert for
                if v.find(":") > 0:
                    self.inInversionMode = 1
                    print("V is : ", v)
            file.close()
            
            # 2) Check preferences
            PreferencesPath = os.path.join(FolderPath,"preferences.txt")
            with open(PreferencesPath) as file:
                for line in file:
                    if "ShowConsole" in line:
                        (key, _, console) = line.split()  
            self.console = int(console)
            
            # 3) Save sample_specific file
            
            if int(self.ParamTable.nbSampleValue.text()) > 0 and self.ParamTable.AgesParam:
                try:
                    DataFolder = os.path.join(self.project_folder,'data',self.ParamTable.input_parameters['data_folder'])
                except KeyError:
                    QErrorMessage(self).showMessage("No folder name provide for the data. Please, provide a name in 'Data folder name")
                    return

                if os.path.exists(DataFolder):
                    files = os.listdir(DataFolder)
                    for f in files:
                        if f.startswith(conf.Variable_names['4He/3He observation file']):
                            pass
                        else:         
                            os.remove(os.path.join(DataFolder,f)) # Remove data
                self.ParamTable.AgesParam.save_obs_file()
            if self.ParamTable.input_parameters[conf.Variable_names['Mode_age_computation']] == '2' and self.ParamTable.AgesParam:
                self.ParamTable.AgesParam.savefile()
                
            # 4) If topo files have been provided from a spm, then write
            # the required files for Pecube (topo, temp, uplift)
            if self.ParamTable.spm and self.ParamTable.spm.spm_Grid.preHistory:
                # Means pre-history has been defined by the user
                # 4.1) Get the number of pre-history steps
                preSteps = int(self.ParamTable.spm.spm_Grid.nPreStep.text())
                
                # has the pre-history already been defined?
                if  self.pre_hist_signal: #yes, 
                    # Ask the user to reload files, this is the easiest way to proceed...
                    win = QMessageBox.question(self, 'PecubeGUI Message',
                                                'The pre-spm history has been changed. Please, reload the spm files',
                                                QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Cancel)
                    if win == QMessageBox.Ok:
                        # First delete all the files in the spm directory of the current project
                        path = os.path.join(self.data_folder,'SPM')
                        shutil.rmtree(path)
                        # Loader
                        options = QFileDialog.Options()
                        options |= QFileDialog.DontUseNativeDialog
                        fileName, _ = QFileDialog.getOpenFileNames(
                            self, "QFileDialog.getOpenFileName()", "", "All Files (*)", options=options)
                        if fileName:
                                # files from spm
                                spmFolder = os.path.join(self.project_folder,'data','SPM')
                                if os.path.exists(spmFolder):
                                    shutil.rmtree(spmFolder)
                                os.mkdir(spmFolder)
                                nfiles = len(fileName)
                                for i in range(nfiles):
                                    shutil.copy(fileName[i],spmFolder)
                    else:
                        return

                # 4.2) Rename files according to the pre-history steps (temporary))
                path = os.path.join(self.data_folder,'SPM/')
                files = os.listdir(path)
                files = sorted(files)
                topoFile = os.path.join(path,"topo000")
                try:
                    elevation = np.loadtxt(topoFile)
                except OSError:
                    QErrorMessage().showMessage("File not found: "+topoFile)
                    return
                InitialTopo = elevation.astype(float)
                upliftFile = os.path.join(path,"uplift000") 
                try:
                    uplift = np.loadtxt(upliftFile)
                except OSError:
                    QErrorMessage().showMessage("File not found: "+upliftFile)
                    return
                uplift_mask = uplift.astype(float)
                uplift_mask[:] = 1 # create mask to handle non-uniform uplift
                
                nsteps = int(len(files)/3) #There is 3 type of files (topo, temp, uplift)
                for i in range(nsteps): #First run to copy the files
                    if i+preSteps < 10:
                        zeroline = '00'
                    elif i+preSteps >= 10 and i+preSteps < 100:
                        zeroline = '0'
                    else:
                        zeroline = ''
                    os.rename(path+files[i],f"{path}ttopo"+zeroline+f"{i+preSteps}")
                    os.rename(path+files[i+nsteps], f"{path}ttemp"+zeroline+f"{i+preSteps}")
                    os.rename(path+files[i+2*nsteps], f"{path}tuplift"+zeroline+f"{i+preSteps}")
                files = os.listdir(path)
                files = sorted(files)
                for i in range(nsteps): #Second run to rename the new files 
                    if i+preSteps < 10:
                        zeroline = '00'
                    elif i+preSteps >= 10 and i+preSteps < 100:
                        zeroline = '0'
                    else:
                        zeroline = ''
                    os.rename(path+files[i], f"{path}topo"+zeroline+f"{i+preSteps}")
                    os.rename(path+files[i+nsteps], f"{path}temp"+zeroline+f"{i+preSteps}")
                    os.rename(path+files[i+2*nsteps], f"{path}uplift"+zeroline+f"{i+preSteps}")
                    
                # 4.3) Now write the pre-history files
                # 4.3a) Check the pre-history topographic evolution (with amplification and offset)
                TimeTable = self.ParamTable.timeTable_dict
                for i in range(preSteps):
                    if i < 10:
                        zeroline = '00'
                    elif i >= 10 and i < 100:
                        zeroline = '0'
                    else:
                        zeroline = ''
                    amplification = float(TimeTable['amplification'+str(i+1)])
                    offset = float(TimeTable['offset'+str(i+1)])
                    #Compute topography
                    topoRef = 0
                    if self.ParamTable.RefElevCombo.currentIndex() == 0: #From sea level
                        topoRef = 0
                    elif self.ParamTable.RefElevCombo.currentIndex() == 1: #From minimum
                        topoRef = np.min(InitialTopo)
                    elif self.ParamTable.RefElevCombo.currentIndex() == 2: #From maximum
                        topoRef = np.max(InitialTopo)
                    elif self.ParamTable.RefElevCombo.currentIndex() == 3: #From mean
                        topoRef = np.mean(InitialTopo)
                        
                    topo = topoRef - (amplification * (topoRef-InitialTopo)) + offset*1e3
                    np.savetxt(f"{path}topo"+zeroline+f"{i}",topo,fmt="%.2f")
                    # 4.3b) Compute surface temperature accordingly
                    Tsl = float(self.ParamTable.spm.spm_Grid.Tsl.text())
                    lapseRate = float(self.ParamTable.spm.spm_Grid.LapseRate.text())/1e3
                    temperature = Tsl - lapseRate * topo
                    np.savetxt(f"{path}temp"+zeroline+f"{i}",temperature,fmt="%.2f")                 
                    # 4.3c) Uplift defined by the user but constrained to the same area as
                    # defined in the spm (non-uniform uplift rates)
                    upliftRate = float(self.ParamTable.spm.spm_Grid.uplift.text())
                    upliftRate = upliftRate * uplift_mask
                    np.savetxt(f"{path}uplift"+zeroline+f"{i}",upliftRate,fmt="%.2f") 
                    
                # We need to write the provided uplift rate into the first file of the spm
                np.savetxt(f"{path}uplift"+zeroline+f"{preSteps}",upliftRate,fmt="%.2f")
                
                # We store the fact that we already went through the process
                self.pre_hist_signal = 1
            
            if not self.parent.runOpen: # if any model is already running
                # Then open a new window
                self.parent.runWindow.show()
                conf.WindowsOpen.append(self.parent.runWindow) 
                self.parent.runOpen = 1
            
            # Disable run button
            self.RunButton.setEnabled(False)
            
            # Then run Pecube if Pecube is allowed to run
            if not self.RunWindow: # Create a new console
                self.RunWindow = U_interface.showMessage(self, self.console, self.parent, self.PNAME)
            
            if self.console == 1:
                self.RunWindow.text.appendPlainText("Pecube is running.")
                
            # Run Pecube
            self.p = QProcess()
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.stateChanged.connect(self.handle_state)
            self.p.finished.connect(self.process_finished)
            # self.p.start('bash',['./bin/compile.sh'])
            os.chdir(self.PecubeFolder)
            if int(self.ParamTable.NbCores.text()) > 0: #"Means we want to run in parallel
                if sys.platform == 'darwin' or sys.platform =='linux':
                    path = os.path.join('bin','PecubeMPI.sh')
                    path = path.replace(os.sep, '/')
                    self.p.start('bash', [path, self.ParamTable.NbCores.text(), str(self.PNAME)])
                    
                elif sys.platform == 'win32' or sys.platform == 'cygwin':
                    # path = os.path.join(conf.FolderPath,'Pecube','bin','mpiexec')
                    self.p.setProgram('mpiexec')
                    self.p.setArguments(['-np',self.ParamTable.NbCores.text(),'bin\PecubeMPI.exe',str(self.PNAME)])
                    self.p.start()
                    
            else:
                if sys.platform == 'darwin' or sys.platform =='linux':
                    if self.ParamTable.InvModeCombo.currentIndex() == 1: # batch mode
                        path = os.path.join('./bin','run_batch.sh')
                        self.p.start('bash', [path, str(self.PNAME)])
                    else: # normal mode
                        path = os.path.join('./bin','run.sh')
                        self.p.start('bash', [path, str(self.PNAME)])
                    print('Path to Pecube executable: ', path)
                elif sys.platform == 'win32' or sys.platform == 'cygwin':
                    # path = os.path.join(conf.FolderPath,'Pecube','bin','Pecube.exe')
                    path = os.path.join(self.PecubeFolder,'bin','Pecube.exe')
                    print('Path to Pecube executable: ', path)
                    self.p.setProgram(path)
                    self.p.setArguments([str(self.PNAME)])
                    self.p.start()
            self.RunWindow.CancelButton.clicked.connect(lambda: self.cancelProcess())
    
    #----------------------------------------------------------
    def run_Ages(self):
        """ To run the production-diffusion model only. """
        
        if self.p is None:  # No process running
            os.chdir(self.PecubeFolder)
            # 1) Write input file
            filePath = self.input_folder
            fileName = os.path.join(filePath, 'Pecube.in')
            if os.path.exists(fileName):
                os.remove(fileName)
            file = open(fileName, 'w')
            inputPar = self.ParamTable.input_parameters
            for k, v in inputPar.items():
                file.write(str(k) + ' = ' + str(v) + '\n')

            file.close()
            
            # 1) Check preferences
            PreferencesPath = os.path.join(FolderPath,"preferences.txt")
            with open(PreferencesPath) as file:
                for line in file:
                    if "ShowConsole" in line:
                        (key, _, console) = line.split()       
            self.console = int(console)
            
            # 2) Save sample_specific file
            if int(self.ParamTable.nbSampleValue.text()) > 0 and self.ParamTable.AgesParam:
                self.ParamTable.AgesParam.save_obs_file()
            if self.ParamTable.input_parameters[conf.Variable_names['Mode_age_computation']] == '2' and self.ParamTable.AgesParam:
                self.ParamTable.AgesParam.savefile()
            
            if not self.parent.runOpen: # if any model is already running
                # Then open a new window
                self.parent.runWindow.show()
                conf.WindowsOpen.append(self.parent.runWindow)
                self.parent.runOpen = 1
                
            if not self.RunWindow: # Create a new console
                self.RunWindow = U_interface.showMessage(self, self.console, self.parent, self.PNAME)
                
            if self.console == 1 and self.RunWindow:
                self.RunWindow.text.appendPlainText("CalcAge is running.")
                
            # 3) run external routine
            self.p2 = QProcess()
            self.p2.readyReadStandardOutput.connect(self.handle_stdout2)
            self.p2.readyReadStandardError.connect(self.handle_stderr2)
            self.p2.stateChanged.connect(self.handle_state)
            self.p2.finished.connect(self.process_finished2)
            # Which model to run?
            print(os.getcwd())
            print(os.path.join(conf.FolderPath,'Pecube','src','CalcAge.exe'))
            if int(self.ParamTable.input_parameters['PDModel_signal']):
                if sys.platform == 'darwin' or sys.platform =='linux':
                    path = os.path.join('bin','run_CalcAge.sh')
                    path = path.replace(os.sep, '/')
                    self.p2.start('bash', [path,str(self.PNAME)])
                elif sys.platform == 'win32' or sys.platform == 'cygwin':
                    print('Run CalcAge')
                    path = os.path.join(conf.FolderPath,'Pecube','bin','CalcAge.exe')
                    self.p2.setProgram(path)
                    self.p2.setArguments([str(self.PNAME)])
                    self.p2.start()
            else:
                path = os.path.join('bin','run_CalcAge.sh')
            
            self.RunWindow.CancelButton.clicked.connect(lambda: self.cancelProcess())
            
    #----------------------------------------------------------
    def cancelProcess(self):
        """ To abort Pecube if 'cancel' is clicked """
        
        # os.chdir('..')
        self.RunButton.setEnabled(True)
        if self.RunWindow:
            self.RunWindow.close()
        # except AttributeError:
        #     print("Handle close runWindow error. Continue.")
    
    #----------------------------------------------------------
    def process_finished(self):
        """ Check for console """
        
        if self.console == 1 and self.RunWindow:
            self.RunWindow.text.appendPlainText("Run is finished.")
            
        self.p = None
        if (sys.platform == 'win32' or sys.platform == 'cygwin') and self.inInversionMode == 0:
            # Run VTK
            # normal mode, don't run VTK in batch mode
            path = os.path.join(conf.FolderPath,'Pecube','bin','VTK.exe')
            self.p = QProcess()
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.stateChanged.connect(self.handle_state)
            self.p.finished.connect(self.process_finished_VTK)
            self.p.setProgram(path)
            self.p.setArguments([str(self.PNAME)])
            self.p.start()
        else:
            # enable run button
            self.RunButton.setEnabled(True)
            

    #----------------------------------------------------------
    def process_finished_VTK(self):
        self.p = None
        inputPar = self.ParamTable.input_parameters
        # enable run button
        self.RunButton.setEnabled(True)
        
        # #When Pecube run is finished, run CalcAge for He ages computation
        # if inputPar['age_computation'] == '2' and self.RunWindow and int(self.ParamTable.NbCores.text()) == 0: #"Means we want to run in parallel:
        #     self.RunWindow.progress.setValue(0)
        #     self.p2 = QProcess()
        #     self.p2.readyReadStandardOutput.connect(self.handle_stdout2)
        #     self.p2.readyReadStandardError.connect(self.handle_stderr2)
        #     self.p2.stateChanged.connect(self.handle_state)
        #     # Which model to run?
        #     print(os.getcwd())
        #     if (sys.platform == 'win32' or sys.platform == 'cygwin') and self.inInversionMode == 0:
        #         path = os.path.join('src','CalcAge.exe')
        #         self.p2.setProgram(path)
        #         self.p2.setArguments([str(self.PNAME)])
        #         self.p2.start()
        #     elif sys.platform == 'darwin' or sys.platform =='linux':
        #         path = os.path.join('bin','run_CalcAge.sh')
        #         path = path.replace(os.sep, '/')
        #         self.p2.start('bash', [path,str(self.PNAME)])
        #         self.p2.finished.connect(lambda: self.process_finished2())
        # else:
        #Enable 'Ok' button to be clicked to close the window
        os.chdir(FolderPath)
        if self.RunWindow:
            self.RunWindow.OkButton.setEnabled(True)
            self.RunWindow.Label.setText('Pecube run is finished!')
        
    #----------------------------------------------------------
    def process_finished2(self):
        #Enable 'Ok' button to be clicked to close the window
        if self.RunWindow:
            self.RunWindow.Label.setText('Computation of ages is finished!')
            self.RunWindow.OkButton.setEnabled(True)
        # enable run button
        self.RunButton.setEnabled(True)
        if not self.aborted: # Means if Pecube has run as expected
            self.RunAges.setEnabled(True)
            # self.RunAges.setEnabled(False)
        os.chdir(FolderPath)
        
        self.p = None
        
    #----------------------------------------------------------       
    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode(conf.encoding_label)
        # self.RunWindow.text.appendPlainText(stdout)
        progress_re = re.compile("Doing step (\d+) of (\d+) - (\d+) percent remaining - (\d+) iterations")
        progress = U_interface.progressbar(stdout,progress_re,1)
        
        if progress:
            try:
                self.RunWindow.progress.setValue(progress)
            except:
                print("RunWindow: 'int' object has no attribute 'progress'")
            
        if self.console == 1 and self.RunWindow:
            self.RunWindow.text.appendPlainText(stdout)

    #----------------------------------------------------------
    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode(conf.encoding_label)
        if self.console == 1 and self.RunWindow:
            self.RunWindow.text.appendPlainText(stderr)

    #----------------------------------------------------------
    def handle_stdout2(self):
        data = self.p2.readAllStandardOutput()
        try:
            stdout = bytes(data).decode(conf.encoding_label)
        except UnicodeDecodeError:
            stdout = bytes(data).decode("latin-1")
            
        progress_re2 = re.compile("Doing grain (\d+) of (\d+)")
        progress = U_interface.progressbar(stdout,progress_re2,2)
        
        if progress:
            self.RunWindow.progress.setValue(progress)
        
        print(self.console, self.RunWindow)
        if self.console == 1 and self.RunWindow:
            self.RunWindow.text.appendPlainText(stdout)

    #----------------------------------------------------------
    def handle_stderr2(self):
        data = self.p2.readAllStandardError()
        stderr = bytes(data).decode(conf.encoding_label)
        
        if self.console == 1 and self.RunWindow:
            self.RunWindow.text.appendPlainText(stderr)

    #----------------------------------------------------------
    def handle_state(self, state):
        states = {
            QProcess.NotRunning: 'Not running',
            QProcess.Starting: 'Starting',
            QProcess.Running: 'Running',
        }
        state_name = states[state]
        
        if self.console == 1 and self.RunWindow:
            self.RunWindow.text.appendPlainText(f"State changed: {state_name}")



#################################################################################
#################################################################################
#################################################################################
def store_input(self, d):
    # This function stores the input parameters provided by the user
    # and update their value from the default value
    self.input_parameters.update(d)
    #print(str(self.input_parameters))


#################################################################################
#################################################################################
#################################################################################
class spm_Window(QWidget):
    """ 
      This class stores and defines the parameters that show when loading from a
      spm, parent = spm_Object
    
    """
    def __init__(self, parent, Param):
        super().__init__()
        
        self.preHistory = 0 # Signal if user choose to set a pre-history
        # second define parameters read
        self.ScrollArea = QScrollArea()
        self.Widget = QWidget()
        if parent.iSOSIA:
            # for model dimensions
            self.nx_spm = QLabel("nx : "+Param['nx_iSOSIA'])
            self.nx_spm.setFont(conf.font10)
            self.ny_spm = QLabel("ny : "+Param['ny_iSOSIA'])
            self.ny_spm.setFont(conf.font10)
            self.L_spm = QLabel("L : "+Param['L_iSOSIA']+' m')
            self.L_spm.setFont(conf.font10)
            self.H_spm = QLabel("H : "+Param['H_iSOSIA'] + ' m')
            self.H_spm.setFont(conf.font10)

        nfiles = str(int(len(parent.fileName)/3))
        self.Nb_output = QLabel("Nb output files : "+nfiles)
        self.Nb_output.setFont(conf.font10)
        # for time evolution
        self.timeInit = QLabel("Initial time (Ma):")
        self.timeInit.setFont(conf.font10)
        self.timeInit.setAlignment(Qt.AlignCenter)
        self.timeInit.setToolTip('When the glacial scenario begins?')
        self.timeInitEdit = QLineEdit()
        self.timeInitEdit.setValidator(DoubleValidator)
        self.timeInitEdit.setText(parent.Param.Param.DParameters["initTime_spm"])
        self.timeEnd = QLabel("End time (Ma):")
        self.timeEnd.setFont(conf.font10)
        self.timeEnd.setAlignment(Qt.AlignCenter)
        self.timeEnd.setToolTip('When the glacial scenario ends?')
        self.timeEndEdit = QLineEdit()
        self.timeEndEdit.setValidator(DoubleValidator)
        self.timeEndEdit.setText(parent.Param.Param.DParameters["endTime_spm"])
        self.NbSteps = QLabel("Number of steps:")
        self.NbSteps.setFont(conf.font10)
        self.NbSteps.setAlignment(Qt.AlignCenter)
        self.NbSteps.setToolTip('how many time steps? (e.g., number output files)')
        self.NbStepsEdit = QLineEdit()
        self.NbStepsEdit.setValidator(QIntValidator())
        self.NbStepsEdit.setText(nfiles)

        # for pre-history
        # it is assumed that sea-level temperature is the same as the first temp file
        # and that uplift rate are uniform and constant over pre-history times
        self.nPreStepLabel = QLabel("Number of pre-steps:")
        self.nPreStepLabel.setFont(conf.font10)
        self.nPreStepLabel.setAlignment(Qt.AlignCenter)
        self.nPreStep = QLineEdit()
        self.nPreStep.setValidator(QIntValidator(0,1000))
        self.nPreStep.setText(parent.Param.Param.DParameters["npreStep"])
        self.upliftLabel = QLabel("Uplift rate (mm/yr):")
        self.upliftLabel.setFont(conf.font10)
        self.upliftLabel.setAlignment(Qt.AlignCenter)
        self.uplift = QLineEdit()
        self.uplift.setValidator(DoubleValidator)
        self.uplift.setText(parent.Param.Param.DParameters["uRateSpm"])
        self.TslLabel = QLabel("sea-level temperature (°C):")
        self.TslLabel.setFont(conf.font10)
        self.TslLabel.setAlignment(Qt.AlignCenter)
        self.Tsl = QLineEdit()
        self.Tsl.setValidator(DoubleValidator)
        self.Tsl.setText(parent.Param.Param.DParameters["Tsl_spm"])
        self.LapseRateLabel = QLabel("Lapse rate (°C/km):")
        self.LapseRateLabel.setFont(conf.font10)
        self.LapseRateLabel.setAlignment(Qt.AlignCenter)
        self.LapseRate = QLineEdit()
        self.LapseRate.setValidator(DoubleValidator)
        self.LapseRate.setText(parent.Param.Param.DParameters["lrate_spm"])
        
        # if text edited, store for input file
        # Check for user edition
        self.timeInitEdit.textEdited.connect(lambda:store_input(
                                parent.Param, {'initTime_spm': self.timeInitEdit.text()}))
        self.timeEndEdit.textEdited.connect(lambda: store_input(
                                            parent.Param,{'endTime_spm': self.timeEndEdit.text()}))
        self.NbStepsEdit.textEdited.connect(lambda: store_input(
                                            parent.Param,{'nstep_spm': self.NbStepsEdit.text()}))
        self.nPreStep.textEdited.connect(lambda: store_input(
                                            parent.Param,{'npreStep': self.nPreStep.text()}))
        self.uplift.textEdited.connect(lambda: store_input(
                                            parent.Param,{'uRateSpm': self.uplift.text()}))
        self.Tsl.textEdited.connect(lambda: store_input(
                                            parent.Param,{'Tsl_spm': self.Tsl.text()}))
        self.LapseRate.textEdited.connect(lambda: store_input(
                                            parent.Param,{'lrate_spm': self.LapseRate.text()}))
        
        # Extra box
        self.VBox = QVBoxLayout()
        #Add boxes for parameters
        Box1 = QGroupBox("Model information")
        Grid1 = QGridLayout()
        if parent.iSOSIA:
            Grid1.addWidget(self.nx_spm, 2, 1)
            Grid1.addWidget(self.ny_spm, 2, 2)
            Grid1.addWidget(self.L_spm, 3, 1)
            Grid1.addWidget(self.H_spm, 3, 2)
        Grid1.addWidget(self.Nb_output, 4, 1)
        Grid1.setSpacing(5)
        Grid1.setColumnStretch(0, 1)
        Grid1.setColumnStretch(1, 1)
        Grid1.setColumnStretch(2, 1)
        Grid1.setColumnStretch(3, 1)
        Box1.setLayout(Grid1)
        Box1.setFlat(True)
        self.VBox.addWidget(Box1)
        Box2 = QGroupBox("Time evolution")
        Grid2 = QGridLayout()
        # Grid2.addWidget(self.Label, 0, 0, 1, 3)
        Grid2.addWidget(self.timeInit, 2, 1)
        Grid2.addWidget(self.timeInitEdit, 3, 1)
        Grid2.addWidget(self.timeEnd, 2, 2)
        Grid2.addWidget(self.timeEndEdit, 3, 2)
        Grid2.addWidget(self.NbSteps, 4, 1)
        Grid2.addWidget(self.NbStepsEdit, 5, 1)
        Grid2.setSpacing(5)
        Grid2.setColumnStretch(0, 1)
        Grid2.setColumnStretch(3, 1)
        Box2.setLayout(Grid2)
        Box3 = QGroupBox("Pre-history")
        Box3.setCheckable(True)
        Box3.setChecked(False)
        Box3.toggled.connect(lambda:self.updateStatus(Box3))
        Grid3 = QGridLayout()
        Grid3.addWidget(self.nPreStepLabel,0,0)
        Grid3.addWidget(self.nPreStep,1,0)
        Grid3.addWidget(self.upliftLabel,0,1)
        Grid3.addWidget(self.uplift,1,1)
        Grid3.addWidget(self.TslLabel,2,0)
        Grid3.addWidget(self.Tsl,3,0)
        Grid3.addWidget(self.LapseRateLabel,2,1)
        Grid3.addWidget(self.LapseRate,3,1)
        Box3.setLayout(Grid3)
        self.VBox.addWidget(Box2)
        self.VBox.addWidget(Box3)
        self.Widget.setLayout(self.VBox)
        self.ScrollArea.setWidget(self.Widget)
    
    #-------------------------------------------------------
    def store_iSOSIAPecube_value(self, d):
        self.input_iSOPec.update(d)

    #-------------------------------------------------------
    def updateStatus(self,box):
        #to update the signal "pre-history" when the check box is checked
        if box.isChecked():
            self.preHistory = 1

        else:
            self.preHistory = 0
            
            
#################################################################################
#################################################################################
#################################################################################
class spm_objects(QWidget):
    """ 
      This class takes care of spm output files to feed Pecube input files
      The files are copied/pasted in the data folder of the current Pecube project
      parent = Windows_Parameters
    
    """
    def __init__(self, parent,filename):
        super().__init__()
        # Get the path to the iSOSIA folder in the Pecube project
        self.Param = parent
        self.iSOSIA = parent.iSOSIA
        self.iSOSIA_Param = {}
        
        # self.open_file()
        self.PecubeFolder = ''
        self.fileName = filename
        if self.fileName:
            # start by create a folder for the spm in data
            if len(self.fileName) > 1:
                self.spm_folder = os.path.join(self.Param.PFolder,"data","SPM")
                self.PecubeFolder = self.spm_folder
                
            # get some information from the iSOSIA model
            if self.iSOSIA:
                self.get_info()
            self.spm_Grid = spm_Window(self,self.iSOSIA_Param)
    
            # Open a window to summarize the model dimension
            self.w = U_interface.WinExtraParameters(text='spm model loader')
            self.w.setFixedWidth(700)
            self.w.setFixedHeight(500)
    
            Grid1 = QGridLayout()
            Grid1.addWidget(self.spm_Grid.ScrollArea, 0, 0, 3, 1)
            if self.iSOSIA:
                Background = QLabel()
                IconIsosia = os.path.join(IconPath,"iSOSIA2.png")
                oImage = QPixmap(IconIsosia)
                OI_resized = oImage.scaled(200, 300)
                Background.setPixmap(OI_resized)
                Background.setScaledContents(True)
                Grid1.addWidget(Background, 0, 1, 3, 2) # Background image
            self.w.layout.addLayout(Grid1, 0, 0)
            Q = QWidget()
            Q.setLayout(self.w.layout)
            self.w.setCentralWidget(Q)
            #add this window in the list of open windows
            conf.WindowsOpen.append(self.w)
            self.w.show()
    
            # Now write files for Pecube
            self.w.OkButton.clicked.connect(lambda: self.write_Pecube_file())

    #----------------------------------------------------------
    def get_info(self):
        """ To create iSOSIA files that can be read by Pecube. """
        self.extracted_Fnames = self.fileName[0].split("/")
        str1 = ''
        strlist = self.extracted_Fnames[-2:]
        str2 = str1.join(strlist)
        Nb_char = len(str2)+2  # count number of character after folder name
        self.folderName = self.fileName[0][:-Nb_char]
        self.inputFolder = os.path.join(self.folderName,'input','SPM.mat')
        # Read the SPM.mat file - data contains 'mesh', 'iprop', 'mprop',
        # 'hwprop','hprop', 'fprop', etc...
        self.iSOSIA_input = loadmat(self.inputFolder)
        # Read mesh to get the model dimensions
        SPM = self.iSOSIA_input['SPM']['mesh']
        self.nx = SPM[0, 0]['ny'][0][0][0][0]
        self.ny = SPM[0, 0]['nx'][0][0][0][0]
        self.iSOSIA_Param['Folder'] = str(self.folderName)
        self.iSOSIA_Param['nx_iSOSIA'] = str(self.nx)
        self.iSOSIA_Param['ny_iSOSIA'] = str(self.ny)
        self.iSOSIA_Param['L_iSOSIA'] = str(SPM[0, 0]['L'][0][0][0][0])
        self.iSOSIA_Param['H_iSOSIA'] = str(SPM[0, 0]['H'][0][0][0][0])
        
    #----------------------------------------------------------
    def write_Pecube_file(self):
        """ 
          This function allows to write spm output to be read by Pecube
          First check if the required text boxes are filled
        
        """
        if self.spm_Grid.timeInitEdit.text() and self.spm_Grid.timeEndEdit.text() and self.spm_Grid.NbStepsEdit.text():
            # Copy paste spm output file to the targeting Pecube folder
            
            # Original path
            OriginalPath = self.fileName
            
            # How many files ?
            nfiles = len(OriginalPath)
            
            if nfiles == 1:
                TargetPath = os.path.join(self.Param.PFolder,'data')
                # Only one file means topo file
                self.Param.topo_file_nameEdit1.setText(self.extracted_Fnames[-1])
                shutil.copy(OriginalPath[0],TargetPath)
                
            else:
                TargetPath = self.PecubeFolder
                # Provide '/iSOSIA' name to topo_file_name
                self.Param.topo_file_nameEdit1.setText('SPM/')
                for i in range(nfiles):
                    shutil.copy(OriginalPath[i],TargetPath)
                
            ###### Fill Table time_topo ######
            # Get number of step
            nb_steps = int(self.spm_Grid.NbStepsEdit.text()) #for spm time evolution
            # Set time steps
            time_start = float(self.spm_Grid.timeInitEdit.text())
            time_end = float(self.spm_Grid.timeEndEdit.text())
            time_topo_array = np.linspace(time_start,time_end,nb_steps)
            time_topo_array = time_topo_array.astype(float)
            
            # include pre-history
            if self.spm_Grid.preHistory:
                npresteps =  int(self.spm_Grid.nPreStep.text())
                nb_steps = int(self.spm_Grid.NbStepsEdit.text()) + npresteps -1
                preHistoryTime = np.zeros((npresteps))
                time_topo_array = np.concatenate((preHistoryTime, time_topo_array))
            
            # Store in TableParameters to be updated in the time table
            time_dict = {'time_topo'+str(i+1): str(round(time_topo_array[i],2)) for i in range(len(time_topo_array))}
            ampl_dict = {'amplification'+str(i+1): '1' for i in range(len(time_topo_array))}
            offset_dict = {'offset'+str(i+1): '0' for i in range(len(time_topo_array))}
            output_dict = {'output'+str(i+1): '0' for i in range(len(time_topo_array))}
            self.Param.Param.TableParameters.update(time_dict)
            self.Param.Param.TableParameters.update(ampl_dict)
            self.Param.Param.TableParameters.update(offset_dict)
            self.Param.Param.TableParameters.update(output_dict)
            # Update 'ntime' so that the table get new values
            self.Param.ntimeEdit1.setText(str(nb_steps))
            self.Param.updateTable(
                self.Param.topoTable, int(self.Param.ntimeEdit1.text()), 1, self.Param.Param)
            if self.iSOSIA:
                self.Param.nxEdit2.setText(str(self.iSOSIA_Param['ny_iSOSIA']))
                self.Param.nyEdit3.setText(str(self.iSOSIA_Param['nx_iSOSIA']))
            self.w.close()
        else:
            QErrorMessage().showMessage("Please, fill all the text boxes.")
            
        
#################################################################################
#################################################################################
#################################################################################
class Geotherm(QWidget):
    """ 
      This class computes the steady-state one-dimensional geotherm
      according to the Thermal parameters provided by the user.
      It uses the equation (8) from Reiners and Brandon (2006)
    
    """
    def __init__(self,Param, MainWin):
        super().__init__()
        
        self.Param = Param
        # Retrieve parameters
        try:
            self.Ts = float(self.Param.sltEdit4.text()) # T° at sea level
            self.TL = float(self.Param.bTemperatureEdit2.text()) # T° at the base
            self.Thickness = float(self.Param.thicknessEdit1.text()) # Thickness of the crust (km)
            self.diff = float(self.Param.TDiffusivityEdit6.text()) # Diffusivity (km²/Myr)
            self.HP =  float(self.Param.HProdEdit7.text()) # Heat production (°C/Myr)
            self.TopoElev = 0 # Elevation of mean topography (km) 
        except ValueError:
            print("Value cannot be read in Geotherm. Perhaps inversion is used.")
            return
        
        # Check for update of these parameters
        Param.sltEdit4.editingFinished.connect(lambda: self.updateGeotherm())
        Param.bTemperatureEdit2.editingFinished.connect(lambda: self.updateGeotherm())
        Param.thicknessEdit1.editingFinished.connect(lambda: self.updateGeotherm())
        Param.TDiffusivityEdit6.editingFinished.connect(lambda: self.updateGeotherm())
        Param.HProdEdit7.textChanged.connect(lambda: self.updateGeotherm())
        self.MainWin = MainWin
        
        # Extra parameters to investigate the geotherm
        self.Label = QLabel('Exploration parameters:')
        self.Label.setFont(conf.fontLine12)
        self.ErLabel = QLabel('Mean erosion rate (km/Myr):')
        self.ErLabel.setAlignment(Qt.AlignCenter)
        self.ErLabel.setFont(conf.font10)
        self.ErEdit = QLineEdit() #Mean erosion rate (km/Myr)
        self.ErEdit.setText('0')
        self.ErEdit.setAlignment(Qt.AlignCenter)
        self.ErEdit.setValidator(DoubleValidator)
        self.SurfHeatFluxLabel = QLabel("Surface heat flux (mW/m2):")
        self.SurfHeatFluxLabel.setFont(conf.font10)
        self.SurfHeatFluxLabel.setAlignment(Qt.AlignCenter)
        self.SurfHeatFlux = QLabel(Param.Param.DParameters['surfHeatFlux'])
        self.SurfHeatFlux.setAlignment(Qt.AlignCenter)
        
        # Heat production parameters
        self.TopLabel = QLabel("Heat production:")
        self.TopLabel.setFont(conf.fontLine12)
        self.SpecificHeatLabel = QLabel("Specific heat capacity (J/kg.K):")
        self.SpecificHeatLabel.setFont(conf.font10)
        self.SpecificHeatLabel.setToolTip("to provide in J/kg.K")
        self.SpecificHeatLabel.setAlignment(Qt.AlignCenter)
        self.SHEdit = QLineEdit(self)
        self.SHEdit.setText(self.Param.Param.DParameters['Specific_heat'])
        self.SHEdit.setValidator(DoubleValidator)
        self.SHEdit.setAlignment(Qt.AlignCenter)
        self.HeatProductionLabel = QLabel("Radioactive heat production (µW/m3):")
        self.HeatProductionLabel.setFont(conf.font10)
        self.HeatProductionLabel.setToolTip("to provide in µW/m3")
        self.HeatProductionLabel.setAlignment(Qt.AlignCenter)
        self.HPEdit = QLineEdit(self)
        self.HPEdit.setText(self.Param.Param.DParameters['Radioactive_heat']) 
        self.HPEdit.setValidator(DoubleValidator)
        self.HPEdit.setAlignment(Qt.AlignCenter)
        self.rhoLabel = QLabel('\u03C1c (kg/m3):')
        self.rhoLabel.setFont(conf.font10)
        self.rhoLabel.setToolTip("Average density of the crust.")
        self.rhoLabel.setAlignment(Qt.AlignCenter)
        self.rhoEdit = QLineEdit(self)
        self.rhoEdit.setText(Param.Param.DParameters['rho_crust'])
        self.rhoEdit.setValidator(DoubleValidator)
        self.rhoEdit.setAlignment(Qt.AlignCenter)
        
        # Define two splitters:
            #One for the parameters
            #One for the plotting area
        self.GeoSplitter = U_interface.QCustomSplitter(Qt.Vertical)
        self.ExploSplitter = QWidget()
        self.PlotSplitter = QWidget()
        
        # Check for updates
        self.ErEdit.editingFinished.connect(lambda: self.updateGeotherm())
        self.SHEdit.editingFinished.connect(lambda: self.updateHeatProd())
        self.HPEdit.editingFinished.connect(lambda: self.updateHeatProd())
        self.rhoEdit.editingFinished.connect(lambda: self.updateHeatProd())
        self.SHEdit.editingFinished.connect(lambda: self.updateGeotherm())
        self.HPEdit.editingFinished.connect(lambda: self.updateGeotherm())
        self.rhoEdit.editingFinished.connect(lambda: self.updateGeotherm())
        
        # Compute thermal profile
        self.Height = self.Thickness +  self.TopoElev #Total thickness of the profile
        nz = 50
        dz = self.Height/(nz-1)
        self.depth = np.arange(0,self.Height,dz) #Depth profile
                
        # Compute and show geotherm
        try:
            Er = float(self.ErEdit.text())
            Pe = Er*self.Height/self.diff # Peclet number
            if Pe > 0.1:
                self.Tz = (self.Ts + (self.TL-self.Ts + (self.HP*self.Height/Er))
                            * (1-np.exp(-Er*self.depth/self.diff))/(1-np.exp(-Er*self.Height/self.diff))
                            - self.HP*self.depth/Er)
            else:
                self.Tz = (self.Ts + (self.TL-self.Ts)*self.depth/self.Height
                            + self.HP*self.depth*(self.Height-self.depth)/(2*self.diff))
        except ValueError:
            print("Value cannot be read in Geotherm. Perhaps inversion is used.")
            return
                
        x = self.Tz
        y = self.depth
        self.plotSpace =  pgu.PlotCanvas()
        self.plotSpace.axes.plot(x,y,'r')
        self.plotSpace.axes.set_xlabel('Temperature (°C)')
        self.plotSpace.axes.set_ylabel('Depth (km)')
        self.plotSpace.axes.grid(True)
        self.plotSpace.axes.invert_yaxis()
        # Add toolbar
        self.toolbar = NavigationToolbar(self.plotSpace,self)
        
        # Grid for the main parameters
        self.VBox1 = QVBoxLayout()
        self.VBox1.addWidget(self.Label)
        self.VBox1.addStretch(1)
        # Grid for the extra parameters
        Gbox = QGridLayout()
        Gbox.addWidget(self.ErLabel, 0, 1)
        Gbox.addWidget(self.ErEdit, 1, 1)
        Gbox.addWidget(self.SurfHeatFluxLabel, 0, 2)
        Gbox.addWidget(self.SurfHeatFlux, 1, 2)
        Gbox.setColumnStretch(0,2)
        Gbox.setColumnStretch(1,1)
        Gbox.setColumnStretch(2,1)
        Gbox.setColumnStretch(3,2)
        self.VBox1.addLayout(Gbox)
        # Set the grid for the heat production parameters
        self.VBox1.addWidget(self.TopLabel)
        Grid = QGridLayout()
        Grid.setSpacing(10)
        Grid.addWidget(self.SpecificHeatLabel,2,1)
        Grid.addWidget(self.HeatProductionLabel,2,2)
        Grid.addWidget(self.rhoLabel,2,3)
        Grid.addWidget(self.SHEdit,3,1)
        Grid.addWidget(self.HPEdit,3,2)
        Grid.addWidget(self.rhoEdit,3,3)
        self.VBox1.addStretch(1)
        self.VBox1.addLayout(Grid)
        self.VBox1.addStretch(1)
        self.VBox2 = QVBoxLayout()
        self.VBox2.addWidget(self.toolbar)
        self.VBox2.addWidget(self.plotSpace)
        # Add the Layouts to the Splitters
        self.ExploSplitter.setLayout(self.VBox1)
        self.PlotSplitter.setLayout(self.VBox2)
        self.GeoSplitter.addWidget(self.ExploSplitter)
        self.GeoSplitter.addWidget(self.PlotSplitter)
        # Add Splitters to the main Tab
        VBox3 = QVBoxLayout()
        VBox3.addWidget(self.GeoSplitter)
        Q = QWidget()
        Q.setLayout(VBox3)
        MainWin.StackedTabs[self.MainWin.CurrentTabs].addTab(Q,"Geotherm")
        
    #-------------------------------------------    
    def check_inversion(self):
        """Check for inversion parameters, and give up updates if inversion."""
        try:
            self.Ts = float(self.Param.sltEdit4.text()) # T° at sea level
            self.TL = float(self.Param.bTemperatureEdit2.text()) # T° at the base
            self.Thickness = float(self.Param.thicknessEdit1.text()) # Thickness of the crust (km)
            self.diff = float(self.Param.TDiffusivityEdit6.text()) # Diffusivity (km²/Myr)
            self.HP =  float(self.Param.HProdEdit7.text()) # Heat production (°C/Myr)
            Er = float(self.ErEdit.text())
            Radiogenic_heat = float(self.HPEdit.text())*1e-6 # in watt/m3
            Heat_capacity = float(self.SHEdit.text())
            Rho_crust = float(self.rhoEdit.text())
            return 0
        except ValueError:
            print("Geotherm cannot be updated. Parameters for inversion have been provided.")
            return 1

    #----------------------------------------------------------
    def updateGeotherm(self):
        """
          To show the resulting 1D geotherm to the user
          Temperature profile array - Reiners and Brandon (2006) eq.(8)
        
        """
        # First, check if inversion parameters (range of values provided)
        Inversion = self.check_inversion()
        if Inversion:
            return
        # Else, retrieve parameters
        self.Ts = float(self.Param.sltEdit4.text()) # T° at sea level
        self.TL = float(self.Param.bTemperatureEdit2.text()) # T° at the base
        self.Thickness = float(self.Param.thicknessEdit1.text()) # Thickness of the crust (km)
        self.diff = float(self.Param.TDiffusivityEdit6.text()) # Diffusivity (km²/Myr)
        self.HP =  float(self.Param.HProdEdit7.text()) # Heat production (°C/Myr)
        self.TopoElev = 0 # Elevation of mean topography (km) 
        
        # Compute thermal profile
        self.Height = self.Thickness +  self.TopoElev # Total thickness of the profile
        nz = 50
        dz = self.Height/(nz-1)
        self.depth = np.arange(0,self.Height,dz) # Depth profile
        
        # Mean erosion rate
        Er = float(self.ErEdit.text())
        
        Pe = Er*self.Height/self.diff # Peclet number
        # E-folding depth ?
        if self.Param.Heat_production_flag == 1:
            self.HP = self.HP*np.exp(-self.depth/float(self.Param.eFolding.text()))
            if Pe > 0.1:
                self.Tz = (self.Ts + (self.TL-self.Ts + (self.HP[-1]*self.Height/Er))
                            * (1-np.exp(-Er*self.depth/self.diff))/(1-np.exp(-Er*self.Height/self.diff))
                            - self.HP*self.depth/Er)
            else:
                self.Tz = (self.Ts + (self.TL-self.Ts)*self.depth/self.Height
                            + self.HP*self.depth*(self.Height-self.depth)/(2*self.diff))
        else:
            if Pe > 0.1:
                self.Tz = (self.Ts + (self.TL-self.Ts + (self.HP*self.Height/Er))
                            * (1-np.exp(-Er*self.depth/self.diff))/(1-np.exp(-Er*self.Height/self.diff))
                            - self.HP*self.depth/Er)
            else:
                self.Tz = (self.Ts + (self.TL-self.Ts)*self.depth/self.Height
                            + self.HP*self.depth*(self.Height-self.depth)/(2*self.diff))
        
        x = self.Tz
        y = self.depth
        self.plotSpace.axes.clear()
        self.plotSpace.axes.plot(x,y,'r')
        self.plotSpace.axes.set_xlabel('Temperature (°C)')
        self.plotSpace.axes.set_ylabel('Depth (km)')
        self.plotSpace.axes.grid(True)
        self.plotSpace.axes.invert_yaxis()
        self.plotSpace.draw()
        
        # Compute resulting surface heat flux
        # get heat conductivity
        k = (float(self.Param.TDiffusivityEdit6.text())*1e6/(3600*24*365*1e6) * float(self.SHEdit.text()) *
                float(self.rhoEdit.text()) * ((self.Tz[1] - self.Tz[0]) / (self.depth[1] - self.depth[0])))
        self.SurfHeatFlux.setText(str(round(k,2)))
        
    #----------------------------------------------------------
    def updateHeatProd(self):
        """
          Update the value of the heat production (°C/Myr)
          according to the extra parameters provided by the user
          in the 'geotherm' Tab
        """
        # First check for inversion
        Inversion = self.check_inversion()
        if Inversion:
            return
        
        # Compute the heat production from detailed parameters
        Radiogenic_heat = float(self.HPEdit.text())*1e-6 #in watt/m3
        Heat_capacity = float(self.SHEdit.text())
        Rho_crust = float(self.rhoEdit.text())
        
        # Heat production
        Myr = 1e6*365.25*24*3600
        Heat_production = round(Radiogenic_heat/(Rho_crust*Heat_capacity)*Myr,2)
        
        # Update the Heat production value within the interface
        self.Param.HProdEdit7.setText(str(Heat_production))
        # Also update the density of the crust for consistency
        self.Param.rhoCrustEdit3.setText(str(Rho_crust))
        
        # Store parameters
        store_input(self.Param, {'Specific_heat': str(self.SHEdit.text())})
        store_input(self.Param, {'Radioactive_heat': str(self.HPEdit.text())})
        store_input(self.Param, {'rho_crust': str(self.rhoEdit.text())})
        store_input(self.Param, {'heat_production': str(self.Param.HProdEdit7.text())})
        

#################################################################################
#################################################################################
################################################################################# 
class build_Topography(QWidget):
    """ 
        To build a sinusoidal topography interactively.
        
    """
    def __init__(self,param):
        super().__init__()
        
        self.canvas = pgu.Canvas3D()
        self.sargs = dict(height=0.25, vertical=True, position_x=0.05,position_y=0.05,interactive=True,color=conf.Colors3Dplot['colorbarLabels'])
        self.specular = 0.3
        self.Param = param
        
        
        #####  Define parameters  ######
        # Labels
        self.ModelDimensions = QLabel("Model dimensions:")
        self.ModelDimensions.setFont(conf.fontLine12)
        self.ModelCharacteristics = QLabel("Model characteristics:")
        self.ModelCharacteristics.setFont(conf.fontLine12)
        self.Length = QLabel("Length (km):")
        self.Length.setFont(conf.font10)
        self.Width = QLabel("Width (km):")
        self.Width.setFont(conf.font10)
        self.resX = QLabel("nx:")
        self.resX.setFont(conf.font10)
        self.resY = QLabel("ny:")
        self.resY.setFont(conf.font10)
        self.AmpLabel = QLabel('Amplitude:')
        self.AmpLabel.setFont(conf.font10)
        self.WaveLabel = QLabel('Wavelength (km):')
        self.WaveLabel.setFont(conf.font10)
        self.MTopoLabel = QLabel("Reference elevation (km):")
        self.MTopoLabel.setFont(conf.font10)
        self.offsetLabel = QLabel("Offset (km):")
        self.offsetLabel.setFont(conf.font10)
        self.PhaseShiftLabel = QLabel('Phase shift (km):')
        self.PhaseShiftLabel.setFont(conf.font10)
        
        # Default values
        L = 50
        W = 50
        nx = 200
        ny = 200
        Amplitude = 1.0
        wavelength = 10
        Offset = 0.0
        Mean_topo = 0.5
        phaseShift = 0.0
        
        # Dimensions of the grid (m)
        self.L = QLineEdit()
        self.L.setText(str(L))
        self.L.setValidator(DoubleValidator)
        self.W = QLineEdit()
        self.W.setText(str(W))
        self.W.setValidator(DoubleValidator)
        self.nxEdit = QLineEdit()
        self.nxEdit.setText(str(nx))
        self.nxEdit.setValidator(QIntValidator())
        self.nyEdit = QLineEdit()
        self.nyEdit.setText(str(ny))
        self.nyEdit.setValidator(QIntValidator())
        self.x = np.linspace(0,L,nx)
        self.y = np.linspace(0,W,ny)
        self.X, self.Y = np.meshgrid(self.x,self.y)

        # Sinusoids characteristics
        self.Amplitude = QLineEdit()
        self.Amplitude.setText(str(Amplitude))
        # self.Amplitude.setValidator(DoubleValidator)
        self.wavelength = QLineEdit() # [km]
        self.wavelength.setText(str(wavelength))
        # self.wavelength.setValidator(DoubleValidator)
        self.Mean_topo = QLineEdit() # [m]
        self.Mean_topo.setText(str(Mean_topo))
        self.Mean_topo.setValidator(DoubleValidator)
        self.Offset =QLineEdit() # [m]
        self.Offset.setText(str(Offset))
        # self.Offset.setValidator(DoubleValidator)
        self.PhaseShift = QLineEdit()
        self.PhaseShift.setText(str(phaseShift))
        
        # Get default topography
        self.update_Topography()
        
        # Build layouts and window
        self.window = U_interface.WinExtraParameters(width=800, height=600,text='Build topography')
        # Define two splitters:
            # One for the parameters
            # One for the plotting area
        self.TopoSplitter = U_interface.QCustomSplitter(Qt.Horizontal)
        self.ParamSplitter = QWidget()
        self.PlotSplitter = QWidget()
        
        # First plot area
        self.VBox = QVBoxLayout()
        self.VBox.addWidget(self.canvas.Frame)
        self.PlotSplitter.setLayout(self.VBox)
        
        # Then, parameters area
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.ModelDimensions, 0, 0, 1, 5)
        self.Grid.addWidget(self.Length, 1, 1, 1, 1)
        self.Grid.addWidget(self.Width, 1, 2, 1, 1)
        self.Grid.addWidget(self.resX, 3, 1, 1, 1)
        self.Grid.addWidget(self.resY, 3, 2, 1, 1)
        self.Grid.addWidget(self.L, 2, 1, 1, 1)
        self.Grid.addWidget(self.W, 2, 2, 1, 1)
        self.Grid.addWidget(self.nxEdit, 4, 1, 1, 1)
        self.Grid.addWidget(self.nyEdit, 4, 2, 1, 1)
        self.Grid.addWidget(self.ModelCharacteristics, 5, 0, 1, 5)
        self.Grid.addWidget(self.AmpLabel, 6, 1, 1, 1)
        self.Grid.addWidget(self.WaveLabel, 6, 2, 1, 1)
        # self.Grid.addWidget(self.MTopoLabel, 5, 3, 1, 1)
        self.Grid.addWidget(self.offsetLabel, 8, 1, 1, 1)
        self.Grid.addWidget(self.Amplitude, 7, 1, 1, 1)
        self.Grid.addWidget(self.wavelength, 7, 2, 1, 1)
        # self.Grid.addWidget(self.Mean_topo, 6, 3, 1, 1)
        self.Grid.addWidget(self.Offset, 9, 1, 1, 1)
        self.Grid.addWidget(self.PhaseShiftLabel, 8, 2, 1, 1)
        self.Grid.addWidget(self.PhaseShift, 9, 2, 1, 1)
        self.Grid.setRowStretch(10,2)
        self.ParamSplitter.setLayout(self.Grid)
        
        # Combine splitters
        self.TopoSplitter.addWidget(self.ParamSplitter)
        self.TopoSplitter.addWidget(self.PlotSplitter)
        
        # Add splitter to the window
        self.window.layout.addWidget(self.TopoSplitter, 0, 0, 1, 3)
        Q = QWidget()
        Q.setLayout(self.window.layout)
        self.window.setCentralWidget(Q)
        conf.WindowsOpen.append(self.window)
        self.window.show()

        # Now write files for Pecube
        self.window.OkButton.clicked.connect(lambda: self.update_Inputs())
        
        # Signals
        self.Amplitude.editingFinished.connect(lambda: self.update_Topography())
        self.wavelength.editingFinished.connect(lambda: self.update_Topography())
        # self.Mean_topo.editingFinished.connect(lambda: self.update_Topography())
        self.Offset.editingFinished.connect(lambda: self.update_Topography())
        self.PhaseShift.editingFinished.connect(lambda: self.update_Topography())
        self.L.editingFinished.connect(lambda: self.update_Topography())
        self.W.editingFinished.connect(lambda: self.update_Topography())
        self.nxEdit.editingFinished.connect(lambda: self.update_Topography())
        self.nyEdit.editingFinished.connect(lambda: self.update_Topography())
        
    #----------------------------------------------
    def update_Topography(self):
        """ To update the topography displayed according to new input parameters
        values.
        """
        # First clear previous topography
        self.canvas.plotter.clear()
        
        # New mesh dimensions
        self.x = np.linspace(0,float(self.L.text()),int(self.nxEdit.text()))
        self.y = np.linspace(0,float(self.W.text()),int(self.nyEdit.text()))
        self.X, self.Y = np.meshgrid(self.x,self.y)
        
        # Topography
        try:
            float(self.wavelength.text())
            float(self.Amplitude.text())
            float(self.Offset.text())
            float(self.PhaseShift.text())
        except:
            QErrorMessage(self).showMessage("Wavelength defined for inversion or batch. Cannot update topography...")
            return

        self.topo = float(self.Amplitude.text())/2 * np.sin(2 * np.pi *
                    (self.X+float(self.PhaseShift.text())) / float(self.wavelength.text())) + float(self.Offset.text())
        
        # Set the grid
        xg = np.reshape(self.X,(self.X.shape[0]*self.X.shape[1]))
        yg = np.reshape(self.Y,(self.X.shape[0]*self.X.shape[1]))
        zg = np.reshape(self.topo,(self.X.shape[0]*self.X.shape[1]))
        self.coordArray = np.transpose(np.asarray([xg,yg,zg]))
        minZ = np.min(zg)
        maxZ = np.max(zg)
        
        grid = pv.StructuredGrid()
        # Set the coordinates from the numpy array
        grid.points = self.coordArray
        # set the dimensions
        grid.dimensions = [self.X.shape[1], self.X.shape[0], 1]
        grid.point_data['Elevation'] = self.coordArray[:,2]
        
        # Add mesh to canva
        self.Mesh = self.canvas.plotter.add_mesh(grid,show_edges=False,cmap='gist_earth',
                            clim=[minZ,maxZ],roughness=1,scalar_bar_args=self.sargs,
                            smooth_shading=True,ambient=0.3)
        self.canvas.plotter.show_bounds(grid='front',location='outer',all_edges=False,color=conf.Colors3Dplot['axes'],
                      xtitle="Distance(km)",ytitle="Distance (km)",ztitle='Elevation (km)')
        #self.canvas.plotter.set_scale(1,1,float(self.ReliefFactor.text()))
        # Remove lights
        self.canvas.plotter.remove_all_lights(True)
        light = pv.Light()
        light.set_direction_angle(30,-20)
        self.canvas.plotter.add_light(light)
        
    #---------------------------------------------------
    def update_Inputs(self):
        """
        To update main inputs once we close the window. Such as nx, ny, dlon and dlat.
        
        """
        # First close plot properly
        self.canvas.plotter.clear()
        self.canvas.plotter.close()
        self.canvas.plotter.deep_clean()
        
        # Then update main parameters
        self.Param.nxEdit2.setText(self.nxEdit.text())
        self.Param.nyEdit3.setText(self.nyEdit.text())
        # Convert to lat long
        dist = float(self.L.text())/(float(self.nxEdit.text())-1)
        self.Param.dLonEdit6.setText(str(GIS.km2lon(dist, 0.0)))
        dist = float(self.W.text())/(float(self.nyEdit.text())-1)
        self.Param.dLatEdit7.setText(str(GIS.km2lat(dist)))
        self.Param.topo_file_nameEdit1.setText(conf.Variable_names['Synthetic topo file name'])
        
        # REconstruct topo according to lon lat
        self.nx = int(self.Param.nxEdit2.text())
        self.ny = int(self.Param.nyEdit3.text())
        self.dx = float(self.Param.dLonEdit6.text())
        self.dy = float(self.Param.dLatEdit7.text())
        self.x0 = float(self.Param.lon0Edit4.text())
        self.y0 = float(self.Param.lat0Edit5.text())
        self.xl = (self.nx-1)*self.dx*111.11*cos((self.y0+self.dy*self.ny/2)*3.141592654/180)
        self.yl = (self.ny-1)*self.dy*111.11
        x = np.linspace(self.x0,self.x0+self.xl,int((self.nx)))
        y = np.linspace(self.y0,self.y0+self.yl,int((self.ny)))
        X,Y = np.meshgrid(x,y)
        try:
            NewTopo = float(self.Amplitude.text())/2 * np.sin(2 * np.pi *
                    (X+float(self.PhaseShift.text())) / float(self.wavelength.text())) + float(self.Offset.text())
            
        except:
            self.Param.input_parameters.update({'topo_wavelength':self.wavelength.text()})
            self.Param.input_parameters.update({'topo_amp':self.Amplitude.text()})
            if float(self.Offset.text()) > 0:
                self.Param.input_parameters.update({'topo_offset':self.Offset.text()})
            print('Problem in updating inputs, assume wavelength is inverted. PecubeGUI.py, line 4662.')
            # Close the window
            self.Param.input_parameters.update({'topo_wavelength':self.wavelength.text()})
            self.Param.input_parameters.update({'topo_amp':self.Amplitude.text()})
            self.Param.input_parameters.update({'topo_phase':self.PhaseShift.text()})
            if float(self.Offset.text()) > 0:
                self.Param.input_parameters.update({'topo_offset':self.Offset.text()})
            self.window.close()
            return
        
        try:
            float(self.wavelength.text())
            path2file = os.path.join(conf.PecubeFolderPath,self.Param.folderName,"data",conf.Variable_names['Synthetic topo file name'])
            with open(path2file, 'w') as file:
                np.savetxt(file,np.reshape(NewTopo*1e3,(len(x)*len(y),1)),delimiter='\n')
                file.seek(0, os.SEEK_END)
                file.seek(file.tell() - conf.NEWLINE_SIZE_IN_BYTES, os.SEEK_SET)
                file.truncate()
        except:
            print('Topography wavelength set for inversion or batch.')
            
        # Close the window
        self.Param.input_parameters.update({'topo_wavelength':self.wavelength.text()})
        self.Param.input_parameters.update({'topo_amp':self.Amplitude.text()})
        self.Param.input_parameters.update({'topo_phase':self.PhaseShift.text()})
        self.Param.wavelength = float(self.wavelength.text())
        self.Param.topo_offset = float(self.Offset.text())
        self.Param.topo_amp = float(self.Amplitude.text())
        self.Param.topo_phase = float(self.PhaseShift.text())
        if float(self.Offset.text()) > 0:
            self.Param.input_parameters.update({'topo_offset':self.Offset.text()})
        self.window.close()
        
        
#################################################################################
#################################################################################
#################################################################################  
class SetTopography(QWidget):
    """ 
      This class contains the parameters and the tools to show
      the topography that will be feed to Pecube.
      It also allows to interactively control the shape of the topography, 
      and its evolution through time.
      Param : all parameters from the WindowsParameters class
      SavedParam : Extra parameters saved in the tables such as topography
    
    """
    def __init__(self, Param):
        super().__init__()
        
        self.Param = Param
        Param.get_topo = None
        self.MainWin = Param.MainWin
        self.cb = 0
        self.update = 0
        self.minElev  = 0
        self.maxElev = 3

        # Define two splitters:
            # One for the parameters
            # One for the plotting area
        self.TopoSplitter = U_interface.QCustomSplitter(Qt.Vertical)
        self.ParamSplitter = QWidget()
        self.PlotSplitter = QWidget()
        self.sargs = dict(height=0.25, vertical=True, position_x=0.05,position_y=0.05,interactive=True,
                          n_labels=0,color=conf.Colors3Dplot['colorbarLabels'])
        self.specular = 0.3
        self.RFstored = 1
        
        # Define plotting area
        self.plotSpace =  pgu.Canvas3D()
        
        # Define toolbar
        # self.toolbar = NavigationToolbar(self.plotSpace,self)
        
        # Define parameters
        self.Ninter = int(self.Param.ntimeEdit1.text())+1
        try:
            self.maxtime = float(self.Param.timeTable_dict['time_topo001'])
        except KeyError: #Means the user did not provide the topographic evolution yet
        # Put a default value
            self.maxtime = 1
        # Get Offset and amplitude values from the Time evolution table
        self.OffsetArray = []
        self.AmplArray = []
        self.TimeArray = []
        self.phaseArray = []
        for i in range(self.Ninter):
            self.OffsetArray.append(self.Param.timeTable_dict['offset'+str(i+1)])
            self.AmplArray.append(self.Param.timeTable_dict['amplification'+str(i+1)])
            self.TimeArray.append(self.Param.timeTable_dict['time_topo'+str(i+1)])
            if self.Param.isSynthTopo:
                self.phaseArray.append(self.Param.timeTable_dict['phase'+str(i+1)])
        # To set the topography - h = offset + amplification * h0
        # Define all labels
        self.offsetLabel = QLabel('Offset (km):')
        self.offsetLabel.setFont(conf.font10)
        self.offsetLabel.setAlignment(Qt.AlignCenter)
        self.ampLabel = QLabel("Amplification:")
        self.ampLabel.setFont(conf.font10)
        self.ampLabel.setAlignment(Qt.AlignCenter)
        self.TitleLabel = QLabel("Time evolution topography:")
        self.TitleLabel.setFont(conf.fontLine12)
        self.TimeEvolution = QLabel("Set time evolution:")
        self.TimeEvolution.setFont(conf.fontLine12)
        self.ReliefFactorLabel = QLabel("Relief factor:")
        self.ReliefFactorLabel.setToolTip("Relief factor for visualization purposes only.")
        self.ReliefFactorLabel.setFont(conf.font10)
        self.ReliefFactorLabel.setAlignment(Qt.AlignCenter)
        self.ReliefFactor = QLineEdit()
        self.ReliefFactor.setText('0.02')
        self.ReliefFactor.setValidator(DoubleValidator)
        self.ReliefFactor.setAlignment(Qt.AlignCenter)
        
        # Define the tools
        self.OffsetEdit = QLabel()
        self.OffsetEdit.setText(self.OffsetArray[0])
        self.OffsetEdit.setAlignment(Qt.AlignCenter)
        self.AmpEdit = QLabel()
        self.AmpEdit.setText(str(self.AmplArray[0]))
        self.AmpEdit.setAlignment(Qt.AlignCenter)
        self.TimeEvolSlider = QSlider(Qt.Horizontal)
        self.TimeEvolSlider.setMinimum(0)
        self.TimeEvolSlider.setMaximum(len(self.OffsetArray))
        self.TimeEvolSlider.setSingleStep(1)
        self.TimeValue = QLabel()
        self.TimeValue.setText('0 Ma')
        self.TimeValue.setAlignment(Qt.AlignCenter)
        self.phaseEdit = QLineEdit()
        if self.Param.isSynthTopo:
            self.phaseEdit.setText(str(self.phaseArray[0]))
        else:
            self.phaseEdit.setText('0.0')
        self.istep = 0
        self.EditButton = QPushButton("edit...")
        self.EditButton.clicked.connect(lambda: self.change_visual_settings())
        
        self.minZ = 0
        self.maxZ = 1
        # If time evolution table is edited
        self.Param.topoTable.itemChanged.connect(lambda: self.storeValues())
        
        # Set the layout for the parameters
        Grid1 = QGridLayout()
        Grid1.addWidget(self.ampLabel,0,1)
        Grid1.addWidget(self.offsetLabel,0,2)
        Grid1.addWidget(self.AmpEdit,1,1)
        Grid1.addWidget(self.OffsetEdit,1,2)
        Grid1.addWidget(self.TimeEvolution,2,0)
        Grid1.addWidget(self.TimeValue, 3,1,1,2)
        Grid1.addWidget(self.TimeEvolSlider,4,1,1,2)
        Grid1.addWidget(self.ReliefFactorLabel,3,3)
        Grid1.addWidget(self.ReliefFactor,4,3)
        Grid1.addWidget(self.EditButton,4,4)
        Grid1.setRowStretch(4,1)
        Grid1.setColumnStretch(0,1)
        Grid1.setColumnStretch(4,1)
        
        self.ParamSplitter.setLayout(Grid1)
        
        # When text boxes are edited by the user
        # self.AmpEdit.textChanged.connect(lambda: self.updateTopo())
        # self.OffsetEdit.textChanged.connect(lambda: self.updateTopo())
        self.TimeEvolSlider.valueChanged.connect(self.updateTime)
        # self.TimeEvolSlider.valueChanged.connect(self.updateTopo)
        self.ReliefFactor.editingFinished.connect(lambda: self.updateTopo())
        
        # Get data 
        self.updateInitTopo()
        
        # Gather the layouts
        self.VBox1 = QVBoxLayout()
        self.VBox2 = QVBoxLayout()
        # self.VBox2.addWidget(self.toolbar)
        self.VBox2.addWidget(self.plotSpace.Frame)
        
        # Add the Layouts to the Splitters
        # self.ParamSplitter.setLayout(self.VBox1)
        self.PlotSplitter.setLayout(self.VBox2)
        self.TopoSplitter.addWidget(self.ParamSplitter)
        self.TopoSplitter.addWidget(self.PlotSplitter)
        # Add Splitters to the main Tab
        self.VBox3 = QVBoxLayout()
        self.VBox3.addWidget(self.TopoSplitter)
        Q = QWidget()
        Q.setLayout(self.VBox3)
        self.MainWin.StackedTabs[self.MainWin.CurrentTabs].addTab(Q,"Topography")
        
    
    #-------------------------------------------  
    def change_visual_settings(self):
        """Change the parameter of 3D visualization"""

        # Make window widget
        self.window = U_interface.WinExtraParameters(200, 100)
        self.window.OkButton.clicked.connect(lambda: self.window.close())
        self.editParam = TopoVizEdit(self)
        self.window.layout.addLayout(self.editParam.Grid, 0, 0, 1, 3)
        Q2 = QWidget()
        Q2.setLayout(self.window.layout)
        self.window.setCentralWidget(Q2)
        conf.WindowsOpen.append(self.window)
        self.window.show()
        self.updateTopo()

    #-------------------------------------------    
    def check_inversion(self):
        """Check for inversion parameters, and give up updates if inversion."""
        try:
            self.nx = int(self.Param.nxEdit2.text())
            self.ny = int(self.Param.nyEdit3.text())
            self.dx = float(self.Param.dLonEdit6.text())
            self.dy = float(self.Param.dLatEdit7.text())
            self.x0 = float(self.Param.lon0Edit4.text())
            self.y0 = float(self.Param.lat0Edit5.text())
            self.nskip = int(self.Param.nskipEdit8.text())
            self.ntime = int(self.Param.ntimeEdit1.text())
            float(self.AmpEdit.text())
            float(self.OffsetEdit.text())
            self.OffsetArray = []
            self.AmplArray = []
            self.TimeArray = []
            self.phaseArray = []
            for i in range(int(self.Param.ntimeEdit1.text())+1):
                self.OffsetArray.append(self.Param.timeTable_dict['offset'+str(i+1)])
                self.AmplArray.append(self.Param.timeTable_dict['amplification'+str(i+1)])
                self.TimeArray.append(self.Param.timeTable_dict['time_topo'+str(i+1)])
                amp = float(self.Param.timeTable_dict['amplification'+str(i+1)])
                offset = float(self.Param.timeTable_dict['offset'+str(i+1)])
                time = float(self.Param.timeTable_dict['time_topo'+str(i+1)])
                if self.Param.isSynthTopo:
                    self.phaseArray.append(self.Param.timeTable_dict['phase'+str(i+1)])
            self.IncisionStart = float(self.Param.IncisionStart.text())
            self.IncisionStop = float(self.Param.IncisionStop.text())
            self.IncisionDepth = float(self.Param.IncisionDepth.text())
            self.Tauh = float(self.Param.Tauh.text())
            return 0
        except ValueError as E:
            print(str(E))
            print("Topography cannot be updated. Parameters for inversion has been provided.")
            return 1
        except KeyError:
            return 0
        
    #-------------------------------------------    
    def updateInitTopo(self):
        """ Update the topo when click on 'Show topo' button. """
        
        # First check if all necessary parameter can be interpreted as value (no inversion)
        Inversion = self.check_inversion()
        print('Update Initial Topo')
        if Inversion:
            return

            
        # Clear plot
        self.plotSpace.plotter.clear()

        # Get the dimension of the grid
        # self.nx = int((int(self.Param.nxEdit2.text())-1)/nskip+1)
        # self.ny = int((int(self.Param.nyEdit3.text())-1)/nskip+1)
        self.nx = int(self.Param.nxEdit2.text())
        self.ny = int(self.Param.nyEdit3.text())
        
        #### Set the plotting space ####
        self.dx = float(self.Param.dLonEdit6.text())
        self.dy = float(self.Param.dLatEdit7.text())
        self.y0 = float(self.Param.lat0Edit5.text())
        self.x0 = float(self.Param.lon0Edit4.text())
        self.nskip = int(self.Param.nskipEdit8.text())
        
        # Check the input topography
        if self.Param.topo_file_nameEdit1.text() == 'SPM/':
            # Get elevation from iSOSIA files
            self.topoFile = os.path.join(self.Param.PFolder,"data","SPM","topo000")
            try:
                elevation = np.loadtxt(self.topoFile)
            except OSError:
                QErrorMessage().showMessage("File not found: "+self.topoFile)
                return
            elevation = elevation.astype(float)
            self.label = ["Y distance (km)", "X distance (km)"]

        # Is it a DEM? (i.e., .tif file)
        elif int(self.Param.Param.DParameters['Raster_flag']) and os.path.exists(os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text()+'.tif')):
            path = os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text()+'.tif')
            try:
                print('Goes to elevation raster')
                elevation,nx,ny,self.dx,self.dy,self.x0,self.y0 = GIS.load_raster_file(path)
                self.label = ["Latitude (°)", "Longitude (°)"]
            except AttributeError:
                U_interface.MessageBox(text1="The file"+self.Param.topo_file_nameEdit1.text()+'.tif has not been found. Please, provide it.')
                return
            
        elif self.Param.topo_file_nameEdit1.text()[-1] != '/' and self.Param.topo_file_nameEdit1.text() != 'Nil':
            # Get elevation from another spm or a DEM
            self.topoFile = os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text())
            try:
                elevation = np.loadtxt(self.topoFile)
            except FileNotFoundError:
                QErrorMessage().showMessage("Topographic file not found, please reload it.")
            elevation = elevation.astype(float)
            self.label = ["Y distance (km)", "X distance (km)"]
            
        else:# Means flat topography
            elevation = np.array(np.zeros((self.ny,self.nx)))
            self.label = ["Y distance (km)", "X distance (km)"]
        # Ensure the dimensions are correct
        try:
            self.Z = np.reshape(elevation,(self.ny,self.nx))
        except ValueError:
            QErrorMessage().showMessage("The dimensions of the grid are wrong. Please, check the values provided.")
            return
        
        # Build the grid
        if int(self.Param.Param.DParameters['Raster_flag']) or os.path.exists(os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text()+'.tif')):
            # Better to see in lat,long for the tif files
            self.xl = (self.nx-1)*self.dx
            self.yl = (self.ny-1)*self.dy
        else:
            self.xl = (self.nx-1)*self.dx*111.11*cos((self.y0+self.dy*self.ny/2)*3.141592654/180)
            self.yl = (self.ny-1)*self.dy*111.11
            self.Z = np.flip(self.Z,1)

        # Compute topoRef
        topoRef = 0
        if self.Param.RefElevCombo.currentIndex() == 0: #From sea level
            topoRef = 0
        elif self.Param.RefElevCombo.currentIndex() == 1: #From minimum
            topoRef = np.min(self.Z)
        elif self.Param.RefElevCombo.currentIndex() == 2: #From maximum
            topoRef = np.max(self.Z)
        elif self.Param.RefElevCombo.currentIndex() == 3: #From mean
            topoRef = np.mean(self.Z)
        elif self.Param.RefElevCombo.currentIndex() == 4: #From custom elevation
            topoRef = float(self.Param.RefCustom.text())
        elif self.Param.RefElevCombo.currentIndex() == 5: #From smooth
            topoRef = float(self.Param.RefCustom.text())

        # Handle nskip
        self.x = np.linspace(self.x0,self.x0+self.xl,int((self.nx-1)/self.nskip+1))
        self.y = np.linspace(self.y0,self.y0+self.yl,int((self.ny-1)/self.nskip+1))
        self.X,self.Y = np.meshgrid(self.x,self.y)
        elevation = self.Z[::self.nskip,::self.nskip]

        # Set the canvas
        self.minZ = np.min(self.Z)/1e3
        self.maxZ = np.max(self.Z)/1e3
        self.minX = np.min(self.x)
        self.maxX = np.max(self.x)
        self.minY = np.min(self.y)
        self.maxY = np.max(self.y)
        self.minElev = self.minZ
        self.maxElev = self.maxZ

        # Set the grid
        xg = np.reshape(self.X,(self.X.shape[0]*self.X.shape[1]))
        yg = np.reshape(self.Y,(self.X.shape[0]*self.X.shape[1]))
        zg = np.reshape(elevation/1e3,(self.X.shape[0]*self.X.shape[1]))
        self.coordArray = np.transpose(np.asarray([xg,yg,zg]))
        minZ = np.min(zg)
        maxZ = np.max(zg)
        
        grid = pv.StructuredGrid()
        # Set the coordinates from the numpy array
        grid.points = self.coordArray
        # set the dimensions
        grid.dimensions = [self.X.shape[1], self.X.shape[0], 1]
        grid.point_data['Elevation'] = self.coordArray[:,2]
        # Add mesh to canva
        rangeElev = self.maxElev - self.minElev
        dz = rangeElev/3
        dz2 = (self.maxElev - self.minElev)/3
        annotations = {self.minElev:round(self.minElev,1),self.minElev+dz:round(self.minElev+dz,1),self.minElev+2*dz:round(self.minElev+2*dz,1),self.maxElev:round(self.maxElev,1)}
        self.Mesh = self.plotSpace.plotter.add_mesh(grid,show_edges=False,cmap='gist_earth',
                            clim=[minZ,maxZ],roughness=1,scalars='Elevation',
                            annotations=annotations, scalar_bar_args=self.sargs,smooth_shading=True,ambient=0.3)
        self.plotSpace.plotter.set_scale(zscale= float(self.ReliefFactor.text()))
        self.plotSpace.plotter.update_bounds_axes()
        self.plotSpace.plotter.show_bounds(mesh=self.Mesh,
                                           location='outer',ticks='outside',color=conf.Colors3Dplot['axes'],
                                           all_edges=True,xtitle=self.label[1],ytitle=self.label[0],ztitle='Elevation (km)',
                                           )
        # Remove lights
        self.plotSpace.plotter.remove_all_lights(True)
        light = pv.Light()
        light.set_direction_angle(30,-20)
        self.plotSpace.plotter.add_light(light)
        
        
        self.update = 1
        
   #----------------------------------------------------------
    def closeEvent(self, event):
        self.deleteWidgets()
        event.accept()
        
   #----------------------------------------------------------
    def deleteWidgets(self):
        self.TitleLabel.setParent(None)
        self.TitleLabel.deleteLater()
        self.ampLabel.deleteLater()
        self.offsetLabel.deleteLater()
        self.AmpEdit.deleteLater()
        self.OffsetEdit.deleteLater()
        self.TimeEvolution.deleteLater()
        self.TimeValue.deleteLater()
        self.TimeEvolSlider.deleteLater()
        # self.toolbar.deleteLater()
        self.plotSpace.deleteLater()
        self.ParamSplitter.deleteLater()
        self.PlotSplitter.deleteLater()
                
    #----------------------------------------------------------
    def updateTopo(self):
        """ 
          Update the topography changed by the user.
          The user can either change the amplitude and/or the offset
          and for each time step.
          We therefore need to define an amplitude and offset value for 
          each time step.
        """

        # check if all necessary parameter can be interpreted as value (no inversion)
        Inversion = self.check_inversion()
        if Inversion:
            print('update_Topo')
            return

        minZinit = self.minZ
        maxZinit = self.maxZ
        nskip = int(self.Param.nskipEdit8.text())
        if not nskip:
            return
        #Headward erosion ?
        self.IncisionStart = float(self.Param.IncisionStart.text())
        self.IncisionStop = float(self.Param.IncisionStop.text())
        self.IncisionDepth = float(self.Param.IncisionDepth.text())
        self.Tauh = float(self.Param.Tauh.text())
        
        # Check if files from iSOSIA
        if self.Param.topo_file_nameEdit1.text() == 'SPM/':
            # Get number of output files
            nfiles = int(len(os.listdir(os.path.join(self.Param.PFolder,"data","SPM")))/3)
            # Compare to the number of time steps defined by the user (ntime)
            ntime = int(self.Param.ntimeEdit1.text())
            
            # if output file number < ntime take the first topo file (fluvial type)
            # spm topo files are always the last one in the list dir
            if nfiles <= ntime: # Means preHistory
                # Get the difference
                ndiff = ntime - nfiles
                # this gives the extra number of steps/topos to create from topo 0
                if self.TimeEvolSlider.value() <= ndiff+1: #take the first topo file
                    self.topoFile = os.path.join(self.Param.PFolder,"data","SPM","topo000")
                elif ndiff+1 < self.TimeEvolSlider.value() <  10+ndiff+1:#take other files
                # Limit pre-history to 9 steps
                    self.topoFile = os.path.join(self.Param.PFolder,"data","SPM","topo00"+str(int(self.TimeEvolSlider.value()-(ndiff+1))))  
                else:
                    self.topoFile = os.path.join(self.Param.PFolder,"data","SPM","topo0"+str(int(self.TimeEvolSlider.value()-(ndiff+1))))  
            else: # nfiles > ntime
                # before the spm files, which topo file to take?
                value = self.TimeEvolSlider.value()
                if value < 10:
                    zeroline = '00'
                elif value >= 10 and value < 100:
                    zeroline = '0'
                else: 
                    zeroline = ''
                self.topoFile = os.path.join(self.Param.PFolder,"data","SPM","topo"+zeroline+str(self.TimeEvolSlider.value()))
            
            try:
                elevation = np.loadtxt(self.topoFile)
            except OSError:
                QErrorMessage().showMessage("File not found: "+self.topoFile)
                return
            topography = elevation.astype(float)
            topography = np.reshape(elevation,(self.ny,self.nx))
            # topography = np.flip(topography,1)
            
        else: # Other than iSOSIA files
            # Get the previous topography
            try:
                topography = self.Z
            except AttributeError:
                print('Error on Z in set_topography.')
                return 

        # Compute the new Topo
        topoRef = 0
        if self.Param.RefElevCombo.currentIndex() == 0: #From sea level
            topoRef = 0
        elif self.Param.RefElevCombo.currentIndex() == 1: #From minimum
            topoRef = np.min(topography)
        elif self.Param.RefElevCombo.currentIndex() == 2: #From maximum
            topoRef = np.max(topography)
        elif self.Param.RefElevCombo.currentIndex() == 3: #From mean
            topoRef = np.mean(topography)
        elif self.Param.RefElevCombo.currentIndex() == 4: #From custom elevation
            topoRef = float(self.Param.RefCustom.text())
        elif self.Param.RefElevCombo.currentIndex() == 5: #From smooth
            topoRef = float(self.Param.RefCustom.text())

        # Handle nskip
        self.x = np.linspace(self.x0,self.x0+self.xl,int((self.nx-1)/nskip+1))
        self.y = np.linspace(self.y0,self.y0+self.yl,int((self.ny-1)/nskip+1))
        self.X,self.Y = np.meshgrid(self.x,self.y)
        elevation = np.asarray(topography[::nskip,::nskip])
        
        # check if smoothing is activated
        Radius = float(self.Param.SmoothingFactor.text())
        if Radius > 0: # -----------------------------------------------------------------
            topoRef = np.ones(np.shape(elevation))*topoMean
            # Compute smoothed topography in range of radius
            dx = self.x[1] - self.x[0] #has to be in km
            RadiusPx = int(Radius/(dx*111.11))
            # SmoothTopo = pgu.moving_average(elevation,RadiusPx)
            # print(topoRef)
            Relief = topoRef - elevation 
            NewTopo = np.zeros(np.shape(elevation))
            indexNeg = elevation > topoRef 
            # print(elevation)
            # print(Relief)
            try:
                NewTopo[indexNeg] = topoRef[indexNeg]  - (1/float(self.AmpEdit.text()) * (topoRef[indexNeg] - elevation[indexNeg])) + float(self.OffsetEdit.text())*1e3
            except ZeroDivisionError:
                NewTopo[indexNeg] = topoRef[indexNeg]  - (float(self.AmpEdit.text()) * (topoRef[indexNeg] -elevation[indexNeg])) + float(self.OffsetEdit.text())*1e3
            indexPos = elevation <= topoRef
            NewTopo[indexPos] = topoRef[indexPos]  - (float(self.AmpEdit.text()) * (topoRef[indexPos] -elevation[indexPos])) + float(self.OffsetEdit.text())*1e3
            
            # form plateau if amplitude too high
            # indexPos = NewTopo > np.max(topography)
            # NewTopo[indexPos] = np.max(topography)
        # Headward propagation of erosion ?
        elif self.Param.DoHeadward.isChecked() == True: # --------------------------------
            # First get topography at time of incision start
            Start_incision = float(self.Param.IncisionStart.text())
            Stop_incision = float(self.Param.IncisionStop.text())
            Depth_incision = float(self.Param.IncisionDepth.text())
            tauh = float(self.Param.Tauh.text())
            [TopoInitIncision, Zmin, Zmax, amp, offset] = update_topo.get_initial_topo(topoRef, elevation,np.float_(self.AmplArray),np.float_(self.OffsetArray)*1e3,np.float_(self.TimeArray),
                                                             Start_incision)
            # then update topography according to time
            if self.istep == 0:
                time_prev = float(self.TimeArray[0])
            else:
                time_prev = float(self.TimeArray[self.istep-1])
            NewTopo = update_topo.Headward_propagation(topoRef, elevation, TopoInitIncision, np.float_(self.AmplArray),np.float_(self.OffsetArray)*1e3, self.istep, np.float_(self.TimeArray),
                                          Start_incision,Stop_incision,Depth_incision,tauh, Zmax, Zmin)
            
        # Simple topographic evolution function
        else: # -----------------------------------------------------------------------
            if self.Param.isSynthTopo: # We use a synthetic sinusoidal topography
                # print('Phase shift: ', self.phaseEdit.text(),self.Param.wavelength)
                amp = float(self.Param.topo_amp)*1e3
                topo_phase = float(self.Param.topo_phase)*1e3
                wavelength = float(self.Param.wavelength)*1e3
                offset = float(self.Param.topo_offset)*1e3
                NewTopo = topoRef - (float(self.AmpEdit.text()) * (topoRef - (amp/2 * np.sin(2 * np.pi *
                            (self.X*1e3+ topo_phase+float(self.phaseEdit.text())) / wavelength)+ offset))) + float(self.OffsetEdit.text())*1e3
                
            elif self.Param.RefElevCombo.currentIndex() == 5:
                NewTopo2 = topoRef - (float(self.AmpEdit.text()) * (topoRef-(elevation))) + float(self.OffsetEdit.text())*1e3
                if not float(self.AmpEdit.text()) == 0:
                    TopoScaled =  topoRef - (1/float(self.AmpEdit.text()) * (topoRef-(elevation))) + float(self.OffsetEdit.text())*1e3
                else:
                    TopoScaled =  topoRef - (float(self.AmpEdit.text()) * (topoRef-(elevation))) + float(self.OffsetEdit.text())*1e3
                NewTopo = np.where(elevation<topoRef,TopoScaled,NewTopo2)
                
            else:
                NewTopo = topoRef - (float(self.AmpEdit.text()) * (topoRef-(elevation))) + float(self.OffsetEdit.text())*1e3
        # ---------------------------------------------------------------------------------
        
        NewTopo = NewTopo/1e3
        minZ = self.minZ
        maxZ = self.maxZ

        # First erase previous plot
        #self.plotSpace.plotter.remove_actor(self.Mesh)
        
        # Set the grid
                # check that DEM size did not change
        try:
            np.reshape(NewTopo,(self.X.shape[0]*self.X.shape[1]))
        except ValueError:
            print('back to initial topo update')
            self.updateInitTopo()   
        xg = np.reshape(self.X,(self.X.shape[0]*self.X.shape[1]))
        yg = np.reshape(self.Y,(self.X.shape[0]*self.X.shape[1]))
        zg = np.reshape(NewTopo,(self.X.shape[0]*self.X.shape[1]))
        self.coordArray = np.transpose(np.asarray([xg,yg,zg]))
        
        grid = pv.StructuredGrid()
        # Set the coordinates from the numpy array
        grid.points = self.coordArray
        # Set the dimensions
        grid.dimensions = [self.X.shape[1], self.X.shape[0], 1]
        grid.point_data['Elevation'] = self.coordArray[:,2]

        # Add mesh to canva
        rangeElev = self.maxElev - self.minElev
        dz = rangeElev/3
        dz2 = (self.maxElev - self.minElev)/3
        annotations = {self.minElev:round(self.minElev,1), self.minElev+dz:round(self.minElev+dz2,1), self.minElev+2*dz:round(self.minElev+2*dz2,1) ,self.maxElev:round(self.maxElev,1)}
        # First clear plotter
        self.camerapos = self.plotSpace.plotter.camera_position
        self.plotSpace.plotter.disable()
        self.plotSpace.plotter.clear()
    
    
        if float(self.ReliefFactor.text()) != self.RFstored:
            self.RFstored = float(self.ReliefFactor.text())
            self.Mesh = self.plotSpace.plotter.add_mesh(grid,show_edges=False,cmap='gist_earth',
                            clim=[self.minElev,self.maxElev],roughness=1,reset_camera=True,scalars='Elevation',
                            annotations=annotations,scalar_bar_args=self.sargs,smooth_shading=True,ambient=0.3)
        else:
            self.Mesh = self.plotSpace.plotter.add_mesh(grid,show_edges=False,cmap='gist_earth',
                        clim=[self.minElev,self.maxElev],roughness=1,reset_camera=False,scalars='Elevation',
                        annotations=annotations,scalar_bar_args=self.sargs,smooth_shading=True,ambient=0.3)

        # self.plotSpace.plotter.show_bounds(mesh=self.Mesh,bounds = [self.minX,self.maxX,self.minY,self.maxY,0,maxZ],
        #                             location='outer',ticks='outside',color=conf.Colors3Dplot['axes'],
        #                             all_edges=True,xtitle=self.label[1],ytitle=self.label[0],ztitle='Elevation (km)',
        #                             )
        self.plotSpace.plotter.set_scale(zscale=float(self.ReliefFactor.text()))
        self.plotSpace.plotter.show_bounds(mesh=self.Mesh,
                                    location='outer',ticks='outside',color=conf.Colors3Dplot['axes'],
                                    all_edges=True,xtitle=self.label[1],ytitle=self.label[0],ztitle='Elevation (km)',
                                    )
        
        self.plotSpace.plotter.enable()
        
   
        # self.plotSpace.plotter.update_bounds_axes()
        # self.plotSpace.plotter.camera_position = self.camerapos
        self.plotSpace.plotter.remove_all_lights(True)
        light = pv.Light()
        light.set_direction_angle(30,-20)
        self.plotSpace.plotter.add_light(light)
        self.plotSpace.plotter.lighting = 'three lights'
        
    #----------------------------------------------------------
    def storeValues(self):
        """ Store the values when the topo table is updated."""

        # check if all necessary parameter can be interpreted as value (no inversion)
        Inversion = self.check_inversion()
        if Inversion:
            print('function storeValues')
            QErrorMessage(self).showMessage("cannot read topographic value, please check that. Topography is not updated.")
            return
        
        # Update the data
        self.Ninter = int(self.Param.ntimeEdit1.text())+1

        # Read the table
        num_rows, num_cols = self.Param.topoTable.rowCount(), self.Param.topoTable.columnCount()
        temp_dict = {}
        try:
            for row in range(num_rows):
                for col in range(num_cols):
                    item = self.Param.topoTable.item(row, col)
                    if item is None:
                        break
                    if col == 0:
                        temp_dict['time_topo'+str(row+1)] = str(item.text())
                    elif col == 1:
                        temp_dict['amplification'+str(row+1)] = float(item.text())
                    elif col == 2:
                        temp_dict['offset'+str(row+1)] = float(item.text())
                    elif col == 3:
                        temp_dict['output'+str(row+1)] = int(item.text())
                    elif col == 4:
                        temp_dict['phase'+str(row+1)] = int(item.text())
                if item is None:
                    break
        except ValueError:
            print("Parameters cannot be read. Perhaps inversion is provided. Topography is not updated.")
            return
        try:
            temp_dict['offset'+str(self.Ninter)]
        except KeyError: #Handle when the size of the topoTable change
            return
        
        # Get Offset and amplitude values from the Time evolution table
        self.OffsetArray.clear()
        self.AmplArray.clear()
        self.TimeArray.clear()
        for i in range(self.Ninter):
            self.OffsetArray.append(temp_dict['offset'+str(i+1)])
            self.AmplArray.append(temp_dict['amplification'+str(i+1)])
            self.TimeArray.append(temp_dict['time_topo'+str(i+1)])
            try: 
                self.phaseArray.append(temp_dict['phase'+str(i+1)])
            except:
                pass
                
        
        # First update the value of the time show in the text boxe
        self.updateTime(0)
    
    #----------------------------------------------------------
    def updateTime(self, value):
        """
          This function is called when the user move the slider and wish
          to define the topography for the next time step
        """

        # check if all necessary parameter can be interpreted as value (no inversion)
        Inversion = self.check_inversion()
        if Inversion:
            print('function updateTime')
            QErrorMessage(self).showMessage("cannot read amplification or offset value, please check that. Topography is not updated.")
            return
        # Then update offset and amplitude values
        try:
            amp = float(self.AmplArray[value])
            offset = float(self.OffsetArray[value]) 

            self.AmpEdit.setText(str(self.AmplArray[value]))
            self.OffsetEdit.setText(str(self.OffsetArray[value]))
        except Exception as E:
            print("cannot read amplification or offset value, please check that.", E)
            QErrorMessage(self).showMessage("cannot read amplification or offset value, please check that. Topography is not updated.")
            return
        
        # First update the value of the time show in the text boxe
        try:
            self.TimeValue.setText(self.TimeArray[value] + ' Ma')
        except IndexError:#the user moved the slider without defining time
            self.TimeValue.setText('x Ma')
            
        
        if self.Param.isSynthTopo:
            try:
                self.phaseEdit.setText(str(self.phaseArray[value]))
            except IndexError:
                print('Handle index error in phase array')
                self.phaseEdit.setText('0.0')
        
        # Update the slider
        self.TimeEvolSlider.setMaximum(len(self.OffsetArray)-1)
        self.istep = value
        
        # Update topo
        # try:
        self.updateTopo()
        # except ValueError as VE:
        #      QErrorMessage(self).showMessage('Fail to update topography:' + str(VE) + "<br>" +
        #                                      'Close topography tab and open again.')
        #      return
        

#############################################################################
class TopoVizEdit(QWidget):
    """ 
      This object store edit parameter for the 3D topographic visualization
    
    """
    def __init__(self,param):
        super().__init__()   

        # Parameters
        self.param = param
        self.minElevLabel = QLabel('min. Elevation (km):') 
        self.minElevLabel.setAlignment(Qt.AlignCenter)
        self.minElevLabel.setFont(conf.font10)
        self.minElev = QLineEdit()
        self.minElev.setText(str(param.minElev))
        self.maxElevLabel = QLabel("max. Elevation (km)")
        self.maxElevLabel.setAlignment(Qt.AlignCenter)
        self.maxElevLabel.setFont(conf.font10)
        self.maxElev = QLineEdit()
        self.maxElev.setText(str(param.maxElev))

        # Grid
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.minElevLabel, 0, 0, 1, 1)
        self.Grid.addWidget(self.maxElevLabel, 0, 1, 1, 1)
        self.Grid.addWidget(self.minElev, 1, 0, 1, 1)
        self.Grid.addWidget(self.maxElev, 1, 1, 1, 1)

        # Signals
        self.minElev.editingFinished.connect(lambda: self.update_Elevation_range())
        self.maxElev.editingFinished.connect(lambda: self.update_Elevation_range())

    # ---------------------------------------------------
    def update_Elevation_range(self):
        """" Update elevation range values for the colorbar. """
        self.param.minElev = float(self.minElev.text())
        self.param.maxElev = float(self.maxElev.text())



#################################################################################
#################################################################################
#################################################################################
class ShowFaultGeometry(QWidget):
    """ 
      This class contains the parameters and tools to show
      the Fault geometry defined by the user.
      It also allows to interactively control the shape of the topography, 
      Param : all parameters from the WindowsParameters class
    """
    def __init__(self, Parent):
        super().__init__()
        
        self.Param = Parent
        self.MainWin = Parent.MainWin
        self.cb = 0
        self.updateFault = 0
        self.MainMesh = 0
        self.sargs = dict(height=0.001, width=0.2, vertical=True, position_x=0.95, position_y=0.05,
                      label_font_size=14,interactive=True,color=conf.Colors3Dplot['colorbarLabels'])

        # Define two splitters:
        # One for the parameters
        # One for the plotting area
        self.FaultSplitter = U_interface.QCustomSplitter(Qt.Vertical)
        self.ParamSplitter = QWidget()
        self.PlotSplitter = QWidget()
        self.ShowVelo = QPushButton("Compute velocities")
        self.TitleLabel = QLabel("Time evolution fault kinematic:")
        self.TitleLabel.setFont(conf.fontLine12)
        self.timeLabel = QLabel('0 Ma')
        self.timeLabel.setAlignment(Qt.AlignCenter)
        self.FaultEvolSlider = QSlider(Qt.Horizontal)
        self.FaultEvolSlider.setMinimum(0)
        self.FaultEvolSlider.setMaximum(int(self.Param.nstepiEdit13.text()))
        self.FaultEvolSlider.setSingleStep(1)
        self.scaleVelo = QSlider(Qt.Horizontal)
        self.scaleVelo.setMinimum(1)
        self.scaleVelo.setMaximum(1000)
        self.scaleVelo.setValue(20)
        self.scaleVelo.setSingleStep(10)
        self.scaleVeloLabel = QLabel("scale velocities:")
        self.scaleVeloLabel.setFont(conf.font10)
        self.scaleVeloLabel.setAlignment(Qt.AlignRight)
        self.scaleVeloLabel.setToolTip("scale the velocity vectors by:")
        self.densityVelo = QSlider(Qt.Horizontal)
        self.densityVelo.setMinimum(1)
        self.densityVelo.setMaximum(200)
        self.densityVelo.setValue(10)
        self.densityVelo.setSingleStep(10)
        self.densityVeloLabel = QLabel("Field density:")
        self.densityVeloLabel.setFont(conf.font10)
        self.densityVeloLabel.setAlignment(Qt.AlignRight)
        self.densityVeloLabel.setToolTip("scale the density of velocity vectors by:")
        
        # Define plotting area
        self.plotSpace =  pgu.Canvas3D()
            
        # Get data 
        self.updateFaults()
        
        # Gather the layouts
        self.Grid1 = QGridLayout()
        hbox = QHBoxLayout()
        self.VBox2 = QVBoxLayout()
        self.Grid1.addWidget(self.ShowVelo,0,1,1,1)
        self.Grid1.addWidget(self.TitleLabel,1,0,1,3)
        self.Grid1.addWidget(self.timeLabel,2,1,1,1)
        hbox.addStretch(1)
        hbox.addWidget(self.FaultEvolSlider)
        hbox.addStretch(1)
        hbox2 = QHBoxLayout()
        hbox2.addStretch(1)
        hbox2.addWidget(self.scaleVeloLabel)
        hbox2.addWidget(self.scaleVelo)
        hbox2.addStretch(1)
        hbox3 = QHBoxLayout()
        hbox3.addStretch(1)
        hbox3.addWidget(self.densityVeloLabel)
        hbox3.addWidget(self.densityVelo)
        hbox3.addStretch(1)
        self.VBox2.addLayout(self.Grid1)
        self.VBox2.addLayout(hbox)
        self.VBox2.addLayout(hbox2)
        self.VBox2.addLayout(hbox3)
        self.VBox2.addWidget(self.plotSpace.Frame)
        
        # Add the Layouts to the Splitters
        # self.ParamSplitter.setLayout(self.VBox1)
        self.PlotSplitter.setLayout(self.VBox2)
        # self.FaultSplitter.addWidget(self.ParamSplitter)
        self.FaultSplitter.addWidget(self.PlotSplitter)
        # Add Splitters to the main Tab
        self.VBox3 = QVBoxLayout()
        self.VBox3.addWidget(self.FaultSplitter)
        Q = QWidget()
        Q.setLayout(self.VBox3)
        self.MainWin.StackedTabs[self.MainWin.CurrentTabs].addTab(Q,"Fault geometry")
        self.ShowVelo.clicked.connect(lambda: self.get_velocities())
        self.scaleVelo.valueChanged.connect(self.update_velocities)
        self.densityVelo.valueChanged.connect(self.update_velocities)
        try:
            self.FaultEvolSlider.valueChanged.connect(self.update_velocities)
        except IndexError:
            QErrorMessage(self).showMessage("Cannot update the velocity field. Please, check the Fault kinematic table...")
            
        
    #-------------------------------------------    
    def get_velocities(self):
        """To compute and load velocities along faults """
        # Check if VTK folder exists
        if not os.path.exists(os.path.join(self.Param.PFolder,'VTK')):
            print('Create VTK directory.')
            os.mkdir(os.path.join(self.Param.PFolder,'VTK'))
        # Save input file
        fileName = os.path.join(self.Param.PFolder,'input', 'Pecube.in')
        if os.path.exists(fileName):
            os.remove(fileName)
        file = open(fileName, 'w')
        inputPar = self.Param.parent.ParamTable.input_parameters
        for k, v in inputPar.items():
            file.write(str(k) + ' = ' + str(v) + '\n')
        file.close()
        
        # Go to Pecube Folder 
        os.chdir(self.Param.parent.PecubeFolder)
        self.p = QProcess()
        self.p.readyReadStandardError.connect(self.handle_stderr)
        self.p.finished.connect(self.process_finished)
        if sys.platform == 'darwin' or sys.platform =='linux':
            path = os.path.join('bin','test.sh')
            path = path.replace(os.sep, '/')
            self.p.start('bash', [path, str(self.Param.parent.PNAME)])
        elif sys.platform == 'win32' or sys.platform == 'cygwin':   
            path = os.path.join(conf.FolderPath,'Pecube','bin','Test.exe')
            self.p.setProgram(path)
            self.p.setArguments([str(self.Param.parent.PNAME)])
            self.p.start()
        # Open message for the user to wait until process finished
        self.Message = QMessageBox()
        self.Message.setText("Please wait until computation is finished.")
        self.Message.setIcon(QMessageBox.Information)
        self.Message.setStandardButtons(QMessageBox.Cancel)
        retval = self.Message.exec_()
        if retval == QMessageBox.Cancel:
            self.p.finished.disconnect()
            self.p.kill()
            
    #-------------------------------------------  
    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode(conf.encoding_label)
        
    #-------------------------------------------    
    def process_finished(self):
        """When files have been produced, load vtk files"""
        path = os.path.join(self.Param.PFolder,'VTK')
        filenames = os.listdir(path)
        target1 = [file for file in filenames if file.startswith('velo')]
        
        sargs = dict(height=0.25, width=0.2, vertical=True, position_x=0.1, position_y=0.75,
                      label_font_size=14,interactive=True,title='Velocity (mm/yr)',color=conf.Colors3Dplot['colorbarLabels'])
        # Read VTK file with Pyvista and store in list
        try:
            # data = pv.get_reader(os.path.join(path,target1[0]))
            self.veloList = [pv.get_reader(os.path.join(path,target1[i])) for i in range(len(target1))]
            data = self.veloList[0].read()
            print("velocity = ", self.veloList)
        except IndexError:
            QErrorMessage(self).showMessage("Cannot read velocity file... Ensure to provide at least one output in"+
                          "tab 'Time evolution' and make another trial.")
            return
        # Plot the first step
        data.points[:,2] = data.points[:,2] - data.points[-1,-1]
        # data2 = data.extract_points(data.points[:,0] > 10)
        nPoints = len(data.points[:,0])
        data2 = data.extract_points(np.arange(0,nPoints,nPoints/10,dtype=int))
        geom = pv.Arrow()
        self.velocities = data2.glyph(scale='velo',geom=geom, factor=float(self.scaleVelo.value()/25))
        # Get min and max values
        DataRange = self.velocities.get_data_range()

        # Set 3D interface
        if self.MainMesh:
            self.plotSpace.plotter.remove_actor(self.MainMesh)
        self.MainMesh= self.plotSpace.plotter.add_mesh(
        self.velocities, cmap='plasma',scalar_bar_args=sargs)
        
        if self.Message:
            self.Message.close()
        # # Remove files
        # for file in filenames:
        #     os.remove(os.path.join(path,file))
        # os.chdir(FolderPath)
    
    #-------------------------------------------    
    def update_velocities(self):
        """To update the velocity field to show """
        
        #Remove current velocity field
        # Set 3D interface
        sargs = dict(height=0.25, width=0.2, vertical=True, position_x=0.1, position_y=0.75,
                      label_font_size=14,interactive=True,title='Velocity (mm/yr)',color=conf.Colors3Dplot['colorbarLabels'])
        if self.MainMesh:
            self.plotSpace.plotter.remove_actor(self.MainMesh)
        # Get new velocity field
        try:
            print("veloslider = ", self.FaultEvolSlider.value())
            data = self.veloList[self.FaultEvolSlider.value()].read()
        except AttributeError:
            print('No velocity files found.')
            return
        data.points[:,2] = data.points[:,2] - data.points[-1,-1]
        # data2 = data.extract_points(data.points[:,0] > 10)
        nPoints = len(data.points[:,0])
        try:
            data2 = data.extract_points(np.arange(0,nPoints,round(nPoints/float(self.densityVelo.value())),dtype=int))
        except ZeroDivisionError:
            data2 = data.extract_points(np.arange(0,nPoints,round(nPoints/10),dtype=int))
        geom = pv.Arrow()
        self.velocities = data2.glyph(scale='velo',geom=geom, factor=float(self.scaleVelo.value()/25))
        # Get min and max values
        DataRange = self.velocities.get_data_range()

        # Set 3D interface
        self.MainMesh= self.plotSpace.plotter.add_mesh(
        self.velocities, cmap='plasma',scalar_bar_args=sargs)
        
        #Update time frame
        self.timeLabel.setText(str(data['TIME'][0])+' Ma')
        
        
    #-------------------------------------------    
    def check_inversion(self):
        """Check for inversion parameters, and give up updates if inversion."""
        try:
            self.nx = int(self.Param.nxEdit2.text())
            self.ny = int(self.Param.nyEdit3.text())
            self.dx = float(self.Param.dLonEdit6.text())
            self.dy = float(self.Param.dLatEdit7.text())
            self.x0 = float(self.Param.lon0Edit4.text())
            self.y0 = float(self.Param.lat0Edit5.text())
            self.nskip = int(self.Param.nskipEdit8.text())
            self.ntime = int(self.Param.ntimeEdit1.text())
            self.Thickness = float(self.Param.thicknessEdit1.text())
            self.nfault = int(self.Param.nfaultEdit1.text())
            self.Lon1 = float(self.Param.lon1Edit2.text())
            self.Lat1 = float(self.Param.lat1Edit3.text())
            self.Lon2 = float(self.Param.lon2Edit4.text())
            self.Lat2 = float(self.Param.lat2Edit5.text())
            self.RefCustom2 = float(self.Param.RefCustom.text())
            return 0
        except ValueError:
            print("Topography cannot be updated. Parameters for inversion has been provided.")
            return 1
        
    #----------------------------------------------------------
    def updateFaults(self):
        """ Update the fault geometry changed by the user """
        # Check parameters
        Inversion = self.check_inversion()
        if Inversion:
            print('Wrong value for fault parameter. Cannot update fault geometry.')
            return
        
        # Then define parameters and plots
        self.nx = int(self.Param.nxEdit2.text())
        self.ny = int(self.Param.nyEdit3.text())
        self.dx = float(self.Param.dLonEdit6.text())
        self.dy = float(self.Param.dLatEdit7.text())
        self.x0 = float(self.Param.lon0Edit4.text())
        self.y0 = float(self.Param.lat0Edit5.text())
        self.Thickness = float(self.Param.thicknessEdit1.text())
        y0 = self.y0
        self.nskip = int(self.Param.nskipEdit8.text())
        
        #########################################
        ######  Check the input topography ######
        #########################################
        textlegend = ['Distance (km)','Distance (km)']
        if self.Param.topo_file_nameEdit1.text() == 'SPM/':
            #Get elevation from iSOSIA files
            self.topoFile = os.path.join(self.Param.PFolder,"data","SPM","topo000")
            try:
                elevation = np.loadtxt(self.topoFile)
            except OSError:
                QErrorMessage().showMessage("File not found: "+self.topoFile)
                return
            elevation = elevation.astype(float)
        #Is it a DEM? (i.e., .tif file)
        elif int(self.Param.Param.DParameters['Raster_flag']) or os.path.exists(os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text()+'.tif')):
            path = os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text()+'.tif')
            # data = gdal.Open(path)
            # Get size of the DEM
            try:
                elevation,nx,ny,self.dx,self.dy,self.x0,self.y0 = GIS.load_raster_file(path)
            except AttributeError:
                QErrorMessage().showMessage("The file"+self.Param.topo_file_nameEdit1.text()+'.tif has not been found. Please, provide it.')
                return
        elif self.Param.topo_file_nameEdit1.text()[-1] != '/' and self.Param.topo_file_nameEdit1.text() != 'Nil':
            # Get elevation from another spm or a DEM
            self.topoFile = os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text())
            elevation = np.loadtxt(self.topoFile)
            elevation = elevation.astype(float)
        else:#means flat topography
            elevation = np.array(np.zeros((self.ny,self.nx)))
        # Ensure the dimensions are correct
        try:
            self.Z = np.reshape(elevation,(self.ny,self.nx))
        except ValueError:
            QErrorMessage().showMessage("The dimensions of the grid are wrong. Please, check the values provided.")
            return
        if int(self.Param.Param.DParameters['Raster_flag']) or os.path.exists(os.path.join(self.Param.PFolder,"data",self.Param.topo_file_nameEdit1.text()+'.tif')):
           self.Z = self.Z
        else:
           self.Z = np.flip(self.Z,1)
           
        # Clear the plot area
        if self.updateFault:
            self.plotSpace.plotter.clear()

        ####################################
        ######### Plot the faults ##########
        ####################################
        # Build the grid
        textlegend = ['Latitude (km)','Longitude (km)']
        self.xl = (self.nx-1)*abs(self.dx)*111.11*cos((self.y0+self.dy*self.ny/2)*3.141592654/180)
        self.yl = (self.ny-1)*abs(self.dy)*111.11
        self.x = np.linspace(0,self.xl,self.nx)
        self.y = np.linspace(0,self.yl,self.ny)
        self.X,self.Y = np.meshgrid(self.x,self.y)
        
        # Read the geometry of the fault
        self.nfault = int(self.Param.nfaultEdit1.text())
        self.Lon1 = float(self.Param.lon1Edit2.text())
        self.Lat1 = float(self.Param.lat1Edit3.text())
        self.Lon2 = float(self.Param.lon2Edit4.text())
        self.Lat2 = float(self.Param.lat2Edit5.text())
        # Ensure coordinate of the fault are correct
        self.npoints = int(self.Param.npointiEdit6.text())
        # Read in Fault table and build faults
        num_rows, num_cols = self.Param.faultTable.rowCount(), self.Param.faultTable.columnCount()
        try:
            self.rcoord = np.zeros((self.npoints, self.nfault))
            self.scoord = np.zeros((self.npoints, self.nfault))
            self.FaultsX1Coordinates = np.zeros((self.npoints,self.nfault))
            self.FaultsY1Coordinates = np.zeros((self.npoints,self.nfault))
            self.FaultsX2Coordinates = np.zeros((self.npoints,self.nfault))
            self.FaultsY2Coordinates = np.zeros((self.npoints,self.nfault))
            self.FaultsZ1Coordinates = np.zeros((self.npoints,self.nfault))
            self.FaultsZ2Coordinates = np.zeros((self.npoints,self.nfault))
        except ValueError:
            QErrorMessage(self).showMessage('Value for fault geometry not allowed, please check the parameters value.')
            return
        fault_id = 0
        row_count = 0
        try:
            for row in range(num_rows):
                if row_count == self.npoints:
                    row_count = 0
                    fault_id += 1
                for col in range(num_cols):
                    item = self.Param.faultTable.item(row, col)
                    if item is None:
                        break
                    if col == 1:
                          self.rcoord[row_count,fault_id] = float(item.text())
                    elif col == 2:
                        self.scoord[row_count,fault_id] = float(item.text())
                if item is None:
                    break
                row_count += 1
        except ValueError:
            print('Wrong value for fault parameter. Cannot update the fault geometry.')
            return
        xlon2 = self.x0 + (self.nx - 1) * abs(self.dx)
        xlat2 = y0 + (self.ny - 1) * abs(self.dy)
        x1 = (self.Lon1-self.x0)/(xlon2-self.x0)*self.xl
        y1 = (self.Lat1-y0)/(xlat2-y0)*self.yl
        x2 = (self.Lon2-self.x0)/(xlon2-self.x0)*self.xl
        y2 = (self.Lat2-y0)/(xlat2-y0)*self.yl
           
        try:
            for i in range(self.nfault):
                xn = y2 - y1
                yn = x1 - x2
                xyn = sqrt(xn**2+yn**2)
                xn = xn/xyn
                yn = yn/xyn
                self.FaultsX1Coordinates[:,i] = x1 + xn * self.rcoord[:,i]
                self.FaultsY1Coordinates[:,i] = y1 + yn * self.rcoord[:,i]
                self.FaultsZ1Coordinates[:,i] = self.scoord[:,i]
                self.FaultsX2Coordinates[:,i] = x2 + xn * self.rcoord[:,i]
                self.FaultsY2Coordinates[:,i] = y2 + yn * self.rcoord[:,i]
                self.FaultsZ2Coordinates[:,i] = self.scoord[:,i]
        except ZeroDivisionError:
            QErrorMessage(self).showMessage("Please define the fault coordinates first.")
            return
        
        #### Set the plotting space ####
        #Set the canvas
        self.minZ = np.min(self.Z)/1e3
        self.maxZ = np.max(self.Z)/1e3
        # Plot fault(s)
        for i in range(self.nfault):
            X = np.asarray([self.FaultsX1Coordinates[:,i], self.FaultsX2Coordinates[:,i]])
            Y = np.asarray([self.FaultsY1Coordinates[:,i], self.FaultsY2Coordinates[:,i]])
            Z = np.asarray([self.FaultsZ1Coordinates[:,i], self.FaultsZ2Coordinates[:,i]])
            #set the grid
            xg = np.reshape(X,(X.shape[0]*X.shape[1]))
            yg = np.reshape(Y,(X.shape[0]*X.shape[1]))
            zg = np.reshape(Z,(X.shape[0]*X.shape[1]))
            self.coordArray = np.transpose(np.asarray([xg,yg,zg]))
            grid = pv.StructuredGrid()
            # Set the coordinates from the numpy array
            grid.points = self.coordArray
            # set the dimensions
            grid.dimensions = [X.shape[1],X.shape[0], 1]
            
            #Add Fault to canva
            self.MeshFault = self.plotSpace.plotter.add_mesh(grid,show_edges=False,
                                roughness=1,color=np.ndarray.tolist(np.random.rand(1,3)[0]),opacity=0.7)
            
        # Set the grid
        xg = np.reshape(self.X,(self.X.shape[0]*self.X.shape[1]))
        yg = np.reshape(self.Y,(self.X.shape[0]*self.X.shape[1]))
        zg = np.reshape(self.Z/1e3,(self.X.shape[0]*self.X.shape[1]))
        self.coordArray = np.transpose(np.asarray([xg,yg,zg]))
        minZ = np.min(zg)
        maxZ = np.max(zg)
        
        grid = pv.StructuredGrid()
        # Set the coordinates from the numpy array
        grid.points = self.coordArray
        # set the dimensions
        grid.dimensions = [self.X.shape[1], self.X.shape[0], 1]
        grid.point_data['Elevation'] = self.coordArray[:,2]
        
        #Add mesh to canva
        self.Mesh = self.plotSpace.plotter.add_mesh(grid,show_edges=False,cmap='gist_earth',
                            clim=[minZ,maxZ],roughness=1,scalar_bar_args=self.sargs,
                            ambient=0.3,color=conf.Colors3Dplot['colorbarLabels'])
        self.plotSpace.plotter.show_bounds(location='outer',all_edges=False,color=conf.Colors3Dplot['axes'],
                   xtitle=textlegend[1],ytitle=textlegend[0],ztitle='Elevation (km)')
        # Remove lights
        self.plotSpace.plotter.remove_all_lights(True)
        light = pv.Light()
        # light.set_direction_angle(30,20)
        light.set_camera_light()
        self.plotSpace.plotter.add_light(light)

        self.updateFault = 1
        

#################################################################################
#################################################################################
#################################################################################
class Set_Fault_geometry(QWidget):
    """ 
      This class defines faults geometry
      Param : all parameters from the WindowsParameters class
     
    """
    def __init__(self, Param):
        super().__init__()
        
        # Define parameters
        self.Param = Param
        self.Label1 = QLabel("Set fault parameters:")
        self.Label1.setFont(conf.fontLine12)
        self.Label2 = QLabel("Set fault(s) geometry:")
        self.Label2.setFont(conf.fontLine12)
        self.Label3 = QLabel()
        self.Label3.setText("To set the fault(s) geometry proceed as follows: \n"+
                            "\n 1) Provide the fault segment coordinates of points X1 and X2 \n"+
                            " The velocity field will be applied on the right-hand side of the fault "+
                            "segment from X1 to X2. \n"
                            "\n 2) Define the geometry points with:\n"+
                            "distance: distance from the fault segment (km)\n"+
                            "depth: depth of the point (depth is negative)\n"+
                            "The velocity field is applied to the right-hand side from point 1 to x.\n"+
                            "\n 3) set the amplitude of the velocity field (Main window - 'Uplift history')\n"+
                            "The direction of the velocity vectors follows the order of the geometry points. "+
                            "Providing a negative value makes the direction opposite to that order."
                            )
        
        # Define layouts
        Grid1 = QGridLayout()
        Grid1.addWidget(self.Label1, 0, 0, 1, 3)
        Grid1.addWidget(self.Param.nfault, 1, 1)
        Grid1.addWidget(self.Param.nfaultEdit1, 2, 1)
        Grid1.addWidget(self.Param.npointi, 1, 2)
        Grid1.addWidget(self.Param.npointiEdit6, 2, 2)
        Grid1.addWidget(self.Param.FaultAdvect, 1, 3)
        Grid1.addWidget(self.Param.FaultAdvectb, 2, 3)
        Grid1.addWidget(self.Param.lon1, 5, 1)
        Grid1.addWidget(self.Param.lon1Edit2, 6, 1)
        Grid1.addWidget(self.Param.lat1, 5, 2)
        Grid1.addWidget(self.Param.lat1Edit3, 6, 2)
        Grid1.addWidget(self.Param.lon2, 5, 3)
        Grid1.addWidget(self.Param.lon2Edit4, 6, 3)
        Grid1.addWidget(self.Param.lat2, 5, 4)
        Grid1.addWidget(self.Param.lat2Edit5, 6, 4)
        Grid1.setColumnStretch(0,1)
        Grid1.setColumnStretch(6,3)
        self.Layout = QVBoxLayout()
        self.Layout.addLayout(Grid1)
        Grid2 = QGridLayout()
        Grid2.addWidget(self.Label2,0,0, 1,2)
        Grid2.addWidget(self.Label3,1,0, 1,2)
        Grid2.addWidget(self.Param.faultTable, 2, 0, 1, 2)
        Grid2.setRowStretch(2,2)
        self.Layout.addLayout(Grid2)
        
        # Define signals
        self.Param.lon1Edit2.textChanged.connect(lambda: store_input(
            self.Param, {'lon1': str(self.Param.lon1Edit2.text())}))
        self.Param.lon1Edit2.editingFinished.connect(
            lambda: pgu.isFloat(str(self.Param.lon1Edit2.text())))
        self.Param.lat1Edit3.textChanged.connect(lambda: store_input(
            self.Param, {'lat1': str(self.Param.lat1Edit3.text())}))
        self.Param.lat1Edit3.editingFinished.connect(
            lambda: pgu.isFloat(str(self.Param.lat1Edit3.text())))
        self.Param.lon2Edit4.textChanged.connect(lambda: store_input(
            self.Param, {'lon2': str(self.Param.lon2Edit4.text())}))
        self.Param.lon2Edit4.editingFinished.connect(
            lambda: pgu.isFloat(str(self.Param.lon2Edit4.text())))
        self.Param.lat2Edit5.textChanged.connect(lambda: store_input(
            self.Param, {'lat2': str(self.Param.lat2Edit5.text())}))
        self.Param.lat2Edit5.editingFinished.connect(
            lambda: pgu.isFloat(str(self.Param.lat2Edit5.text())))
        self.Param.faultTable.cellChanged.connect(lambda: self.Param.updateTableFault(
            self.Param.faultTable, int(self.Param.nfaultEdit1.text()), self.Param.npointiEdit6.text(), 1, Param))
 
#################################################################################
#################################################################################
#################################################################################
class eFolding_HP(QWidget):
    """Exponential decrease of heat production rate with depth."""
    
    def __init__(self, parent, Param):
        super().__init__()
        
        self.parent = parent
        self.Param = Param
        # Parameters
        self.HeatProduction_Label = QLabel("Heat production rate (°C/Myr):")
        self.HeatProduction_Label.setFont(conf.font12)
        self.HeatProduction_Label.setAlignment(Qt.AlignCenter)
        self.HeatProduction = QLineEdit(self.parent.HProdEdit7.text())
        self.eFolding_Label = QLabel('e-folding depth (km):')
        self.eFolding_Label.setFont(conf.font12)
        self.eFolding_Label.setAlignment(Qt.AlignCenter)
        
        # Grid
        self.GBox = QGridLayout()
        self.GBox.addWidget(self.HeatProduction_Label,0,0)
        self.GBox.addWidget(self.HeatProduction,1,0)
        self.GBox.addWidget(self.eFolding_Label,0,1)
        self.GBox.addWidget(self.parent.eFolding,1,1)
        
        # Signals
        self.HeatProduction.editingFinished.connect(lambda: self.update_HP())
        self.parent.eFolding.editingFinished.connect(lambda: self.updateGeotherm())
        
    def update_HP(self):
        """ Update value of heat production rate in the interface."""
        self.parent.HProdEdit7.setText(str(self.HeatProduction.text()))
        self.updateGeotherm()
    
    def updateGeotherm(self):
        """ to update geotherm if shown."""
        if self.parent.GeothermTab:
            self.parent.GeothermTab.updateGeotherm()

#################################################################################
#################################################################################
#################################################################################
class WindowParameters(QWidget, object):
    """
      This class defines all the input parameters position in the main window
      with the associated edit lines, or Tables.
      
    """
    def __init__(self, parent, MainWin, Param, IParam, PFolder):
        super(QWidget, self).__init__(parent)

        self.input_parameters = IParam.storeInput  # list to store new input parameters
        self.Param = Param
        self.PFolder = PFolder
        self.MainWin = MainWin
        self.parent = parent
        self.folderName = parent.PNAME
        self.AgesParam = 0
        self.spm = 0
        self.iSOSIA = 0 # Flag if iSOSIA files are loaded
        self.setTopo = 0 # Signal for opening of topo tab
        self.setFault = 0 # Signal for opening of fault tab
        self.GeothermTab = 0 # Signal for Geotherm Tab
        self.HeTab = 0 # Signal for Helium parameters tab
        store_input(self, {'PDModel_signal': '0'})
        # Initialize tab screen
        self.tabs = QTabWidget()
        # self.tabs.setStyleSheet(TabStyleSheet)
        self.tabs.setTabBar(U_interface.HorizontalTabWidget())
        self.tab1 = QWidget(self.tabs)  # Topography
        self.tab2 = QWidget()  # Time evolution
        self.tab3 = QWidget()  # Thermal
        self.tab4 = QWidget()  # Data
        self.tab5 = QWidget()  # Tectonic
        self.tab6 = QWidget()  # Output
        self.tab7 = QWidget()  # Isostasy
        self.tab8 = QWidget()  # Inversion
        self.tab9 = QWidget()  # Output
        # self.tabs.resize(300,200)
        self.tabs.setTabPosition(self.tabs.West)

        # Add tabs
        self.tabs.addTab(self.tab1, "Topography")
        self.tabs.addTab(self.tab2, "Time evolution")
        self.tabs.addTab(self.tab3, "Thermal")
        self.tabs.addTab(self.tab4, "Data")
        self.tabs.addTab(self.tab5, "Tectonic")
        #self.tabs.addTab(self.tab6, "Ages")
        self.tabs.addTab(self.tab7, "Isostasy")
        self.tabs.addTab(self.tab8, "Inversion")
        self.tabs.addTab(self.tab9, "Output")

        ############ Create Topography tab ###########
        # Create text
        self.tab1.layout = QGridLayout(self)
        self.spmFile = QPushButton('load file(s)...')
        self.spmFile.adjustSize()
        self.spmFile.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.spmFile.clicked.connect(self.file_open_spm)
        self.BuildTopoButton = QPushButton('build topo...')
        self.BuildTopoButton.adjustSize()
        self.BuildTopoButton.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.BuildTopoButton.clicked.connect(self.load_BuildTopo)
        self.GetTopoButton = QPushButton('Get topo...')
        self.GetTopoButton.setEnabled(True)
        self.GetTopoButton.adjustSize()
        self.GetTopoButton.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.GetTopoButton.clicked.connect(self.load_from_map)
        self.nx = QLabel('nx:')
        self.nx.setFont(conf.font10)
        self.nx.setAlignment(Qt.AlignCenter)
        self.nx.setToolTip(
            "Number of points describing the topography in the x-(or longitude-) direction.")
        self.ny = QLabel('ny: ')
        self.ny.setFont(conf.font10)
        self.ny.setAlignment(Qt.AlignCenter)
        self.ny.setToolTip(
            "Number of points describing the topography in the y-(or longitude-) direction.")
        self.lon0 = QLabel('Longitude 0: ')
        self.lon0.setFont(conf.font10)
        self.lon0.setAlignment(Qt.AlignCenter)
        self.lon0.setToolTip(
            "Longitude of the bottom left corner of the topographic data/file/grid(s).")
        self.lat0 = QLabel('Latitude 0: ')
        self.lat0.setFont(conf.font10)
        self.lat0.setAlignment(Qt.AlignCenter)
        self.lat0.setToolTip(
            "Latitude of the bottom left corner of the topographic data/file/grid(s).")
        self.deltaLon = QLabel('\u0394lon: ')
        self.deltaLon.setFont(conf.font10)
        self.deltaLon.setAlignment(Qt.AlignCenter)
        self.deltaLon.setToolTip(
            "Distance in decimal degrees between two points on the topographic grid in the longitude direction.")
        self.deltaLat = QLabel('\u0394lat:')
        self.deltaLat.setFont(conf.font10)
        self.deltaLat.setAlignment(Qt.AlignCenter)
        self.deltaLat.setToolTip(
            "Distance in decimal degrees between two points on the topographic grid in the latitude direction.")
        self.nskip = QLabel('nskip: ')
        self.nskip.setFont(conf.font10)
        self.nskip.setAlignment(Qt.AlignCenter)
        self.nskip.setToolTip("Number of topographic point to skip. ! If the 3D topography looks more cubic as you decrease nskip (increase resolution)," +
                              " this is purely artefact of visualization and this does not impact the input resolution of the DEM into Pecube. Double check DEM resolution in tab 'Data->Check sample locations'.")
        self.GridResol = QLabel('Grid resolution: ')
        self.GridResol.setFont(conf.fontLine12)
        self.GridResol.setAlignment(Qt.AlignLeft)
        self.topo_file_name = QLabel('Topography file name:')
        self.topo_file_name.setFont(conf.fontLine12)
        self.topo_file_name.setAlignment(Qt.AlignLeft)
        self.topo_file_name.setToolTip(
            'Name of the file containing the topography as a grid of elevation points.')
        self.ConvertLatLon = QPushButton("Convert to lat/long")
        self.ConvertLatLon.setFont(conf.font8)
        self.ConvertLatLon.adjustSize()
        self.ShowTopo = QPushButton("Show topography")
        self.ShowTopo.setFont(conf.font8)
        self.ShowTopo.adjustSize()
        self.topo_offset = 0
        self.topo_amp = 1
        try:
            self.wavelength = float(Param.DParameters['topo_wavelength'])
            self.topo_offset = float(Param.DParameters['topo_offset'])
            self.topo_amp = float(Param.DParameters['topo_amp'])
            self.topo_phase = float(Param.DParameters['topo_phase'])
        except:
            print('topo wavelength is inverted, take a default value of 10 km.')
            self.wavelength = 10
        
        
        # Create text boxes
        self.topo_file_nameEdit1 = QLineEdit(self)
        self.topo_file_nameEdit1.setText(Param.DParameters['topo_file_name'])
        self.topo_file_nameEdit1.setAlignment(Qt.AlignCenter)
        self.nxEdit2 = QLineEdit(self)
        self.nxEdit2.setValidator(QIntValidator(0,10000))
        self.nxEdit2.setText(Param.DParameters['nx'])
        self.nxEdit2.setAlignment(Qt.AlignCenter)
        self.nyEdit3 = QLineEdit(self)
        self.nyEdit3.setValidator(QIntValidator(0,10000))
        self.nyEdit3.setText(Param.DParameters['ny'])
        self.nyEdit3.setAlignment(Qt.AlignCenter)
        self.lon0Edit4 = QLineEdit(self)
        self.lon0Edit4.setText(Param.DParameters['lon0'])
        self.lon0Edit4.setAlignment(Qt.AlignCenter)
        self.lon0Edit4.setValidator(DoubleValidator)
        self.lat0Edit5 = QLineEdit(self)
        self.lat0Edit5.setText(Param.DParameters['lat0'])
        self.lat0Edit5.setValidator(DoubleValidator)
        self.lat0Edit5.setAlignment(Qt.AlignCenter)
        self.dLonEdit6 = QLineEdit(self)
        self.dLonEdit6.setText(Param.DParameters['dlon'])
        self.dLonEdit6.setAlignment(Qt.AlignCenter)
        self.dLonEdit6.setValidator(DoubleValidator)
        self.dLatEdit7 = QLineEdit(self)
        self.dLatEdit7.setText(Param.DParameters['dlat'])
        self.dLatEdit7.setAlignment(Qt.AlignCenter)
        self.dLatEdit7.setValidator(DoubleValidator)
        self.nskipEdit8 = QLineEdit(self)
        self.nskipEdit8.setText(Param.DParameters['nskip'])
        self.nskipEdit8.setAlignment(Qt.AlignCenter)

        # check input changes
        self.topo_file_nameEdit1.textChanged.connect(lambda: self.update_filename())
        self.nxEdit2.textChanged.connect(lambda: store_input(
            self, {'nx': str(self.nxEdit2.text())}))
        self.nyEdit3.textChanged.connect(lambda: store_input(
            self, {'ny': str(self.nyEdit3.text())}))
        self.lon0Edit4.textChanged.connect(lambda: store_input(
            self, {'lon0': str(self.lon0Edit4.text())}))
        self.lon0Edit4.editingFinished.connect(
            lambda: pgu.isFloat(str(self.lon0Edit4.text())))
        self.lat0Edit5.textChanged.connect(lambda: store_input(
            self, {'lat0': str(self.lat0Edit5.text())}))
        self.lat0Edit5.editingFinished.connect(
            lambda: pgu.isFloat(str(self.lat0Edit5.text())))
        self.dLonEdit6.textChanged.connect(lambda: store_input(
            self, {'dlon': str(self.dLonEdit6.text())}))
        self.dLonEdit6.textChanged.connect(
            lambda: pgu.isFloat(str(self.dLonEdit6.text())))
        self.dLatEdit7.textChanged.connect(lambda: store_input(
            self, {'dlat': str(self.dLatEdit7.text())}))
        self.dLatEdit7.textChanged.connect(
            lambda: pgu.isFloat(str(self.dLatEdit7.text())))
        self.nskipEdit8.textEdited.connect(lambda: store_input(
            self, {'nskip': str(self.nskipEdit8.text())}))
        self.ConvertLatLon.clicked.connect(lambda: self.Convert2LatLon())
        self.ShowTopo.clicked.connect(lambda: self.set_Topography())
        
        # Add text boxes
        self.tab1.layout.setSpacing(10)
        self.tab1.layout.addWidget(self.topo_file_name, 0, 0)
        self.tab1.layout.addWidget(self.topo_file_nameEdit1, 1, 1)
        self.tab1.layout.addWidget(self.spmFile, 1, 3)
        self.tab1.layout.addWidget(self.BuildTopoButton,1,4)
        self.tab1.layout.addWidget(self.GetTopoButton,1,2)

        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        self.tab1.layout.addWidget(Separator, 3, 0, 1, 8)

        # Add text
        self.tab1.layout.addWidget(self.GridResol, 4, 0, 1, 2)
        self.tab1.layout.addWidget(self.nx, 6, 1)
        self.tab1.layout.addWidget(self.ny, 6, 2)
        self.tab1.layout.addWidget(self.lon0, 8, 1)
        self.tab1.layout.addWidget(self.lat0, 8, 2)
        self.tab1.layout.addWidget(self.deltaLon, 8, 3)
        self.tab1.layout.addWidget(self.deltaLat, 8, 4)
        self.tab1.layout.addWidget(self.ConvertLatLon, 9, 5)
        self.tab1.layout.addWidget(self.nskip, 10, 1)
        self.tab1.layout.addWidget(self.nxEdit2, 7, 1)
        self.tab1.layout.addWidget(self.nyEdit3, 7, 2)
        self.tab1.layout.addWidget(self.lon0Edit4, 9, 1)
        self.tab1.layout.addWidget(self.lat0Edit5, 9, 2)
        self.tab1.layout.addWidget(self.dLonEdit6, 9, 3)
        self.tab1.layout.addWidget(self.dLatEdit7, 9, 4)
        self.tab1.layout.addWidget(self.nskipEdit8, 11, 1)
        self.tab1.layout.addWidget(self.ShowTopo, 12, 1)

        self.tab1.layout.setRowStretch(2, 1)
        self.tab1.layout.setRowStretch(13, 2)
        self.tab1.layout.setColumnStretch(7,1)

        self.tab1.setLayout(self.tab1.layout)
        ################# End of Tab 1 ##########################

        ############ Create Time evolution tab ###########
        self.tab2.layout = QVBoxLayout(self)
        self.isSynthTopo = 0
        self.ntime = QLabel('ntime:')
        self.ntime.setFont(conf.font10)
        self.ntime.setAlignment(Qt.AlignCenter)
        self.ntime.setToolTip(
            "Number of time steps needed to describe the scenario for the evolution of the topography.")
        if self.topo_file_nameEdit1.text() ==  conf.Variable_names['Synthetic topo file name']:
            data = {'time topo (Ma)': ['1']*int(self.Param.DParameters['ntime']),
                    "amplification": ['1']*int(self.Param.DParameters['ntime']),
                    "offset (km)": ["0"]*int(self.Param.DParameters['ntime']),
                    "output": ["0"]*int(self.Param.DParameters['ntime']),
                    "phase shift (km)": ["0"]*int(self.Param.DParameters['ntime'])}
            self.isSynthTopo = 1
            self.topoTable = U_interface.WinTable(data, int(Param.DParameters['ntime']), 5)
        else:
            self.isSynthTopo = 0
            data = {'time topo (Ma)': ['1']*int(Param.DParameters['ntime']),
                    "amplification": ['1']*int(Param.DParameters['ntime']),
                    "offset (km)": ["0"]*int(Param.DParameters['ntime']),
                    "output": ["0"]*int(Param.DParameters['ntime'])}
            self.topoTable = U_interface.WinTable(data, int(Param.DParameters['ntime']), 4)
        if int(Param.DParameters['ntime']) > 0:
            self.updateTable(self.topoTable, int(
                Param.DParameters['ntime']), 1, Param)
        self.EroTimeScale = QLabel('Erosional time scale: ')
        self.EroTimeScale.setFont(conf.font10)
        self.EroTimeScale.setAlignment(Qt.AlignCenter)
        self.EroTimeScale.setToolTip(
            "Time scale (Myr) that determines how the topography is interpolated between each time steps.")
        self.TipText = QLabel('Provide the times of topography evolution in the table below ' +
                              'from past to present-day. Example: time 1 = 20 Ma, time 2 = 5 Ma, time 3 = 0 Ma')
        self.TipText.setFont(conf.font11)
        self.TipText.setAlignment(Qt.AlignLeft)
        self.TopoTableName = QLabel('Topography evolution table')
        self.TopoTableName.setFont(conf.fontLine12)
        self.TopoTableName.setAlignment(Qt.AlignLeft)
        self.TopoTableName.setToolTip("time topo: time in Myr in the past" +
                                      "\namplification: topographic amplification factor used at time(i)." +
                                      "\noffset: tvertical offset factor used at time(i)." +
                                      "\noutput: output thermochronological ages for the current time(i) (output = 1).")
        self.TimeParam = QLabel('Time parameters:')
        self.TimeParam.setFont(conf.fontLine12)
        self.TimeParam.setAlignment(Qt.AlignLeft)
        self.RefElevationLabel = QLabel("Reference elevation:")
        self.RefElevationLabel.setFont(conf.font10)
        self.RefElevationLabel.setAlignment(Qt.AlignLeft)
        self.RefElevationLabel.setToolTip("Elevation from which to apply the amplication parameter. (see documentation)")
        self.RefElevCombo = QComboBox()
        self.RefElevCombo.addItems(['seal level','minimum','maximum','mean','custom'])
        self.RefElevCombo.setCurrentIndex(int(Param.DParameters['topo_ref']))
        self.RefCustom = QLineEdit(Param.DParameters[conf.Variable_names['Custom ref topo']])
        self.RefCustom.setEnabled(False)
        try:
            if float(self.RefCustom.text()) > 0:
                self.RefCustom.setEnabled(True)
        except:
            self.RefCustom.setEnabled(True)
            print('Reference elevation cannot be set as a float. Must be iversion format')
        self.SmoothingFactorLabel = QLabel("Radius smoothing (km):")
        self.SmoothingFactorLabel.setFont(conf.font10)
        self.SmoothingFactor = QLineEdit('0')
        self.SmoothingFactor.setEnabled(True)
        self.DoHeadward = QCheckBox("Include headward erosion: ")
        self.DoHeadward.setToolTip("Set an incision that propagate headward")
        if Param.DParameters[conf.Variable_names['Headward erosion flag']] == '1':
            self.DoHeadward.setChecked(True)
        else:
            self.DoHeadward.setChecked(False)
        self.IncisionStart = QLineEdit(self)
        self.IncisionStart.setText(Param.DParameters[conf.Variable_names['Starting incision']])
        self.IncisionStop = QLineEdit(self)
        self.IncisionStop.setText(Param.DParameters[conf.Variable_names['Ending incision']])
        self.IncisionDepth = QLineEdit(self)
        self.IncisionDepth.setText(Param.DParameters[conf.Variable_names['Incision depth']])
        self.Tauh = QLineEdit(self)
        self.Tauh.setText(Param.DParameters[conf.Variable_names['Scaling_propagation_rate']])
        
        # Create text boxes
        self.ntimeEdit1 = QLineEdit(self)
        self.ntimeEdit1.setValidator(QIntValidator(1,10000))
        self.ntimeEdit1.setAlignment(Qt.AlignCenter)
        self.ntimeEdit1.setText(Param.DParameters['ntime'])
        self.EroTimeScaleEdit5 = QLineEdit(self)
        self.EroTimeScaleEdit5.setText(Param.DParameters['erosional_time_scale'])
        self.EroTimeScaleEdit5.setAlignment(Qt.AlignCenter)

        # Add text boxes
        Grid1 = QGridLayout()
        Grid1.setSpacing(5)
        Grid1.addWidget(self.TimeParam, 0, 0, 1, 4)
        Grid1.addWidget(self.ntime, 2, 1, 1, 1)
        Grid1.addWidget(self.ntimeEdit1, 3, 1, 1, 1)
        Grid1.addWidget(self.EroTimeScale, 2, 2, 1, 1)
        Grid1.addWidget(self.EroTimeScaleEdit5, 3, 2, 1, 1)
        Grid1.setColumnStretch(0, 1)
        Grid1.setColumnStretch(1, 1)
        Grid1.setColumnStretch(2, 1)
        Grid1.setColumnStretch(3, 1)
        Grid1.setColumnStretch(4,1)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        self.tab2.layout.addLayout(Grid1)
        self.tab2.layout.addWidget(Separator)
        Grid2 = QGridLayout()
        Grid2.setSpacing(5)
        Grid2.addWidget(self.TopoTableName, 0, 0, 1, 4)
        Grid2.addWidget(self.TipText, 1, 0, 1, 4)
        Hbox = QHBoxLayout()
        Hbox.addWidget(self.RefElevationLabel)
        Hbox.addWidget(self.RefElevCombo)
        Hbox.addWidget(self.RefCustom)
        # Hbox.addWidget(self.DoHeadward)
        Hbox.addStretch(1)
        Hbox2 = QHBoxLayout()
        # Grid2.addWidget(self.SmoothingFactorLabel,3,0,1,1)
        # Grid2.addWidget(self.SmoothingFactor,3,1,1,1)
        Hbox2.addWidget(self.topoTable)
        self.tab2.layout.addLayout(Hbox,1)
        self.tab2.layout.addLayout(Hbox2,1)

        # check input changes
        self.ntimeEdit1.textChanged.connect(lambda: store_input(
            self, {'ntime': str(self.ntimeEdit1.text())}))
        self.ntimeEdit1.editingFinished.connect(lambda: self.updateTable(
            self.topoTable, int(self.ntimeEdit1.text()), 1, Param))
        self.EroTimeScaleEdit5.editingFinished.connect(lambda: store_input(
            self, {conf.Variable_names['Scale erosion time']: str(self.EroTimeScaleEdit5.text())}))
        self.EroTimeScaleEdit5.editingFinished.connect(
            lambda: pgu.isFloat(str(self.EroTimeScaleEdit5.text())))
        self.topoTable.cellChanged.connect(lambda: self.updateTable(
            self.topoTable, int(self.ntimeEdit1.text()), 0, Param))
        self.RefElevCombo.currentIndexChanged.connect(lambda: self.updateRefElev())
        self.DoHeadward.stateChanged.connect(lambda: self.set_Headward_Erosion())
        self.IncisionStart.editingFinished.connect(lambda:store_input(
            self, {conf.Variable_names['Starting incision']: str(self.IncisionStart.text())}))
        self.IncisionStop.editingFinished.connect(lambda:store_input(
            self, {conf.Variable_names['Ending incision']: str(self.IncisionStop.text())}))
        self.RefCustom.editingFinished.connect(lambda:store_input(
            self, {conf.Variable_names['Custom ref topo']: str(self.RefCustom.text())}))
        self.IncisionDepth.editingFinished.connect(lambda:store_input(
            self, {conf.Variable_names['Incision depth']: str(self.IncisionDepth.text())}))
        self.Tauh.editingFinished.connect(lambda:store_input(
            self, {conf.Variable_names['Scaling_propagation_rate']: str(self.Tauh.text())}))
        # self.SmoothingFactor.editingFinished.connect(lambda: self.updateSmoothingFactor())

        self.tab2.setLayout(self.tab2.layout)
        ################# End of Tab 2 ##########################

        ############ Create Thermal tab ###########
        self.tab3.layout = QGridLayout(self)
        self.tab3.layout.setSpacing(10)
        self.Heat_production_flag = int(Param.DParameters[conf.Variable_names['Heat_production_flag']])
        self.thickness = QLabel('Thickness (km):')
        self.thickness.setFont(conf.font10)
        self.thickness.setAlignment(Qt.AlignCenter)
        self.thickness.setToolTip("Crustal thickness.")
        self.bTemperature = QLabel('Basal Temperature (°C):')
        self.bTemperature.setFont(conf.font10)
        self.bTemperature.setAlignment(Qt.AlignCenter)
        self.bTemperature.setToolTip(
            "Temperature imposed at the base of the model.")
        self.nz = QLabel('nz: ')
        self.nz.setFont(conf.font10)
        self.nz.setAlignment(Qt.AlignCenter)
        self.nz.setToolTip(
            "Number of points used to discretize the crust in the z-direction.")
        self.slt = QLabel('Sea level Temperature (°C): ')
        self.slt.setFont(conf.font10)
        self.slt.setAlignment(Qt.AlignCenter)
        self.slt.setToolTip("Temperature at sea level.")
        self.lapseRate = QLabel('Lapse Rate (°C/km): ')
        self.lapseRate.setFont(conf.font10)
        self.lapseRate.setAlignment(Qt.AlignCenter)
        self.lapseRate.setToolTip(
            "Rate of change of temperature with elevation in the atmosphere.")
        self.TDiffusivity = QLabel('Thermal diffusivity (km²/Myr): ')
        self.TDiffusivity.setFont(conf.font10)
        self.TDiffusivity.setAlignment(Qt.AlignCenter)
        self.TDiffusivity.setToolTip("Thermal diffusivity (\u039A).")
        self.HProd = QLabel('Heat production (°C/Myr): ')
        self.HProd.setFont(conf.font10)
        self.HProd.setAlignment(Qt.AlignCenter)
        self.HProd.setToolTip(
            "Heat production used to solve the heat equation.")
        self.EfoldingHeatProd = QCheckBox("use e-folding HP")
        self.EfoldingHeatProd.setToolTip("Use an exponential decrease of heat production rate.")
        if self.Heat_production_flag:
            self.EfoldingHeatProd.setChecked(True)
        else:
            self.EfoldingHeatProd.setChecked(False)  
        self.Crust = QLabel("Crust:")
        self.Crust.setFont(conf.fontLine12)
        self.Crust.setAlignment(Qt.AlignLeft)
        self.Atmosphere = QLabel("Atmosphere:")
        self.Atmosphere.setFont(conf.fontLine12)
        self.Atmosphere.setAlignment(Qt.AlignLeft)
        self.VerticalResol = QLabel("Vertical resolution:")
        self.VerticalResol.setFont(conf.fontLine12)
        self.VerticalResol.setAlignment(Qt.AlignLeft)
        self.Geotherm = QPushButton('Show geotherm')
        self.Geotherm.adjustSize()
        self.Geotherm.setFont(conf.font8)
        self.shearHeating = QLabel('Shear heating:')
        self.shearHeating.setFont(conf.font10)
        self.shearHeating.setAlignment(Qt.AlignCenter)
        self.shearHeating.setToolTip(
            "Friction coefficient used to multiply an arbitrary stress value of 100 MPa to comute the heat produced by friction.")
        self.shearHeatingEdit = QLineEdit(self)
        self.shearHeatingEdit.setText(Param.DParameters['shear_heating'])
        self.shearHeatingEdit.setAlignment(Qt.AlignCenter)
        self.minimum_dt_label = QLabel('Minimum dt (Myr):')
        self.minimum_dt_label.setFont(conf.font10)
        self.minimum_dt_label.setAlignment(Qt.AlignCenter)
        self.minimum_dt_label.setToolTip(
            "Minimum time step to use by Pecube. Useful if internal dt = dt(time_topo_i-1 - time_topo_i).")

        # Create text boxes
        self.thicknessEdit1 = QLineEdit(self)
        self.thicknessEdit1.setText(Param.DParameters['thickness'])
        self.thicknessEdit1.setAlignment(Qt.AlignCenter)
        self.bTemperatureEdit2 = QLineEdit(self)
        self.bTemperatureEdit2.setText(Param.DParameters['basal_temperature'])
        self.bTemperatureEdit2.setAlignment(Qt.AlignCenter)
        self.nzEdit3 = QLineEdit(self)
        self.nzEdit3.setText(Param.DParameters['nz'])
        self.nzEdit3.setAlignment(Qt.AlignCenter)
        self.nzEdit3.setValidator(QIntValidator(0,1000))
        self.sltEdit4 = QLineEdit(self)
        self.sltEdit4.setText(Param.DParameters['sea_level_temperature'])
        self.sltEdit4.setAlignment(Qt.AlignCenter)
        self.lapseRateEdit5 = QLineEdit(self)
        self.lapseRateEdit5.setText(Param.DParameters['lapse_rate'])
        self.lapseRateEdit5.setAlignment(Qt.AlignCenter)
        self.TDiffusivityEdit6 = QLineEdit(self)
        self.TDiffusivityEdit6.setText(Param.DParameters['thermal_diffusivity'])
        self.TDiffusivityEdit6.setAlignment(Qt.AlignCenter)
        self.HProdEdit7 = QLineEdit(self)
        self.HProdEdit7.setText(str(Param.DParameters['heat_production']))
        self.HProdEdit7.setAlignment(Qt.AlignCenter)
        self.eFolding = QLineEdit(Param.DParameters[conf.Variable_names['eFoldingDepth']])
        self.minimum_dt = QLineEdit(Param.DParameters[conf.Variable_names['minimum time step']])
        self.minimum_dt.setAlignment(Qt.AlignCenter)

        # Add text
        self.tab3.layout.addWidget(self.Crust, 0, 0, 1, 3)
        self.tab3.layout.addWidget(self.thickness, 2, 1)
        self.tab3.layout.addWidget(self.thicknessEdit1, 3, 1)
        self.tab3.layout.addWidget(self.bTemperature, 2, 2)
        self.tab3.layout.addWidget(self.bTemperatureEdit2, 3, 2)
        self.tab3.layout.addWidget(self.Geotherm, 3, 3)
        self.tab3.layout.addWidget(self.TDiffusivity, 4, 1)
        self.tab3.layout.addWidget(self.TDiffusivityEdit6, 5, 1)
        self.tab3.layout.addWidget(self.HProd, 4, 2)
        self.tab3.layout.addWidget(self.HProdEdit7, 5, 2)
        self.tab3.layout.addWidget(self.EfoldingHeatProd, 5, 3)
        self.tab3.layout.addWidget(self.shearHeating, 6, 1)
        self.tab3.layout.addWidget(self.shearHeatingEdit, 7, 1)
        # self.tab3.layout.addWidget(self.minimum_dt_label, 6, 2)
        # self.tab3.layout.addWidget(self.minimum_dt, 7, 2)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        self.tab3.layout.addWidget(Separator, 6, 0, 1, 4)
        self.tab3.layout.addWidget(self.VerticalResol, 8, 0, 1, 3)
        self.tab3.layout.addWidget(self.nz, 10, 1)
        self.tab3.layout.addWidget(self.nzEdit3, 11, 1)
        self.tab3.layout.addWidget(Separator, 12, 0, 1, 4)
        self.tab3.layout.addWidget(self.Atmosphere, 13, 0, 1, 3)
        self.tab3.layout.addWidget(self.slt, 15, 1)
        self.tab3.layout.addWidget(self.lapseRate, 15, 2)
        self.tab3.layout.addWidget(self.sltEdit4, 16, 1)
        self.tab3.layout.addWidget(self.lapseRateEdit5, 16, 2)

        # Size of boxes
        self.tab3.layout.setColumnStretch(0, 1)
        self.tab3.layout.setColumnStretch(1, 1)
        self.tab3.layout.setColumnStretch(2, 1)
        self.tab3.layout.setColumnStretch(4, 3)
        self.tab3.layout.setRowStretch(17,1)

        #### check updates ####
        self.thicknessEdit1.textEdited.connect(lambda: store_input(
            self, {'thickness': str(self.thicknessEdit1.text())}))
        self.bTemperatureEdit2.textEdited.connect(lambda: store_input(
            self, {'basal_temperature': str(self.bTemperatureEdit2.text())}))
        self.nzEdit3.textEdited.connect(lambda: store_input(
            self, {'nz': str(self.nzEdit3.text())}))
        self.sltEdit4.textEdited.connect(lambda: store_input(
            self, {'sea_level_temperature': str(self.sltEdit4.text())}))
        self.lapseRateEdit5.textEdited.connect(lambda: store_input(
            self, {'lapse_rate': str(self.lapseRateEdit5.text())}))
        self.TDiffusivityEdit6.textEdited.connect(lambda: store_input(
            self, {'thermal_diffusivity': str(self.TDiffusivityEdit6.text())}))
        self.HProdEdit7.textChanged.connect(lambda: self.store_HeatProd())
        self.shearHeatingEdit.textEdited.connect(lambda: store_input(
            self, {'shear_heating': str(self.shearHeatingEdit.text())}))
        #Geotherm box is updated
        self.Geotherm.clicked.connect(lambda: self.set_Geotherm())
        self.EfoldingHeatProd.stateChanged.connect(lambda: self.get_EFolding())
        self.eFolding.editingFinished.connect(lambda: store_input(
            self, {conf.Variable_names['eFoldingDepth']: str(self.eFolding.text())}))
        self.minimum_dt.editingFinished.connect(lambda:store_input(
            self, {conf.Variable_names['minimum time step']: str(self.minimum_dt.text())}) )

        self.tab3.setLayout(self.tab3.layout)
        ################# End of Tab 3 ##########################

        # ############ Create Data tab ###########
        self.tab4.layout = QGridLayout(self)
        self.tab4.layout.setSpacing(10)
        # Select computation style (for all nodes or sample specific)
        self.ComputeAge = QLabel("Compute ages:")
        self.ComputeAge.setFont(conf.fontLine12)
        self.ComputeAge.setAlignment(Qt.AlignLeft)
        self.AgeCombo = QComboBox(self)
        self.AgeCombo.setGeometry(0,0,50,50)
        ItemList = ['none','for all nodes','sample specific']
        self.AgeCombo.addItems(ItemList)
        self.AgeCombo.setCurrentIndex(
            int(Param.DParameters[conf.Variable_names['Mode_age_computation']]))
        self.AgeCombo.resize(self.AgeCombo.sizeHint())
        self.ShowAgesButton = QPushButton("Show/update ages tab")
        self.ShowAgesButton.setFont(conf.font8)
        # self.ShowAgesButton.adjustSize()
        self.dataFolder = QLabel('Data folder name:')
        self.dataFolder.setFont(conf.font10)
        self.dataFolder.setToolTip(
            "Name of a folder stored in the 'data/' directory.")
        # Create text boxes
        self.dataFolderEdit1 = QLineEdit(self)
        self.dataFolderEdit1.setText(Param.DParameters['data_folder'])
        self.ModelButton = QPushButton('Check sample locations')
        self.ModelButton.setFont(conf.font8)
        self.CheckElevButton = QPushButton('Check sample elevations')
        self.CheckElevButton.setFont(conf.font8)
        self.groupBoxData = QGroupBox("Provide sample location(s)")
        self.groupBoxData.setCheckable(False)
        self.Samples = {}
        self.Samples['ngrains1'] = Param.DParameters['Nb_Samples']
        self.Samples['SampleName1'] = Param.DParameters['SampleName'] 
        self.nbSample = QLabel('Number of samples:')
        self.nbSample.setFont(conf.font8)
        self.nbSample.setAlignment(Qt.AlignRight)
        self.nbSampleValue = QLineEdit()
        self.nbSampleValue.setText(Param.DParameters['Nb_Samples'])
        self.nbSampleValue.setValidator(QIntValidator())
        data = {'Sample': ['1']*int(Param.DParameters['ngrains']),
                "Longitude": ['1']*int(Param.DParameters['ngrains']),
                "Latitude": ["0"]*int(Param.DParameters['ngrains']),
                "Elevation (m)": ["0"]*int(Param.DParameters['ngrains']),
                "Th. L.": ["0"]*int(Param.DParameters['ngrains']),
                "OSL": ["0"]*int(Param.DParameters['ngrains']),
                "ESR": ["0"]*int(Param.DParameters['ngrains']),
                "AHe": ["0"]*int(Param.DParameters['ngrains']),
                "ZHe": ["0"]*int(Param.DParameters['ngrains']),
                "AFT": ["0"]*int(Param.DParameters['ngrains']),
                "ZFT": ["0"]*int(Param.DParameters['ngrains']),
                "KAr": ["0"]*int(Param.DParameters['ngrains']),
                "BAr": ["0"]*int(Param.DParameters['ngrains']),
                "MAr": ["0"]*int(Param.DParameters['ngrains']),
                "HAr": ["0"]*int(Param.DParameters['ngrains']),}
        self.SampleTable = U_interface.WinTable(data, int(Param.DParameters['ngrains']), 15)
        
        # Widgets for "for all nodes" computation style
        self.Label0 = QLabel('Thermochronological systems: ')
        self.Label0.setFont(conf.fontLine12)
        self.Label0.setAlignment(Qt.AlignLeft)
        # if Param.DParameters['oldInput'] == '1':

        self.Label1 = QLabel('Trapped-charge: ')
        self.Label1.setFont(conf.fontBold10)
        self.Label1.setAlignment(Qt.AlignRight)

        self.Label2 = QLabel('(U-Th)/He: ')
        self.Label2.setFont(conf.fontBold10)
        self.Label2.setAlignment(Qt.AlignRight)

        self.Label3 = QLabel('Fission tracks: ')
        self.Label3.setFont(conf.fontBold10)
        self.Label3.setAlignment(Qt.AlignRight)

        self.Label4 = QLabel('Argon: ')
        self.Label4.setFont(conf.fontBold10)
        self.Label4.setAlignment(Qt.AlignRight)

        self.ageESRb = QCheckBox()
        self.ageESRb.setEnabled(False)
        if self.Param.DParameters['age_ESR_flag'] == '1' or self.Param.DParameters['age_ESR_flag'] == 1:
            self.ageESRb.setEnabled(False)
            self.ageESRb.setChecked(False)
        self.ageESR = QLabel('ESR (disabled): ')
        self.ageESR.setFont(conf.font10)
        self.ageESR.setAlignment(Qt.AlignRight)
        self.ageESR.setToolTip("Compute ESR ages. (disabled for now, test)")
		
        self.ageTLb = QCheckBox()
        self.ageTLb.setEnabled(False)
        if self.Param.DParameters['age_TL_flag'] == '1' or self.Param.DParameters['age_TL_flag'] == 1:
            self.ageTLb.setEnabled(False)
            self.ageTLb.setChecked(False)
        self.ageTL = QLabel('Th.L (disabled): ')
        self.ageTL.setFont(conf.font10)
        self.ageTL.setAlignment(Qt.AlignRight)
        self.ageTL.setToolTip("Thermoluminescence (TL, disabled for now, test).")
		
        self.ageOSLb = QCheckBox()
        self.ageOSLb.setEnabled(False)
        if self.Param.DParameters['age_OSL_flag'] == '1' or self.Param.DParameters['age_OSL_flag'] == 1:
            self.ageOSLb.setEnabled(False)
            self.ageOSLb.setChecked(False)
        self.ageOSL = QLabel('OSL (disabled): ')
        self.ageOSL.setFont(conf.font10)
        self.ageOSL.setAlignment(Qt.AlignRight)
        self.ageOSL.setToolTip("Optical lumiscence (OSL, disabled for now, test).")
		
        self.ageAHeb = QCheckBox()
        self.ageAHeb.setEnabled(False)
        if self.Param.DParameters['age_AHe_flag'] == '1' or self.Param.DParameters['age_AHe_flag'] == 1:
            self.ageAHeb.setEnabled(True)
            self.ageAHeb.setChecked(True)

        self.ageAHe = QLabel('AHe: ')
        self.ageAHe.setFont(conf.font10)
        self.ageAHe.setAlignment(Qt.AlignRight)
        self.ageAHe.setToolTip("Compute helium ages in Apatite.")
        self.ageZHeb = QCheckBox()
        self.ageZHeb.setEnabled(False)
        if self.Param.DParameters['age_ZHe_flag'] == '1' or self.Param.DParameters['age_ZHe_flag'] == 1:
            self.ageZHeb.setEnabled(True)
            self.ageZHeb.setChecked(True)

        self.ageZHe = QLabel('ZHe: ')
        self.ageZHe.setFont(conf.font10)
        self.ageZHe.setAlignment(Qt.AlignRight)
        self.ageZHe.setToolTip("Compute helium ages in Zircon.")
        self.ageAFTb = QCheckBox()
        self.ageAFTb.setEnabled(False)
        if self.Param.DParameters['age_AFT_flag'] == '1' or self.Param.DParameters['age_AFT_flag'] == 1:
            self.ageAFTb.setEnabled(True)
            self.ageAFTb.setChecked(True)

        self.ageAFT = QLabel('AFT: ')
        self.ageAFT.setFont(conf.font10)
        self.ageAFT.setAlignment(Qt.AlignRight)
        self.ageAFT.setToolTip("Compute fission track ages in Apatite.")
        self.ageZFTb = QCheckBox()
        self.ageZFTb.setEnabled(False)
        if self.Param.DParameters['age_ZFT_flag'] == '1' or self.Param.DParameters['age_ZFT_flag'] == 1:
            self.ageZFTb.setEnabled(True)
            self.ageZFTb.setChecked(True)

        self.ageZFT = QLabel('ZFT: ')
        self.ageZFT.setFont(conf.font10)
        self.ageZFT.setAlignment(Qt.AlignRight)
        self.ageZFT.setToolTip("Compute fission track ages in Zircon.")

        self.ageKAr = QLabel('KAr: ')
        self.ageKAr.setFont(conf.font10)
        self.ageKAr.setAlignment(Qt.AlignCenter)
        self.ageKArb = QCheckBox()
        self.ageKArb.setEnabled(False)
        if self.Param.DParameters['age_KAr_flag'] == '1' or self.Param.DParameters['age_KAr_flag'] == 1:
            self.ageKArb.setEnabled(True)
            self.ageKArb.setChecked(True)
        self.ageKArb.setToolTip("Compute Argon ages on K-Felsdpar.")

        self.ageBArb = QCheckBox()
        self.ageBArb.setEnabled(False)
        if self.Param.DParameters['age_BAr_flag'] == '1' or self.Param.DParameters['age_BAr_flag'] == 1:
            self.ageBArb.setEnabled(True)
            self.ageBArb.setChecked(True)
        self.ageBAr = QLabel('BAr')
        self.ageBAr.setFont(conf.font10)
        self.ageBAr.setAlignment(Qt.AlignRight)
        self.ageBAr.setToolTip("Compute Argon ages in Biotite.")

        self.ageMArb = QCheckBox()
        self.ageMArb.setEnabled(False)
        if self.Param.DParameters['age_MAr_flag'] == '1' or self.Param.DParameters['age_MAr_flag'] == 1:
            self.ageMArb.setEnabled(True)
            self.ageMArb.setChecked(True)
        self.ageMAr = QLabel('MAr: ')
        self.ageMAr.setFont(conf.font10)
        self.ageMAr.setAlignment(Qt.AlignRight)
        self.ageMAr.setToolTip("Compute Argon ages in Muscovite.")

        self.ageHArb = QCheckBox()
        self.ageHArb.setEnabled(False)
        if self.Param.DParameters['age_HAr_flag'] == '1' or self.Param.DParameters['age_HAr_flag'] == 1:
            self.ageHArb.setEnabled(True)
            self.ageHArb.setChecked(True)
        self.ageHAr = QLabel('HAr: ')
        self.ageHAr.setFont(conf.font10)
        self.ageHAr.setAlignment(Qt.AlignRight)
        self.ageHAr.setToolTip("Compute Argon ages in Hornblende.")

        self.defaultAge = QLabel('Default age (Ma):')
        self.defaultAge.setFont(conf.font10)
        self.defaultAge.setAlignment(Qt.AlignRight)
        self.defaultAge.setToolTip("Default age given to rock particles that have not" +
                                    " been reset (i.e. not reached the the 'closure temperature'" +
                                    "for the system considered).")
        self.defaultAgeEdit2 = QLineEdit(self)
        self.defaultAgeEdit2.setText(Param.DParameters[conf.Variable_names['Default_age']])

        # Add text
        VBox = QVBoxLayout()
        # VBox.addWidget(self.Label0,0,0)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        VBox.addWidget(Separator)
        Grid1 = QGridLayout()
        Grid1.setSpacing(10)
        Grid1.addWidget(self.Label0, 0, 0, 1, 3)
        Grid1.addWidget(self.Label1, 3, 1)
        Grid1.addWidget(self.Label2, 4, 1)
        Grid1.addWidget(self.Label3, 5, 1)
        Grid1.addWidget(self.ageTLb, 3, 3)
        Grid1.addWidget(self.ageTL, 3, 2)
        Grid1.addWidget(self.ageOSL, 3, 4)
        Grid1.addWidget(self.ageOSLb, 3, 5)
        Grid1.addWidget(self.ageESR, 3, 6)
        Grid1.addWidget(self.ageESRb, 3, 7)
        Grid1.addWidget(self.ageAHe, 4, 2)
        Grid1.addWidget(self.ageAHeb, 4, 3)
        Grid1.addWidget(self.ageZHe, 4, 4)
        Grid1.addWidget(self.ageZHeb, 4, 5)
        Grid1.addWidget(self.ageAFT, 5, 2)
        Grid1.addWidget(self.ageAFTb, 5, 3)
        Grid1.addWidget(self.ageZFT, 5, 4)
        Grid1.addWidget(self.ageZFTb, 5, 5)
        Grid1.setColumnStretch(8, 1)
        Grid1.addWidget(self.Label4, 7, 1)
        Grid1.addWidget(self.ageKAr, 7, 2)
        Grid1.addWidget(self.ageKArb, 7, 3)
        Grid1.addWidget(self.ageBAr, 7, 4)
        Grid1.addWidget(self.ageBArb, 7, 5)
        Grid1.addWidget(self.ageMAr, 8, 2)
        Grid1.addWidget(self.ageMArb, 8, 3)
        Grid1.addWidget(self.ageHAr, 8, 4)
        Grid1.addWidget(self.ageHArb, 8, 5)
        Grid1.addWidget(self.defaultAge,9,1)
        Grid1.addWidget(self.defaultAgeEdit2,9,2,1,2)
        
        VBox.addLayout(Grid1)
        VBox.addStretch(4)
        
        #### SET LAYOUTS ###
        # grid for Computation style combo box
        Hbox = QHBoxLayout()
        # Hbox.addWidget(self.ComputeAge)
        Hbox.addWidget(self.AgeCombo)
        Hbox.addWidget(self.ShowAgesButton)
        Hbox.addStretch(4)
        # Widget 1 for "all nodes" computation
        self.widget1 = QWidget()
        self.widget1.setLayout(VBox)
        # Widget 2 for "sample specific" computation
        self.widget2 = QWidget()
        HBox2 = QHBoxLayout()
        HBox2.addWidget(self.nbSample)
        HBox2.addWidget(self.nbSampleValue)
        HBox2.addWidget(self.ModelButton)
        HBox2.addWidget(self.CheckElevButton)
        HBox2.addStretch(1)
        HBox3 = QHBoxLayout()
        HBox3.addWidget(self.dataFolder)
        HBox3.addWidget(self.dataFolderEdit1)
        HBox3.addStretch(1)
        HBox4 = QHBoxLayout()
        HBox4.addWidget(self.defaultAge)
        HBox4.addWidget(self.defaultAgeEdit2)
        HBox4.addStretch(1)
        GBLayout = QGridLayout()
        GBLayout.addLayout(HBox3,0,0)
        GBLayout.addLayout(HBox4,1,0)
        GBLayout.addLayout(HBox2,2,0)
        GBLayout.addWidget(QWidget(),2,3,1,2)
        GBLayout.addWidget(self.SampleTable,3,0,3,3)
        self.groupBoxData.setLayout(GBLayout)
        self.widget2.setLayout(GBLayout)

        # Add text
        self.tab4.layout.addLayout(Hbox, 0,0)
        self.tab4.layout.addWidget(self.widget1,1,0,3,3)
        self.tab4.layout.addWidget(self.widget2,2,0,3,3)
        
        # Size of boxes
        self.tab4.layout.setRowStretch(3, 1)
        
        # Update table sample
        self.updateTableSample(int(self.nbSampleValue.text()), 1, Param,0)
        if Param.DParameters['oldInput'] == '1':
            try:
                for i in range(int(Param.DParameters['Nb_Samples'])):
                    self.Samples['ngrains'+str(i+1)] = Param.TableParameters['nb_grains'][i][0]
                    self.Samples['SampleName'+str(i+1)]  = Param.TableParameters['SampleName'+str(i+1)]
            except KeyError:
                self.Samples['SampleName'+str(i+1)] = Param.DParameters['SampleName']
        
        # Initialize ages
        self.AgesParam = th.set_Thermochron_Parameters(self, self.Param, self.PFolder)
        
        # check input changes
        self.ShowAgesButton.clicked.connect(lambda:self.open_ages_tab())
        self.dataFolderEdit1.textEdited.connect(lambda: store_input(
            self, {'data_folder': str(self.dataFolderEdit1.text())}))
        self.dataFolderEdit1.editingFinished.connect(
            lambda: pgu.isString(self.dataFolderEdit1.text()))
        self.nbSampleValue.editingFinished.connect(lambda: self.updateTableSample(
            int(self.nbSampleValue.text()), 1, Param, 0))
        self.SampleTable.cellChanged.connect(lambda: self.updateTableSample(
            int(self.nbSampleValue.text()), 0, Param, 1))
        self.nbSampleValue.editingFinished.connect(lambda: store_input(
            self, {'Nb_Samples': str(self.nbSampleValue.text())}))
        self.ModelButton.clicked.connect(lambda: self.showModel(self.ModelButton))
        self.CheckElevButton.clicked.connect(lambda: self.showModel(self.CheckElevButton))

        self.tab4.setLayout(self.tab4.layout)
        
        ################# End of Tab 4 ##########################

        ############ Create Tectonic tab ###########
        self.tab5.layout = QVBoxLayout(self)
        self.faultParam = 0
        self.UpliftButton1 = QRadioButton("no uplift")
        self.UpliftButton2 = QRadioButton("bloc uplift")
        self.UpliftButton3 = QRadioButton("faulting")
        self.nfaultEdit1 = QLineEdit(self)
        self.nfaultEdit1.setValidator(QIntValidator(0,100))
        self.nfaultEdit1.setText(Param.DParameters['nfault'])
        self.nfaultEdit1.setAlignment(Qt.AlignCenter)
        self.nfault = QLabel('Number of faults:')
        self.nfault.setFont(conf.font10)
        self.nfault.setAlignment(Qt.AlignCenter)
        self.nfault.setToolTip(
            "Number of faults used to describe the velocity field.")
        self.lon1 = QLabel('Longitude X1:')
        self.lon1.setFont(conf.font10)
        self.lon1.setAlignment(Qt.AlignCenter)
        self.lon1.setToolTip("Longitude of X1, the first of two points used to" +
                              " define the fault strike and width")
        self.lat1 = QLabel('Latitude X1: ')
        self.lat1.setFont(conf.font10)
        self.lat1.setAlignment(Qt.AlignCenter)
        self.lat1.setToolTip("Latitude of X1, the first of two points used to" +
                              " define the fault strike and width")
        self.lon2 = QLabel('Longitude X2: ')
        self.lon2.setFont(conf.font10)
        self.lon2.setAlignment(Qt.AlignCenter)
        self.lon2.setToolTip("Longitude of X2, the second of two points used to" +
                              " define the fault strike and width")
        self.lat2 = QLabel('Latitude X2: ')
        self.lat2.setFont(conf.font10)
        self.lat2.setAlignment(Qt.AlignCenter)
        self.lat2.setToolTip("Latitude of X2, the second of two points used to" +
                              " define the fault strike and width")
        self.faultTable = U_interface.createTable(int(Param.DParameters['nfault']), 3)
        HLabels = ["fault(s)", "distance (km)", "depth (km)"]
        self.faultTable.setHorizontalHeaderLabels(HLabels)
        if int(Param.DParameters['nfault']) > 1 or int(Param.DParameters['uplift']) == 2:
            self.updateTableFault(self.faultTable, int(Param.DParameters['nfault']), int(
                Param.TableParameters['npoint1']), 2, Param)
        self.nstepKinTable = U_interface.createTable(
            int(Param.DParameters['nstep1']), 4)
        HLabels = ["ID", "time start (Ma)", "time end (Ma)", "velo (km/Myr)"]
        self.nstepKinTable.setHorizontalHeaderLabels(HLabels)
        self.ThrustCombo = QComboBox()
        self.ThrustCombo.addItem("Thrust")
        self.ThrustCombo.addItem("Normal")
        self.npointi = QLabel('Number of points: ')
        self.npointi.setFont(conf.font10)
        self.npointi.setAlignment(Qt.AlignCenter)
        self.npointi.setToolTip(
            "Number of points used to define the fault segments in the vertical direction for fault i.")
        self.bottomLeft = QLabel('Bottom left: ')
        self.bottomLeft.setFont(conf.font10)
        self.bottomLeft.setAlignment(Qt.AlignCenter)
        self.bottomLeft.setToolTip("Amplification factor for uplift velocity imposed across the entire depth range at the bottom left corner" +
                                    "of the Pecube domain.")
        self.bottomRight = QLabel('Bottom right: ')
        self.bottomRight.setFont(conf.font10)
        self.bottomRight.setAlignment(Qt.AlignCenter)
        self.bottomRight.setToolTip("Amplification factor for uplift velocity imposed across the entire depth range at the bottom right corner" +
                                    "of the Pecube domain.")
        self.upperRight = QLabel('Upper right: ')
        self.upperRight.setFont(conf.font10)
        self.upperRight.setAlignment(Qt.AlignCenter)
        self.upperRight.setToolTip("Amplification factor for uplift velocity imposed across the entire depth range at the upper right corner" +
                                   "of the Pecube domain.")
        self.upperLeft = QLabel('Upper left: ')
        self.upperLeft.setFont(conf.font10)
        self.upperLeft.setAlignment(Qt.AlignCenter)
        self.upperLeft.setToolTip("Amplification factor for uplift velocity imposed across the entire depth range at the upper left corner" +
                                  "of the Pecube domain.")
        self.nstepi = QLabel('nstep: ')
        self.nstepi.setFont(conf.font10)
        self.nstepi.setAlignment(Qt.AlignCenter)
        self.nstepi.setToolTip(
            "Number of time intervals needed to define the velocity history for fault/bloc i.")
        self.logVelo = QLabel('Logarithmic velocity: ')
        self.logVelo.setFont(conf.font10)
        self.logVelo.setToolTip(
            "Checking this box gives velocity in logarithmic values.")
        self.logVelo.setAlignment(Qt.AlignCenter)
        self.logVelob = QCheckBox()
        self.logVelob.setEnabled(False)
        self.UpliftLabel = QLabel("Define uplift:")
        self.UpliftLabel.setFont(conf.fontLine12)
        self.UpliftLabel.setAlignment(Qt.AlignLeft)
        self.FaultGeometry = QLabel("Fault geometry:")
        self.FaultGeometry.setFont(conf.fontLine12)
        self.FaultGeometry.setAlignment(Qt.AlignLeft)
        self.FaultKinematics = QLabel("Model kinematics:")
        self.FaultKinematics.setFont(conf.fontLine12)
        self.FaultKinematics.setAlignment(Qt.AlignLeft)
        self.UpliftHistory = QLabel("Uplift history:")
        self.UpliftHistory.setFont(conf.fontLine12)
        self.UpliftHistory.setAlignment(Qt.AlignLeft)
        self.UpliftHistory.setToolTip("Set the time evolution of uplift rates.")
        self.ShowFault = QPushButton("show faults")
        self.ShowFault.adjustSize()
        self.ShowFault.setFont(conf.font8)
        self.ShowFault.setEnabled(False)
        self.setFaultButton = QPushButton("set faults")
        self.setFaultButton.adjustSize()
        self.setFaultButton.setFont(conf.font8)
        self.setFaultButton.setEnabled(False)
        
        #---------------------------------
        #See "Set_Fault_geometry" 
        self.FaultAdvectb = QCheckBox()
        self.FaultAdvectb.setText('yes/no')
        self.FaultAdvect = QLabel('Enable fault advection:')
        self.FaultAdvect.setFont(conf.font10)
        self.FaultAdvect.setAlignment(Qt.AlignRight)
        self.FaultAdvect.setToolTip(
            "Checking this box enable the motion on any given fault to affect the others. Warning: need to define the fault geometry with"+
            " a large number of points.")
        #----------------------------------

        # Create text boxes
        self.lon1Edit2 = QLineEdit(self)
        self.lon1Edit2.setText(Param.DParameters['lon1'])
        self.lon1Edit2.setAlignment(Qt.AlignCenter)
        self.lon1Edit2.setValidator(DoubleValidator)
        self.lat1Edit3 = QLineEdit(self)
        self.lat1Edit3.setText(Param.DParameters['lat1'])
        self.lat1Edit3.setAlignment(Qt.AlignCenter)
        self.lat1Edit3.setValidator(DoubleValidator)
        self.lon2Edit4 = QLineEdit(self)
        self.lon2Edit4.setText(Param.DParameters['lon2'])
        self.lon2Edit4.setAlignment(Qt.AlignCenter)
        self.lon2Edit4.setValidator(DoubleValidator)
        self.lat2Edit5 = QLineEdit(self)
        self.lat2Edit5.setText(Param.DParameters['lat2'])
        self.lat2Edit5.setAlignment(Qt.AlignCenter)
        self.lat2Edit5.setValidator(DoubleValidator)
        self.npointiEdit6 = QLineEdit(self)
        self.npointiEdit6.setValidator(QIntValidator())
        self.npointiEdit6.setText(Param.DParameters['npoint1'])
        self.npointiEdit6.setAlignment(Qt.AlignCenter)
        # if not self.UpliftButton1.isChecked():
        #     store_input(self, {'npoint1': '-1'})
        self.bottomLeftEdit9 = QLineEdit(self)
        self.bottomLeftEdit9.setText(Param.DParameters['bottom_left'])
        self.bottomLeftEdit9.setEnabled(False)
        self.bottomLeftEdit9.setAlignment(Qt.AlignCenter)
        self.bottomRightEdit10 = QLineEdit(self)
        self.bottomRightEdit10.setText(Param.DParameters['bottom_right'])
        self.bottomRightEdit10.setEnabled(False)
        self.bottomRightEdit10.setAlignment(Qt.AlignCenter)
        self.upperRightEdit11 = QLineEdit(self)
        self.upperRightEdit11.setText(Param.DParameters['top_right'])
        self.upperRightEdit11.setEnabled(False)
        self.upperRightEdit11.setAlignment(Qt.AlignCenter)
        self.upperLeftEdit12 = QLineEdit(self)
        self.upperLeftEdit12.setText(Param.DParameters['top_left'])
        self.upperLeftEdit12.setEnabled(False)
        self.upperLeftEdit12.setAlignment(Qt.AlignCenter)
        self.nstepiEdit13 = QLineEdit(self)
        self.nstepiEdit13.setValidator(QIntValidator(0,1000))
        self.nstepiEdit13.setText(Param.DParameters['nstep1'])
        self.nstepiEdit13.setEnabled(False)
        self.nstepiEdit13.setAlignment(Qt.AlignCenter)
        
        self.faultParam = Set_Fault_geometry(self) #A fault has been defined in old input file

        # Add text
        Grid1 = QGridLayout()
        Grid1.setSpacing(10)
        Grid1.addWidget(self.UpliftLabel, 0, 0)
        Grid1.addWidget(self.UpliftButton1, 0, 1)
        Grid1.addWidget(self.UpliftButton2, 0, 2)
        Grid1.addWidget(self.UpliftButton3, 0, 3)
        Grid1.addWidget(self.FaultKinematics, 1, 0, 1, 2)
        Grid1.addWidget(self.nstepi, 2, 1)
        Grid1.addWidget(self.nstepiEdit13, 3, 1)
        Grid1.addWidget(self.logVelo, 2, 2)
        Grid1.addWidget(self.logVelob, 3, 2)
        Grid1.addWidget(self.bottomLeft, 4, 1)
        Grid1.addWidget(self.bottomLeftEdit9, 5, 1)
        Grid1.addWidget(self.bottomRight, 4, 2)
        Grid1.addWidget(self.bottomRightEdit10, 5, 2)
        Grid1.addWidget(self.upperRight, 4, 3)
        Grid1.addWidget(self.upperRightEdit11, 5, 3)
        Grid1.addWidget(self.upperLeft, 4, 4)
        Grid1.addWidget(self.upperLeftEdit12, 5, 4)
        Grid1.addWidget(self.ShowFault, 6, 2)
        Grid1.addWidget(self.setFaultButton, 6, 1)
        Grid1.setColumnStretch(0,1)
        Grid1.setColumnStretch(6,3)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        self.tab5.layout.addLayout(Grid1)
        self.tab5.layout.addWidget(Separator)
        Grid2 = QGridLayout()
        Grid2.addWidget(self.UpliftHistory, 0, 0, 1, 3)
        Grid2.addWidget(self.nstepKinTable, 1, 0, 1, 3)
        Grid2.setSpacing(10)
        self.tab5.layout.addLayout(Grid2)
        
        # check input changes 
        self.UpliftButton1.toggled.connect(lambda: self.updateUplift(self.UpliftButton1))
        self.UpliftButton2.toggled.connect(lambda: self.updateUplift(self.UpliftButton2))
        self.UpliftButton3.toggled.connect(lambda: self.updateUplift(self.UpliftButton3))
        self.nfaultEdit1.textChanged.connect(lambda: store_input(
            self, {'nfault': str(self.nfaultEdit1.text())}))
        self.nfaultEdit1.editingFinished.connect(lambda: self.updateKinTable(
            self.nstepKinTable, int(self.nfaultEdit1.text()), int(self.nstepiEdit13.text()), 0, Param))
        self.npointiEdit6.editingFinished.connect(lambda: self.updateTableFault(
            self.faultTable, int(self.nfaultEdit1.text()), self.npointiEdit6.text(), 0, Param))
        self.nstepiEdit13.editingFinished.connect(lambda: self.updateKinTable(
            self.nstepKinTable, int(self.nfaultEdit1.text()), int(self.nstepiEdit13.text()), 0, Param))
        self.nstepKinTable.cellChanged.connect(lambda: self.updateKinTable(
            self.nstepKinTable, int(self.nfaultEdit1.text()), int(self.nstepiEdit13.text()), 1, Param))
        self.bottomLeftEdit9.editingFinished.connect(lambda: store_input(
            self, {'bottom_left': str(self.bottomLeftEdit9.text())}))
        self.bottomRightEdit10.textEdited.connect(lambda: store_input(
            self, {'bottom_right': str(self.bottomRightEdit10.text())}))
        self.upperRightEdit11.editingFinished.connect(lambda: store_input(
            self, {'top_right': str(self.upperRightEdit11.text())}))
        self.upperLeftEdit12.editingFinished.connect(lambda: store_input(
            self, {'top_left': str(self.upperLeftEdit12.text())}))
        self.logVelob.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.logVelob, 'logarithmic_velocity'))
        self.ShowFault.clicked.connect(lambda: self.set_Faults())
        self.setFaultButton.clicked.connect(lambda: self.define_faults())
        
        if int(Param.DParameters['uplift']) == 0:
            self.UpliftButton1.setChecked(True)
            self.UpliftButton2.setChecked(False)
            self.UpliftButton3.setChecked(False)
        elif int(Param.DParameters['uplift']) == 1:
            self.UpliftButton1.setChecked(False)
            self.UpliftButton2.setChecked(True)
            self.UpliftButton3.setChecked(False)
        elif int(Param.DParameters['uplift']) == 2:
            self.UpliftButton1.setChecked(False)
            self.UpliftButton2.setChecked(False)
            self.UpliftButton3.setChecked(True)
            self.ShowFault.setEnabled(True)
            
        if int(Param.DParameters['nstep1']) > 0:
            self.updateKinTable(self.nstepKinTable, int(
                Param.DParameters['nfault']), int(Param.TableParameters['nstep1']), 2, Param)

        self.tab5.setLayout(self.tab5.layout)
        ################# End of Tab 5 ##########################

        ############ Create Ages tab ###########
        
        #------------------------
        
        #self.tab6.layout = VBox

        # Check input changes
        self.ageESRb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageESRb, 'age_ESR_flag'))
        self.ageTLb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageTLb, 'age_TL_flag'))
        self.ageOSLb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageOSLb, 'age_OSL_flag'))
        self.ageAHeb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageAHeb, 'age_AHe_flag'))
        self.ageZHeb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageZHeb, 'age_ZHe_flag'))
        self.ageAFTb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageAFTb, 'age_AFT_flag'))
        self.ageZFTb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageZFTb, 'age_ZFT_flag'))
        self.ageKArb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageKArb, 'age_KAr_flag'))
        self.ageBArb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageBArb, 'age_BAr_flag'))
        self.ageMArb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageMArb, 'age_MAr_flag'))
        self.ageHArb.stateChanged.connect(
            lambda: self.store_input_CheckBox(self.ageHArb, 'age_HAr_flag'))
        self.defaultAgeEdit2.textEdited.connect(lambda: store_input(
            self, {conf.Variable_names['Default_age']: str(self.defaultAgeEdit2.text())}))
        self.defaultAgeEdit2.editingFinished.connect(
            lambda: pgu.isString(self.defaultAgeEdit2.text()))
        
        ################# End of Tab 6 ##########################

        ############ Create Isostasy tab ###########
        # Create text
        self.isostasyb = QCheckBox()
        self.isostasyb.setText('yes/no')
        if int(Param.DParameters['isostasy']) == 1:
            self.isostasyb.setChecked(True)
        else: 
            self.isostasyb.setChecked(False)
        self.EET = QLabel('EET (km):')
        self.EET.setFont(conf.font10)
        self.EET.setAlignment(Qt.AlignRight)
        self.EET.setToolTip(
            "Effective Elastic thickness used to compute isostatic deflection.")
        self.rhoCrust = QLabel('\u03C1c (kg/m\u00B3):')
        self.rhoCrust.setFont(conf.font10)
        self.rhoCrust.setAlignment(Qt.AlignRight)
        self.rhoCrust.setToolTip(
            "Crustal density used to compute isostatic deflection.")
        self.rhoAstenosphere = QLabel('\u03C1a (kg/m\u00B3): ')
        self.rhoAstenosphere.setFont(conf.font10)
        self.rhoAstenosphere.setAlignment(Qt.AlignRight)
        self.rhoAstenosphere.setToolTip(
            "Asthenospheric density used to compute isostatic deflection.")
        self.youngModulus = QLabel('\u0395 (Pa): ')
        self.youngModulus.setFont(conf.font10)
        self.youngModulus.setAlignment(Qt.AlignRight)
        self.youngModulus.setToolTip(
            "Young modulus used to compute isostatic deflection.")
        self.poissonRatio = QLabel('\u03BD: ')
        self.poissonRatio.setFont(conf.font10)
        self.poissonRatio.setAlignment(Qt.AlignRight)
        self.poissonRatio.setToolTip(
            "Poisson's ratio used to compute isostatic deflection.")
        self.nxIsostasy = QLabel('nx (isostasy):')
        self.nxIsostasy.setFont(conf.font10)
        self.nxIsostasy.setAlignment(Qt.AlignCenter)
        self.nxIsostasy.setToolTip("Resolution of the FFT (Fast Fourier Transform) in the x-direction used to solve the flexure equation." +
                                    "Must be a power of 2.")
        self.nyIsostasy = QLabel('ny (isostasy): ')
        self.nyIsostasy.setFont(conf.font10)
        self.nyIsostasy.setAlignment(Qt.AlignCenter)
        self.nyIsostasy.setToolTip("Resolution of the FFT in the y-direction used to solve the flexure equation." +
                                    "Must be a power of 2.")
        self.GridResolution = QLabel('Grid resolution:')
        self.GridResolution.setFont(conf.fontLine12)
        self.GridResolution.setAlignment(Qt.AlignLeft)
        self.Lithosphere = QLabel('Lithosphere characteristics:')
        self.Lithosphere.setFont(conf.fontLine12)
        self.Lithosphere.setAlignment(Qt.AlignLeft)
        self.isostasy = QLabel('Enable isostasy?')
        self.isostasy.setFont(conf.fontLine12)
        self.isostasy.setAlignment(Qt.AlignLeft)
        self.isostasy.setToolTip(
            "Checking this box triggers flexural isostatic adjustment to Pecube.")

        # Create text boxes
        self.EETEdit2 = QLineEdit(self)
        self.EETEdit2.setText(Param.DParameters['EET'])
        self.EETEdit2.setValidator(DoubleValidator)
        self.rhoCrustEdit3 = QLineEdit(self)
        self.rhoCrustEdit3.setText(Param.DParameters['rho_crust'])
        self.rhoAstenosphereEdit4 = QLineEdit(self)
        self.rhoAstenosphereEdit4.setText(Param.DParameters['rho_astenosphere'])
        self.youngModulusEdit5 = QLineEdit(self)
        self.youngModulusEdit5.setText(Param.DParameters['young_modulus'])
        self.poissonRatioEdit6 = QLineEdit(self)
        self.poissonRatioEdit6.setText(Param.DParameters['poisson_ratio'])
        self.nxIsostasyEdit7 = QLineEdit(self)
        self.nxIsostasyEdit7.setText(Param.DParameters['nx_isostasy'])
        self.nxIsostasyEdit7.setValidator(QIntValidator(0,100000))
        self.nxIsostasyEdit7.setAlignment(Qt.AlignCenter)
        self.nyIsostasyEdit8 = QLineEdit(self)
        self.nyIsostasyEdit8.setText(Param.DParameters['ny_isostasy'])
        self.nyIsostasyEdit8.setValidator(QIntValidator(0,100000))
        self.nyIsostasyEdit8.setAlignment(Qt.AlignCenter)

        # Add text
        Grid = QGridLayout()
        Grid.setSpacing(10)
        VBox = QVBoxLayout()
        HBox = QHBoxLayout()
        HBox.addWidget(self.isostasy)
        HBox.addWidget(self.isostasyb)
        HBox.addStretch(1)
        VBox.addLayout(HBox)
        VBox.addStretch(1)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        VBox.addWidget(Separator)

        Grid.addWidget(self.Lithosphere, 0, 0)
        Grid.addWidget(self.EET, 1, 0)
        Grid.addWidget(self.rhoCrust, 2, 0)
        Grid.addWidget(self.rhoAstenosphere, 3, 0)
        Grid.addWidget(self.youngModulus, 4, 0)
        Grid.addWidget(self.poissonRatio, 5, 0)
        # Add text boxes
        Grid.addWidget(self.EETEdit2, 1, 1)
        Grid.addWidget(self.rhoCrustEdit3, 2, 1)
        Grid.addWidget(self.rhoAstenosphereEdit4, 3, 1)
        Grid.addWidget(self.youngModulusEdit5, 4, 1)
        Grid.addWidget(self.poissonRatioEdit6, 5, 1)
        Grid.setColumnStretch(3, 2)
        VBox.addLayout(Grid)
        VBox.addStretch(1)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        Grid2 = QGridLayout()
        Grid2.addWidget(Separator, 13, 0, 1, 4)
        Grid2.addWidget(self.GridResolution, 14, 0)
        Grid2.addWidget(self.nxIsostasy, 16, 1)
        Grid2.addWidget(self.nyIsostasy, 16, 2)
        Grid2.addWidget(self.nxIsostasyEdit7, 17, 1)
        Grid2.addWidget(self.nyIsostasyEdit8, 17, 2)
        Grid2.setColumnStretch(3, 2)
        VBox.addLayout(Grid2)
        VBox.addStretch(2)
        self.tab7.setLayout(VBox)
        
        # Check input parameter changes
        self.isostasyb.stateChanged.connect(
            lambda: self.store_input_CheckBox_savePTT(self.isostasyb, 'isostasy'))
        self.EETEdit2.textEdited.connect(lambda: store_input(
            self, {'EET': str(self.EETEdit2.text())}))
        self.rhoCrustEdit3.textEdited.connect(lambda: store_input(
            self, {'rho_crust': str(self.rhoCrustEdit3.text())}))
        self.rhoAstenosphereEdit4.textEdited.connect(lambda: store_input(
            self, {'rho_astenosphere': str(self.rhoAstenosphereEdit4.text())}))
        self.youngModulusEdit5.textEdited.connect(lambda: store_input(
            self, {'young_modulus': str(self.youngModulusEdit5.text())}))
        self.poissonRatioEdit6.textEdited.connect(lambda: store_input(
            self, {'poisson_ratio': str(self.poissonRatioEdit6.text())}))
        self.nxIsostasyEdit7.textEdited.connect(lambda: store_input(
            self, {'nx_isostasy': str(self.nxIsostasyEdit7.text())}))
        self.nyIsostasyEdit8.textEdited.connect(lambda: store_input(
            self, {'ny_isostasy': str(self.nyIsostasyEdit8.text())}))
        
        ################# End of Tab 7 ##########################

        ############ Create Inversion tab ###########
        # Create text
        self.tab8.layout = QVBoxLayout(self)
        self.tab8.layout.setSpacing(10)
        self.RunParallel = QLabel('Number of cores:')
        self.RunParallel.setFont(conf.font10)
        self.RunParallel.setAlignment(Qt.AlignCenter)
        self.RunParallel.setToolTip("Number of cores to be used inversion.")
        self.NbCores = QLineEdit()
        self.NbCores.setText(Param.DParameters['NbCores'])
        self.NbCores.setAlignment(Qt.AlignCenter)
        self.NbCores.setValidator(QIntValidator(0,1000))
        self.maxNIt = QLabel('Maximum number of iterations:')
        self.maxNIt.setFont(conf.font10)
        self.maxNIt.setAlignment(Qt.AlignCenter)
        self.maxNIt.setToolTip(
            "Total number of NA iterations to be performed (not includint he first one)")
        self.sampleSizeFirst = QLabel('Sample size for the first iteration:')
        self.sampleSizeFirst.setFont(conf.font10)
        self.sampleSizeFirst.setAlignment(Qt.AlignCenter)
        self.sampleSizeFirst.setToolTip(
            "Number of model runs (and mode parameter sets) performed during the first, initial iteration.")
        self.sampleSizeOther = QLabel('Sample size for all other iterations:')
        self.sampleSizeOther.setFont(conf.font10)
        self.sampleSizeOther.setAlignment(Qt.AlignCenter)
        self.sampleSizeOther.setToolTip(
            "Number of model runs performed during every subsequent iteration.")
        self.nbCellResample = QLabel('Number of cells to resample: ')
        self.nbCellResample.setFont(conf.font10)
        self.nbCellResample.setAlignment(Qt.AlignCenter)
        self.nbCellResample.setToolTip(
            "Number of Voronoi cells that will be resampled uniformly during any given iteration (except the first one)")
        self.misfitWAge = QLabel('Misfit weight Age: ')
        self.misfitWAge.setFont(conf.font10)
        self.misfitWAge.setAlignment(Qt.AlignCenter)
        self.misfitWAge.setToolTip(
            "Factor or weight mutliplying the part of the misfit function that measures the difference between observed and predicted ages.")
        self.misfitWFTLD = QLabel('Misfit weight FTLD: ')
        self.misfitWFTLD.setFont(conf.font10)
        self.misfitWFTLD.setAlignment(Qt.AlignCenter)
        self.misfitWFTLD.setToolTip(
            "Factor or weight mutliplying the part of the misfit function that measures the difference between observed and track length distributions.")
        self.misfitWTH = QLabel('Misfit weight TH:')
        self.misfitWTH.setFont(conf.font10)
        self.misfitWTH.setAlignment(Qt.AlignCenter)
        self.misfitWTH.setToolTip(
            "Factor or weight mutliplying the part of the misfit function that measures the difference between observed and predicted thermal histories.")
        self.misfitW43He = QLabel('Misfit weight 43He: ')
        self.misfitW43He.setFont(conf.font10)
        self.misfitW43He.setAlignment(Qt.AlignCenter)
        self.misfitW43He.setToolTip(
            "Factor or weight mutliplying the part of the misfit function that measures the difference between observed and predicted 43He ages.")
        self.misfitWTL = QLabel('Misfit weight TL: ')
        self.misfitWTL.setFont(conf.font10)
        self.misfitWTL.setAlignment(Qt.AlignCenter)
        self.misfitWTL.setToolTip(
            "Factor or weight mutliplying the part of the misfit function that measures the difference between observed and predicted Thermoluminescence data.")
        self.misfitSlopeb = QCheckBox()
        self.misfitSlopeb.setText('yes/no')
        self.misfitSlope = QLabel("Misfit slope")
        self.misfitSlope.setFont(conf.font10)
        self.misfitSlope.setAlignment(Qt.AlignCenter)
        self.misfitSlope.setToolTip(
            "To set whether the ages themselves are to be used to define the age part of the misfit function (unchecked) or if the slope and intercept of the age-elevation relationship derived form the predictions and observations are to be used (checked).")
        self.misfitCorrb = QCheckBox()
        self.misfitCorrb.setText('yes/no')
        self.misfitCorr = QLabel("Misfit corrected")
        self.misfitCorr.setFont(conf.font10)
        self.misfitCorr.setAlignment(Qt.AlignCenter)
        self.misfitCorr.setToolTip(
            "To correct the misfit function by the number of observations minus the number of parameters +1 (checked).")
        self.IterationLabel = QLabel('Iteration performance:')
        self.IterationLabel.setFont(conf.fontLine12)
        self.MisfitFunctionLabel = QLabel("Misfit function controls:")
        self.MisfitFunctionLabel.setFont(conf.fontLine12)
        self.MisfitTypeLabel = QLabel("Misfit function: ")
        self.MisfitTypeLabel.setFont(conf.font10)
        self.MisfitTypeLabel.setToolTip("Choose the misfit function to be computed by Pecube.")
        self.MisfitTypeLabel.setAlignment(Qt.AlignCenter)
        self.MisfitType = QComboBox()
        self.MisfitType.setMinimumContentsLength(15)
        self.MisfitType.addItems(['Chi-squared','Reduced Chi-squared','L2-norm','L1-norm','Reduced L1-norm','log-scale'])
        self.MisfitType.setCurrentIndex(int(Param.DParameters['misfit_function'])-1)
        self.MisfitOn43Label = QLabel("Misfit 43He on:")
        self.MisfitOn43Label.setFont(conf.font10)
        self.MisfitOn43Label.setAlignment(Qt.AlignCenter)
        self.MisfitOn43 = QComboBox()
        # self.MisfitOn43.addItems(['4He/3He spectrum',"Normalized step age"])
        self.MisfitOn43.addItems(['4He/3He spectrum'])
        self.WeightAHE = QLabel('AHE weight: ')
        self.WeightAHE.setFont(conf.font10)
        self.WeightAHE.setToolTip("Multiplication factor on the AHE misfit value.")
        self.WeightAFT = QLabel('AFT weight: ')
        self.WeightAFT.setFont(conf.font10)
        self.WeightAFT.setToolTip("Multiplication factor on the AFT misfit value.")
        self.WeightZHE = QLabel('ZHE weight: ')
        self.WeightZHE.setFont(conf.font10)
        self.WeightZHE.setToolTip("Multiplication factor on the ZHE misfit value.")
        self.WeightZFT = QLabel('ZFT weight: ')
        self.WeightZFT.setFont(conf.font10)
        self.WeightZFT.setToolTip("Multiplication factor on the ZFT misfit value.")
        self.WeightBAR = QLabel('BAR weight: ')
        self.WeightBAR.setFont(conf.font10)
        self.WeightBAR.setToolTip("Multiplication factor on the BAR misfit value.")
        self.WeightMAR = QLabel('MAR weight: ')
        self.WeightMAR.setFont(conf.font10)
        self.WeightMAR.setToolTip("Multiplication factor on the MAR misfit value.")
        self.WeightHAR = QLabel('HAR weight: ')
        self.WeightHAR.setFont(conf.font10)
        self.WeightHAR.setToolTip("Multiplication factor on the HAR misfit value.")
        self.WeightKAR = QLabel('KAR weight: ')
        self.WeightKAR.setFont(conf.font10)
        self.WeightKAR.setToolTip("Multiplication factor on the KAR misfit value.")
        self.ErrorOnPredLabel = QLabel('Error on predictions (%):')
        self.ErrorOnPredLabel.setToolTip("Use error on prediction in %.This will be applied in the misfit calculation.")
        self.ErrorOnPredLabel.setFont(conf.font10)
        self.ErrorOnPredEdit = QLineEdit(Param.DParameters['error_predictions'])
        self.ErrorOnPredEdit .setAlignment(Qt.AlignCenter)
        # Inversion mode
        self.InvModeTitle = QLabel('Choose inversion mode:')
        self.InvModeTitle.setFont(conf.fontLine12)
        self.InvModeLabel = QLabel('Inversion mode:')
        self.InvModeLabel.setFont(conf.font10)
        self.InvModeLabel.setAlignment(Qt.AlignRight)
        self.InvModeLabel.setToolTip("Mode of inversion: 'inversion' uses the NA algorithm, 'batch' run \n"+ 
                                     "Pecube models with in batch mode and provide prediction for a set of \n"+
                                     "pre-defined parameter values.")
        items = ['inversion','inversion from previous file','batch']
        self.InvModeCombo = QComboBox()
        self.InvModeCombo.addItems(items)
        self.InvModeCombo.setMinimumContentsLength(15)
        self.InvModeCombo.setCurrentIndex(int(Param.DParameters[conf.Variable_names['use inversion or batch']]))
        self.nIntervalsLabel = QLabel("Number of intervals: ")
        self.nIntervalsLabel.setFont(conf.font10)
        self.nIntervalsLabel.setToolTip("Set the number of intervals for the range of values for the parameter to batch.")
        self.nIntervalsEdit = QLineEdit()
        self.nIntervalsEdit.setText(str(int(Param.DParameters['sample_size_for_first_iteration'])-1))
        
        # Inversion mode grid
        Grid0 = QGridLayout()
        Grid0.addWidget(self.InvModeTitle, 0, 0, 1, 3)
        Grid0.addWidget(self.InvModeLabel, 2, 1)
        Grid0.addWidget(self.InvModeCombo, 2, 2)
        Grid0.addWidget(self.nIntervalsLabel, 3, 1)
        Grid0.addWidget(self.nIntervalsEdit, 3, 2)
        Grid0.setColumnStretch(4, 1)
        # Add text
        Grid1 = QGridLayout()
        Grid1.addWidget(self.IterationLabel,0,0,1,3)
        Grid1.addWidget(self.RunParallel,1,1)
        Grid1.addWidget(self.NbCores,2,1)
        Grid1.addWidget(self.maxNIt, 3, 1)
        Grid1.addWidget(self.sampleSizeFirst, 3, 2)
        Grid1.addWidget(self.sampleSizeOther, 5, 1)
        Grid1.addWidget(self.nbCellResample, 5, 2)
        Grid2 = QGridLayout()
        Grid2.addWidget(self.MisfitFunctionLabel,0,0,1,3)
        Grid2.addWidget(self.MisfitTypeLabel, 1, 1)
        Grid2.addWidget(self.MisfitType,1,2)
        Grid2.addWidget(self.misfitWAge, 2, 1)
        Grid2.addWidget(self.misfitWFTLD, 3, 1)
        Grid2.addWidget(self.misfitWTH, 3, 2)
        Grid2.addWidget(self.misfitW43He, 5, 1)
        Grid2.addWidget(self.MisfitOn43Label, 5, 2)
        Grid2.addWidget(self.ErrorOnPredLabel, 5, 3)
        Grid2.addWidget(self.misfitWTL, 3,3)
        # Grid2.addWidget(self.misfitSlopeb, 7,2)
        # Grid2.addWidget(self.misfitSlope, 7, 1)

        # Size of boxes
        Grid1.setColumnStretch(0, 1)
        Grid1.setColumnStretch(7, 2)
        Grid1.setRowStretch(7, 1)
        Grid2.setColumnStretch(0, 1)
        Grid2.setColumnStretch(4, 2)
        Grid2.setRowStretch(6, 1)

        # Create text boxes
        self.maxNItEdit1 = QLineEdit(self)
        self.maxNItEdit1.setValidator(QIntValidator(0,1000))
        self.maxNItEdit1.setText(Param.DParameters['maximum_number_of_iterations'])
        self.maxNItEdit1.setAlignment(Qt.AlignCenter)
        self.sampleSizeFirstEdit2 = QLineEdit(self)
        self.sampleSizeFirstEdit2.setValidator(QIntValidator(0,1000))
        self.sampleSizeFirstEdit2.setText(Param.DParameters['sample_size_for_first_iteration'])
        self.sampleSizeFirstEdit2.setAlignment(Qt.AlignCenter)
        self.sampleSizeOtherEdit3 = QLineEdit(self)
        self.sampleSizeOtherEdit3.setValidator(QIntValidator(0,1000))
        self.sampleSizeOtherEdit3.setText(Param.DParameters['sample_size_for_all_other_iterations'])
        self.sampleSizeOtherEdit3.setAlignment(Qt.AlignCenter)
        self.nbCellResampleEdit4 = QLineEdit(self)
        self.nbCellResampleEdit4.setValidator(QIntValidator(0,1000))
        self.nbCellResampleEdit4.setText(Param.DParameters['number_of_cells_to_resample'])
        self.nbCellResampleEdit4.setAlignment(Qt.AlignCenter)
        self.misfitWAgeBox= QCheckBox()
        self.misfitWAgeBox.setChecked(False)
        self.misfitWFTLDEdit6 = QLineEdit(self)
        self.misfitWFTLDEdit6.setText(Param.DParameters['misfit_weight_FTLD'])
        self.misfitWFTLDEdit6.setValidator(DoubleValidator)
        self.misfitWFTLDEdit6.setAlignment(Qt.AlignCenter)
        self.misfitWTHEdit7 = QLineEdit(self)
        self.misfitWTHEdit7.setText(Param.DParameters['misfit_weight_TH'])
        self.misfitWTHEdit7.setValidator(DoubleValidator)
        self.misfitWTHEdit7.setAlignment(Qt.AlignCenter)
        self.misfitW43HeEdit8 = QLineEdit(self)
        self.misfitW43HeEdit8.setText(Param.DParameters['misfit_weight_43He'])
        self.misfitW43HeEdit8.setValidator(DoubleValidator)
        self.misfitW43HeEdit8.setAlignment(Qt.AlignCenter)
        self.misfitWTLEdit9 = QLineEdit(self)
        self.misfitWTLEdit9.setText(Param.DParameters['misfit_weight_TL'])
        self.misfitWTLEdit9.setValidator(DoubleValidator)
        self.misfitWTLEdit9.setAlignment(Qt.AlignCenter)
        self.weightAHEedit = QLineEdit()
        self.weightAHEedit.setText(Param.DParameters['misfit_weight_AHE'])
        self.weightAHEedit.setValidator(DoubleValidator)
        self.weightAHEedit.setAlignment(Qt.AlignCenter)
        self.weightAFTedit = QLineEdit("1.0")
        self.weightAFTedit.setValidator(DoubleValidator)
        self.weightAFTedit.setAlignment(Qt.AlignCenter)
        self.weightZHEedit = QLineEdit("1.0")
        self.weightZHEedit.setValidator(DoubleValidator)
        self.weightZHEedit.setAlignment(Qt.AlignCenter)
        self.weightZFTedit = QLineEdit("1.0")
        self.weightZFTedit.setValidator(DoubleValidator)
        self.weightZFTedit.setAlignment(Qt.AlignCenter)
        self.weightBARedit = QLineEdit("1.0")
        self.weightBARedit.setValidator(DoubleValidator)
        self.weightBARedit.setAlignment(Qt.AlignCenter)
        self.weightMARedit = QLineEdit("1.0")
        self.weightMARedit.setValidator(DoubleValidator)
        self.weightMARedit.setAlignment(Qt.AlignCenter)
        self.weightHARedit = QLineEdit("1.0")
        self.weightHARedit.setValidator(DoubleValidator)
        self.weightHARedit.setAlignment(Qt.AlignCenter)
        self.weightKARedit = QLineEdit("1.0")
        self.weightKARedit.setValidator(DoubleValidator)
        self.weightKARedit.setAlignment(Qt.AlignCenter)
        
        self.update_InvMode()
        # Add text boxes
        Grid1.addWidget(self.maxNItEdit1, 4, 1)
        Grid1.addWidget(self.sampleSizeFirstEdit2, 4, 2)
        Grid1.addWidget(self.sampleSizeOtherEdit3, 6, 1)
        Grid1.addWidget(self.nbCellResampleEdit4, 6, 2)
        Grid2.addWidget(self.misfitWAgeBox, 2, 2)
        Grid2.addWidget(self.misfitWFTLDEdit6, 4, 1)
        Grid2.addWidget(self.misfitWTHEdit7, 4, 2)
        Grid2.addWidget(self.misfitW43HeEdit8, 6, 1)
        Grid2.addWidget(self.MisfitOn43, 6, 2)
        Grid2.addWidget(self.ErrorOnPredEdit, 6, 3)
        Grid2.addWidget(self.misfitWTLEdit9, 4, 3)
        self.tab8.layout.addLayout(Grid0)
        self.tab8.layout.addStretch(1)
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        self.tab8.layout.addWidget(Separator)
        self.tab8.layout.addLayout(Grid1)
        self.tab8.layout.addStretch(1)
        self.tab8.layout.addWidget(Separator)
        self.tab8.layout.addLayout(Grid2)
        self.tab8.layout.addStretch(1)

        # Check input parameter changes
        self.NbCores.textEdited.connect(lambda: store_input(
            self, {'NbCores': str(self.NbCores.text())}))
        self.maxNItEdit1.textChanged.connect(lambda: store_input(
            self, {'maximum_number_of_iterations': str(self.maxNItEdit1.text())}))
        self.sampleSizeFirstEdit2.textChanged.connect(lambda: store_input(
            self, {'sample_size_for_first_iteration': str(self.sampleSizeFirstEdit2.text())}))
        self.sampleSizeOtherEdit3.textChanged.connect(lambda: store_input(
            self, {'sample_size_for_all_other_iterations': str(self.sampleSizeOtherEdit3.text())}))
        self.nbCellResampleEdit4.textChanged.connect(lambda: store_input(
            self, {'number_of_cells_to_resample': str(self.nbCellResampleEdit4.text())}))
        self.misfitWAgeBox.stateChanged.connect(lambda: self.open_AgeWeight())
        self.weightAHEedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_AHE': str(self.weightAHEedit.text())}))
        self.weightAFTedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_AFT': str(self.weightAFTedit.text())}))
        self.weightZHEedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_ZHE': str(self.weightZHEedit.text())}))
        self.weightZFTedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_ZFT': str(self.weightZFTedit.text())}))
        self.weightBARedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_BAR': str(self.weightBARedit.text())}))
        self.weightMARedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_MAR': str(self.weightMARedit.text())}))
        self.weightHARedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_HAR': str(self.weightHARedit.text())}))
        self.weightKARedit.editingFinished.connect(lambda: store_input(
            self, {'misfit_weight_KAR': str(self.weightKARedit.text())}))
        self.misfitWFTLDEdit6.textEdited.connect(lambda: store_input(
            self, {'misfit_weight_FTLD': str(self.misfitWFTLDEdit6.text())}))
        self.misfitWFTLDEdit6.editingFinished.connect(
            lambda: pgu.isFloat(str(self.misfitWFTLDEdit6.text())))
        self.misfitWTHEdit7.textEdited.connect(lambda: store_input(
            self, {'misfit_weight_TH': str(self.misfitWTHEdit7.text())}))
        self.misfitWTHEdit7.editingFinished.connect(
            lambda: pgu.isFloat(str(self.misfitWTHEdit7.text())))
        self.misfitW43HeEdit8.textEdited.connect(lambda: store_input(
            self, {'misfit_weight_43He': str(self.misfitW43HeEdit8.text())}))
        self.misfitW43HeEdit8.editingFinished.connect(
            lambda: pgu.isFloat(str(self.misfitW43HeEdit8.text())))
        self.misfitWTLEdit9.textEdited.connect(lambda: store_input(
            self, {'misfit_weight_TL': str(self.misfitWTLEdit9.text())}))
        self.misfitWTLEdit9.editingFinished.connect(
            lambda: pgu.isFloat(str(self.misfitWTLEdit9.text())))
        self.ErrorOnPredEdit.textEdited.connect(lambda: store_input(
            self, {'error_predictions': str(self.ErrorOnPredEdit.text())}))
        self.misfitSlopeb.stateChanged.connect(lambda: self.store_input_CheckBox_savePTT(self.misfitSlopeb, 'misfit_slope'))
        self.MisfitType.currentIndexChanged.connect(lambda: store_input(self, {'misfit_function': str(self.MisfitType.currentIndex()+1)}))
        self.MisfitOn43.currentIndexChanged.connect(lambda: store_input(self, {conf.Variable_names["Inversion_43"]: str(self.MisfitOn43.currentIndex())}))
        self.InvModeCombo.currentIndexChanged.connect(lambda: self.update_InvMode())
        self.nIntervalsEdit.editingFinished.connect(lambda: self.update_ninterval())
        
        self.tab8.setLayout(self.tab8.layout)
        ################# End of Tab 8 ##########################

        # ############ Create Output tab ###########
        self.tab9.layout = QGridLayout(self)
        self.tab9.layout.setSpacing(10)
        self.Label17c = QLabel('Save file:')
        self.Label17c.setFont(conf.fontLine12)
        self.Label17c.setAlignment(Qt.AlignLeft)
        self.echoInputFile = QLabel('Echo input file:')
        self.echoInputFile.setFont(conf.fontLine12)
        self.debugb = QCheckBox()
        self.debug = QLabel('Debug: ')
        self.debug.setFont(conf.font10)
        self.debug.setAlignment(Qt.AlignRight)
        self.debug.setToolTip("Debug information to be send to a log file.")
        self.savePTTb = QCheckBox()
        self.savePTT = QLabel('save PTT paths: ')
        self.savePTT.setFont(conf.font10)
        self.savePTT.setAlignment(Qt.AlignRight)
        self.savePTT.setToolTip(
            "Save thermal histories of observation points.")
        self.saveEroVolb = QCheckBox()
        self.saveEroVol = QLabel('save eroded volume:')
        self.saveEroVol.setFont(conf.font10)
        self.saveEroVol.setAlignment(Qt.AlignRight)
        self.saveEroVol.setToolTip("Save eroded volume through time.")
        self.saveCoolingRates = QLabel('Cooling rates:')
        self.saveCoolingRates.setFont(conf.font10)
        self.saveCoolingRates.setAlignment(Qt.AlignRight)
        self.saveCoolingRates.setToolTip("Save Tt paths for all nodes to compute cooling rates.")
        self.saveCoolingRatesBox = QCheckBox()
        if self.Param.DParameters['save_cooling_rates'] == '1':
            self.saveCoolingRatesBox.setChecked(True)
        self.Label22 = QComboBox()
        self.Label22.addItem('no echo')
        self.Label22.addItem('with parameter and value')
        self.Label22.addItem('with default value')
        self.Label22.addItem('with short description')
        self.Label22.addItem('list parameters only')
        # Get previous value for echo_input_file
        self.Label22.setCurrentIndex(int(Param.DParameters['echo_input_file']))
        
        # Set grids
        VBox = QVBoxLayout()
        Separator = QFrame()
        Separator.setFrameShape(QFrame.HLine)
        Separator.setFrameShadow(QFrame.Raised)
        Separator.setLineWidth(1)
        Separator.setMidLineWidth(1)
        Grid = QGridLayout()
        Grid.setSpacing(10)
        Grid.addWidget(self.Label17c, 0, 0, 1, 2)
        Grid.addWidget(self.debug, 3, 2)
        Grid.addWidget(self.debugb, 3, 3)
        Grid.addWidget(self.savePTT, 4, 2)
        Grid.addWidget(self.savePTTb, 4, 3)
        Grid.addWidget(self.saveEroVol, 5, 2)
        Grid.addWidget(self.saveEroVolb, 5, 3)
        Grid.addWidget(self.saveCoolingRates, 6, 2)
        Grid.addWidget(self.saveCoolingRatesBox, 6, 3)
        Grid.setColumnStretch(3, 1)
        VBox.addLayout(Grid)
        VBox.addStretch(1)
        VBox.addWidget(Separator)
        Gbox = QHBoxLayout()
        Gbox.addWidget(self.echoInputFile)
        Gbox.addWidget(self.Label22)
        Gbox.addStretch(1)
        VBox.addLayout(Gbox)
        VBox.addStretch(3)
        
        self.tab9.layout = VBox

        # Check input parameter changes
        self.debugb.stateChanged.connect(
            lambda: self.store_input_CheckBox_savePTT(self.debugb, 'debug'))
        self.savePTTb.stateChanged.connect(
            lambda: self.store_input_CheckBox_savePTT(self.savePTTb, 'save_PTT_paths'))
        self.saveEroVolb.stateChanged.connect(
            lambda: self.store_input_CheckBox_savePTT(self.saveEroVolb, 'save_eroded_volume'))
        self.saveCoolingRatesBox.stateChanged.connect(
            lambda: self.store_input_CheckBox_savePTT(self.saveCoolingRatesBox, 'save_cooling_rates'))
        self.Label22.currentIndexChanged.connect(
            lambda: self.store_input_ComboBoxEcho(self.Label22, 'echo_input_file'))
        #self.store_input_CheckBox(self.savePTTb, 'save_PTT_paths')
        self.FaultAdvectb.stateChanged.connect(
            lambda: self.store_input_CheckBox_savePTT(self.FaultAdvectb, 'fault_advect_flag'))

        self.tab9.setLayout(self.tab9.layout)
        ################# End of Tab 9 ##########################
        self.AgeCombo.currentIndexChanged.connect(
            lambda: self.store_input_ComboBoxAge(self.AgeCombo, conf.Variable_names['Mode_age_computation']))
        # Check old input
        self.store_input_ComboBoxAge(self.AgeCombo, conf.Variable_names['Mode_age_computation'])
        if self.AgesParam:
            self.AgesParam.updateComboBox()


    #----------------------------------------------------------
    def Convert2LatLon(self):
        """ To convert the model dimension given in km to lat lon. """
        
        # First get the distance steps
        dx = float(self.dLonEdit6.text())
        dy = float(self.dLatEdit7.text())
        # Then convert to lat/long
        dx = dx*0.008333
        dy = dy*0.008333
        self.dLonEdit6.setText(str(round(dx,6)))
        self.dLatEdit7.setText(str(round(dy,6)))
            
    #----------------------------------------------------------    
    def define_faults(self):
        """ This function shows the window to define fault geometry."""
        
        self.wfault = U_interface.WinExtraParameters(1000, 800)
        self.wfault.OkButton.clicked.connect(lambda: self.wfault.close())
        self.faultParam = Set_Fault_geometry(self)
        self.wfault.layout.addLayout(self.faultParam.Layout, 0, 0, 1, 3)
        self.wfault.layout.setRowStretch(0,3)
        Q2 = QWidget()
        Q2.setLayout(self.wfault.layout)
        self.wfault.setCentralWidget(Q2)
        conf.WindowsOpen.append(self.wfault)
        self.ShowFault.setEnabled(True)
        self.wfault.show()

    #----------------------------------------------------------
    def file_open_spm(self):
        """This function allows to load files from DEM or spm into Pecube."""
        
        #Get the file(s) to load
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dialog = QFileDialog()
        dialog.setDirectory(os.path.join(conf.PecubeFolderPath,self.folderName))
        self.fileName, _ = dialog.getOpenFileNames(
            self, "QFileDialog.getOpenFileName()", "", "All Files (*);;Pecube Files (*.in)", options=options)
        if self.fileName:
            # Check if it comes from iSOSIA 
            self.extracted_Fnames = self.fileName[0].split("/")
            str1 = ''
            strlist = self.extracted_Fnames[-2:]
            str2 = str1.join(strlist)
            Nb_char = len(str2)+2  # count number of character after folder name
            self.folderName = self.fileName[0][:-Nb_char]
            nfiles = len(self.fileName)

            if nfiles == 1: # a single topo file, can be single column file or .tif
                self.Raster_flag = 0
                store_input(self, {'Raster_flag': str(self.Raster_flag)})
                TargetPath = os.path.join(self.PFolder,'data')
                # a single file means topo file
                # Copy the file to the working Pecube directory 
                if not os.path.exists(os.path.join(TargetPath,os.path.basename(self.fileName[0]))):
                    shutil.copy(self.fileName[0],TargetPath)
                # If DEM (.tif), try to get some information
                topo_file_name = self.extracted_Fnames[-1]
                if topo_file_name[-3:] == 'tif':
                    self.Raster_flag = 1
                    store_input(self, {'Raster_flag': str(self.Raster_flag)})
                    self.topo_file_nameEdit1.setText(self.extracted_Fnames[-1][:-4])
                    path = os.path.join(self.PFolder,"data",self.topo_file_nameEdit1.text()+'.tif')
                    elevation,nx,ny,dlon,dlat,lon0,lat0 = GIS.load_raster_file(path)
                    self.nxEdit2.setText(str(nx))
                    self.nyEdit3.setText(str(ny))
                    self.lon0Edit4.setText(str(round(lon0,6)))
                    self.lat0Edit5.setText(str(round(lat0,6)))
                    self.dLonEdit6.setText(str(round(abs(dlon),6)))
                    self.dLatEdit7.setText(str(round(abs(dlat),6)))
                    # Create readable topo file for Pecube
                    path = os.path.join(self.PFolder,"data",self.topo_file_nameEdit1.text())
                    np.savetxt(f"{path}",elevation,fmt="%.2f")
                    # close Raster
                    data = None
                else:
                    self.topo_file_nameEdit1.setText(self.extracted_Fnames[-1])

            else: # files from spm
                self.Raster_flag = 0
                store_input(self, {'Raster_flag': str(self.Raster_flag)})
                self.spmFolder = os.path.join(self.PFolder,'data','SPM')
                if os.path.exists(self.spmFolder):
                    win = QMessageBox.question(self, 'Pecube Message',
                                                'This topographic folder already exists, do you want to erase it?',
                                                QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Cancel)
                    if win == QMessageBox.Ok:
                        shutil.rmtree(self.spmFolder)
                    else:
                        return
                os.mkdir(self.spmFolder)
                # Provide '/spm' name to topo_file_name
                self.topo_file_nameEdit1.setText('SPM/')
                for i in range(nfiles):
                    shutil.copy(self.fileName[i],self.spmFolder)
                # Check whether files come from iSOSIA
                self.inputFolder = os.path.join(self.folderName, 'input','SPM.mat') #iSOSIA files
                if os.path.exists(self.inputFolder):
                    self.iSOSIA = 1
                self.spm = spm_objects(self,self.fileName)
                # Update nfault
                self.nfaultEdit1.setText("1")
            if self.setTopo: #if topography visualization already open
                # Update initial topography
                self.setTopo.updateInitTopo()
            
            elevation = nx = ny =dlon = dlat= lon0=lat0 = None
    
    #----------------------------------------------------------
    def  load_BuildTopo(self):
        """ To load window to set sinusoidal topography."""
        self.build_topo = build_Topography(self)
    
    #----------------------------------------------------------
    def  load_from_map(self):
        """ To load topo from interactive map."""
        # Get and save raster file
        try:
            self.get_topo = MapWindow(self)
        except Exception as E:
            raise
            QErrorMessage(self).showMessage("Get Topo function encountered an error when trying to open: <br>"+
                                            str(type(E)) + "<br>" +
                                            str(E.args))
            return
        
    #----------------------------------------------------------
    def  set_AgesParam(self):
        """ This function call the class set_Thermochron_Parameters and keep track of it. """
        # if set_Thermochron_Parameters is already opened, then make sure the tab is visible
        try:
            if not self.HeTab: #If not already opened
                self.AgesParam = th.set_Thermochron_Parameters(self, self.Param, self.PFolder)
                Q = QWidget()
                Q.setLayout(self.AgesParam.Layout)
                self.HeTab = 1
                getTab = 0
                for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
                    if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Ages":
                        getTab = 1
                if not getTab:
                    self.MainWin.StackedTabs[self.MainWin.CurrentTabs].addTab(Q,"Ages")
                
            #set tab visible
            for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
                if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Ages":
                    self.MainWin.StackedTabs[self.MainWin.CurrentTabs].setTabVisible(i,True)
        except AttributeError:
            print('\n Attribute Error in "set_AgesParam"')
            pass

    #----------------------------------------------------------
    def get_EFolding(self):
        """Call the classe EFolding to get the heat production rate and e-folding depth
        from the user"""
        
        if self.EfoldingHeatProd.isChecked() == True:
            self.Heat_production_flag = 1
            self.wEF = U_interface.WinExtraParameters(200, 100)
            self.wEF.OkButton.clicked.connect(lambda: self.wEF.close())
            self.EFParam = eFolding_HP(self,self.Param)
            self.wEF.layout.addLayout(self.EFParam.GBox, 0, 0, 1, 3)
            Q2 = QWidget()
            Q2.setLayout(self.wEF.layout)
            self.wEF.setCentralWidget(Q2)
            conf.WindowsOpen.append(self.wEF)
            self.wEF.show()
            store_input(self, {conf.Variable_names['Heat_production_flag']: '1'})
        else:
            self.Heat_production_flag = 0
            if self.GeothermTab:
                self.GeothermTab.updateGeotherm()
            store_input(self, {conf.Variable_names['Heat_production_flag']: '0'})
            
    #----------------------------------------------------------
    def set_Faults(self):
        """ This function call the class SetFaultGeometry and keep track of it."""
        
        # if setTopography is already opened, then make sure the tab is visible
        if not self.setFault: #If not already opened
            self.setFault = ShowFaultGeometry(self)
            for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
                if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Fault geometry":
                    self.MainWin.StackedTabs[self.MainWin.CurrentTabs].setTabVisible(i,True)
            
        else: #set tab visible
            for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
                if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Fault geometry":
                    self.setFault.update()
                    self.setFault.updateFaults()
                    
                    self.MainWin.StackedTabs[self.MainWin.CurrentTabs].setTabVisible(i,True)

    #----------------------------------------------------------
    def set_Geotherm(self):
        """ This function call the class Geotherm and keep track of it."""
        
        # if Geotherm is already opened, then make sure the tab is visible
        if not self.GeothermTab: #If not already opened
            self.GeothermTab = Geotherm(self,self.MainWin)
            
        #set tab visible
        for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
            if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Geotherm":
                self.MainWin.StackedTabs[self.MainWin.CurrentTabs].setTabVisible(i,True)
                
    #----------------------------------------------------------
    def set_Headward_Erosion(self):
        """ Call *** class to get the input parameter for headward propagation of erosion
        from the user"""
        
        if self.DoHeadward.isChecked() == True:
            self.windowHeadwardProp = U_interface.WinExtraParameters(200, 100)
            self.windowHeadwardProp.OkButton.clicked.connect(lambda: self.windowHeadwardProp.close())
            self.HeadwardPropParam = Headward_erosion_param(self)
            self.windowHeadwardProp.layout.addLayout(self.HeadwardPropParam.grid, 0, 0, 1, 3)
            Q2 = QWidget()
            Q2.setLayout(self.windowHeadwardProp.layout)
            self.windowHeadwardProp.setCentralWidget(Q2)
            conf.WindowsOpen.append(self.windowHeadwardProp)
            self.windowHeadwardProp.show()
            store_input(self, {conf.Variable_names['Headward erosion flag']: '1'})
        else:
            store_input(self, {conf.Variable_names['Headward erosion flag']: '0'})
        
    #----------------------------------------------------------
    def set_Topography(self):
        """ This function call the class setTopography and keep track of it."""
        
        # if setTopography is already opened, then make sure the tab is visible
        try:
            if not self.setTopo: #If not already opened
                self.setTopo = SetTopography(self)
                for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
                    if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Topography":
                        self.MainWin.StackedTabs[self.MainWin.CurrentTabs].setTabVisible(i,True)
            else: 
                # set tab visible
                for i in range(self.MainWin.StackedTabs[self.MainWin.CurrentTabs].count()):
                    if self.MainWin.StackedTabs[self.MainWin.CurrentTabs].tabText(i) == "Topography":
                        self.MainWin.StackedTabs[self.MainWin.CurrentTabs].setTabVisible(i,True)
        except Exception as E:
             QErrorMessage(self).showMessage('Something wrong happen when trying to update topography:' + str(E))
             return
        
    #----------------------------------------------------------
    def showModel(self,button):
        """ To show the surface model with sample location(s).
        Show in subplot predicted height agains obs height"""
        
        # Which button has been clicked
        if button.text() == 'Check sample locations':
            WinName = 'Sample locations'
        else:
            WinName = 'Sample elevations'
        
        self.ModelWindow = U_interface.WinExtraParameters(text=WinName,width=600, height=600)
        self.ModelWindow.OkButton.clicked.connect(lambda: self.ModelWindow.close())
        try:
            nskip = int(self.nskipEdit8.text())
            nx = int(self.nxEdit2.text())
            ny = int(self.nyEdit3.text())
            dx = float(self.dLonEdit6.text())
            dy = float(self.dLatEdit7.text())
            y0 = float(self.lat0Edit5.text())
            x0 = float(self.lon0Edit4.text())
        except ValueError:
            print("Wrong parameter value when trying to update topography for sample location.")
            return
        if self.topo_file_nameEdit1.text() == 'SPM/':
            # Get topo files
            filenames = os.listdir(os.path.join(self.PFolder,"data","SPM"))
            # Get last topo file
            topoFiles = [file for file in filenames if file.startswith("topo")][-1]
            # Get elevation from iSOSIA file
            self.topoFile = os.path.join(self.PFolder,"data","SPM",topoFiles)
            try:
                elevation = np.loadtxt(self.topoFile)
            except OSError:
                QErrorMessage().showMessage("File not found: "+self.topoFile)
                return
            elevation = elevation.astype(float)
            Z = np.flip(np.reshape(elevation,(int(self.nyEdit3.text()),int(self.nxEdit2.text()))),1)
            label = ["Latitude (km)", "Longitude (km)"]
        elif int(self.Param.DParameters['Raster_flag'])or os.path.exists(os.path.join(self.PFolder,"data",self.topo_file_nameEdit1.text()+'.tif')):
            path = os.path.join(self.PFolder,"data",self.topo_file_nameEdit1.text()+'.tif')
            try:
                elevation,nx,ny,dx,dy,x0,y0 = GIS.load_raster_file(path)
            except AttributeError:
                QErrorMessage().showMessage("The file"+self.topo_file_nameEdit1.text()+'.tif has not been found. Please, provide it.')
                return
            
            Z = np.reshape(elevation,(ny,nx))
            Z = np.fliplr(Z)
            label = ["Latitude (°)", "Longitude (°)"]
        elif self.topo_file_nameEdit1.text()[-1] != '/' and self.topo_file_nameEdit1.text() != 'Nil':
            self.topoFile = os.path.join(self.PFolder,"data",self.topo_file_nameEdit1.text())
            elevation = np.loadtxt(self.topoFile)
            elevation = elevation.astype(float)
            Z = np.reshape(elevation,(ny,nx))
            label = ["Latitude (°)", "Longitude (°)"]
        # elif topo30
        else:#means flat topography
            Z = np.array(np.zeros((ny,nx)))
            label = ["Latitude (°)", "Longitude (°)"]
        Zinit = np.flip(Z,1)
        # Handle nskip
        nx0 = nx
        ny0 = ny
        nx = int((nx-1)/nskip+1)
        ny = int((ny-1)/nskip+1)
        Z = Zinit[::nskip,::nskip]
        xl = (nx-1)*dx*nskip
        yl = (ny-1)*dy*nskip
        x = np.linspace(x0,x0+xl,nx)
        y = np.linspace(y0,y0+yl,ny)
        X,Y = np.meshgrid(x,y)
        # Get sample location(s)
        num_rows, num_cols = self.SampleTable.rowCount(), self.SampleTable.columnCount()
        lat_grains = np.array(np.zeros(int(self.nbSampleValue.text())))
        lon_grains = np.array(np.zeros(int(self.nbSampleValue.text())))
        Elev_grains = np.array(np.zeros(int(self.nbSampleValue.text())))
        for row in range(num_rows):
            for col in range(num_cols):
                item = self.SampleTable.item(row, col)
                if item is None:
                    break
                if col == 1:
                    lon_grains[row] = float(item.text())
                elif col == 2:
                    lat_grains[row] = float(item.text())
                elif col == 3:
                    Elev_grains[row] = float(item.text())
            if item is None:
                break  
        
        self.plotSpace =  pgu.PlotCanvas()
        if WinName == 'Sample locations':
            #First check if data are in DEM
            latg2 = []
            long2 = []
            xl = x0+(nx0-1)*dx 
            yl = y0+(ny0-1)*dy
            Nb_discard = 0
            for i in range(len(lon_grains)):
                xg = lon_grains[i]
                yg = lat_grains[i]
                if (xg-x0)*(xg-xl) > 0 or (yg-y0)*(yg-yl) > 0:
                    latg2.append(np.nan)
                    long2.append(np.nan)
                    Nb_discard += 1
                else:
                    latg2.append(lat_grains[i])
                    long2.append(lon_grains[i])
                
            try:
                self.plotSpace.axes.pcolormesh(X,Y,Z)
            except TypeError:
                QErrorMessage().showMessage("The dimensions of Z are incompatible with X and/or Y.")
                print('Z dimensions:', np.shape(Z))
                print('X dimensions:', np.shape(X))
                print('Y dimensions:', np.shape(Y))
                return
            self.plotSpace.axes.plot(long2,latg2,'w.')
            self.plotSpace.axes.set_aspect('equal')
            self.plotSpace.axes.set_xlabel(label[1])
            self.plotSpace.axes.set_ylabel(label[0])
            if Nb_discard:
                self.plotSpace.fig.suptitle("\nNumber of data out of DEM = "+str(Nb_discard) +"/"+str(len(latg2)))
        else:
            #Find the 4th closest neighbors for each data coordinates
            Zr = np.reshape(Z,(1,nx*ny))
            Elevations_predicted = []
            xl = x0+(nx0-1)*dx 
            yl = y0+(ny0-1)*dy
            try:
                Nb_discard = 0
                for i in range(len(lon_grains)):
                    xg = lon_grains[i]
                    yg = lat_grains[i]
                    #First check if data are in DEM
                    if (xg-x0)*(xg-xl) > 0 or (yg-y0)*(yg-yl) > 0:
                        Elevations_predicted.append(np.nan)
                        Nb_discard += 1
                    else:
                        # Get euclidean distance
                        xtemp = np.reshape(X,(1,nx*ny))
                        ytemp = np.reshape(Y,(1,nx*ny))
                        # Get nxx
                        nxx = int((xg-x0)/(dx*nskip))
                        # Get nyy
                        nyy = int((yg-y0)/(dy*nskip))
                        # indexes of the 4 neighboors, left bot,left top, right top, right bot
                        neighboors_nx = np.asarray([nxx,nxx,nxx+1,nxx+1],dtype=int)
                        neighboors_ny = np.asarray([nyy,nyy+1,nyy+1,nyy],dtype=int)
                        # Coordinates
                        xt = X[neighboors_ny[:],neighboors_nx[:]]
                        yt = Y[neighboors_ny[:],neighboors_nx[:]]
                        zt = Z[neighboors_ny[:],neighboors_nx[:]]
                        # Bilinear interpolation
                        xg = x0 + abs(xl - x0) * ((xg - x0) / (xl - x0))
                        yg = y0 + abs(yl - y0) * ((yg - y0) / (yl - y0))
                        i1=int((xg-x0)/(dx*nskip))
                        if (i1 == nx):
                            i1=nx-1
                        j1=int((yg-y0)/(dy*nskip))
                        if (j1 == ny):
                            j1=ny-1
                        r=(xg-(i1)*dx*nskip-x0)/(dx*nskip) 
                        r=-1.+2.*r
                        s=(yg-(j1)*dy*nskip-y0)/(dy*nskip) 
                        s=-1.+2.*s
                        wobs1=(1.-r)*(1.-s)/4.
                        wobs2=(1.+r)*(1.-s)/4.
                        wobs3=(1.+r)*(1.+s)/4.
                        wobs4=(1.-r)*(1.+s)/4.
                        Elevations_predicted.append(wobs1*zt[0]+wobs4*zt[1]+wobs3*zt[2]+wobs2*zt[3])
                    
                # Plot 1:1 line
                self.plotSpace.axes.plot([0,8000], [0,8000], marker = 'None', linestyle = '-', 
                             label = '1:1 line', color = 'lightgrey', alpha = 1)
                # Plot data obs vs pred
                self.plotSpace.axes.plot(Elev_grains,Elevations_predicted,marker = 'o', linestyle = 'None',color='g')
                self.plotSpace.axes.set_xlim(left = min(min(Elev_grains), min(Elevations_predicted)), 
                                             right = max(max(Elev_grains), max(Elevations_predicted))+20)
                self.plotSpace.axes.set_ylim(bottom = min(min(Elev_grains), min(Elevations_predicted)), 
                                             top = max(max(Elev_grains), max(Elevations_predicted)))
                self.plotSpace.axes.set_aspect('equal')
                self.plotSpace.axes.set_xlabel(u'Observed elevation (m)')
                self.plotSpace.axes.set_ylabel(u'Predicted elevation (m)')
                # Compute mean difference in elevation
                MeanDiffElev = np.mean(abs(Elevations_predicted-Elev_grains))
                StdDiffElev = np.std(abs(Elevations_predicted-Elev_grains))
                if Nb_discard:
                    self.plotSpace.fig.suptitle("nskip = " + str(nskip) +
                                                "\nNumber of data out of DEM = "+str(Nb_discard) +"/"+str(len(Elevations_predicted)))
                else: 
                    self.plotSpace.fig.suptitle("nskip = " + str(nskip) +
                                                "\nAverage elevation difference = "+str(round(MeanDiffElev,2)) +"+/-"+str(round(StdDiffElev,2)))
                
            except:
                raise
                QErrorMessage(self).showMessage('Issue when trying to plot obs vs pred elevations.')
                return
        # Add toolbar
        self.toolbar = NavigationToolbar(self.plotSpace,self)
        Vbox = QVBoxLayout()
        Vbox.addWidget(self.toolbar)
        Vbox.addWidget(self.plotSpace)
        self.ModelWindow.layout.addLayout(Vbox,0,0,1,3)
        Q = QWidget()
        Q.setLayout(self.ModelWindow.layout)
        self.ModelWindow.setCentralWidget(Q)
        conf.WindowsOpen.append(self.ModelWindow)
        self.ModelWindow.show()
    
    #---------------------------------------------------------
    def store_HeatProd(self):
        """Store Heat production value to write in input file.
        if e-folding is selected then heat production is negative."""
        try:
            if self.EfoldingHeatProd.isChecked() == True:
                store_input(
                    self, {'heat_production': str(float(self.HProdEdit7.text()))})
            else:
                store_input(
                    self, {'heat_production': str(float(self.HProdEdit7.text()))})
        except ValueError:
            print("Error in value for heat production rate !! Assume range is provided.")
            store_input(
                self, {'heat_production': str(self.HProdEdit7.text())})    
            
    #----------------------------------------------------------
    def store_input_CheckBox_savePTT(self, d, s):
        """ 
          Update the value of the check box each time it is changed
          d is the check box
          s is the parameter name
        """
        if d.isChecked() == True:
            store_input(self, {s: '1'})
        else:
            store_input(self, {s: '0'})
     
    #----------------------------------------------------------
    def open_ages_tab(self):
        """
          Open the age tab when click on button.

        """
        self.set_AgesParam()
        self.AgesParam.updateComboBox()
        self.AgesParam.updateWidgets()
        self.AgesParam.updateTables()
        
    #----------------------------------------------------------
    def store_input_CheckBox(self, d, s):
        """ 
          Update the value of the check box each time it is changed
          d is the check box
          s is the parameter name
        """
        if d.isChecked() == True:
            store_input(self, {s: '1'})
        else:
            store_input(self, {s: '0'})
            # if self.AgeCombo.currentIndex() == 2:

    #----------------------------------------------------------
    def store_input_ComboBoxEcho(self, d, s):
        """ Update the value of the Combo box each time it is changed."""
        if d.currentText() == 'no echo':
            flag = 0
        elif d.currentText() == 'with parameter and value':
            flag = 1
        elif d.currentText() == 'with default value':
            flag = 2
        elif d.currentText() == 'with short description':
            flag = 3
        elif d.currentText() == 'with all parameters':
            flag = 4
        store_input(self, {s: str(flag)})

    #----------------------------------------------------------
    def store_input_ComboBoxAge(self, d, s):
        """ Update the value of the Combo box each time it is changed."""
        if d.currentText() == 'none':
            flag = 0
            self.widget1.setVisible(False)
            self.widget2.setVisible(False)
            #First uncheck all check boxes
            self.ageESRb.setChecked(False)
            self.ageTLb.setChecked(False)
            self.ageAHeb.setChecked(False)
            self.ageOSLb.setChecked(False)
            self.ageZHeb.setChecked(False)
            self.ageAFTb.setChecked(False)
            self.ageZFTb.setChecked(False)
            self.ageKArb.setChecked(False)
            self.ageBArb.setChecked(False)
            self.ageMArb.setChecked(False)
            self.ageHArb.setChecked(False)
            #Then disable them
            self.ageESRb.setEnabled(False)
            self.ageTLb.setEnabled(False)
            self.ageAHeb.setEnabled(False)
            self.ageOSLb.setEnabled(False)
            self.ageZHeb.setEnabled(False)
            self.ageAFTb.setEnabled(False)
            self.ageZFTb.setEnabled(False)
            self.ageKArb.setEnabled(False)
            self.ageBArb.setEnabled(False)
            self.ageMArb.setEnabled(False)
            self.ageHArb.setEnabled(False)
            
        elif d.currentText() == 'for all nodes':
            flag = 1
            self.widget1.setVisible(True)
            self.widget2.setVisible(False)
            self.ageESRb.setEnabled(True)
            self.ageTLb.setEnabled(False) #↕ Thermoluminescence disabled
            self.ageAHeb.setEnabled(True)
            self.ageOSLb.setEnabled(True)
            self.ageZHeb.setEnabled(True)
            self.ageAFTb.setEnabled(True)
            self.ageZFTb.setEnabled(True)
            self.ageKArb.setEnabled(True)
            self.ageBArb.setEnabled(True)
            self.ageMArb.setEnabled(True)
            self.ageHArb.setEnabled(True)
            
        elif d.currentText() == 'sample specific':
            flag = 2
            self.widget1.setVisible(False)
            self.widget2.setVisible(True)
            self.ageESRb.setEnabled(True)
            self.ageTLb.setEnabled(False)
            self.ageAHeb.setEnabled(True)
            self.ageOSLb.setEnabled(True)
            self.ageZHeb.setEnabled(True)
            self.ageAFTb.setEnabled(True)
            self.ageZFTb.setEnabled(True)
            self.ageKArb.setEnabled(True)
            self.ageBArb.setEnabled(True)
            self.ageMArb.setEnabled(True)
            self.ageHArb.setEnabled(True)
            self.savePTTb.setChecked(True)
        if self.AgesParam:
            self.AgesParam.updateWidgets()

        store_input(self, {s: str(flag)})
        
    # #----------------------------------------------------------
    # def store_input_ComboBoxFT(self, d, s):
    #     """ Update the value of the Combo box each time it is changed."""
    #     if d.currentText() == 'van der Beek (1995)':
    #         flag = 0
    #     elif d.currentText() == 'Ketcham et al. (1999)':
    #         flag = 1
    #     store_input(self, {s: str(flag)})
    
    #----------------------------------------------------------
    def open_AgeWeight(self):
        """ When the user click on the Age weight check box, this function is called.
        It open a window with the weight parameters on each thermochronometer."""
        
        self.windowWeightAge = U_interface.WinExtraParameters(200, 100)
        self.windowWeightAge.OkButton.clicked.connect(lambda: self.windowWeightAge.close())
        self.AgeWeightParam = Weight_on_ages(self)
        self.windowWeightAge.layout.addLayout(self.AgeWeightParam.grid, 0, 0, 1, 3)
        Q2 = QWidget()
        Q2.setLayout(self.windowWeightAge.layout)
        self.windowWeightAge.setCentralWidget(Q2)
        conf.WindowsOpen.append(self.windowWeightAge)
        self.windowWeightAge.show()
    
    #----------------------------------------------------------
    def update_InvMode(self):
        """ Update the inversion mode. In batch mode set the parameter for inversion
        to 1. """
        
        # Store parameter value
        store_input(self, {conf.Variable_names["use inversion or batch"]: str(self.InvModeCombo.currentIndex())})
        
        # Set inversion parameters
        if self.InvModeCombo.currentIndex() == 2: # Batch mode
            self.sampleSizeOtherEdit3.setText('1')
            self.sampleSizeFirstEdit2.setText(str(int(self.nIntervalsEdit.text())+1))
            self.maxNItEdit1.setText('1')
            self.nbCellResampleEdit4.setText('1')
            
            self.sampleSizeOtherEdit3.setEnabled(False)
            self.sampleSizeFirstEdit2.setEnabled(False)
            self.maxNItEdit1.setEnabled(False)
            self.nbCellResampleEdit4.setEnabled(False)
            self.nIntervalsEdit.setEnabled(True)
        
        else: # Inversion mode
            self.sampleSizeOtherEdit3.setEnabled(True)
            self.sampleSizeFirstEdit2.setEnabled(True)
            self.maxNItEdit1.setEnabled(True)
            self.nbCellResampleEdit4.setEnabled(True)
            self.nIntervalsEdit.setEnabled(False)
        
    #----------------------------------------------------------    
    def update_filename(self):
        """ update some parameters according to the topography file name."""
        
        # First, store the new topography file name
        store_input(
            self, {'topo_file_name': str(self.topo_file_nameEdit1.text())})
        
        # If use synthetic sinusoidal topo then update time topo table
        if self.topo_file_nameEdit1.text() ==  conf.Variable_names['Synthetic topo file name']:
            data = {'time topo (Ma)': ['1']*int(self.Param.DParameters['ntime']),
                    "amplification": ['1']*int(self.Param.DParameters['ntime']),
                    "offset (km)": ["0"]*int(self.Param.DParameters['ntime']),
                    "output": ["0"]*int(self.Param.DParameters['ntime']),
                    "phase shift (km)": ["0"]*int(self.Param.DParameters['ntime'])}
            self.isSynthTopo = 1
            # self.topoTable = U_interface.WinTable(data, int(self.Param.DParameters['ntime']), 5)
            self.topoTable.insertColumn(4)
            self.topoTable.setHorizontalHeaderLabels(['time topo (Ma)',"amplification","offset (km)","output","phase shift (km)"])
        else: 
            data = {'time topo (Ma)': ['1']*int(self.Param.DParameters['ntime']),
                    "amplification": ['1']*int(self.Param.DParameters['ntime']),
                    "offset (km)": ["0"]*int(self.Param.DParameters['ntime']),
                    "output": ["0"]*int(self.Param.DParameters['ntime'])}
            self.topoTable.removeColumn(4)
            self.topoTable.setHorizontalHeaderLabels(['time topo (Ma)',"amplification","offset (km)","output"])
            # self.topoTable = U_interface.WinTable(data, int(self.Param.DParameters['ntime']), 4)
            self.isSynthTopo = 0
        # update table
        self.updateTable(self.topoTable, int(self.Param.DParameters['ntime']), 1, self.Param)
            
    #----------------------------------------------------------    
    def update_ninterval(self):
        """ Update number of intervals for the parameter to batch. We use the Pecube
        parameter sample_size_for_first_iteration to set the subdivision of the 
        parameter range values."""
        
        # Set the value for the sample_size_for_first_iteration
        self.sampleSizeFirstEdit2.setText(str(int(self.nIntervalsEdit.text())+1))
        
    #----------------------------------------------------------
    def updateKinTable(self, d, f, g, h, p):
        """ 
          To store input value for the nstep fault table
          d is the table, f is the number of fault, g is the number of nstep
          h = 0 means we just want to set the number of rows of the table
          h = 1 means we want to store the value of the table
          h = 2 means we load an old Pecube input file
        """
        uplift =''
        #Check if we do bloc uplift or faulting
        if self.UpliftButton1.isChecked():
            #No uplift, erase table
            d.clearContents()
            uplift = 'Bloc'
            h = 3
        elif self.UpliftButton2.isChecked():
            uplift = "Bloc"
        elif self.UpliftButton3.isChecked():
            uplift = "Fault"
        
        if h == 0:
            d.clearContents()
            d.setRowCount(g*f)
            count = 0
                
            for k in range(f):
                for j in range(g):
                    d.setItem(j+count, 0, QTableWidgetItem(uplift+str(k+1)))
                count += g
            for k in range(f):
                store_input(self, {'nstep'+str(k+1): str(g)})
        elif h == 1:
            num_rows, num_cols = d.rowCount(), d.columnCount()
            nrow = 0
            for row in range(num_rows):
                if nrow > g-1:
                    nrow = 1
                else:
                    nrow += 1
                for col in range(num_cols):
                    item = d.item(row, col)
                    if item is None:
                        break
                    if col == 1:
                        store_input(
                            self, {'time_start'+str(floor((row+1e-3)/g)+1)+'_'+str(nrow): str(item.text())})
                    elif col == 2:
                        store_input(
                            self, {'time_end'+str(floor((row+1e-3)/g)+1)+'_'+str(nrow): str(item.text())})
                    elif col == 3:
                        store_input(
                            self, {'velo'+str(floor((row+1e-3)/g)+1)+'_'+str(nrow): str(item.text())})
                if item is None:
                    break
        elif h == 2:
            d.clearContents()
            d.setRowCount(g*f)
            count = 0
            for j in range(f):
                for i in range(g):
                    d.setItem(i+count, 0, QTableWidgetItem(uplift+str(j+1)))
                    d.setItem(
                        i+count, 1, QTableWidgetItem(p.TableParameters['time_start'+str(j+1)+'_'+str(i+1)]))
                    d.setItem(
                        i+count, 2, QTableWidgetItem(p.TableParameters['time_end'+str(j+1)+'_'+str(i+1)]))
                    d.setItem(
                        i+count, 3, QTableWidgetItem(p.TableParameters['velo'+str(j+1)+'_'+str(i+1)]))
                count += g
                    
    #---------------------------------------------------------- 
    def updateRefElev(self):
        """ To update the reference elevation from which to compute topographic evolution. """
        store_input(self,{'topo_ref':str(self.RefElevCombo.currentIndex())})
        if self.RefElevCombo.currentIndex() == 4 or self.RefElevCombo.currentIndex() == 5: # smooth topography
            # self.SmoothingFactor.setEnabled(True)
            self.RefCustom.setEnabled(True)
        else:
            self.RefCustom.setEnabled(False)
    
    #---------------------------------------------------------- 
    def updateSmoothingFactor(self):
        """ To update the reference elevation from which to compute topographic evolution. """
        if self.setTopo:
            factor = float(self.SmoothingFactor.text())
            self.setTopo.updateTopo()
       
    #----------------------------------------------------------
    def updateTable(self, d, f, h, p):
        """ To store input value for the time table. """
        
        d.setRowCount(f+1)
        num_rows, num_cols = d.rowCount(), d.columnCount()
        if h:
            dParam = conf.DefaultParameterValues()
            for j in range(f):
                # from old input file
                try:
                    d.setItem(j, 0, QTableWidgetItem(
                            p.TableParameters['time_topo'+str(j+1)]))
                except KeyError:
                    # The user loaded an old input file and update the table
                    d.setItem(j, 0, QTableWidgetItem(
                        dParam.DParameters['time_topo1']))
    
                # from old input file
                try:
                    d.setItem(j, 1, QTableWidgetItem(
                        p.TableParameters['amplification'+str(j+1)]))
                except KeyError:
                    # The user loaded an old input file and update the table
                    d.setItem(j, 1, QTableWidgetItem(
                        dParam.DParameters['amplification1']))
                    
                # from old input file
                try:
                    d.setItem(j, 2, QTableWidgetItem(
                        p.TableParameters['offset'+str(j+1)]))
                except KeyError:
                    # The user loaded an old input file and update the table
                    #or is a new input file
                    d.setItem(j, 2, QTableWidgetItem(
                        dParam.DParameters['offset1']))
                try:
                    # from old input file
                    d.setItem(j, 3, QTableWidgetItem(
                        p.TableParameters['output'+str(j+1)]))
                except KeyError:
                    # The user loaded an old input file and update the table
                    d.setItem(j, 3, QTableWidgetItem(
                        dParam.DParameters['output1']))
                if self.isSynthTopo: #is a synthetic sinusoidal topo, then read phase shift
                    try:
                        # from old input file
                        d.setItem(j, 4, QTableWidgetItem(
                            p.TableParameters['phase'+str(j+1)]))
                    except KeyError:
                        # The user loaded an old input file and update the table
                        d.setItem(j, 4, QTableWidgetItem(
                            dParam.DParameters['phase1']))
                    
    
            # The last row is alway at time 0 Ma        
            d.setItem(f,0,QTableWidgetItem('0'))
            try:
                # from old input file
                d.setItem(j+1, 1, QTableWidgetItem(
                    p.TableParameters['amplification'+str(j+2)]))
            except KeyError:
                # The user loaded an old input file and update the table
                d.setItem(j+1, 1, QTableWidgetItem(
                    dParam.DParameters['amplification1']))
            try:
                d.setItem(j+1, 2, QTableWidgetItem(
                     p.TableParameters['offset'+str(j+2)]))
            except KeyError:
                # The user loaded an old input file and update the table
                # or is a new input file
                d.setItem(j+1, 2, QTableWidgetItem(
                     dParam.DParameters['offset1']))   
            try:
                # from old input file
                d.setItem(j+1, 3, QTableWidgetItem(
                     p.TableParameters['output'+str(j+2)]))
            except KeyError:
                # The user loaded an old input file and update the table
                d.setItem(j+1, 3, QTableWidgetItem('1'))
                
            if self.isSynthTopo: #is a synthetic sinusoidal topo, then read phase shift
                try:
                    # from old input file
                    d.setItem(j+1, 4, QTableWidgetItem(
                         p.TableParameters['phase'+str(j+2)]))
                except KeyError:
                    # The user loaded an old input file and update the table
                    d.setItem(j+1, 4, QTableWidgetItem(dParam.DParameters['phase1']))
                    
        self.timeTable_dict = {}
        # Read cells from the table
        for row in range(num_rows):
            for col in range(num_cols):
                item = d.item(row, col)
                if item is None:
                    break
                if col == 0:
                    self.timeTable_dict['time_topo'+str(row+1)] = str(item.text())
                    p.TableParameters['time_topo'+str(row+1)] = str(item.text())
                elif col == 1:
                    self.timeTable_dict['amplification'+str(row+1)] = str(item.text())
                    p.TableParameters['amplification'+str(row+1)] = str(item.text())
                elif col == 2:
                    self.timeTable_dict['offset'+str(row+1)] = str(item.text())
                    p.TableParameters['offset'+str(row+1)] = str(item.text())
                elif col == 3:
                    cond = int(item.text()) == 0 or int(item.text()) == 1
                    if not cond:
                        QErrorMessage().showMessage("Wrong value for 'output'. Accepted values are '0' or '1'.")
                    self.timeTable_dict['output'+str(row+1)] = str(item.text())
                    p.TableParameters['output'+str(row+1)] = str(item.text())
                elif col == 4:
                    self.timeTable_dict['phase'+str(row+1)] = str(item.text())
                    p.TableParameters['phase'+str(row+1)] = str(item.text())
                
            if item is None:
                break
        store_input(self, self.timeTable_dict)
        
    #----------------------------------------------------------
    def updateTableFault(self, d, f, g, h, p):
        """ 
          To store input value for the fault table
          d is the table, f is the number of fault, g is the number of npoint
          h = 0 means we just want to set the number of rows of the table
          h =1 means we want to store the value of the empty table
          h =2 means we load an old Pecube input file
        
        """
        if str(g) == '-1':
            store_input(self, {'npoint1': str(g)})
            d.clearContents()
            d.setRowCount(0)
        elif int(g) >= 0:
            g = int(g)
            if h == 0:
                d.clearContents()
                d.setRowCount(g*f)
                count = 0
                for k in range(f):
                    for j in range(g):
                        d.setItem(
                            j+count, 0, QTableWidgetItem("Fault "+str(k+1)))
                    count += g
                for k in range(f):
                    store_input(self, {'npoint'+str(k+1): str(g)})

            elif h == 1:
                num_rows, num_cols = d.rowCount(), d.columnCount()
                nrow = 0
                for row in range(num_rows):
                    if nrow > g-1:
                        nrow = 1
                    else:
                        nrow += 1
                    try:
                        for col in range(num_cols):
                            item = d.item(row, col)
                            if item is None:
                                break
                            if col == 1:
                                float(item.text())
                                store_input(
                                    self, {'r'+str(floor((row+1e-3)/g)+1)+'_'+str(nrow): str(item.text())})
                            elif col == 2:
                                float(item.text())
                                store_input(
                                    self, {'s'+str(floor((row+1e-3)/g)+1)+'_'+str(nrow): str(item.text())})
                        if item is None:
                            break
                    except ValueError:
                        QErrorMessage().showMessage("In Fault table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
                        return
            elif h == 2:
                d.clearContents()
                d.setRowCount(g*f)
                count = 0
                for j in range(f):
                    for i in range(g):
                        d.setItem(
                            i+count, 0, QTableWidgetItem('Fault'+str(j+1)))
                        d.setItem(
                            i+count, 1, QTableWidgetItem(p.TableParameters['r'+str(j+1)+'_'+str(i+1)]))
                        d.setItem(
                            i+count, 2, QTableWidgetItem(p.TableParameters['s'+str(j+1)+'_'+str(i+1)]))
                    count += g
        else:
            QErrorMessage().showMessage("Authorized values for npoint '-1' or values >= 0")

    #----------------------------------------------------------
    def updateTableSample(self, f, h, p, c):
        """ 
          To store input value for the Grain table
          f : nb sample
          h : parameter
          c: if cell changed
          p : parameters
        """  

        if h:
            try:
                self.SampleTable.cellChanged.disconnect() #Shut down update from the user
            except TypeError:
                pass
            self.SampleTable.setRowCount(f)
            # dParam = conf.DefaultParameterValues()
            count = 0
            for j in range(f):
                if "SampleName"+str(j+1) in p.TableParameters:
                    sampleName = p.TableParameters['SampleName'+str(j+1)]
                    self.SampleTable.setItem(j, 0, QTableWidgetItem(sampleName))
                else:
                    sampleName = p.DParameters['SampleName']+str(j+1)
                    self.SampleTable.setItem(j, 0, QTableWidgetItem(sampleName))
                if sampleName+"_lon" in p.TableParameters:
                    self.SampleTable.setItem(j, 1, QTableWidgetItem(
                        str(p.TableParameters[sampleName+'_lon'])))
                else:
                    p.DParameters[sampleName+'_lon'] = '0.0'
                    self.SampleTable.setItem(j, 1, QTableWidgetItem(
                        p.DParameters[sampleName+'_lon']))
                if sampleName+"_lat" in p.TableParameters:
                    self.SampleTable.setItem(j, 2, QTableWidgetItem(
                        str(p.TableParameters[sampleName+'_lat'])))
                else:
                    p.DParameters[sampleName+'_lat'] = '0.0'
                    self.SampleTable.setItem(j, 2, QTableWidgetItem(
                        p.DParameters[sampleName+'_lat']))
                if sampleName+"_Elev" in p.TableParameters:
                    self.SampleTable.setItem(j, 3, QTableWidgetItem(
                        str(p.TableParameters[sampleName+'_Elev'])))
                else:
                    p.DParameters[sampleName+'_Elev'] = '0.0'
                    self.SampleTable.setItem(j,3, QTableWidgetItem(
                        p.DParameters[sampleName+'_Elev']))
                if "TL_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith("TL_"+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(str, temp.values()))))
                        self.SampleTable.setItem(j, 4, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 4, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 4, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "OSL_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('OSL_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(str, temp.values()))))
                        self.SampleTable.setItem(j, 5, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 5, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 5, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                    
                if "ESR_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('ESR_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(str, temp.values()))))
                        self.SampleTable.setItem(j, 6, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 6, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 6, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                    
                if "AHE_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('AHE_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 7, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 7, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 7, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "ZHE_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('ZHE_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 8, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 8, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 8, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "AFT_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('AFT_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 9, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 9, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 9, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "ZFT_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('ZFT_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 10, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 10, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 10, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "KAR_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('KAR_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 11, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 11, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 11, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "BAR_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('BAR_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 12, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 12, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 12, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "MAR_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('MAR_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 13, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 13, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 13, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                if "HAR_"+sampleName+"_1" in p.TableParameters: # at least one grain, start counting how much
                    try:
                        # Extract information from old input storage (dictionary)
                        temp = {k:v for k,v in p.TableParameters.items() if k.startswith('HAR_'+sampleName+'_')}
                        # remove zero value (mean no data)
                        countg = np.count_nonzero(np.array(list(map(float, temp.values()))))
                        self.SampleTable.setItem(j, 14, QTableWidgetItem(
                        str(countg)))
                    except IndexError:
                        self.SampleTable.setItem(j, 14, QTableWidgetItem(
                            p.DParameters['ngrains1']))
                else:
                    self.SampleTable.setItem(j, 14, QTableWidgetItem(
                        p.DParameters['ngrains1']))
                
        num_rows, num_cols = self.SampleTable.rowCount(), self.SampleTable.columnCount()
        temp_dict = {}
        temp_dict['TL_TotGrain'] = temp_dict['OSL_TotGrain'] = temp_dict['AHe_TotGrain']  = 0
        temp_dict['ZHe_TotGrain']  = temp_dict['AFT_TotGrain']  = temp_dict['ZFT_TotGrain'] = 0
        temp_dict['KAr_TotGrain']  = temp_dict['BAr_TotGrain']  = temp_dict['MAr_TotGrain']  = 0
        temp_dict['HAr_TotGrain']  = temp_dict['ESR_TotGrain'] = 0
        try:
            count = 0
            countgrains = 0
            sampleName_prev = 'blank'
            for row in range(num_rows):
                maxnbgrain = []
                for col in range(num_cols):
                    item = self.SampleTable.item(row, col)
                    if item is None:
                        break
                    if col == 0:
                        sampleName = str(item.text())
                        if not sampleName == sampleName_prev:
                            count = 0
                        count += 1
                        temp_dict['SampleName'+str(row+1)] = str(item.text())
                        sampleName_temp = temp_dict['SampleName'+str(row+1)] 
                        if c:
                            p.TableParameters['SampleName'+str(row+1)] = str(item.text())
                    if col == 1:
                        float(item.text())
                        temp_dict[sampleName_temp+'_lon'] = str(item.text())
                        if c:
                            p.TableParameters[sampleName_temp+'_lon'] = str(item.text())
                    elif col == 2:
                        float(item.text())
                        temp_dict[sampleName_temp+'_lat'] = str(item.text())
                        if c:
                            p.TableParameters[sampleName_temp+'_lat'] = str(item.text())
                    elif col == 3:
                        float(item.text())
                        temp_dict[sampleName_temp+'_Elev'] = str(item.text())
                        if c:
                            p.TableParameters[sampleName_temp+'_Elev'] = str(item.text())
                    elif col == 4:
                        int(item.text())
                        temp_dict['TL_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['TL_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['TL_ngrains'+str(row+1)] = str(item.text())
                    elif col == 5:
                        int(item.text())
                        temp_dict['OSL_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['OSL_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['OSL_ngrains'+str(row+1)] = str(item.text())
                    elif col == 6:
                        int(item.text())
                        temp_dict['ESR_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['ESR_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['ESR_ngrains'+str(row+1)] = str(item.text())
                    elif col == 7:
                        int(item.text())
                        temp_dict['AHe_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['AHe_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['AHe_ngrains'+str(row+1)] = str(item.text())
                    elif col == 8:
                        int(item.text())
                        temp_dict['ZHe_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['ZHe_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['ZHe_ngrains'+str(row+1)] = str(item.text())
                    elif col == 9:
                        int(item.text())
                        temp_dict['AFT_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['AFT_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['AFT_ngrains'+str(row+1)] = str(item.text())
                    elif col == 10:
                        int(item.text())
                        temp_dict['ZFT_TotGrain'] += int(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['ZFT_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['ZFT_ngrains'+str(row+1)] = str(item.text())
                    elif col == 11:
                        int(item.text())
                        temp_dict['KAr_TotGrain'] += float(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['KAr_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['KAr_ngrains'+str(row+1)] = str(item.text())
                    elif col == 12:
                        int(item.text())
                        temp_dict['BAr_TotGrain'] += float(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['BAr_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['BAr_ngrains'+str(row+1)] = str(item.text())
                    elif col == 13:
                        int(item.text())
                        temp_dict['MAr_TotGrain'] += float(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['MAr_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['MAr_ngrains'+str(row+1)] = str(item.text())
                    elif col == 14:
                        int(item.text())
                        temp_dict['HAr_TotGrain'] += float(item.text())
                        maxnbgrain.append(int(item.text()))
                        temp_dict['HAr_ngrains'+str(row+1)] = str(item.text())
                        if c:
                            p.TableParameters['HAr_ngrains'+str(row+1)] = str(item.text())
                        # Record max nb grain for single thermochronometer
                        temp_dict['MaxNbGrain_'+str(row+1)] = max(maxnbgrain)
                if maxnbgrain:
                    countgrains += max(maxnbgrain)
        
                sampleName_prev = sampleName
                if item is None:
                    break
        except ValueError:
            # raise
            print("Value error in update table sample: col = "+str(col) + ", row = " +str(row))
            QErrorMessage().showMessage("In Sample table: \n'"+item.text()+"' is not a float/int number. Please, provide a float/int number.")
            

        # # Check age_XXX flags
        if  int(temp_dict['HAr_TotGrain']) > 0: self.ageHArb.setChecked(True), self.store_input_CheckBox(self.ageHArb, 'age_HAr_flag')
        if  int(temp_dict['MAr_TotGrain']) > 0: self.ageMArb.setChecked(True), self.store_input_CheckBox(self.ageMArb, 'age_MAr_flag')
        if  int(temp_dict['BAr_TotGrain']) > 0: self.ageBArb.setChecked(True), self.store_input_CheckBox(self.ageBArb, 'age_BAr_flag')
        if  int(temp_dict['KAr_TotGrain']) > 0: self.ageKArb.setChecked(True), self.store_input_CheckBox(self.ageKArb, 'age_KAr_flag')
        if  int(temp_dict['ZFT_TotGrain']) > 0: self.ageZFTb.setChecked(True), self.store_input_CheckBox(self.ageZFTb, 'age_ZFT_flag')
        if  int(temp_dict['ZHe_TotGrain']) > 0: self.ageZHeb.setChecked(True), self.store_input_CheckBox(self.ageZHeb, 'age_ZHe_flag')
        if  int(temp_dict['AFT_TotGrain']) > 0: self.ageAFTb.setChecked(True), self.store_input_CheckBox(self.ageAFTb, 'age_AFT_flag')
        if  int(temp_dict['AHe_TotGrain']) > 0: self.ageAHeb.setChecked(True), self.store_input_CheckBox(self.ageAHeb, 'age_AHe_flag')
        if  int(temp_dict['TL_TotGrain']) > 0: self.ageTLb.setChecked(True), self.store_input_CheckBox(self.ageTLb, 'age_TL_flag')
        if  int(temp_dict['OSL_TotGrain']) > 0: self.ageOSLb.setChecked(True), self.store_input_CheckBox(self.ageOSLb, 'age_OSL_flag')
        if  int(temp_dict['ESR_TotGrain']) > 0: self.ageESRb.setChecked(True), self.store_input_CheckBox(self.ageESRb, 'age_ESR_flag')
        
        if h:
            self.SampleTable.cellChanged.connect(lambda: self.updateTableSample(
                int(self.nbSampleValue.text()), 0, self.Param, 1))
            
        temp_dict['Total_maxGrain'] = countgrains # = total of samples name with grain id
        self.Samples.update(temp_dict)
        p.DParameters['Nb_Samples'] = self.nbSampleValue.text()
        
       
        
    #---------------------------------------------------------- 
    def updateUplift(self,uplift):
         """
         This function updates what parameters to display to define the
         the tectonic uplift pattern.
         """
         # Combo is the comboBox
         if uplift.text() == "no uplift": #Means no uplift  
             self.bottomLeftEdit9.setEnabled(False)
             self.bottomRightEdit10.setEnabled(False)
             self.upperRightEdit11.setEnabled(False)
             self.upperLeftEdit12.setEnabled(False)
             self.nfaultEdit1.setText('0')
             self.npointiEdit6.setText('-1')
             store_input(self, {'npoint1':'-1'})
             store_input(self, {'nfault':'-1'})
             self.setFaultButton.setEnabled(False)
             self.ShowFault.setEnabled(False)
             self.nstepiEdit13.setEnabled(False)
             self.logVelob.setEnabled(False)
             store_input(self,{'uplift':'0'})
             self.UpliftButton2.setChecked(False)
             self.UpliftButton2.setChecked(False)
             
         elif uplift.text() == "bloc uplift": #Means bloc uplift
             self.bottomLeftEdit9.setEnabled(True)
             self.bottomRightEdit10.setEnabled(True)
             self.upperRightEdit11.setEnabled(True)
             self.upperLeftEdit12.setEnabled(True)
             self.nfaultEdit1.setText('1')
             self.npointiEdit6.setText('-1')
             store_input(self, {'npoint1':'-1'})
             self.setFaultButton.setEnabled(False)
             self.ShowFault.setEnabled(False)
             self.nstepiEdit13.setEnabled(True)
             self.logVelob.setEnabled(True)
             store_input(self,{'uplift':'1'})
             self.UpliftButton1.setChecked(False)
             self.UpliftButton2.setChecked(False)
             
         elif uplift.text() == "faulting": #Means faulting
             self.bottomLeftEdit9.setEnabled(False)
             self.bottomRightEdit10.setEnabled(False)
             self.upperRightEdit11.setEnabled(False)
             self.upperLeftEdit12.setEnabled(False)
             self.setFaultButton.setEnabled(True)
             if self.faultParam and int(self.npointiEdit6.text()) >0:
                 self.ShowFault.setEnabled(True)
             self.nstepiEdit13.setEnabled(True)
             self.logVelob.setEnabled(True)
             store_input(self,{'uplift':'2'})
             self.UpliftButton1.setChecked(False)
             self.UpliftButton2.setChecked(False)
         # Update Kinematic Table
         try:
             self.updateKinTable(self.nstepKinTable, int(self.nfaultEdit1.text()), int(self.nstepiEdit13.text()), 0, self.Param)
         except Exception as E:
             QErrorMessage(self).showMessage('Something wrong happen when trying to update the kinematic table:' + str(E))
             return
         

class MapWindow(QWidget):
    """
      This class enables to load a dem from an interactive map.
      It uses the Folium library to show the interactive map and the
      'elevation' library to load the srtm dem. 
      The coordinates of the rectangle area are used to extract the dem and save
      the raster file (.tif) in the 'data' directory of the Pecube project, then
      this raster can be loaded.
      
    """
    def __init__(self,parent):
        super().__init__()
    
        self.ProjectName = parent.folderName
        conf.ProjectName = parent.folderName
        self.nx = int(parent.nxEdit2.text())
        self.ny = int(parent.nyEdit3.text())
        self.lon0 = float(parent.lon0Edit4.text())
        self.lat0 = float(parent.lat0Edit5.text())
        self.dlon = float(parent.dLonEdit6.text())
        self.dlat = float(parent.dLatEdit7.text())
        
        # Read API key
        if os.path.exists(os.path.join(conf.FolderPath,'OpenTopoKey.txt')):
            with open(os.path.join(conf.FolderPath,'OpenTopoKey.txt')) as f:
                      self.apikey = f.readlines()
        else:
            QErrorMessage(self).showMessage("API key not found in root directory. <br>"+
                                            "Please, first request an API key from OpenTopography.org and follow instruction in the PecubeGUI documentation. <br>"+
                                            "If you have an API key, then make sure to save the key in a text file called OpenTopoKey.txt (in the root directory)")
            return
        
        lower_left = self.lat0
        upper_left = self.lat0 + self.dlat*self.ny
        Easternmost = self.lon0 + self.dlon*self.nx
        Westernmost = self.lon0
        params = Topography.DEFAULT.copy()
        params["south"] = lower_left
        params["north"] = upper_left
        params["west"] = Westernmost
        params["east"] = Easternmost
        # params["dem_type"] = 'SRTMGL1' # 30m
        params["dem_type"] = 'SRTMGL3' # 90m
        params["cache_dir"] = os.path.join(conf.PecubeFolderPath,conf.ProjectName,"data")
        params["api_key"] = self.apikey[0]
        params["output_format"] = "GTiff"
        # Save the dem in a raster file
        try:
            # srtmelev.clip(bounds=tuple(roi), product='SRTM3', output=os.path.join(conf.PecubeFolderPath,conf.ProjectName,"data","topoRaster.tif"))
            # srtmelev.clean()
            # if os.path.exists(os.path.join(conf.PecubeFolderPath,conf.ProjectName,"data")):
            #     shutil.rmtree(os.path.join(conf.PecubeFolderPath,conf.ProjectName,"data"))
            topo = Topography(**params)
            topo.fetch()
        except Exception as inst:
            print("Error in downloading dem: " + str(inst))
                              
            
#############################################################
#############################################################
#############################################################
class Weight_on_ages(QWidget):
    """ 
      This class defines weight for each thermochronometer for the
      calculation of the misfit value (inversion mode).
      parent : all parameters from the WindowsParameters class
     
    """
    def __init__(self, parent):
        
        self.parent = parent
        
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.grid.addWidget(self.parent.WeightAHE,0,1)
        self.grid.addWidget(self.parent.weightAHEedit,1,1)
        self.grid.addWidget(self.parent.WeightAFT,0,2)
        self.grid.addWidget(self.parent.weightAFTedit,1,2)
        self.grid.addWidget(self.parent.WeightZHE,0,3)
        self.grid.addWidget(self.parent.weightZHEedit,1,3)
        self.grid.addWidget(self.parent.WeightZFT,2,1)
        self.grid.addWidget(self.parent.weightZFTedit,3,1)
        self.grid.addWidget(self.parent.WeightBAR,2,2)
        self.grid.addWidget(self.parent.weightBARedit,3,2)
        self.grid.addWidget(self.parent.WeightMAR,2,3)
        self.grid.addWidget(self.parent.weightMARedit,3,3)
        self.grid.addWidget(self.parent.WeightKAR,4,1)
        self.grid.addWidget(self.parent.weightKARedit,5,1)
        self.grid.addWidget(self.parent.WeightHAR,4,2)
        self.grid.addWidget(self.parent.weightHARedit,5,2)
        
#############################################################
#############################################################
#############################################################
class Headward_erosion_param(QWidget):
    """ 
      This class defines parameters for headward propagation of 
      erosion.
      parent : all parameters from the WindowsParameters class
     
    """
    def __init__(self, parent):
        
        self.parent = parent
        
        # Parameters label
        self.IncisionStartLabel = QLabel('Starting time of incision (Ma):')
        self.IncisionStartLabel.setFont(conf.font10)
        self.IncisionStartLabel.setToolTip("Choose time in past when the incision start.")
        self.IncisionStartLabel.setAlignment(Qt.AlignCenter)
        self.IncisionStopLabel = QLabel('Ending time of incision (Ma):')
        self.IncisionStopLabel.setFont(conf.font10)
        self.IncisionStopLabel.setToolTip("Choose time in past when the incision stop.")
        self.IncisionStopLabel.setAlignment(Qt.AlignCenter)
        self.IncisionDepthLabel = QLabel('Depth of incision (%):')
        self.IncisionDepthLabel.setFont(conf.font10)
        self.IncisionDepthLabel.setAlignment(Qt.AlignCenter)
        self.IncisionDepthLabel.setToolTip("Relative depth of incision. Computed relative to the incision depth remaining from elevation at time of incision and present-day elevation.")
        self.TauhLabel = QLabel('Scale propagation rate:')
        self.TauhLabel.setFont(conf.font10)
        self.TauhLabel.setAlignment(Qt.AlignCenter)
        self.TauhLabel.setToolTip('Scaling factor to set a time-varying non-linear propation rate. A positive value set an acceleration of propagation rate at the start of incision, a negative value set an acceleration toward the end of incision propagation.')
        
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.grid.addWidget(self.IncisionStartLabel,0,1)
        self.grid.addWidget(self.parent.IncisionStart,1,1)
        self.grid.addWidget(self.IncisionStopLabel,0,2)
        self.grid.addWidget(self.parent.IncisionStop,1,2)
        self.grid.addWidget(self.IncisionDepthLabel,2,1)
        self.grid.addWidget(self.parent.IncisionDepth,3,1)
        self.grid.addWidget(self.TauhLabel,2,2)
        self.grid.addWidget(self.parent.Tauh,3,2)
        
        