# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:47:40 2021


 This is the main file of the user interface. It creates the application and
 the main tools to navigate through the interface. The package includes the
 following modules:
     - configs.py : this is where we define the default parameters
     - PGUI_utils.py: contains functions and classes for the interface
     - Outputs.py: contains classes used to plot results 
     - Thermochronometers.py: contains classes related to the thermochronometers
     settings.
     
@author: Maxime Bernard.

"""
############################################################################
############################# Import modules ###############################
############################################################################
import sys
import os
# if hasattr(sys, '_MEIPASS'):
#       os.environ['QTWEBENGINEPROCESS_PATH'] = os.path.normpath(os.path.join(
#                   sys._MEIPASS, 'PyQt5', 'Qt5', 'lib',
#                   'QtWebEngineCore.framework', 'Versions','5','Helpers','QtWebEngineProcess'
#               ))
import configs as conf
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import (QWidget, QMainWindow, QApplication, QPushButton,
                             QVBoxLayout, QHBoxLayout, QGridLayout,
                             QAction, qApp, QStyle, QTabWidget, QLabel,
                             QLineEdit, QPlainTextEdit, QMessageBox, QCheckBox,
                             QSizePolicy, QMenu, QFrame,QComboBox,QTableWidgetItem, 
                             QFileDialog, QDockWidget,QTreeView, QStackedWidget,
                             QListView,QGroupBox,QSlider, QScrollArea, QRadioButton,
                             QErrorMessage)
from PyQt5.QtGui import (QPixmap, QIcon, QFont,QIntValidator)
from PyQt5.QtCore import QProcess, Qt, QSize
from PyQt5.Qt import QStandardItemModel, QStandardItem
import re
import stat
import shutil
from math import *
import pandas as pd
import io
from scipy.interpolate import interp1d
from scipy.io import loadmat
from scipy.optimize import curve_fit
import numpy as np
import matplotlib
import matplotlib.backends.backend_svg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import pyvista as pv
import pyvistaqt as pvqt
import folium, json
from folium.plugins import Draw
import subprocess
import webbrowser
import xarray as xr
import Utils.PGUI_utils as pgu
from Utils import misfits
import Utils.interface as U_interface
import Utils.GIS as GIS
import Utils.graphics.set_Thermochronometers as plot_thermo
import Utils.graphics.batch_results as batch
import Utils.graphics.plot_NA as plot_NA
import Utils.graphics.plot_NA as load_results
import Utils.graphics.plot_2Dmaps as plot_2D
import Utils.interface_graphics as GUI_graph
import Utils.interface_inputs as GUI_inputs
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

os.environ["QT_API"] = "pyqt5"  
    
##############################################################################
################### Check for the operating system ###########################
##############################################################################

#for Debug
#FolderPath = '.'

# Path of executable
FolderPath = conf.FolderPath

# Force to go to the directory of the executable
print("Working Directory: ", FolderPath)
os.chdir(FolderPath)
    
#Locate Icones directory
IconPath = conf.IconPath

 
# ##############################################################################
# ############################ Global parameters ##############################
# ##############################################################################

#For the format of doubles
#Force dot notation for the decimals
DoubleValidator = conf.DoubleValidator

#All the fonts to used in the interface
font6 = conf.font6
font8 = conf.font8
fontBold8 =  conf.fontBold8
font10 = conf.font10
fontBold10 = conf.fontBold10
fontItalic10 = conf.fontItalic10
font11 =  conf.font11
fontBold11 =  conf.fontBold11
font12 = conf.font12
fontLine12 = conf.fontLine12
fontBold12 =  conf.fontBold12


##############################################################################
############################ Functions #######################################
##############################################################################

#----------------------------------------------------------
def store_input(self, d):
    # This function stores the input parameters provided by the user
    # and update their value from the default value
    self.input_parameters.update(d)
    #print(str(self.input_parameters))



##############################################################################
################## Classes for Input Parameters area #######################
##############################################################################

#############################################################
class InputToWrite:
    """ This class stores the input parameters either provided by the user or by
      an old Pecube input file
    """
    def __init__(self):
        super().__init__()
        self.storeInput = {}



#############################################################
class MainWindow(QMainWindow, object):
    """
       This class handle the main widow of PecubeGUI.
       
    """
    def __init__(self, app, parent=None):
        super(MainWindow, self).__init__(parent)
        
        self.app = app            
        self.PrefList = conf.PrefList 
        self.FolderPath = FolderPath
        self.PecubePath = os.path.abspath(self.PrefList['PecubePath'])
        self.UI()
        self.oldInput = 0 #signal for old input file
        self.InputParamSignal = 0 #Signal for input parameters provided
        self.IconePlanxy = os.path.join(IconPath,"View_xy.png")
        self.IconePlanxz = os.path.join(IconPath,"View_xz.png")
        self.IconePlanyz = os.path.join(IconPath,"View_yz.png")
        self.runOpen = 0 # Signal for a runing process active
        self.runWindow = U_interface.run_Window(self)
        self.UsePredictedElevation = QComboBox()
        self.UsePredictedElevation.addItem('yes')
        self.UsePredictedElevation.addItem('no')
        if conf.PrefList['UsePredictedElevation'] == 'yes':
            self.UsePredictedElevation.setCurrentIndex(0)
        else:
            self.UsePredictedElevation.setCurrentIndex(1)
        self.UsePredictedElevation.currentIndexChanged.connect(lambda: self.update_PredElevation())

    #----------------------------------------------------------    
    def UI(self):
    
        # initiate Graphic window signal
        self.PlotWin = 0 #To check if the output
        self.ParamTable = 0
        self.GUITheme = 0
        

        self.projectTabs = QStackedWidget() # to store all opened tabs of projects
        
        #Define the Stacked layouts that will contain the main windows:
        #   - welcome message
        #   - Input parameters
        #   - Plot results
        self.centralWidget = QStackedWidget()
        self.setCentralWidget(self.centralWidget)
        
        #Welcome message
        self.wtext = QPlainTextEdit()
        self.wtext.setReadOnly(True)
        self.wtext.setFont(font10)
        welcomeText = "Welcome to PecubeGUI beta version !"
        descriptionText = """ \n\nto create a new input file : crtl+N
                            \nto load a previous input file: crtl+O"""
        self.wtext.appendPlainText(welcomeText)
        self.wtext.appendPlainText(descriptionText)
        self.centralWidget.addWidget(self.wtext)
        
        # Input Parameters
        self.ParamWidget = QWidget()
        self.centralWidget.addWidget(self.ParamWidget)
        
        # Select loaded projects
        self.ProjectLabel = QLabel('Active project:')
        self.ProjectLabel.setFont(fontBold12)
        self.ProjectLabel.adjustSize()
        self.ProjectSelection = QComboBox()
        self.ProjectSelection.setMinimumContentsLength(8)
        self.ProjectSelection.currentIndexChanged.connect(lambda: 
                    self.select_Project(self.ProjectSelection))
        self.indexProject = 0
        self.HBoxProject = QHBoxLayout()
        self.HBoxProject.addWidget(self.ProjectLabel)
        self.HBoxProject.addWidget(self.ProjectSelection)
        self.HBoxProject.addStretch(1)
        self.HBoxProject.setSpacing(20)
        self.GridProject = QVBoxLayout()
        self.GridProject.addLayout(self.HBoxProject)
        self.projectParam = QStackedWidget() # to store all the opened projects
        self.StackedTabs = [] # to stack tabs (right hand side of the interface) for each project
        self.CurrentTabs = 0 # Index of the active Tabs
        self.GridProject.addWidget(self.projectParam)
        self.widgetProject = QWidget()
        
        # Plotting area
        self.OutputWidget = QWidget()
        self.centralWidget.addWidget(self.OutputWidget)
        
        # ------------------ Create Menu bar ---------------------------------
        # What I need...
        # New input file button
        NewInput = QAction("New Input file       ", self)
        NewInput.setShortcut('Ctrl+N')
        NewInput.setStatusTip('New Input file')
        NewInput.triggered.connect(self.MainParameters)

        # Open Pecube Input file
        OpenFile = QAction("Open...", self)
        OpenFile.setShortcut('Ctrl+O')
        OpenFile.setStatusTip('Open input file')
        OpenFile.triggered.connect(self.file_open)

        # Exit button
        exitAction = QAction(self.style().standardIcon(QStyle.SP_DialogCancelButton),
                              '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit PecubeGUI')
        exitAction.triggered.connect(qApp.quit)

        self.statusBar()
        # MenuBar - add buttons
        mainMenu = self.menuBar()
        fileMenu = QMenu("&File", self)
        fileMenu.addAction(NewInput)
        fileMenu.addSeparator()
        fileMenu.addAction(OpenFile)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)
        mainMenu.addMenu(fileMenu)

        # Edit menu
        # Read saved file
        Preferences = QAction("Preferences", self)
        Preferences.setShortcut('Ctrl+P')
        EditMenu = mainMenu.addMenu("&Edit")
        EditMenu.addAction(Preferences)
        Preferences.triggered.connect(self.openPref)
        mainMenu.addMenu(EditMenu)
        helpMenu = mainMenu.addMenu("&Help")
        Documentation = QAction("PecubeGUI documentation...",self)
        Documentation.triggered.connect(lambda:self.openDoc())
        helpMenu.addAction(Documentation)
        PecubeDoc = QAction("Pecube documentation...",self)
        PecubeDoc.triggered.connect(lambda:self.openPecubeDoc())
        helpMenu.addAction(PecubeDoc)

        # ------------------- Create toolbar ----------------------------------
        self.toolbar = self.addToolBar('New')
        # #Input file button
        self.FolderPath = FolderPath
        IconNewInput = os.path.join(IconPath,"New_Input.png")
        NewInputt = QAction(QIcon(IconNewInput), "New Input file", self)
        NewInputt.setStatusTip('New input file')
        NewInputt.triggered.connect(self.MainParameters)

        # Open Input file button
        IconPecubeLoad = os.path.join(IconPath,"Pecube_load.png")
        openButton = QAction(QIcon(IconPecubeLoad),
                              "open Pecube file", self)
        openButton.setStatusTip('Open input file')
        openButton.triggered.connect(self.file_open)

        # Graphic output button
        IconKG1 = os.path.join(IconPath,"KG1.jpg")
        GraphButton = QAction(QIcon(IconKG1), "show output", self)
        GraphButton.setStatusTip('Show output')
        GraphButton.setCheckable(True)
        GraphButton.toggled.connect(lambda: self.plotOutput(GraphButton))
        
        # Add tools to toolbar
        self.toolbar.addAction(NewInputt)  # Add input file button
        self.toolbar.addAction(openButton)  # Add open Pecube button
        self.toolbar.addAction(GraphButton)  # Add Graph button

        # ------------------- Set Central Widget ------------------------------
        self.statusBar().showMessage("PecubeGUI")
        #Add welcome message to the central window
        self.centralWidget.setCurrentWidget(self.wtext)

        # ------------------ Show the Window ---------------------------------
        # Window settings
        self.setGeometry(200, 500, 1200, 800)
        self.setMinimumSize(1000,800)
        self.IconPecube = os.path.join(IconPath,"Pecube_icon.ico")
        self.setWindowIcon(QtGui.QIcon(self.IconPecube))
        self.setWindowTitle('PecubeGUI')
        self.showMaximized()

    #----------------------------------------------------------
    def select_Project(self, project):
        """ To navigate through opened projects """
        self.indexProject = project.currentIndex()
        self.projectParam.setCurrentIndex(self.indexProject)
        self.CurrentTabs = self.indexProject
        
    #----------------------------------------------------------
    def openPref(self):
        """ Open the preference window. """
        
        self.PrefWin = U_interface.WinExtraParameters(width=500,height=500,text='PecubeGUI preferences',fixed=True)
        self.PrefWin.OkButton.clicked.connect(lambda: self.PrefWin.close())
        self.Title = QLabel('Preferences')
        self.Title.setFont(fontBold12)
        # Theme Label
        self.Label = QLabel('GUI theme:')
        self.Label.setScaledContents(True)
        self.Label.setFont(font11)
        self.Label.setAlignment(Qt.AlignLeft)
        # Theme Combo box
        self.Themes = QComboBox()
        self.Themes.addItem('White')
        self.Themes.addItem('Dark')
        if self.PrefList['Theme'] == 'Dark':
            self.Themes.setItemText(0, 'Dark  ')
            self.Themes.setItemText(1, 'White  ')
        else:
            self.Themes.setItemText(0, 'White  ')
            self.Themes.setItemText(1, 'Dark  ')
        # Show Console
        self.ShowConsole = QLabel('Show Console: ')
        self.ShowConsole.setFont(font11)
        self.ShowConsole.setAlignment(Qt.AlignLeft)
        self.ShowConsole.setToolTip('Show Pecube run state in a console: yes/no')
        self.ShowConsoleb = QCheckBox()
        # Pecube's folder path
        self.PathLabel = QLabel("Pecube's folder path: " + self.PecubePath)
        self.PathLabel.setScaledContents(True)
        self.PathLabel.setFont(font11)
        # Predicted elevations or observed elevations
        self.UsePredictedElevationLabel = QLabel('Use predicted elevations:')
        self.UsePredictedElevationLabel.setFont(font11)
        # Save button
        self.SavePref = QPushButton('save preferences...')
        self.SavePref.clicked.connect(lambda: self.saveFile())
        
        #Check for updates
        self.Themes.currentIndexChanged.connect(lambda: self.ComboChanged())
        self.ShowConsoleb.stateChanged.connect(lambda: self.store_CheckBox(self.ShowConsoleb,'ShowConsole'))
        
        # Set up the grid
        Grid1 = QGridLayout()
        Grid1.addWidget(self.Title, 0, 0, 2, 2)
        Grid1.addWidget(self.PathLabel, 3, 0, 1, 2)
        Grid1.addWidget(self.Label, 4, 0)
        Grid1.addWidget(self.Themes, 4, 1)
        Grid1.addWidget(self.ShowConsole,5,0)
        Grid1.addWidget(self.ShowConsoleb,5,1)
        Grid1.addWidget(self.UsePredictedElevationLabel,6,0)
        Grid1.addWidget(self.UsePredictedElevation,6,1)
        Grid1.addWidget(self.SavePref, 7, 0)
        Grid1.setRowStretch(7, 2)
        Grid1.setColumnStretch(3, 1)
        Grid1.setSpacing(20)
        self.PrefWin.layout.addLayout(Grid1, 0, 0)
        Q = QWidget()
        Q.setLayout(self.PrefWin.layout)
        self.PrefWin.setCentralWidget(Q)
        conf.WindowsOpen.append(self.PrefWin)
        self.PrefWin.show()
    
    #----------------------------------------------------------
    def update_PredElevation(self):
        """ update whether the user wants to plot data with their predicted elevation or from
        the observed elevation of samples."""
        
        print('update predicted elevation setting.')
        if self.UsePredictedElevation.currentIndex() == 0:
            self.PrefList['UsePredictedElevation'] = 'yes'
        else:
            self.PrefList['UsePredictedElevation'] = 'no'
        
    #----------------------------------------------------------
    def store_CheckBox(self, d, s):
        """
          Update the value of the check box each time it is changed
          d is the check box
          s is the parameter name
        """
        if d.isChecked() == True:
            self.PrefList[s] = '1'
        else:
            self.PrefList[s] = '0'
    
    #----------------------------------------------------------
    def openDoc(self):
        """ Open the online documentation. """
        webbrowser.open('https://pecubegui-doc.readthedocs.io/en/latest/')
    
    #----------------------------------------------------------
    def openPecubeDoc(self):
        """ Open the Pecube documentation (pdf). """
        try:
            subprocess.call([os.path.join(conf.PecubeFolderPath,'docs'),'Pecube.pdf'])
        except PermissionError:
            QErrorMessage(self).showMessage("Permission issue when trying to open the Pecube documentation.")
            return
        
    #----------------------------------------------------------
    def saveFile(self):
        """ To save the file "preferences.txt". """
        
        name = os.path.join(self.FolderPath,"preferences.txt")
        try:
            file = open(name, 'w')
        except PermissionError:
            file = os.open(name, os.O_RDWR|os.O_CREAT)
            os.close(file)
            os.chmod(name,stat.S_IWUSR)
            file = open(name, 'w')
        for k, v in self.PrefList.items():
            file.write(str(k) + ' = ' + str(v) + '  \n')
        file.close()

    #----------------------------------------------------------
    def ComboChanged(self):
        """ To set the theme of PecubeGUI. """
        
        if self.Themes.currentText() == 'Dark  ':
            self.GUITheme = 1
            self.app.setStyleSheet(conf.styleSheet)
            self.PrefList['Theme'] = 'Dark  '
            self.app.setStyle('Fusion')
        elif self.Themes.currentText() == 'White  ':
            self.GUITheme = 0
            self.app.setStyleSheet("")
            self.app.setStyle('Fusion')
            self.PrefList['Theme'] = 'White  '

    #----------------------------------------------------------
    def MainParameters(self):
        """ 
          This function asks the user to define a name of the new project
          and shows the input parameters to be set
          Get New project Name set by the user
          
        """
        
        EmptyName = U_interface.emptyText()
        PName = U_interface.WindowProject(self, EmptyName)
        conf.ProjectName = PName.ProjectName
        
        if PName.ProjectName and len(PName.ProjectName) == 5:
            DParam = conf.DefaultParameterValues()
            # set Table for input parameters
            IParam = InputToWrite
            IParam.storeInput = {}
            self.SplitTabs = U_interface.QCustomTabWidget()
            ProjectPath = os.path.join(conf.PecubeFolderPath,PName.ProjectName)
            if os.path.exists(ProjectPath):
                win = U_interface.MessageBox(text1="This Project already exists. If you want to work on this project click 'yes'. If you want to make a copy click 'no'. ")
                if win.QMessage == QMessageBox.Yes:
                    if os.path.exists(os.path.join(ProjectPath,'NA')):
                        shutil.rmtree(os.path.join(ProjectPath,'NA'))
                    else: 
                        os.mkdir(os.path.join(ProjectPath,'NA'))
                    if os.path.exists(os.path.join(ProjectPath,'VTK')):
                        shutil.rmtree(os.path.join(ProjectPath,'VTK'))
                    else:
                        os.mkdir(os.path.join(ProjectPath,'VTK'))
                    if os.path.exists(os.path.join(ProjectPath,'LOG')):
                        shutil.rmtree(os.path.join(ProjectPath,'LOG'))
                    else: 
                        os.mkdir(os.path.join(ProjectPath,'LOG'))
                    if os.path.exists(os.path.join(ProjectPath,'output')):
                        shutil.rmtree(os.path.join(ProjectPath,'output'))
                    else:
                        os.mkdir(os.path.join(ProjectPath,'output')) 
                        
                    ProjectNameNew = PName.ProjectName
                elif win.QMessage == QMessageBox.No:
                    EmptyName = U_interface.emptyText()
                    PName = U_interface.WindowProject(self, EmptyName)
                    if PName.ProjectName == 'thisisempty':
                        return
                    conf.ProjectName = PName.ProjectName
                    ProjectNameNew = PName.ProjectName
                    ProjectPath = os.path.join(conf.PecubeFolderPath,PName.ProjectName)
                    os.mkdir(ProjectPath)
                    os.mkdir(os.path.join(ProjectPath,'input'))
                    os.mkdir(os.path.join(ProjectPath,'output'))
                    os.mkdir(os.path.join(ProjectPath,'NA'))
                    os.mkdir(os.path.join(ProjectPath,'VTK'))
                    os.mkdir(os.path.join(ProjectPath,'data'))
                    os.mkdir(os.path.join(ProjectPath,'LOG'))
                elif win.QMessage == QMessageBox.Cancel:
                    return
            else:
                os.mkdir(ProjectPath)
                os.mkdir(os.path.join(ProjectPath,'input'))
                os.mkdir(os.path.join(ProjectPath,'output'))
                os.mkdir(os.path.join(ProjectPath,'NA'))
                os.mkdir(os.path.join(ProjectPath,'VTK'))
                os.mkdir(os.path.join(ProjectPath,'data'))
                os.mkdir(os.path.join(ProjectPath,'LOG'))
            ProjectNameNew = PName.ProjectName
                
            # Get the input Parameters
            self.ParamTable = GUI_inputs.ParamWin(self,ProjectNameNew, DParam, IParam, 0)
            
            #Define the splitters for the Input parameter area
            self.SplitterParam = U_interface.QCustomSplitter(Qt.Horizontal) #QSplitter(Qt.Horizontal)
            
            # Store the project
            # self.storeProjects.append(self.ParamTable)
            if self.ParamTable.success == 0:
                return
            self.InputParamSignal = 1
            # Check whether QSplitter already contains widgets
            # Add the input parameters area into the first splitter (left)
            self.ProjectSelection.addItem(PName.ProjectName)
            # self.projectParam.setCurrentIndex()
            self.SplitterParam.addWidget(self.ParamTable.widget)
            # Same check with the second splitter
            self.SplitterParam.addWidget(self.SplitTabs)
            self.SplitterParam.setSizes([500,800])
            
            # set the vertical box to embed in the main Window
            self.projectParam.addWidget(self.SplitterParam)
            self.StackedTabs.append(self.SplitTabs)
            self.widgetProject.setLayout(self.GridProject)
            self.ParamWidget.setLayout(self.GridProject)
            # Embed to the main Window
            self.centralWidget.setCurrentWidget(self.ParamWidget)
            self.ProjectSelection.setCurrentIndex(self.ProjectSelection.count() -1)
         
    #----------------------------------------------------------        
    def file_open(self):
        """ To open and load parameters from an old input file. """
        
        OpenFile = U_interface.oldProjectName(self)
        if OpenFile.fileName:
            # Check whether project already exists
            ProjectPath = os.path.join(conf.PecubeFolderPath,OpenFile.ProjectName)
            if os.path.exists(ProjectPath):
                win = U_interface.MessageBox(text1="This Project already exists. If you want to work on this project click 'yes'. If you want to make a copy click 'no'. ")
                if win.QMessage == QMessageBox.Yes:
                    # if os.path.exists(os.path.join(ProjectPath,'NA')):
                    #     shutil.rmtree(os.path.join(ProjectPath,'VTK'))
                    #     shutil.rmtree(os.path.join(ProjectPath,'NA'))
                    #     shutil.rmtree(os.path.join(ProjectPath,'LOG'))
                    #     shutil.rmtree(os.path.join(ProjectPath,'output'))
                    #     os.mkdir(os.path.join(ProjectPath,'VTK'))
                    #     os.mkdir(os.path.join(ProjectPath,'output'))
                    # else:
                    #     os.mkdir(ProjectPath)
                    #     os.mkdir(os.path.join(ProjectPath,'input'))
                    #     os.mkdir(os.path.join(ProjectPath,'output'))
                    #     os.mkdir(os.path.join(ProjectPath,'NA'))
                    #     os.mkdir(os.path.join(ProjectPath,'VTK'))
                    #     os.mkdir(os.path.join(ProjectPath,'data'))
                    #     os.mkdir(os.path.join(ProjectPath,'LOG'))
                    ProjectNameNew = OpenFile.ProjectName
                elif win.QMessage == QMessageBox.No:
                    EmptyName = U_interface.emptyText()
                    PName = U_interface.WindowProject(self, EmptyName)
                    if PName.ProjectName == 'thisisempty':
                        return
                    conf.ProjectName = PName.ProjectName
                    ProjectNameNew = PName.ProjectName
                    # Copy old input file into new folder
                    # Is there any topo files that should be copied?
                    old_folder_path = os.path.join(self.PecubePath,OpenFile.ProjectName,"input")
                    new_folder_path = os.path.join(self.PecubePath,ProjectNameNew,"input")
                    print(old_folder_path,new_folder_path)
                    if os.path.exists(old_folder_path) and old_folder_path != new_folder_path:
                        try:
                            shutil.copytree(old_folder_path,new_folder_path)
                        except FileExistsError:
                            #Erase first and copy again
                            shutil.rmtree(new_folder_path)
                            shutil.copytree(old_folder_path,new_folder_path)
                        except:
                            print("An error occured when trying to copy the input files...")
                    # Is there any topo files that should be copied?
                    old_folder_path = os.path.join(self.PecubePath,OpenFile.ProjectName,"data","SPM")
                    new_folder_path = os.path.join(self.PecubePath,ProjectNameNew,"data","SPM")
                    if os.path.exists(old_folder_path):
                        try:
                            shutil.copytree(old_folder_path,new_folder_path)
                        except FileExistsError:
                            print("Spm files already exist in the targeted directory, the files are not copied.")
                        except:
                            print("An error occured when trying to copy the spm files...")
                    # Is there any output files that should be copied?
                    old_folder_path = os.path.join(self.PecubePath,OpenFile.ProjectName,"output")
                    new_folder_path = os.path.join(self.PecubePath,ProjectNameNew,"output")
                    if os.path.exists(old_folder_path):
                        try:
                            shutil.copytree(old_folder_path,new_folder_path)
                        except FileExistsError:
                            print("Output files already exist in the targeted directory, the files are not copied.")
                        except:
                            print("An error occured when trying to copy the output files...")
                    # Is there any NA files that should be copied?
                    old_folder_path = os.path.join(self.PecubePath,OpenFile.ProjectName,"NA")
                    new_folder_path = os.path.join(self.PecubePath,ProjectNameNew,"NA")
                    if os.path.exists(old_folder_path):
                        try:
                            shutil.copytree(old_folder_path,new_folder_path)
                        except FileExistsError:
                            print("NA files already exist in the targeted directory, the files are not copied.")
                        except:
                            print("An error occured when trying to copy the output files...")
                    # Is there any topo files that should be copied?
                    old_folder_path = os.path.join(self.PecubePath,OpenFile.ProjectName,"data")
                    new_folder_path = os.path.join(self.PecubePath,ProjectNameNew,"data")
                    if os.path.exists(old_folder_path):
                        try:
                            shutil.copytree(old_folder_path,new_folder_path)
                        except FileExistsError:
                            print("data files already exist in the targeted directory, the files are not copied.")
                        except:
                            print("An error occured when trying to copy the data files...")
                    new_folder_path = os.path.join(self.PecubePath,ProjectNameNew,"VTK")
                    if os.path.exists(new_folder_path):
                        pass
                    else:
                        os.mkdir(new_folder_path)
                elif win.QMessage == QMessageBox.Cancel:
                    return
                
            self.oldInput = 1
            try:
                self.oldInput = pgu.readOldInput(OpenFile.ProjectName)
            except Exception as E:
                raise
                QErrorMessage(self).showMessage('Error with loading data:' + str(E))
                return
            self.SplitTabs = U_interface.QCustomTabWidget()
            self.ParamTable = GUI_inputs.ParamWin(self,ProjectNameNew, self.oldInput.ParametersInput, self.oldInput, 1)
            
            if self.ParamTable.success:
                self.InputParamSignal = 1
                
                # Add the name of the project
                self.ProjectSelection.addItem(ProjectNameNew)
                
                # Splitter
                self.SplitterParam = U_interface.QCustomSplitter(Qt.Horizontal) #QSplitter(Qt.Horizontal)
                
                # Add to splitter
                self.SplitterParam.addWidget(self.ParamTable.widget)
                self.SplitterParam.addWidget(self.SplitTabs)
                self.SplitterParam.setSizes([500,800])
        
                # set the vertical box to embed in the main Window
                self.projectParam.addWidget(self.SplitterParam)
                self.StackedTabs.append(self.SplitTabs)
                self.widgetProject.setLayout(self.GridProject)
                self.ParamWidget.setLayout(self.GridProject)
                # Embed to the main Window
                self.centralWidget.setCurrentWidget(self.ParamWidget)
                self.ProjectSelection.setCurrentIndex(self.ProjectSelection.count() -1)
                
                #Check whether the plotting area is not open
                if self.PlotWin:
                    self.removeDockWidget(self.PlotWin.DockData)
                    self.removeDockWidget(self.PlotWin.Dock2dProperties)
                    try:
                        self.removeDockWidget(self.PlotWin.Dock3dProperties)
                    except AttributeError:
                        print("no Dock3dProperties")

    #-----------------------------------------------------------
    def plotOutput(self, b):
        """ 
          To handle the window that show within the main window of PecubeGUI
          Either input parameters or plot area
          
        """
        if b.isChecked() == True:
            if self.PlotWin and self.ParamTable: #Already defined
                #Show the plotting area window
                self.centralWidget.setCurrentWidget(self.OutputWidget)
                #Restore the Dockable windows of the plotting area
                self.restoreDockWidget(self.PlotWin.DockData)
                self.restoreDockWidget(self.PlotWin.Dock2dProperties)
                for i in range(len(self.PlotWin.List3DProp)):
                    self.restoreDockWidget(self.PlotWin.List3DProp[i])
            else:
                self.PlotWin = GUI_graph.GraphWin(self)
                self.centralWidget.setCurrentWidget(self.OutputWidget)

        else:
            if self.PlotWin and self.InputParamSignal:
                #remove the Dockable windows of the plotting area
                self.removeDockWidget(self.PlotWin.DockData)
                self.removeDockWidget(self.PlotWin.Dock2dProperties)
                for i in range(len(self.PlotWin.List3DProp)):
                    self.removeDockWidget(self.PlotWin.List3DProp[i])
                #Show the Input parameters
                self.centralWidget.setCurrentWidget(self.ParamWidget)
            else: #Means no input parameters, return to welcome message
                self.removeDockWidget(self.PlotWin.DockData)
                self.removeDockWidget(self.PlotWin.Dock2dProperties)
                for i in range(len(self.PlotWin.List3DProp)):
                    self.removeDockWidget(self.PlotWin.List3DProp[i])
                self.centralWidget.setCurrentWidget(self.wtext)
            
    #----------------------------------------------------------
    def closeEvent(self, event):
        """ To ask the user to confirm the closure of PecubeGUI """
        reply = QMessageBox.question(self, "Window Close", "Are you sure you want to close PecubeGUI?",
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        
        if reply == QMessageBox.Yes:
            if self.PlotWin:
                if self.PlotWin.Object3D:
                    self.PlotWin.Object3D.plotter.close()
            try:
                if self.ParamTable.ParamTable.setTopo:
                    self.ParamTable.ParamTable.setTopo.plotSpace.plotter.close()
            except AttributeError:
                pass
            
            try:
                if self.ParamTable.ParamTable.setFault:
                    self.ParamTable.ParamTable.setFault.plotSpace.plotter.close()
            except AttributeError:
                pass
                
            #check if other windows are opened
            for i in conf.WindowsOpen:
                if i:
                    i.close()
            conf.WindowsOpen = []
            self.app.quit()
            event.accept()
    
        else:
            event.ignore()
            
            
            
##############################################################################            
#############################    Main    #####################################
############################################################################## 
def main():
    app = QApplication(sys.argv)  # initializer
    app.setStyleSheet(conf.styleSheet)
    app.setStyle(conf.appStyle)
        
    # QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)

    # Try to locate the Pecube folder with the default path
    # if no success, ask the user to provide the path    
    try:
        TestFile = open(os.path.join(conf.PecubeFolderPath,"bin","compile.sh"),'r')
        TestFile.close()
    except FileNotFoundError:
        try:
            path = os.path.join(conf.FolderPath,"Pecube","bin","compile.sh")
            TestFile = open(path,'r')
            TestFile.close()
            path = os.path.abspath(os.path.join(conf.FolderPath,"Pecube"))
            conf.PrefList['PecubePath'] = path
            conf.PecubeFolderPath = path
            file = open(conf.PreferencesPath, 'w+')
            for k, v in conf.PrefList.items():
                file.write(str(k) + ' = ' + str(v) + ' \n')
            file.close()
        except FileNotFoundError:
            reply = QMessageBox()
            reply.setText("The Pecube's folder has not been found.\n"+
                          "Please, provide the current path of Pecube.")
            reply.setIcon(QMessageBox.Information)
            reply.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            retval = reply.exec_()
            if retval == QMessageBox.Ok:
                PecubeFolderPath = QFileDialog.getExistingDirectory()
                conf.PecubeFolderPath = os.path.abspath(PecubeFolderPath)
                print(conf.PecubeFolderPath)
                conf.PrefList['PecubePath'] = PecubeFolderPath
                # Save new path
                conf.PreferencesPath = os.path.join(FolderPath,"preferences.txt")
                # if sys.platform == 'win32' or sys.platform =='cygwin':
                    # change_access_rights(conf.PreferencesPath, AccessRight.Full)
                file = open(conf.PreferencesPath, 'w+')
                for k, v in conf.PrefList.items():
                    file.write(str(k) + ' = ' + str(v) + ' \n')
                file.close()
                
    # Show the Pecube's interface
    window = MainWindow(app)  # create the window
    conf.WindowsOpen.append(window)
    appEx =app.exec_()
    
    try:
        sys.exit(appEx)
    except SystemExit:
        print( "Closing window...")
            
    
if __name__ == '__main__':
    main()
    
