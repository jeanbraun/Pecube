#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This module contains functions and classes used in the interface to perform specific tasks.

@author: Maxime Bernard

"""


# Import modules
import configs as conf
from PyQt5.QtWidgets import (QErrorMessage,QSplitter,QTabWidget,QMainWindow,
                             QMdiArea,QDesktopWidget,QMessageBox,QMdiSubWindow,
                             QProgressBar,QLabel,QHBoxLayout,QPushButton,
                             QVBoxLayout,QPlainTextEdit,QWidget, QFileDialog,
                             QInputDialog,QLineEdit,QToolBar,QAction,QDockWidget,
                             QCheckBox,QComboBox,QGroupBox,QGridLayout,QTableWidget,
                             QTableWidgetItem,QApplication,qApp,QTabBar,QStylePainter,
                             QStyleOptionTab,QStyle,QFrame)
from PyQt5.QtGui import (QIcon,QKeySequence,QFont)
from PyQt5.QtCore import (Qt)
from PyQt5.Qt import QStandardItemModel, QStandardItem
from PyQt5 import QtCore
import matplotlib
import matplotlib.backends.backend_svg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import os
import pandas as pd
import io
import numpy as np
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import pyvistaqt as pvqt
import shutil
import rioxarray
import xarray as xr
import math

##############################################################################
############################ Functions #######################################
##############################################################################

################### General ####################
#----------------------------------------------------------
def message(self, s):
    # Text to show in the console when a Pecube model runs
    self.text.appendPlainText(s)

#----------------------------------------------------------
def errorMessage(text):
    # This function show error message
    QErrorMessage().showMessage(text)
    # print(str(self.input_parameters))

#----------------------------------------------------------
def findOccurrences(text,ch):
    # This function find occurence of a character in strings
    return [i for i, letter in enumerate(text) if letter == ch]

############## To check imputs formats ####################
#----------------------------------------------------------
def isFloat(text):
    """ To check if a value is a float """
    try:
        float(text)
    except:
        QErrorMessage().showMessage("Only float number allowed.")

#----------------------------------------------------------
def isString(text):
    """ To check if a value is a string """
    if type(text) != str:
        QErrorMessage().showMessage("Only characters allowed.")
        
        
        
############## Actions on data reading ####################
#----------------------------------------------------------   
def readOldInput(folderName):
    """
    To read some parameters from the pecube project we 
    want to plot results
    
    """
    input_file = os.path.join(conf.PecubeFolderPath,folderName,"input","Pecube.in")
    oldParameters = old_Input(folderName,input_file)
    
    return oldParameters
    
#----------------------------------------------------------      
def removeZeros(a, n):
    """ To remove zeros values at the end of an array
    index to store the first non-zero number
    """
    ind = -1
 
    # go throughout the array backward
    # and find the first non-zero number
    for i in reversed(range(n)):
        if (a[i] != float(0)):
            ind = i
            break
 
    # if no non-zero
    # number is there
    if (ind == -1):
        return

    # removing zeros
    b = a[0:ind+1]
    
    return b


#----------------------------------------------------------
def isNaN(num):
    return num != num

#----------------------------------------------------------
def isNeg(num):
    return num < 0


########### To return mathematical functions ##############
#----------------------------------------------------------
def objective_poly(x, a, b, c):
    """ Define a polynomial function of second order
    Used in the class "GraphWin" """
    return a + b * x + c * x**2 

#----------------------------------------------------------
def objective_poly3(x, a, b, c, d):
    """ Define a polynomial function of third order
    Used in the class "GraphWin" """
    return a + (b * x) + (c * x**2) + (d * x**3) 

#----------------------------------------------------------
def objective_poly4(x, a, b, c, d, e):
    """ Define a polynomial function of fourth order
    Used in the class "GraphWin" """
    return a + (b * x) + (c * x**2) + (d * x**3) + (e * x**4)

#----------------------------------------------------------
def objective_poly5(x, a, b, c, d, e, f):
    """ Define a polynomial function of fifth order
    Used in the class "GraphWin" """
    return a + (b * x) + (c * x**2) + (d * x**3) + (e * x**4) + (f * x**5)

#----------------------------------------------------------
def objective_power(x, a, b):
    """ Define a power-law function
    Used in the class "GraphWin" """
    return a * x**b 

#----------------------------------------------------------
def moving_average(a,n):
    """ Compute average in radius n in a matrix a."""
    nx = np.size(a,0)
    ny = np.size(a,1)
    R = a
    for i in np.arange(1+n,nx-n):
        for j in np.arange(1+n,ny-n):
            M = a[i-n:i+n,j-n:j+n]
            R[i,j] = np.mean(M[:])
            
    return R
    



#----------------------------------------------------------
def clean_input_file(oldParam,DefaultParam):
    """
    This function clean the input file. If the user loaded a previous Pecube
    project and changed some parameter value, we want to remove parameters that
    are not used anymore.
    
    """
    # start from oldParam dictionary
    # if old parameter value == default value then remove from old dictionary
    defaultKeys = DefaultParam.DParameters.keys()
    new_dic = {key:value for key, value in oldParam.items() if key in defaultKeys and DefaultParam.DParameters[key] != value}
    
    # Handle time evolution topographic step
    if 'ntime' in oldParam:
        ntime = int(oldParam['ntime'])
        for i in range(ntime):
            new_dic['time_topo'+str(i+1)] = oldParam['time_topo'+str(i+1)]
            new_dic['amplification'+str(i+1)] = oldParam['amplification'+str(i+1)] 
            new_dic['offset'+str(i+1)] = oldParam['offset'+str(i+1)] 
            new_dic['output'+str(i+1)] = oldParam['output'+str(i+1)]
            try:
                new_dic['phase'+str(i+1)] = oldParam['phase'+str(i+1)] 
            except:
                continue
            
    # Handle uplift kinetic table
    if 'nstep1' in oldParam and 'nfault' in oldParam:
        nstep1 = int(oldParam['nstep1'])
        nfault = int(oldParam['nfault'])
        for i in range(nfault):
            new_dic['nstep'+str(i+1)] = oldParam['nstep'+str(i+1)]
            for j in range(nstep1):
                new_dic['time_start'+str(i+1)+'_'+str(j+1)] = oldParam['time_start'+str(i+1)+'_'+str(j+1)]
                new_dic['time_end'+str(i+1)+'_'+str(j+1)] = oldParam['time_end'+str(i+1)+'_'+str(j+1)]
                new_dic['velo'+str(i+1)+'_'+str(j+1)] = oldParam['velo'+str(i+1)+'_'+str(j+1)]
                
    # Handle Fault geometry table
    if 'nfault' in oldParam and 'npoint1' in oldParam:
        npoint1 = int(oldParam['npoint1'])
        nfault = int(oldParam['nfault'])
        for i in range(nfault):
            new_dic['npoint'+str(i+1)] = oldParam['npoint'+str(i+1)]
            for j in range(npoint1):
                new_dic['r'+str(i+1)+'_'+str(j+1)] = oldParam['r'+str(i+1)+'_'+str(j+1)]
                new_dic['s'+str(i+1)+'_'+str(j+1)] = oldParam['s'+str(i+1)+'_'+str(j+1)]
                
    ## We should not erase e.g., nstep2 if load from previous project but still want to use it
    return new_dic

    
##############################################################################
############################## Classes #######################################
##############################################################################


#############################################################
class old_Input:
    """ 
    This class load and store old input parameters.
    
    """
    
    def __init__(self,folderName,input_file):
        Read_parameters = {}
        try:
            open(input_file)
        except:
            print("No input file found. The user has to load a project to plot data.")
            return
        with open(input_file) as file:
            for line in file:
                if "=" in line:
                    (key, _, par) = line.split()
                    Read_parameters[str(key)] = par
        self.ParametersInput = conf.DefaultParameterValues()
        for k, v in Read_parameters.items():
            for i in self.ParametersInput.DParameters:
                if k == i:
                    self.ParametersInput.DParameters[k] = v
                    break
                else:  # Means parameters inside table such as topo1...
                    self.ParametersInput.TableParameters[k] = v
        self.storeInput = Read_parameters
        #Check for sample specific files
        if int(self.ParametersInput.DParameters[conf.Variable_names['Mode_age_computation']]) == 2: # Means data have been provided
            data_folder_name = self.ParametersInput.DParameters['data_folder']
            try:
                with open(os.path.join(conf.PecubeFolderPath,folderName,"data","Samples_settings.txt")) as file:
                    parameters_list = []
                    for line in file:
                        line = line.rstrip('\n')
                        parameters_list.append(line.split())
                    #store parameters     
                    self.ParametersInput.DParameters['Nb_Samples'] = parameters_list[1][0]
                    nbSamples = int(parameters_list[1][0])
                    self.ParametersInput.DParameters['Nb_TotGrains'] = parameters_list[2][0]
                    self.ParametersInput.DParameters['age_AHe_flag'] = parameters_list[3][0]
                    self.ParametersInput.DParameters['age_ZHe_flag'] = parameters_list[3][1]
                    self.ParametersInput.DParameters['age_AFT_flag'] = parameters_list[3][2]
                    self.ParametersInput.DParameters['43He_flag'] = parameters_list[4][0]
                    self.ParametersInput.DParameters['zonation_flag'] = parameters_list[5][0]
                    for z in range(int(parameters_list[1][0])): 
                        self.ParametersInput.TableParameters['SampleName'+str(z+1)] = parameters_list[6+z][0]
                        sampleName_temp = parameters_list[6+z][0]
                        self.ParametersInput.TableParameters[sampleName_temp+'_lon'] = parameters_list[6+z][1]
                        self.ParametersInput.TableParameters[sampleName_temp+'_lat'] = parameters_list[6+z][2]
                        self.ParametersInput.TableParameters[sampleName_temp+'_Elev'] = parameters_list[6+z][3]
                    self.ParametersInput.TableParameters['nb_grains'] = parameters_list[6+nbSamples:6+nbSamples*2]
                    idnumber = 7+z
                    # count = 0
                    # for i in range(len(self.ParametersInput.TableParameters['nb_grains'])):
                    #     nb_grains = self.ParametersInput.TableParameters['nb_grains'][i][0]
                    #     for ii in range(int(nb_grains)):
                    #         self.ParametersInput.TableParameters['grainsizeApatite'+str(i+1)+"_"+str(ii+1)] = parameters_list[count+idnumber+nbSamples][0]
                    #         self.ParametersInput.TableParameters['U_grainApatite'+str(i+1)+"_"+str(ii+1)] = parameters_list[count+idnumber+nbSamples][1]
                    #         self.ParametersInput.TableParameters['Th_grainApatite'+str(i+1)+"_"+str(ii+1)] = parameters_list[count+idnumber+nbSamples][2]
                    #         self.ParametersInput.TableParameters['grainsizeZircon'+str(i+1)+"_"+str(ii+1)] = parameters_list[count+idnumber+nbSamples][3]
                    #         self.ParametersInput.TableParameters['U_grainZircon'+str(i+1)+"_"+str(ii+1)] = parameters_list[count+idnumber+nbSamples][4]
                    #         self.ParametersInput.TableParameters['Th_grainZircon'+str(i+1)+"_"+str(ii+1)] = parameters_list[count+idnumber+nbSamples][5]
                    #         count+=1
                    # i = int(self.ParametersInput.DParameters['Nb_TotGrains'])-1
                    i=0
                    #AHe parameters
                    self.ParametersInput.DParameters['Alpha_Ejec_Mod_Flag'] = parameters_list[i+idnumber+nbSamples][0]
                    self.ParametersInput.DParameters['Alpha_Ejec_Flag'] = parameters_list[i+idnumber+nbSamples][1]
                    self.ParametersInput.DParameters['diff_model'] = parameters_list[i+idnumber+nbSamples][2]
                    self.ParametersInput.DParameters[conf.Variable_names['D0_Apatite']] = parameters_list[i+idnumber+nbSamples][3]
                    self.ParametersInput.DParameters[conf.Variable_names['Ea_Apatite']] = parameters_list[i+idnumber+nbSamples][4]
                    #ZHe parameters
                    idnumber += 1
                    self.ParametersInput.DParameters['DiffzModel'] = parameters_list[i+idnumber+nbSamples][0]
                    self.ParametersInput.DParameters[conf.Variable_names['D0_Zircon']] = parameters_list[i+idnumber+nbSamples][1]
                    self.ParametersInput.DParameters[conf.Variable_names['Ea_Zircon']] = parameters_list[i+idnumber+nbSamples][2]
                    #AFT parameters
                    idnumber += 1
                    self.ParametersInput.DParameters['AnnModel'] = parameters_list[i+idnumber+nbSamples][0]
                    self.ParametersInput.DParameters['rhoST'] = parameters_list[i+idnumber+nbSamples][1]
                    idnumber += 1
                    # self.ParametersInput.DParameters['NstepsHe43'] = parameters_list[i+idnumber+nbSamples][0]
                    # self.ParametersInput.DParameters['nColHe43'] = parameters_list[i+idnumber+nbSamples][1]
                    # self.ParametersInput.DParameters['cRowHe43'] = parameters_list[i+idnumber+nbSamples][2]
                    # idnumber += 1
                    # for j in range(int(self.ParametersInput.DParameters['NstepsHe43'])):
                    #     self.ParametersInput.TableParameters['Heating'+str(j+1)] = parameters_list[i+idnumber+j+nbSamples][0]
                    #     self.ParametersInput.TableParameters['Duration'+str(j+1)] = parameters_list[i+idnumber+j+nbSamples][1]
                    j=0
                    idnumber += 1
                    if self.ParametersInput.DParameters['zonation_flag'] == '1': #Means old zonation   
                        self.ParametersInput.DParameters['nLayers'] = parameters_list[i+idnumber+j+nbSamples][0]
                        counter = 0
                        idnumber += 1
                        for k in range(int(self.ParametersInput.DParameters['Nb_TotGrains'])):
                            for l in range(int(self.ParametersInput.DParameters['nLayers'])):
                                self.ParametersInput.TableParameters['layer'+str(k+1)+'_'+str(l+1)] = parameters_list[i+idnumber+counter+j+nbSamples][0]
                                self.ParametersInput.TableParameters['zWeight'+str(k+1)+'_'+str(l+1)] = parameters_list[i+idnumber+counter+j+nbSamples][1]
                                counter += 1   
            except FileNotFoundError:
                print("Cannot find 'Sample_settings.txt' file in ./data folder")
            except AttributeError:
                print("Cannot find 'Sample_settings.txt' file in ./data folder")
                
            # Read predictions
            try:
                fileName = os.path.join(conf.PecubeFolderPath,folderName,"data",data_folder_name,data_folder_name+".csv")
                dataObs = pd.DataFrame(pd.read_csv(fileName,skipinitialspace=True,encoding=conf.encoding_label))
                #Get column names, aka thermochronometers
                Headers = dataObs.columns
                #First get sample names, and coordinates
                for i in range(len(Headers)):
                    if Headers[i] == "SAMPLE" or Headers[i] == "LAT" or Headers[i] == "LON" or Headers[i] == "HEIGHT":
                        continue
                    else:
                        for j in range(len(dataObs[Headers[0]])):
                            sampleName = dataObs["SAMPLE"][j]
                            data = dataObs[Headers[i]][j].astype(float)
                            if np.isnan(data):
                                data = 0.0
                            self.ParametersInput.TableParameters[Headers[i]+"_"+sampleName] = str(data)
                    # self.ParametersInput.TableParameters['lon_grain'+str(j+1)] = dataObs["LON"][j].astype(float)
                    # self.ParametersInput.TableParameters['lat_grain'+str(j+1)] = dataObs["LAT"][j].astype(float)
                    # self.ParametersInput.TableParameters['Elev_grain'+str(j+1)] = dataObs["HEIGHT"][j].astype(float)
                # #Get index column of observed ages
                # indexCol = [Headers.get_loc(Headers[i]) for i in range(len(Headers)) 
                #             if Headers[i].startswith('AHE')
                #             or Headers[i].startswith('AFT')
                #             or Headers[i].startswith('ZHE')
                #             or Headers[i].startswith('ZFT')
                #             or Headers[i].startswith('KAR')
                #             or Headers[i].startswith('BAR')
                #             or Headers[i].startswith('MAR')
                #             or Headers[i].startswith('HAR')
                #             or Headers[i].startswith('SIZE')
                #             or Headers[i].startswith('UPPM')
                #             or Headers[i].startswith('THPPM')]

                # #Get observed ages
                # for i in range(len(indexCol)):
                #     for j in range(len(dataObs[Headers[i]])):
                #         value = dataObs[Headers[indexCol[i]]][j]
                #         sampleName = dataObs["SAMPLE"][j]
                #         if pd.isna(value):
                #             value = 0.0
                #         #Observed age
                #         self.ParametersInput.TableParameters['Age'+Headers[indexCol[i]]+'_'+sampleName] = str(value)
                #         #Observed Error
                #         value = dataObs[Headers[indexCol[i]+1]][j]
                #         if pd.isna(value):
                #             value = 0.0
                #         self.ParametersInput.TableParameters['ErrorAge'+Headers[indexCol[i]+1]+'_'+sampleName] = str(value)
                       
                       # self.ParametersInput.TableParameters['grainsize'+Headers[indexCol[i]+1]+'_'+sampleName] = str(value)
                       # self.ParametersInput.TableParameters['U_grainApatite'+Headers[indexCol[i]+1]+'_'+sampleName] = str(value)
                       # self.ParametersInput.TableParameters['Th_grainApatite'+Headers[indexCol[i]+1]+'_'+sampleName]= str(value)
            except FileNotFoundError:
                QErrorMessage(self).showMessage("Cannot find: './data/"+os.path.join(conf.PecubeFolderPath,folderName,"data",data_folder_name,data_folder_name+".csv"))
                exit()
                
            # Is there any Thermoluminescence data ? File is ThL_data.csv
            try:
                fileName = os.path.join(conf.PecubeFolderPath,folderName,"data",data_folder_name,"ThL_data.csv")
                dataObs = pd.DataFrame(pd.read_csv(fileName,skipinitialspace=True,encoding=conf.encoding_label))
                #Get column names, aka thermochronometers
                Headers = dataObs.columns
                #First get sample names, and coordinates
                for i in range(len(Headers)):
                    for j in range(len(dataObs[Headers[0]])):
                        sampleName = dataObs["SAMPLETL"][j]
                        if Headers[i] == "SAMPLETL":
                            self.ParametersInput.TableParameters["TL_"+sampleName] = dataObs["SAMPLETL"][j]
                        else:
                            data = dataObs[Headers[i]][j].astype(float)
                            if np.isnan(data):
                                data = 0.0
                            self.ParametersInput.TableParameters[sampleName+"_"+Headers[i]+"_TL"] = str(data)
            except FileNotFoundError:
                print("No ThL_data.csv file")
           
            # Is there any OSL data ? File is OSL_data.csv
            try:
                fileName = os.path.join(conf.PecubeFolderPath,folderName,"data",data_folder_name,"OSL_data.csv")
                dataObs = pd.DataFrame(pd.read_csv(fileName,skipinitialspace=True,encoding=conf.encoding_label))
                #Get column names, aka thermochronometers
                Headers = dataObs.columns
                #First get sample names, and coordinates
                for i in range(len(Headers)):
                    for j in range(len(dataObs[Headers[0]])):
                        sampleName = dataObs["SAMPLEOSL"][j]
                        if Headers[i] == "SAMPLEOSL":
                            self.ParametersInput.TableParameters["OSL_"+sampleName] = dataObs["SAMPLEOSL"][j]
                        else:
                            data = dataObs[Headers[i]][j].astype(float)
                            if np.isnan(data):
                                data = 0.0
                            self.ParametersInput.TableParameters[sampleName+"_"+Headers[i]+"_OSL"] = str(data)
            except FileNotFoundError:
                print("No OSL_data.csv file")
            
            # Is there any ESR data ? File is ESR_data.csv
            try:
                fileName = os.path.join(conf.PecubeFolderPath,folderName,"data",data_folder_name,"ESR_data.csv")
                dataObs = pd.DataFrame(pd.read_csv(fileName,skipinitialspace=True,encoding=conf.encoding_label))
                #Get column names, aka thermochronometers
                Headers = dataObs.columns
                #First get sample names, and coordinates
                for i in range(len(Headers)):
                    for j in range(len(dataObs[Headers[0]])):
                        sampleName = dataObs["SAMPLEESR"][j]
                        if Headers[i] == "SAMPLEESR":
                            self.ParametersInput.TableParameters["ESR_"+sampleName] = dataObs["SAMPLEESR"][j]
                        else:
                            data = dataObs[Headers[i]][j].astype(float)
                            if np.isnan(data):
                                data = 0.0
                            self.ParametersInput.TableParameters[sampleName+"_"+Headers[i]+"_ESR"] = str(data)
            except FileNotFoundError:
                print("No ESR_data.csv file")
        self.ParametersInput.DParameters['oldInput'] = '1'
  
 
###################################################################
class TableRead:
    """
      This table will show up when loading csv file to plot 2D data
      It updates according to the pandas dataframe provided.
      Used when reading Cooling_rates.csv file. This file is big so better
      to load only when user asks.
      
    """
    def __init__(self,df,HeadersLabel):
        super().__init__()
        self.df = df
        self.Table = QTableWidget()
        #set table dimension
        nRows, nColumns = self.df.shape
        self.Table.setColumnCount(nColumns)
        self.Table.setRowCount(nRows)
        
        self.Table.setHorizontalHeaderLabels(HeadersLabel)
        #Insert data
        for i in range(self.Table.rowCount()):
            for j in range(self.Table.columnCount()):
                self.Table.setItem(i,j, QTableWidgetItem(str(self.df.iloc[i,j])))
                
        self.Table.cellChanged[int, int].connect(self.updateDF)
    
    #----------------------------------------------------------
    def updateDF(self, row, column):
        text = self.Table.item(row,column).text()
        self.df.iloc[row, column] = text
        

#################################################################
class StandardItem(QStandardItem):
    """ 
      Defines fonts for the items of the tree View.
      
    """
    def __init__(self, txt='', font_size=12, set_bold=False,set_checked=True):
        super().__init__()

        fnt = QFont('Open Sans', font_size)
        fnt.setBold(set_bold)

        self.setEditable(False)
        self.setFont(fnt)
        self.setText(txt)
        if set_checked == True:
            self.setCheckState(Qt.Checked)
        else:
            self.setCheckState(Qt.Unchecked)
        self.setCheckable(True)
        
#############################################################
# Classes related to graphical outputs section
#############################################################

###################################################################
class Canvas3D:
    """ To create the area where to draw the 3D model. """
    
    def __init__(self):
        super().__init__()
        self.Frame = QFrame()
        vbox = QVBoxLayout()
        self.plotter = pvqt.QtInteractor(self.Frame)
        self.plotter.enable_3_lights()
        self.plotter.set_background(conf.Colors3Dplot['background'])
        vbox.addWidget(self.plotter.interactor)
        
        # pvqt.MainWindow.signal_close.connect(self.plotter.close)
        
        self.Frame.setLayout(vbox)

###################################################################
class toolbar_3Dplot(QWidget):
    """
    This class defines a custom toolbar for the 3D plots
    parent = GraphWin
    
    """
    def __init__(self, parent, MainWin, Object3D):
        super().__init__()
        
        # ------ Create toolbar ---------------
        self.toolbar = QToolBar(self)
        #View plan xy
        planXY = QAction(QIcon(MainWin.IconePlanxy),"View xy",parent)
        planXY.setShortcut('Ctrl+Z')
        planXY.triggered.connect(lambda: Object3D.plotter.view_xy())
        #View plan xz
        planXZ = QAction(QIcon(MainWin.IconePlanxz),"View xz",parent)
        planXZ.setShortcut('Ctrl+X')
        planXZ.triggered.connect(lambda: Object3D.plotter.view_xz())
        #View plan xy
        planYZ = QAction(QIcon(MainWin.IconePlanyz),"View yz",parent)
        planYZ.setShortcut('Ctrl+Y')
        planYZ.triggered.connect(lambda: Object3D.plotter.view_yz())
        # Add to Toolbar
        self.toolbar.addAction(planXY)
        self.toolbar.addAction(planXZ)
        self.toolbar.addAction(planYZ)
        
        

#############################################################
class Properties3D(QWidget):
    """
    This class set the properties tab of the 3D plots area.
    Parent is GraphWin.
    
    """
    
    def __init__(self,parent,MainWin,Object3D,MainMesh,SampleLoc,grid,activePlot,sargs,poly):
        super().__init__()
        
        self.Object3D = Object3D
        self.Samples_loc = SampleLoc
        self.grid = grid
        self.activePlot = activePlot
        self.MainMesh = MainMesh
        self.sargs = sargs
        self.poly = poly
        self.parent = parent
        
        # Create widget
        self.Dock3dProperties = QDockWidget("Properties", self)
        
        #Colorbar range
        CbRange = QLabel("data range:")
        CbRange.setFont(conf.fontBold8)
        CbRange.setAlignment(Qt.AlignLeft)
        #Reset camera position
        ResCamButton = QPushButton("Reset camera position")
        ResCamButton.clicked.connect(
            lambda: self.Object3D.plotter.reset_camera(True))
        #Clear plots
        # ClearPlotButton = QPushButton("Clear plots")
        # ClearPlotButton.clicked.connect(
        #     lambda: self.Object3D.plotter.clear())
        #Save a screenshot
        ScreenshotButton = QPushButton("Export image...")
        ScreenshotButton.clicked.connect(
            lambda: self.SaveImage())
        # Add colormap selection
        self.ColormapLabel = QLabel("Colormap:")
        self.ColormapLabel.setFont(conf.font11)
        self.ColormapCombo = QComboBox()
        self.ColormapCombo.addItems(matplotlib.colormaps())
        self.ColormapCombo.setCurrentText('RdBu')
        #Show bounds
        self.BoundsCheck = QCheckBox("show box")
        self.BoundsCheck.setChecked(False)
        self.BoundsCheck.stateChanged.connect(self.Show_bounds)
        #Show Edges
        self.EdgesCheck = QCheckBox('show edges')
        self.EdgesCheck.setChecked(True)
        self.EdgesCheck.stateChanged.connect(self.Show_edges)
        #Show Sample locations
        self.SLocCheck = QCheckBox('show sample location(s)')
        self.SLocCheck.setEnabled(False)
        self.SLocCheck.stateChanged.connect(lambda:self.Show_SampleLoc(self.SLocCheck))
        #Set value range for colorbar
        self.MinVal3D = QLineEdit('0')
        self.MaxVal3D = QLineEdit('1')
        self.MinVal3D.editingFinished.connect(self.SetValueRange)
        self.MaxVal3D.editingFinished.connect(self.SetValueRange)
        #Select data from vtk file
        self.Select3dData = QComboBox()
        #Set scales of the model 
        self.XScale = QLabel("X scale:")
        self.XScale.setFont(conf.fontBold8)
        self.XScale.setAlignment(Qt.AlignCenter)
        self.YScale = QLabel("Y scale:")
        self.YScale.setFont(conf.fontBold8)
        self.YScale.setAlignment(Qt.AlignCenter)
        self.ZScale = QLabel("Z scale:")
        self.ZScale.setFont(conf.fontBold8)
        self.ZScale.setAlignment(Qt.AlignCenter)
        self.XScaleEdit = QLineEdit()
        self.XScaleEdit.setValidator(conf.DoubleValidator)
        self.XScaleEdit.setText('1')
        self.XScaleEdit.editingFinished.connect(lambda:self.editScale())
        self.YScaleEdit = QLineEdit()
        self.YScaleEdit.setValidator(conf.DoubleValidator)
        self.YScaleEdit.setText('1')
        self.YScaleEdit.editingFinished.connect(lambda:self.editScale())
        self.ZScaleEdit = QLineEdit()
        self.ZScaleEdit.setValidator(conf.DoubleValidator)
        self.ZScaleEdit.setText('1')
        self.ZScaleEdit.editingFinished.connect(lambda:self.editScale())

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
        #Current data
        CurData = QLabel("current data:")
        CurData.setFont(conf.fontBold8)
        CurData.setAlignment(Qt.AlignLeft)
        
        ####### Build the grid for 3D properties #######
        VBoxMain = QVBoxLayout()
        #Data Properties
        BoxData = QGroupBox("Data properties")
        GboxDataProp = QGridLayout()
        GboxDataProp.setRowStretch(8,1)
        # GboxDataProp.addWidget(DataProp,0,0)
        GboxDataProp.addWidget(CbRange,2,0)
        GboxDataProp.addWidget(self.MinVal3D,3,0)
        GboxDataProp.addWidget(self.MaxVal3D,3,1)
        GboxDataProp.addWidget(CurData,5,0)
        GboxDataProp.addWidget(self.Select3dData,6,0)
        hbox = QGridLayout()
        hbox.addWidget(self.XScale,0,0)
        hbox.addWidget(self.YScale,0,1)
        hbox.addWidget(self.ZScale,0,2)
        hbox.addWidget(self.XScaleEdit,1,0)
        hbox.addWidget(self.YScaleEdit,1,1)
        hbox.addWidget(self.ZScaleEdit,1,2)
        hbox.setSpacing(10)
        GboxDataProp.addLayout(hbox,7,0,1,3)
        GboxDataProp.setSpacing(10)
        BoxData.setLayout(GboxDataProp)
        VBoxMain.addWidget(BoxData)
        #Plot properties
        BoxPlot = QGroupBox("Plot properties")
        GboxPlotProp = QGridLayout()
        GboxPlotProp.addWidget(ResCamButton,1,0,1,2)
        # GboxPlotProp.addWidget(ClearPlotButton,2,0)
        GboxPlotProp.addWidget(ScreenshotButton,3,0,1,2)
        GboxPlotProp.addWidget(self.ColormapLabel,4,0)
        GboxPlotProp.addWidget(self.ColormapCombo,4,1)
        GboxPlotProp.addWidget(self.BoundsCheck,5,0,1,2)
        GboxPlotProp.addWidget(self.EdgesCheck,6,0,1,2)
        GboxPlotProp.addWidget(self.SLocCheck,7,0,1,2)
        GboxPlotProp.setRowStretch(8,1)
        GboxPlotProp.setSpacing(10)
        BoxPlot.setLayout(GboxPlotProp)
        # VBoxMain.addStretch(1)
        VBoxMain.addWidget(BoxPlot)
        widgetDockProp = QWidget()
        widgetDockProp.setLayout(VBoxMain)
        self.Dock3dProperties.setWidget(widgetDockProp)
        self.Dock3dProperties.setFixedWidth(300)
        
        # Signals
        self.ColormapCombo.currentIndexChanged.connect(lambda: self.update_colormap())
    
    #----------------------------------------------------------
    def update_colormap(self):
        """ To update colormap of the 3D model """
        
        if self.activePlot and self.Object3D.plotter and self.MainMesh:
            # Remove the mesh
            self.Object3D.plotter.remove_actor(self.MainMesh)
            # Rebuild the mesh with updated colormap
            self.MainMesh= self.Object3D.plotter.add_mesh(self.grid,
                                                          show_edges=False, cmap=self.ColormapCombo.currentText(),
                                                          flip_scalars=True, scalar_bar_args=self.sargs,
                                                          clim=[float(self.MinVal3D.text()),float(self.MaxVal3D.text())])
        
        
    #----------------------------------------------------------
    def Show_bounds(self):
        """ To show grid axes on the 3D model """
        if self.BoundsCheck.checkState() == Qt.Checked:
            self.Object3D.plotter.show_bounds(grid='front',location='outer',all_edges=True)
        else:
            self.Object3D.plotter.remove_bounds_axes()
            self.Object3D.plotter.remove_bounding_box()
    
    #----------------------------------------------------------
    def Show_edges(self):
        """ To show edges of the box on the 3D model """
        if self.EdgesCheck.checkState() == Qt.Checked:
            self.Object3D.plotter.show_edges = True
        else:
            self.Object3D.plotter.show_edges = False
            
    #----------------------------------------------------------
    def SetValueRange(self):
        """ 
          To update color bar range
          get min and max value set by the user`
          
        """
        minval = float(self.MinVal3D.text())
        maxval = float(self.MaxVal3D.text())
        if self.activePlot and self.Object3D.plotter and self.MainMesh:
            # Remove the mesh
            try: 
                self.Object3D.plotter.remove_actor(self.MainMesh)
                # Rebuild the mesh with updated color bar
                self.MainMesh= self.Object3D.plotter.add_mesh(
                self.grid,show_edges=False, cmap=self.ColormapCombo.currentText(), flip_scalars=True, scalar_bar_args=self.sargs, clim=[minval,maxval])
            except:
                print("Fail when try to remove actor. PGUI_utils/SetValueRange")
            
            # info: on the pyvista version 0.31 the maximum value of scalar bar can
            # only increases. To circumvent the issue, need to go to .../pyvista/plotting/scalar_bars.py
            # and change line 330 the operator from > to !=...
    
    #----------------------------------------------------------
    def Show_SampleLoc(self,box):
        """ To show/hide samples locations """
        if box.isChecked() and not self.Samples_loc:
            # self.Samples_loc.SetVisibility(True)
            self.Samples_loc = self.Object3D.plotter.add_point_labels(self.poly,"My Labels",point_size=10,
                                                    font_size=10,name='Sample_location')
        elif self.Samples_loc and not box.isChecked():
            self.Samples_loc.SetVisibility(False)
            self.Object3D.plotter.remove_actor('Sample_location')
            self.Samples_loc = 0
    
    #----------------------------------------------------------
    def editScale(self):
        """ To scale the 3D model """
        xs = float(self.XScaleEdit.text())
        ys = float(self.YScaleEdit.text())
        zs = float(self.ZScaleEdit.text())
        self.Object3D.plotter.set_scale(xs,ys,zs)
    
    #----------------------------------------------------------
    def DataChanged(self):
        """ To change the data to show """
        dataToPlot = self.Select3dData.currentText()
        active_scalar = self.grid.active_scalars_name
        # Remove mesh from plottter
        self.Object3D.plotter.remove_actor(self.MainMesh)
        # Remove the current color bar
        self.Object3D.plotter.remove_scalar_bar(title=active_scalar)
        # Plot the new mesh
        self.MainMesh= self.Object3D.plotter.add_mesh(
        self.grid,show_edges=False, cmap=self.ColormapCombo.currentText(), scalars=dataToPlot, flip_scalars=True, scalar_bar_args=self.sargs)
        # Get min and max values
        DataRange = self.grid.get_data_range()
        # Set min and max value in properties tab
        self.MinVal3D.setText(str(DataRange[0]))
        self.MaxVal3D.setText(str(DataRange[1]))
        
    #----------------------------------------------------------
    def SaveImage(self):
        """ To export screenshot of the 3D model """
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save Image",
                          r"H:\Image.pdf","Image files (*.svg *.eps *.ps *.pdf *.tex)",options=options)
        if fileName:
            self.Object3D.plotter.save_graphic(fileName,painter=True)
     
            
###################################################################
class PlotCanvas(FigureCanvas):
    """
      Set plot area object (for 1D and 2D plots).
      
    """
    def __init__(self, parent=None, sp = 111, width=5, height=4, dpi=100,Axe3=None):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        if Axe3:
            self.axes = self.fig.add_axes([0,0,1,1],projection="3d")
        else:
            self.axes = self.fig.add_subplot(sp)
        plt.ion()
        #plt.rcParams['toolbar'] = 'toolmanager'
        # FigureCanvas.__init__(self, fig)
        # # self.setParent(parent)

        # FigureCanvas.setSizePolicy(self,
        #         QSizePolicy.Expanding,
        #         QSizePolicy.Expanding)
        # FigureCanvas.updateGeometry(self)
        super(PlotCanvas, self).__init__(self.fig)
    
    #----------------------------------------------------------
    def resizeEvent(self,event):
        """
          To avoid the plotting area to reach zero size when the user is
          dragging the edge of the splitter 
        """
        if event.size().width() <= 0 or event.size().height() <= 0:
            return
        super(PlotCanvas, self).resizeEvent(event)

###################################################################

    