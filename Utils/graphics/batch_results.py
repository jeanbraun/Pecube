#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to plot output results from batch mode within the interface. 

@author: maxime Bernard

"""

import os
from PyQt5.QtWidgets import (QWidget,QGridLayout,
                             QVBoxLayout,QLabel,QComboBox,
                             QCheckBox)
from PyQt5.QtCore import Qt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import configs as conf
import Thermochronology.thermochronometers as th
import Thermochronology.settings as Thermo_settings
import Utils.PGUI_utils as pgu



#############################################################      
class AgeBatch(QWidget):
    """ This class handle the plot when Pecube is used in batch mode.
    Batch mode means that the user provides a range of value for one parameter as inversion mode.
    But instead of using NA algorithm the range of value is divided in N (sample_size_for_first_iteration)
    intervals regularly spaced. This is useful to explore the influence of one parameters on e.g. age-elevation 
    relationships."""
    
    def __init__(self,ProjectFolder,GraphWin):
        super().__init__()
        # GraphinWin is all parameters defined in GraphWin class (see PecubeGUI.py)
        # ProjectFolder is path to Project directory of the current Pecube project
        
        #### Parameters ####
        self.ProjectFolder = ProjectFolder
        self.GraphWin = GraphWin
        self.NAPath = os.path.join(self.ProjectFolder,"NA")
        self.dataFile = self.GraphWin.ParametersInput.DParameters['data_folder']+'.csv'
        self.DataFolder = os.path.join(self.ProjectFolder,"data",self.GraphWin.ParametersInput.DParameters['data_folder'])
        
        # Set parameters to ask user what to plot
        # Select chart type
        self.chartTypeLabel = QLabel("chart type: ")
        self.chartTypeLabel.setFont(conf.font10)
        self.chartTypeLabel.setAlignment(Qt.AlignRight)
        self.chartTypeLabel.setToolTip("Choose the type of chart you want to plot.")
        items = ["age-elevation"]
        self.chartCombo = QComboBox()
        self.chartCombo.addItems(items)
        self.chartCombo.setMinimumContentsLength(15)
        self.chartType = self.chartCombo.currentText()
        # Select which thermochronometer
        self.SelectThermochronometer = QLabel('Select thermochronometer:')
        self.SelectThermochronometer.setFont(conf.font10)
        self.SelectThermochronometer.setAlignment(Qt.AlignRight)
        self.ThermochronometerCombo = QComboBox()
        self.ThermochronometerCombo.setMinimumContentsLength(8)
        self.set_ThermoCombo()
        self.ThermochronometerFile = []
        self.ThermochronometerID = []
        self.update_ThermoCombo()
        self.Normalize_y_label = QLabel('Normalize y: ')
        self.Normalize_y_label.setFont(conf.font10)
        self.Normalize_y = QCheckBox()
       
        # Grid
        self.Layout = QVBoxLayout()
        Grid1 = QGridLayout()
        Grid1.addWidget(self.chartTypeLabel,0, 1)
        Grid1.addWidget(self.chartCombo,0, 2)
        Grid2 = QGridLayout()
        Grid2.addWidget(self.SelectThermochronometer,0, 1)
        Grid2.addWidget(self.ThermochronometerCombo,0, 2)
        Grid2.addWidget(self.Normalize_y_label,1,1)
        Grid2.addWidget(self.Normalize_y,1,2)
        
        self.Layout.addLayout(Grid1)
        self.Layout.addLayout(Grid2)
        
        # Signals
        self.ThermochronometerCombo.currentIndexChanged.connect(lambda: self.update_ThermoCombo())
        

    #------------------------------------------------------------------
    def plot_Batch_results(self):
        """ Plot the results from batch models according to the choice of the user."""
        
        ###### First read the data #######
        OutputFileName = os.path.join(self.NAPath,self.ThermochronometerFile)
        # Read file where ages are stored
        outputData = pd.DataFrame(pd.read_csv(OutputFileName,sep=',',header=None,index_col=False,
                                    encoding=conf.encoding_label))
        Elevation_file = os.path.join(self.NAPath,'ElevationPredicted.csv')
        # Read input File
        if self.ThermochronometerID == 'ESR':
            inputFileName = os.path.join(self.DataFolder,self.DataFolder,conf.Variable_names["ESR_file_input"])
        elif self.ThermochronometerID == 'TL':
            inputFileName = os.path.join(self.DataFolder,self.DataFolder,conf.Variable_names["ThL_file_input"])
        elif self.ThermochronometerID == 'OSL':
            inputFileName = os.path.join(self.DataFolder,self.DataFolder,conf.Variable_names["OSL_file_input"])
        else:
            inputFileName = os.path.join(self.DataFolder,self.DataFolder,self.dataFile)
            
        inputData = pd.DataFrame(pd.read_csv(inputFileName,sep=',',index_col=False,
                                    encoding=conf.encoding_label))
        # ElevationData = pd.DataFrame(pd.read_csv(Elevation_file,sep=',',header=None,index_col=False,
        #                             encoding=conf.encoding_label))
        # Read model parameter values in NA_int_results.csv file
        NAFileName = os.path.join(self.NAPath,conf.Variable_names['File int inversion'])
        ModelParamValues = pd.DataFrame(pd.read_csv(NAFileName,sep=',',index_col=False,
                                    encoding=conf.encoding_label))
        
        agecol, errname, agename, predname, colores, profdict, He43dict,coloramp, ageMarker  =  Thermo_settings.dict_pecube()
        
        ###### Second, sort the data #######
        # Elevations of data (prediction)
        ObservedElevations = inputData['HEIGHT']
        Observed_Data = inputData[self.ThermochronometerID]
        Observed_error = inputData['D'+self.ThermochronometerID]
        # inputElevations = ElevationData.iloc[0:-1,1:]
        Misfit = outputData.iloc[0:-1,0]
        # Output Ages - alternate line elevation and age
        outputAges = outputData.iloc[0:-1,1:] # Remove last parameter value (random value)
        # Number of models
        nmodels = int(outputAges.shape[0]/2)

        # Parameter values (assume a single parameter)
        paramValues = ModelParamValues['Param0001']
        ###### Third, plot the data #######
        # Create the canva
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        # Create a customed toolbar to show/hide ages from specific thermochronometer
        self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
        # plot age-elevation profiles
        if self.chartType == 'age-elevation':
            
            cmap = plt.get_cmap('magma')
            print('nmodels : ', nmodels)
            models_IDs = np.linspace(0,1,nmodels)
            colors = cmap(models_IDs)
            for i in range(nmodels):
                if self.Normalize_y.isChecked() == True: # Normalize elevation
                    minElev = np.min(outputAges.iloc[i*2])
                    maxElev = np.max(outputAges.iloc[i*2])
                    Elevations = outputAges.iloc[i*2].values
                    Elevations = (Elevations + abs(minElev)) / (maxElev+abs(minElev))
                else:
                    Elevations = outputAges.iloc[i*2].values
                # sort data by elevations
                indices = np.argsort(Elevations)
                elev = Elevations[indices[:]]
                
                ages = outputAges.iloc[i*2+1].values; ages = ages[indices]
                self.plotSpace.axes.plot(ages,elev/1e3,'-',marker = ageMarker[self.ThermochronometerID],
                                         color = colors[i], label='Param = '+str(round(paramValues[i],3)) + ', Misfit = '+str(round(Misfit[i],2)),
                                         markeredgecolor='black')
            # Plot observations
            self.plotSpace.axes.errorbar(Observed_Data, ObservedElevations/1e3,xerr = Observed_error,
                      fmt = ageMarker[self.ThermochronometerID], label = 'observed_data', color = 'black',
                      markeredgecolor='black',markerfacecolor = 'white',zorder=9)
    
            self.plotSpace.axes.set_xlabel('Age (Ma)')
            self.plotSpace.axes.set_ylabel('Elevation (km)')
            self.plotSpace.axes.legend(loc = 'best')
            self.plotSpace.axes.set_box_aspect(1)
            # Add plot to interface
            VBox = QVBoxLayout()
            VBox.addWidget(self.toolbar)
            VBox.addWidget(self.plotSpace)
            Q = QWidget()
            Q.setLayout(VBox)
            self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'Batch_Age_Elevation')
        
        
    #------------------------------------------------------------------
    def set_ThermoCombo(self):
        """Set list of thermochronometer to plot. This depends on the presence of files
        AgeThermo.csv in NA directory."""
        
        filenames = os.listdir(self.NAPath)
        self.listFiles = [file for file in filenames if file.startswith('Age')]
        items = []
        if conf.Variable_names["File Age ESR batch"] in self.listFiles:
            items.append('ESR')
        if conf.Variable_names["File Age AHe batch"] in self.listFiles:
            items.append('AHE')
        if conf.Variable_names["File Age ZHe batch"] in self.listFiles:
            items.append('ZHE')
        if conf.Variable_names["File Age AFT batch"] in self.listFiles:
            items.append('AFT')
        if conf.Variable_names["File Age ZFT batch"] in self.listFiles:
            items.append('ZFT')
        if conf.Variable_names["File Age KAr batch"] in self.listFiles:
            items.append('KAR')
        if conf.Variable_names["File Age MAr batch"] in self.listFiles:
            items.append('MAR')
        if conf.Variable_names["File Age BAr batch"] in self.listFiles:
            items.append('BAR')
        if conf.Variable_names["File Age HAr batch"] in self.listFiles:
            items.append('HAR')
        
        self.ThermochronometerCombo.addItems(items)
    
    #------------------------------------------------------------------
    def update_ThermoCombo(self):
        """Update which thermochronometer to plot."""
        
        self.ThermochronometerID = self.ThermochronometerCombo.currentText()
        if self.ThermochronometerCombo.currentText() == 'ESR':
            self.ThermochronometerFile = conf.Variable_names["File Age ESR batch"]
        elif self.ThermochronometerCombo.currentText() == 'AHE':
            self.ThermochronometerFile = conf.Variable_names["File Age AHe batch"]
        elif self.ThermochronometerCombo.currentText() == 'AFT':
            self.ThermochronometerFile = conf.Variable_names["File Age AFT batch"]
        elif self.ThermochronometerCombo.currentText() == 'ZFT':
            self.ThermochronometerFile = conf.Variable_names["File Age ZFT batch"]
        elif self.ThermochronometerCombo.currentText() == 'ZHE':
            self.ThermochronometerFile = conf.Variable_names["File Age ZHe batch"]
        elif self.ThermochronometerCombo.currentText() == 'KAR':
            self.ThermochronometerFile = conf.Variable_names["File Age KAr batch"]
        elif self.ThermochronometerCombo.currentText() == 'BAR':
            self.ThermochronometerFile = conf.Variable_names["File Age BAr batch"]
        elif self.ThermochronometerCombo.currentText() == 'MAR':
            self.ThermochronometerFile = conf.Variable_names["File Age MAr batch"]
        elif self.ThermochronometerCombo.currentText() == 'HAR':
            self.ThermochronometerFile = conf.Variable_names["File Age HAr batch"]