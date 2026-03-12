#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to plot output results in the interface, related to thermochronometers 
predictions.

@author: maxime Bernard

"""

import os
from PyQt5.QtWidgets import (QErrorMessage, QWidget,QGridLayout,
                             QVBoxLayout,QHBoxLayout,QGroupBox,QLabel,QComboBox,
                             QCheckBox,QLineEdit)
from PyQt5.QtCore import Qt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.collections import LineCollection
import configs as conf
import Thermochronology.thermochronometers as th
import Thermochronology.settings as Thermo_settings
import Utils.interface as U_interface
import Utils.PGUI_utils as pgu
import Utils.graphics.set_plot as set_plot
import Utils.GIS as GIS



################################## GENERAL ###################################
#-----------------------------------------------------------------------------
#--------------------------------- Classes -----------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

class transect(QWidget):
    """ This class handle parameters to plot data along transect. It set
    a window that ask the user to provide information according to the choice 
    of transect : 'Latitude','Longitude','projected'. If 'projected' the coordinates
    of the 2 points that define the transect are asked."""
    
    def __init__(self,inputData,thermochronometers,CompareAge,inputc,agecol,errname,agename,colores,predname,ageMarker,profdict,GraphWin):
        super().__init__()
        
        self.dataplot = thermochronometers
        self.datac = CompareAge
        self.inputc = inputc
        self.agecol = agecol
        self.errname = errname
        self.colores = colores
        self.predname = predname
        self.profdict = profdict
        self.inputdata = inputData
        self.agename = agename
        self.ageMarker = ageMarker
        self.GraphWin = GraphWin
        
        # Load topographic file to plot line
        topofile = self.GraphWin.ParametersInput.DParameters['topo_file_name']
        path = os.path.join(self.GraphWin.PecubePath,self.GraphWin.folderName,"data",topofile+'.tif')
        try:
            elevation,nx,ny,dx,dy,x0,y0 = GIS.load_raster_file(path)
        except:
            QErrorMessage(self).showMessage(topofile+'.tif not found. This works only with tif file (so far). Cannot plot transect.')
            return
        self.Z = np.reshape(elevation,(ny,nx))
        xl = (nx-1)*dx
        yl = (ny-1)*dy
        x = np.linspace(x0,x0+xl,nx)
        y = np.linspace(y0,y0+yl,ny)
        self.X,self.Y = np.meshgrid(x,y)
        self.plotDEM =  pgu.PlotCanvas()
        #plot DEM
        self.plotDEM.axes.pcolormesh(self.X,self.Y,self.Z)
        # plot Data if any
        if self.datac:
            self.plotDEM.axes.scatter(self.datac['LONOBS'],self.datac['LATOBS'],
                                  s = 10, c = 'w')
        self.plotDEM.axes.set_aspect('equal')
        self.toolbarDEM = NavigationToolbar(self.plotDEM,self.GraphWin)
        
        #Set parameters
        self.profileComboLabel = QLabel("Choose transect: ")
        self.profileComboLabel.setFont(conf.font11)
        self.profileCombo = QComboBox()
        self.profileCombo.addItems(['Latitude','Longitude','Projected'])
        self.profileCombo.setMinimumContentsLength(15)
        self.text1 = QLabel("Provide coordinates of the 2 points of the transect: ")
        self.text1.setFont(conf.fontBold11)
        self.text1.setAlignment(Qt.AlignLeft)
        self.pointCoords = []
        self.PointLatLabel = QLabel("Latitude:")
        self.PointLatLabel.setFont(conf.font10)
        self.PointLatLabel.setAlignment(Qt.AlignCenter)
        self.PointLonLabel = QLabel("Longitude:")
        self.PointLonLabel.setFont(conf.font10)
        self.PointLonLabel.setAlignment(Qt.AlignCenter)
        self.Point1Label = QLabel("Point 1:")
        self.Point1Label.setFont(conf.font10)
        self.Point1Label.setAlignment(Qt.AlignCenter)
        self.Point2Label = QLabel("Point 2:")
        self.Point2Label.setFont(conf.font10)
        self.Point2Label.setAlignment(Qt.AlignCenter)
        self.Lat1 = QLineEdit("0.0")
        self.Lat1.setAlignment(Qt.AlignCenter)
        self.Lat2 = QLineEdit("0.0")
        self.Lat2.setAlignment(Qt.AlignCenter)
        self.Lon1 = QLineEdit("0.0")
        self.Lon1.setAlignment(Qt.AlignCenter)
        self.Lon2 = QLineEdit("0.0")
        self.Lon2.setAlignment(Qt.AlignCenter)
        self.Lat1.setEnabled(False)
        self.Lat2.setEnabled(False)
        self.Lon1.setEnabled(False)
        self.Lon2.setEnabled(False)
        
        # Set Layouts
        self.MainLayout = QVBoxLayout()
        hbox = QHBoxLayout()
        hbox.addWidget(self.profileComboLabel)
        hbox.addWidget(self.profileCombo)
        hbox.addStretch(1)
        
        self.Box1 = QGroupBox()
        self.Box1.setCheckable(False)
        # Widget for projected
        self.widget1 = QWidget()
        self.Grid1 = QGridLayout()
        self.Grid1.setSpacing(10)
        self.Grid1.addWidget(self.text1,0,0,1,3)
        self.Grid1.addWidget(self.PointLatLabel,1,1)
        self.Grid1.addWidget(self.PointLonLabel,1,2)
        self.Grid1.addWidget(self.Lat1,2,1)
        self.Grid1.addWidget(self.Lon1,2,2)
        self.Grid1.addWidget(self.Lat2,3,1)
        self.Grid1.addWidget(self.Lon2,3,2)
        self.Grid1.addWidget(self.Point1Label,2,0)
        self.Grid1.addWidget(self.Point2Label,3,0)
        self.Grid1.setColumnStretch(4,1)
        self.Grid1.setRowStretch(6,1)
        Grid2 = QVBoxLayout()
        Grid2.addWidget(self.toolbarDEM)
        Grid2.addWidget(self.plotDEM)
        self.widget1.setLayout(Grid2)
        self.widget1.setVisible(True)
        VBox = QVBoxLayout()
        VBox.addWidget(self.widget1)
        self.Box1.setLayout(VBox)
        self.MainLayout.addLayout(hbox)
        self.MainLayout.addStretch(1)
        self.MainLayout.addLayout(self.Grid1)
        self.MainLayout.addWidget(self.Box1)
        self.MainLayout.addStretch(1)
        
        self.clicked = 0
        # Signal and connections
        self.profileCombo.currentIndexChanged.connect(lambda: self.select_transect())
        self.plotDEM.fig.canvas.mpl_connect('button_press_event',self.onclickDEM)

    #-------------------------------------------------------------------------
    def onclickDEM(self,event):
        """ on click return coordinates of the point a show it on plot.
            'clicked' records number of time the user click to remove points in the plot"""
            
        print(self.clicked)
        if self.clicked > 2: # delete first node in self.pointCoords
            self.pointCoords = []
            self.plotDEM.axes.clear()
            self.plotDEM.axes.pcolormesh(self.X,self.Y,self.Z)
            if self.datac:
                self.plotDEM.axes.scatter(self.datac['LONOBS'],self.datac['LATOBS'],
                                      s = 5, c = 'w')
            self.clicked = 0
        self.pointCoords.append([event.xdata,event.ydata])
        self.clicked += 1
        if len(self.pointCoords) == 2:
            self.scatter = self.plotDEM.axes.scatter(event.xdata,event.ydata,s=20,c='r',edgecolors='k')
            self.line = self.plotDEM.axes.plot([self.pointCoords[0][0],self.pointCoords[1][0]],
                                               [self.pointCoords[0][1],self.pointCoords[1][1]],
                                               'r-')
        else:
            self.scatter = self.plotDEM.axes.scatter(event.xdata,event.ydata,s=20,c='r',edgecolors='k')
        self.plotDEM.fig.canvas.draw()
        self.clicked += 1
        self.Lat1.setText(str(round(self.pointCoords[0][1],4)))
        self.Lon1.setText(str(round(self.pointCoords[0][0],4)))
        self.Lat2.setText(str(round(self.pointCoords[1][1],4)))
        self.Lon2.setText(str(round(self.pointCoords[1][0],4)))
        
        return event.xdata, event.ydata
    
    #-------------------------------------------------------------------------
    def select_transect(self):
        """To select the transect and set the parameters to provide accordingly"""
        
        if self.profileCombo.currentIndex() == 0: # Latitude
            self.widget1.setVisible(True)
            self.Lat1.setEnabled(False)
            self.Lat2.setEnabled(False)
            self.Lon1.setEnabled(False)
            self.Lon2.setEnabled(False)
        elif self.profileCombo.currentIndex() == 1: # Longitude
            self.widget1.setVisible(True)
            self.Lat1.setEnabled(False)
            self.Lat2.setEnabled(False)
            self.Lon1.setEnabled(False)
            self.Lon2.setEnabled(False)
        elif self.profileCombo.currentIndex() == 2: # Projected
            self.widget1.setVisible(True)
            self.Lat1.setEnabled(True)
            self.Lat2.setEnabled(True)
            self.Lon1.setEnabled(True)
            self.Lon2.setEnabled(True)
            
    #-------------------------------------------------------------------------
    def plotTransect(self):
        """
        Author: Xavier Robert, university of Grenoble
        Adapted by: Maxime Bernard, university of Potsdam (12/05/2023)
        
    	Args:
    		inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars
    
    		profiletype (list, optional): type of profile = ['Latitude', 'Longitude', 'Altitude', 'Projected']
    						              If [], no age profile is plotted.
    
    		dataplot (list): List of data to plot ; By default, the altitude will be plotted
    						 ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
    						 Rem: For the moment, MTL and TTp not implemented.
    
    		datac (string): results of Pecube forward modeling
    
    		inputc (array, optional): Array of data to plot.
    		                          Defaults to None.
    
    		coordproj (array, optional): Predictions projected along A-B.
    		                             Defaults to None.
    		
    		coordprojinputc (array, optional): Observations projected along A-B.
    		                                   Defaults to None.
    
    		agerange (2*1 array of floats, optional): range of the ages to plot on the profiles
    												  [min, max]. Defaults to None.
    
    		agecol (dict, optional): respective column number of the data in the file comparison.txt.
    		                         Defaults to None.
    		
    		errname (dict, optional): respective column number of the error on data in the data input file. 
    		                          Defaults to None.
    		
    		colores (dict, optional): Colors used for the different age system. 
    		                          Defaults to None.
    				
    		predname (dict, optional): legend of each predicted system. 
    		                           Defaults to None.
    
    		agename (dict, optional): legend of each data system. 
    		                          Defaults to None.
    
    		profdict (dict, optional): Columns of profile types. 
    		                           Defaults to None.
    
    		graphpath (str, optional): name of the folder where the plot will be written
    						            Usually you do not have to change it. 
    									Defaults to 'Graphs'.
    
    		graphtitle (str, optional): title to write on the graph.
    		                            Defaults to None.
    
    	"""
    
    	# Dictionnary to build axis' legends
        profname = {'Latitude'  : 'Latitude (°)', 
    	            'Longitude' : 'Longitude (°)', 
    				'Altitude'  : 'Elevation (m)', 
    				'Projected' : 'Distance along profile (km)'}
        # Get profile type
        profiletype = [self.profileCombo.currentText()]
        if 'Projected' in profiletype:
            A = [float(self.Lon1.text()), float(self.Lat1.text())]
            B = [float(self.Lon2.text()), float(self.Lat2.text())]
        # if two coordinates are provided
        if 'Projected' in profiletype and A and B:
			# plot Projected transect
			# Compute projected coordinates
            coordproj = self.project(A, B, self.datac,1)
            if self.inputdata: coordprojinputc = self.project(A, B, self.datac,1)
            
    	# fig1 = plt.figure()
        # Create the canva
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        # Create a customed toolbar to show/hide ages from specific thermochronometer
        self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
        
        for profile in profiletype:
            print(u'\tPlotting %s transect' %(profile))
            for item in self.dataplot:
                #Filter age predictions, we want ages > -9999
                sampleName = self.datac["SID"][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                sampleNamePred = self.inputc['SAMPLE'][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                # sort data by sample (alphabetically)
                sampleNametemp = [s for s in sampleName.values if pgu.isNaN(s) == False]
                sampleNamePredtemp = [s for s in sampleNamePred.values if pgu.isNaN(s) == False]
                sampleNameIndexes = np.argsort(sampleNametemp)
                sampleNamePredIndexes = np.argsort(sampleNamePredtemp)
                errortemp = self.inputc[self.errname[item]][np.logical_not(self.inputc[self.errname[item]] == -9999)]
                errortemp2 = np.asarray([s for s in errortemp.values if pgu.isNaN(s) == False])
                errortemp = errortemp2[sampleNamePredIndexes]
                obs = self.datac[self.agecol[item]+'OBS'][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                obs = obs[sampleNameIndexes]
                pred = self.datac[self.agecol[item]+'PRED'][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                pred = pred[sampleNameIndexes]
                # Now we have the right uncertainties we need to sort the data, meaning we drop nan values at the end
                # of the array
                error = np.empty(pred.shape)
                error[:] = np.nan
                count = 0
                for n in range(len(pred)):
                    if not pgu.isNaN(errortemp[n]):
                        error[count] = errortemp[n]
                        count += 1
                
                if item.casefold() != 'TTp'.casefold() and item.upper() != 'MTL':
                    item = item.upper()
                    if self.inputdata:
    					# Plot input data with Error bars
                        if profile != 'Projected' and profile != 'Altitude':
                            location = self.datac[self.profdict[profile]+'OBS'][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                            location = location[sampleNameIndexes]
                            # plot error bars
                            segments = [
                                [(xi, yi - err), (xi, yi + err)]
                                for xi, yi, err in zip(location, obs, error)
                            ]

                            # Create LineCollection for all error bars
                            error_lines = LineCollection(segments, colors=self.colores[item], linewidths=1)
                            self.plotSpace.axes.add_collection(error_lines)
                            self.plotSpace.axes.plot(location, obs, 
    		    			     marker = self.ageMarker[item], linestyle = 'None', 
    							 label = self.agename[item], color = self.colores[item],zorder = 9)
                        else:
                            coordprojinputc2 = coordprojinputc.values[np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                            coordprojinputc2 = coordprojinputc2[sampleNameIndexes]
                            # plot error bars
                            segments = [
                                [(xi, yi - err), (xi, yi + err)]
                                for xi, yi, err in zip(coordprojinputc2, obs, error)
                            ]

                            # Create LineCollection for all error bars
                            error_lines = LineCollection(segments, colors=self.colores[item], linewidths=1)
                            self.plotSpace.axes.add_collection(error_lines)
                            self.plotSpace.axes.plot(coordprojinputc2, obs, marker = self.ageMarker[item],
    			    		          linestyle = 'None',  label = self.agename[item], color = self.colores[item],zorder=9)
                    else:
    					# Plot input data with NO error bars
                        if profile != 'Projected' and profile != 'Altitude':
                            location = self.datac[self.profdict[profile]+'OBS'][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                            location = location[sampleNameIndexes]
                            self.plotSpace.axes.plot(self.datac[self.profdict[profile]], self.datac[self.agecol[item]+'OBS'], 
    			    			     marker = self.ageMarker[item], linestyle = 'None', 
    								 label = self.agename[item], color = self.colores[item],zorder=9)
                        else:
                            coordproj2 = coordproj.values[np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                            coordprojinputc2 = coordprojinputc2[sampleNameIndexes]
                            self.plotSpace.axes.plot(coordproj2, obs, 
    			    			     marker = self.ageMarker[item], linestyle = 'None', 
    								 label = self.agename[item], color = self.colores[item],zorder=9)
    				# Plot predictions
                    if profile != 'Projected' and profile != 'Altitude':
                        location = self.datac[self.profdict[profile]+'OBS'][np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                        location = location[sampleNameIndexes]
                        self.plotSpace.axes.plot(location, pred, 
    		    			     marker = self.ageMarker[item], linestyle = 'None', 
    							 label = self.predname[item], color = self.colores[item], alpha = 0.3,zorder=10)

                    else:
                        coordproj2 = coordproj.values[np.logical_not(self.datac[self.agecol[item]+'OBS'] == -9999)]
                        coordproj2 = coordproj2[sampleNameIndexes]
                        self.plotSpace.axes.plot(coordproj2, pred, 
    		    			     marker = self.ageMarker[item], linestyle = 'None', 
    							 label = self.predname[item], color = self.colores[item], alpha = 0.3,zorder=10)
  
    
    		# set legend
            self.plotSpace.axes.legend(loc = 'best')
    		#plt.xlabel(u'Longitude (°)')
            if profile != 'Altitude':
                self.plotSpace.axes.set_xlabel(profname[profile])
                self.plotSpace.axes.set_ylabel(u'Age (Ma)')
            else:
                self.plotSpace.axes.set_ylabel(profname[profile])
                self.plotSpace.axes.set_xlabel(u'Age (Ma)')
            # plt.title(graphtitle)
            
            # Add plot to interface
            VBox = QVBoxLayout()
            VBox.addWidget(self.toolbar)
            VBox.addWidget(self.plotSpace)
            Q = QWidget()
            Q.setLayout(VBox)
            self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'Age-transect')
    
        return
    
    #-------------------------------------------------------------------------
    def calc_length(self,XA,YA,XB,YB):
        """
        function to calcule the distance on a sphere
    	between two points A and B given by their long/lat coordinates,   

    	Read the page http://gis.stackexchange.com/questions/44064/how-to-calculate-distances-in-a-point-sequence
        
        """

        # Earth radius in km
        R = 6371
    	
        calc_lengthd = R * np.arccos(np.sin(YA * np.pi / 180) * np.sin(YB * np.pi / 180) + \
    	                          np.cos(YA * np.pi / 180) * np.cos(YB * np.pi / 180) * \
    	                          np.cos(-XA * np.pi / 180 + XB * np.pi / 180))
        return calc_lengthd      
    	

    
    #-------------------------------------------------------------------------
    def bearing(self,XA,YA,XB,YB):
        """
    	function to calcule the bearing on a sphere
    	between two points A and B given by their long/lat coordinates,   
        """   
               
        cosdeltaAB = np.sin(YA * np.pi / 180) * np.sin(YB * np.pi / 180) + \
    	                    np.cos(YA * np.pi / 180) * np.cos(YB * np.pi / 180) * \
    	                    np.cos(-XA * np.pi / 180 + XB * np.pi / 180)

        bearingd = np.arccos((np.sin(YB * np.pi / 180) - cosdeltaAB * np.sin(YA * np.pi / 180)) \
    	             / (np.sqrt(1 - cosdeltaAB * cosdeltaAB) * np.cos(YA * np.pi / 180)))

        return bearingd 
 
    	
    
    #-------------------------------------------------------------------------
    def project(self,A, B, datac,flag):
        """
    	Function to project lat/long data along a line defined by the lat/long coordinates
    	   of the two points defining the line ends.
    
    	Args:
    		A (2*1 np array of floats): beginning of the line on which to project
    		B (2*1 np array of floats): end of the line on which to project
    		datac (np array): data to project
    	
    	(c) licence CCby-nc-sa : http://creativecommons.org/licenses/by-nc-sa/4.0/ 2021
    
    	"""
    
    	# Earth radius in km
        R = 6371
    
    	# Begin stepping on data
        if flag == 1:
            lon = 'LONOBS'
            lat = 'LATOBS'
        else:
            lon = 'LON'
            lat = 'LAT'
            
        longM = datac[lon]
        # longM.where(longM2.Longitude > 180,longM2.Longitude - 360,longM2.Longitude)
        latM = datac[lat]

    	# If lat long         
    	# compute on great circles
        AM = self.calc_length(A[0], A[1], longM, latM)
    	# Calcul azimuth AB
        AzAB = self.bearing(A[0], A[1], B[0], B[1])
    	# Calcul azimuth MM'
        AzMM = AzAB + np.pi / 2 # Verifier sens du signe
    	# Calcul azimuth AM
        AzAM = self.bearing(A[0], A[1], longM, latM)
    	# calcul Bheta = angle AM-MM from bearings of AM and MM'
        Bheta = (np.pi - AzMM + AzAM)            
    	# calcul de la longueur
        coordproj = R * np.arcsin(np.sin(AM / R) * np.sin(Bheta))
    
        return coordproj


##############################################################################
############################ TRAPPED-CHARGE ##################################
##############################################################################
#-----------------------------------------------------------------------------
#--------------------------------- Classes -----------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
class AgeTime(QWidget):
    """ This class load and plot the age (or n/N) evolution with time. """
    
    def __init__(self,outputPath,GraphWin):
        super().__init__()
        
        # Parameters - choose the thermochronometers to plot
        self.MainLabel = QLabel("Choose the thermochronometers to plot:")
        self.MainLabel.setFont(conf.fontBold12)
        self.TLbox = QCheckBox("Th.L.")
        self.TLbox.setEnabled(False)
        self.OSLbox = QCheckBox("OSL")
        self.OSLbox.setEnabled(False)
        self.ESRbox = QCheckBox("ESR")
        self.ESRbox.setEnabled(False)
        self.storeSystems = {}
        self.outputPath = outputPath
        agecol, errname, agename, predname, self.colores, profdict,He43dict,self.coloramp,self.ageMarker = Thermo_settings.dict_pecube()
        self.GraphWin = GraphWin
        filenames = os.listdir(outputPath)
        self.listFiles = [file for file in filenames if file.startswith('ESRTime') or file.startswith('OSLTime') or file.startswith('TLTime')]
        
        # Layout
        self.MainLayout = QVBoxLayout()
        self.MainLayout.addWidget(self.MainLabel)
        self.MainLayout.addStretch(1)
        self.Grid1 = QGridLayout()
        self.Grid1.setSpacing(10)
        self.Grid1.addWidget(self.TLbox,1,1)
        self.Grid1.addWidget(self.OSLbox,1,2)
        self.Grid1.addWidget(self.ESRbox,1,3)
        self.Grid1.setColumnStretch(4,1)
        self.Grid1.setRowStretch(6,1)
        
        self.MainLayout.addLayout(self.Grid1)
        self.MainLayout.addStretch(4)
        
        self.updateThermEnabled(self.listFiles)
        
        # Signals
        self.TLbox.stateChanged.connect(lambda: self.updateTherm2plot(self.TLbox,1,conf.Variable_names['ThL vs time file'],'ThL'))
        self.OSLbox.stateChanged.connect(lambda: self.updateTherm2plot(self.OSLbox,2,conf.Variable_names['OSL vs time file'],'OSL'))
        self.ESRbox.stateChanged.connect(lambda: self.updateTherm2plot(self.ESRbox,3,conf.Variable_names['ESR vs time file'],'ESR'))
    
    
    #-------------------------------------------------------------------------
    def updateThermEnabled(self,files):
        """Check which thermochronometer the user can select."""
        for i in range(len(files)):
            if files[i].startswith(conf.Variable_names['ThL vs time file']):
                self.TLbox.setEnabled(True)
            if files[i].startswith(conf.Variable_names['OSL vs time file']):
                self.OSLbox.setEnabled(True)
            if files[i].startswith(conf.Variable_names['ESR vs time file']):
                self.ESRbox.setEnabled(True)
        
    #-------------------------------------------------------------------------
    def updateTherm2plot(self,sys,flag,file,label):
        if sys.isChecked():
            self.storeSystems[file] = label
        else:
            self.storeSystems[file] = 0
            
    #-------------------------------------------------------------------------       
    def plot_ageTime(self):
        """ Plot the age vs time for each thermochronometers"""
        
        # Create the canva
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        # Create a customed toolbar to show/hide ages from specific thermochronometer
        self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
        self.plot_list = []
        self.categories = []
        
        # Iterate through thermochronometers to plot
        for k, v in self.storeSystems.items():
            if v:
                # load data
                filePath = os.path.join(self.outputPath,k)
                try:
                    dataset = pd.read_csv(filePath)
                except FileNotFoundError as E:
                    QErrorMessage(self).showMessage(filePath + " file not found. Please check if it exists: <br>" +
                                                    str(E))
                    break
                # Get data
                dataFrame = pd.DataFrame(dataset)
                # Get sample names
                HeadersLabel = list(dataFrame.columns.values)
                nsamples = len(HeadersLabel) - 1
                data = dataFrame.to_numpy().astype(np.float32)
                # Time is first column
                Time = data[:,0]

                # iterate through samples
                if nsamples > 1:
                    c1 = self.coloramp[self.storeSystems[k]][0]
                    c2 = self.coloramp[self.storeSystems[k]][1]
                    for p in range(nsamples):
                        c = set_plot.get_gradient_color(c1,c2,p/nsamples)
                        p1 = self.plotSpace.axes.plot(Time, data[:,p+1],linestyle="-",
                                                 label = self.storeSystems[k]+"_"+HeadersLabel[p+1], color = c)
                        # Add plot to list
                        self.plot_list.append(p1)
                else:
                    p1 = self.plotSpace.axes.plot(Time, data[:,1],linestyle="-",
                                             label = self.storeSystems[k]+"_"+HeadersLabel[1], color = self.colores[self.storeSystems[k]])
                    # Add plot to list
                    self.plot_list.append(p1)
                self.categories.append(self.storeSystems[k])
                
        # set legend
        self.plotSpace.axes.legend(loc = 'best')
        if (int(self.GraphWin.ParametersInput.DParameters[conf.Variable_names['ESR_misfit_target']]) or
        int(self.GraphWin.ParametersInput.DParameters[conf.Variable_names['OSL_misfit_target']])):
              self.plotSpace.axes.set_ylabel(u'n/N')      
        else:
            self.plotSpace.axes.set_ylabel(u'Age (Ma)')
        self.plotSpace.axes.set_xlabel(u'Time (Myr)')
        
        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'Age-Time')




##############################################################################
############################ 4He/3He #########################################
##############################################################################
#-----------------------------------------------------------------------------
#--------------------------------- Classes -----------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
class He43_plot(QWidget):
    """ 
      This class contains the parameters to plot 4He/3He data.
      Table = Table with observed 4He/3He spectra.
    
    """
    def __init__(self,inputData,Compare43,inputc,GraphWin):
        super().__init__()
        
        self.agecol, self.errname, self.agename, self.predname, self.colores, self.profdict, self.He43dict, coloramp, self.ageMarker =  Thermo_settings.dict_pecube()
        
        self.datac = Compare43
        self.inputc = inputc
        self.inputdata = inputData
        self.GraphWin = GraphWin
        self.dataplot = ['Observed','Predicted']

        #Define parameters
        self.Title = QLabel("Plot 4He/3He data")
        self.Title.setFont(conf.fontBold12)
        self.DataComboLabel = QLabel("Data to plot:")
        self.DataComboLabel.setFont(conf.font12)
        self.DataCombo = QComboBox()
        self.DataCombo.addItems(["4He/3He spectrum", "Step age","Normalized step age"])
        self.DataCombo.setMinimumContentsLength(15)
        self.SampleComboLabel = QLabel("Grain to plot:")
        self.SampleComboLabel.setFont(conf.font12)
        self.SampleCombo = QComboBox()
        self.SampleCombo.setMinimumContentsLength(15)
        # Get Sample list
        Samples = set(self.inputc["Sample43"].to_numpy().tolist())
        Samples = [i for i in Samples if isinstance(i, str)]
        self.SampleCombo.addItems(Samples)
        self.plotStepAge = 0 #plot 4He/3He by default, 1 = plot step age
        self.HelpLabel = QLabel("Select a range to plot, or press 'Ok' to plot all data.")
        self.HelpLabel.setFont(conf.font12)
        self.observedData = 0
        self.data2Plot = 0
        
        #Define splitters
        self.Splitter = U_interface.QCustomSplitter(Qt.Vertical)
        self.ParamSplitter = QWidget()
        self.TableSplitter = QWidget()
        
        #Define layout
        self.MainLayout = QVBoxLayout()
        VBox = QVBoxLayout()
        Grid = QHBoxLayout()
        self.MainLayout.addWidget(self.Title)
        Grid.addWidget(self.DataComboLabel)
        Grid.addWidget(self.DataCombo)
        Grid.addStretch(2)
        hbox = QHBoxLayout()
        hbox.addWidget(self.SampleComboLabel)
        hbox.addWidget(self.SampleCombo)
        hbox.addStretch(2)
        self.MainLayout.addLayout(Grid)
        self.MainLayout.addLayout(hbox)
        self.MainLayout.addStretch(1)
        self.MainLayout.setSpacing(10)
        # Grid2 = QVBoxLayout()
        # Grid2.addWidget(self.HelpLabel)
        # Grid2.addWidget(Table)
        
        #Build splitters
        self.ParamSplitter.setLayout(VBox)
        # self.TableSplitter.setLayout(Grid2)
        self.Splitter.addWidget(self.ParamSplitter)
        # self.Splitter.addWidget(self.TableSplitter)
        
        #Check updates
        self.DataCombo.currentIndexChanged.connect(lambda: self.updateCombo())
    
    #------------------------------------------------------------------------
    def updateCombo(self):
        
        self.plotStepAge = self.DataCombo.currentIndex() 
   
    
    #-------------------------------------------------------------------------
    def plotHe43(self):
        """
        Author: Xavier Robert, university of Grenoble
        Adapted by: Maxime Bernard, university of Potsdam (12/05/2023)
        
    	Args:
    		inputdata (string): input data declared in Pecube.in; This is used to plot the errorbars
    
    		profiletype (list, optional): type of profile = ['Latitude', 'Longitude', 'Altitude', 'Projected']
    						              If [], no age profile is plotted.
    
    		dataplot (list): List of data to plot ; By default, the altitude will be plotted
    						 ['AHe', 'AFT', 'ZHe', 'ZFT', 'KAr', 'MAr', 'BAr', 'MTL', 'TTp']
    						 Rem: For the moment, MTL and TTp not implemented.
    
    		datac (string): results of Pecube forward modeling
    
    		inputc (array, optional): Array of data to plot.
    		                          Defaults to None.
    
    		coordproj (array, optional): Predictions projected along A-B.
    		                             Defaults to None.
    		
    		coordprojinputc (array, optional): Observations projected along A-B.
    		                                   Defaults to None.
    
    		agerange (2*1 array of floats, optional): range of the ages to plot on the profiles
    												  [min, max]. Defaults to None.
    
    		agecol (dict, optional): respective column number of the data in the file comparison.txt.
    		                         Defaults to None.
    		
    		errname (dict, optional): respective column number of the error on data in the data input file. 
    		                          Defaults to None.
    		
    		colores (dict, optional): Colors used for the different age system. 
    		                          Defaults to None.
    				
    		predname (dict, optional): legend of each predicted system. 
    		                           Defaults to None.
    
    		agename (dict, optional): legend of each data system. 
    		                          Defaults to None.
    
    		profdict (dict, optional): Columns of profile types. 
    		                           Defaults to None.
    
    		graphpath (str, optional): name of the folder where the plot will be written
    						            Usually you do not have to change it. 
    									Defaults to 'Graphs'.
    
    		graphtitle (str, optional): title to write on the graph.
    		                            Defaults to None.
    
    	"""
    

        # Get profile type
        profiletype = self.DataCombo.currentText()
        SampleID = self.SampleCombo.currentText()
        
        # Get data from the sample in the data and predictions
        IndexObs = np.where(self.inputc["Sample43"]==SampleID)
        IndexPred = np.where(self.datac["Sample43Pred"]==SampleID)
        # If step age, compute step age
        if profiletype == 'Step age' or profiletype == 'Normalized step age':
            # Ejection only
            EjecProfile = self.datac[self.He43dict['Ejection']][IndexPred]
            # Get Age
            Age = self.datac[self.agename['AHE43']+'PRED'][IndexPred]
            AgePred = Age[0]
            Age = self.datac[self.agename['AHE43']+'OBS'][IndexPred]
            AgeObs = Age[0]
    	# fig1 = plt.figure()
        # Create the canva
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        # Create a customed toolbar to show/hide ages from specific thermochronometer
        self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
        
        print(u'\tPlotting %s' %(profiletype))
        if self.inputdata:
			    # Plot input data with Error bars (observations)
            if profiletype == '4He/3He spectrum' : # Plot Rs/Rb / Rbulk
                # The reference point of the box is left bottom corner
                xdata = self.inputc[self.He43dict["S3He"]][IndexObs]
                ydata = self.inputc[self.He43dict["4HE"]][IndexObs]
                set_plot.make_error_boxes(self.plotSpace.axes,xdata,ydata,self.inputc[self.errname["S3HE"]][IndexObs],self.inputc[self.errname["HE43"]][IndexObs],
                                 facecolor='g',edgecolor='none',alpha=0.5,label='Portion extracted')
                self.plotSpace.axes.errorbar(self.inputc[self.He43dict["S3He"]][IndexObs], self.inputc[self.He43dict["4HE"]][IndexObs], yerr = self.inputc[self.errname["HE43"]][IndexObs], 
			    			         xerr = self.inputc[self.errname["S3HE"]][IndexObs],fmt = 'o', label = 'Observed', color = self.colores['He43obs'])
                
            elif profiletype == 'Step age':
                  stepAgeobs = self.inputc[self.He43dict["4HE"]][IndexObs]/EjecProfile*AgeObs
                  error = stepAgeobs * (self.inputc[self.errname["HE43"]][IndexObs]/self.inputc[self.He43dict["4HE"]][IndexObs])
                  self.plotSpace.axes.errorbar(self.inputc[self.He43dict["S3He"]][IndexObs],stepAgeobs , yerr = error, 
			    		         fmt = 'o', label = 'Observed', color = self.colores['He43obs'])
            elif profiletype == 'Normalized step age':
                  stepAge = self.inputc[self.He43dict["4HE"]][IndexObs]/EjecProfile
                  # stepAge = stepAge/max(stepAge[:])
                  error = stepAge * (self.inputc[self.errname["HE43"]][IndexObs]/self.inputc[self.He43dict["4HE"]][IndexObs])
                  self.plotSpace.axes.errorbar(self.inputc[self.He43dict["S3He"]][IndexObs],stepAge , yerr = error, 
			    		         fmt = 'o', label = 'Observed', color = self.colores['He43obs'])

			# Plot predictions
        if profiletype == '4He/3He spectrum':
            # prediction are plotted at Sum3He observations
            xdata = self.inputc[self.He43dict["S3He"]][IndexObs]
            self.plotSpace.axes.plot(xdata, self.datac[self.He43dict["4HePred"]][IndexPred], 
		    			     marker = 's', linestyle = '-', 
							 label = 'Predicted', color = self.colores['He43pred'], alpha = 0.3)
            # self.plotSpace.axes.plot(self.datac[self.He43dict["S3HePred"]][IndexPred], self.datac[self.He43dict["4HePred"]][IndexPred], 
		    # 			     marker = 's', linestyle = '-', 
			# 				 label = 'Predicted', color = self.colores['He43pred'], alpha = 0.3)
            self.plotSpace.axes.set_xlim(0,1.2)
            self.plotSpace.axes.set_ylim(0,1.8)
            
        elif profiletype == 'Step age':
            # prediction are plotted at Sum3He observations
            xdata = self.inputc[self.He43dict["S3He"]][IndexObs]
            stepAgepred = self.datac[self.He43dict["4HePred"]][IndexPred]/EjecProfile*AgePred
            self.plotSpace.axes.plot(xdata, stepAgepred, 
		    	    		 marker = 's', linestyle = '-',
							 label = 'Predicted', color = 'r', alpha = 0.3)	
            self.plotSpace.axes.set_xlim(0,1.2)
            self.plotSpace.axes.set_ylim(0,max(max(stepAgeobs),max(stepAgepred))+1)
            
        elif profiletype == 'Normalized step age':
            # prediction are plotted at Sum3He observations
            xdata = self.inputc[self.He43dict["S3He"]][IndexObs]
            stepAge = self.datac[self.He43dict["4HePred"]][IndexPred]/EjecProfile
            # stepAge = stepAge/max(stepAge[:])
            self.plotSpace.axes.plot(xdata, stepAge, 
		    	    		 marker = 's', linestyle = '-',
							 label = 'Predicted', color = (1,69/255,0), alpha = 0.3)	
            self.plotSpace.axes.set_xlim(0,1.2)
            self.plotSpace.axes.set_ylim(0,1.2)

		# set legend
        self.plotSpace.axes.legend(loc = 'best')
		#plt.xlabel(u'Longitude (°)')
        if profiletype == '4He/3He spectrum':
            self.plotSpace.axes.set_xlabel(u'$SF^3He$')
            self.plotSpace.axes.set_ylabel(u'$R_{step}/R_{bulk}$')
        elif profiletype == 'Step age':
            self.plotSpace.axes.set_xlabel(u'$SF^3He$')
            self.plotSpace.axes.set_ylabel(u'Step age (Ma)')
        elif profiletype == 'Normalized step age':
            self.plotSpace.axes.set_xlabel(u'$SF^3He$')
            self.plotSpace.axes.set_ylabel(u'Normalized step age (Ma)')
            
        # plt.title(graphtitle)
        
        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'4He/3He')
   
        return
    

##############################################################################
############################ Fission track ###################################
##############################################################################
#-----------------------------------------------------------------------------
#--------------------------------- Classes -----------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

class ChooseMFTL(QWidget):
    """This class ask the user whether he/she wants to plot MFTL or Std MFTL data."""
    def __init__(self):
        super().__init__()
        
        self.Data2Plot_label = QLabel('Choose data to plot: ')
        self.Data2Plot_label.setAlignment(Qt.AlignRight)
        self.Data2Plot_label.setFont(conf.font10)
        items = ['Mean fission track length','std of MFTL']
        self.Data2Plot_list = QComboBox()
        self.Data2Plot_list.setMinimumContentsLength(15)
        self.Data2Plot_list.addItems(items)
        self.Data2Plot = 0 # Default is MFTL
        
        self.Layout = QHBoxLayout()
        self.Layout.addWidget(self.Data2Plot_label)
        self.Layout.addWidget(self.Data2Plot_list)
        
        self.Data2Plot_list.currentIndexChanged.connect(lambda: self.update_Data2plot())
    
    def update_Data2plot(self):
        """ update data2plot"""
        self.Data2Plot = self.Data2Plot_list.currentIndex()