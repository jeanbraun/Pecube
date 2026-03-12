#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to plot output results in the interface. 
It contains classes that handle inversion results from Pecube.

@author: maxime - maxime.bernard@uni-potsdam.de

"""


import os
from PyQt5.QtWidgets import (QErrorMessage, QWidget,QGridLayout,
                             QVBoxLayout,QHBoxLayout,QLabel,QComboBox,
                             QCheckBox,QLineEdit,QErrorMessage)

from PyQt5.QtGui import (QIntValidator)
from PyQt5.QtCore import Qt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from scipy.stats import gmean
import peakutils    # To find number of peaks
from peakutils.plot import plot as pplot
from scipy.optimize import curve_fit    # To find gaussian parameters


import configs as conf
import Utils.PGUI_utils as pgu

#-----------------------------------------------------------------------------
#-------------------------------- Functions ----------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#########################################################################################
class plotPecubeNA(QWidget):
    """
    Author: Xavier Robert, University of Grenoble
    Adapted by: Maxime bernard, University of Potsdam (11/05/2023)
    
    This class is designed to make the plotPecubeNA function from Xavier Robert works with PecubeGUI.
    It reads and plot the results from inversion.
    
    Main code of NA Plot for Pecube inversions with NA

    INPUTS:
        param (list of str): Define as many variable as you have, with their unit
                             Check the order in na.sum (open with text editor) or in NA_Results
                             If you use the later, first column always is the misfit
                             The others paramaters are from topo_parameters.txt then from fault_parameters.txt
                             Exemple : param=['Offset (km)','Basal Temperature (°C)','Slip rate (km∕Ma)']
                             /!\ IF YOU WANT TO USE A SLASH USE THIS ONE --> '∕' <-- , IT'S A UNICODE DIVISION SYMBOL
                            WINDOWS AND OSX DON'T ALLOW THE USE OF THE REGULAR SLASH 

        dataplot (list of couple of integers): Set the couple of variables to plot against each other
                                               Exemple : Offset vs Slip rates
                                               dataplot = [(1,3)] -->  plot=(1,3)
                                               If you plot 2D pdfs (contours), please, CHECK that the couple of parameters to plot
                                               are the same and in the same order than in the nab.in file.
                                               This python script checks it and will insult you if this is not compatible !!!

        tick_space (array of floats): Set the space between ticks for x and y axes for each parameters
                                      (same order than the list param).
                                      If the tick format does not fit your variables, 
                                      you may need to modify the dictionnary tick_order in the function multiplot

        inv_results (str, optional): Name of the NA file with the inversion results; Usually NA_Results.
                                     Defaults to 'NA/NA_results.csv'.

        data_nab (str, optional): Name of the NAB file with the inversion results. Usually nab.out.
                                  Defaults to 'NA/NAB/nab.out'.

        graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                     Do not forget the '/' at the end.
                                     Defaults to 'NA/Graphs/'.

        PDF_1D (bool, optional): Choose if you want the 1-PDFs (Probability Density Function) 
                                 Defaults to True.

        PDF_2D (bool, optional): Choose if you want the 2-PDFs (Probability Density Function) 
                                 Defaults to False.

        pdf1d_results (string, optional): Print the 1-pdfs in a text file.
                                          Defaults to None.

        pdf2d_results (string, optional): Print the 2-pdfs in a text file.
                                          Defaults to None.

        size_x (int, optional): Size of the font for the x axes label. 
                                Defaults to 15.

        size_y (int, optional): Size of the font for the y axes label. 
                                Defaults to 15.

        size_m (int, optional): Size of the font for the markers label. 
                                Defaults to 15.

        size_mis (int, optional): Size of the markers of the misfits. 
                                    Defaults to 50.

        peak_thres (float, optional): Threshold to find peaks; between 0. and 1.
                                      Default = 0.05.

        peak_min_dist (interger, optional): Minimum distance between the peaks
                                            Default = 30.
    
    Raises:
        NameError  : Input file names does not exists
        ValueError : Problems with number of values or with NaN
        ImportError: File format not supported, please, change the format
     
    
    """
    
    def __init__(self,fileResults,fileNAB,FolderNA,GraphWin):
        super().__init__()
    
        ###### Define input parameters #######
        # First initiate the name of the parameters inverted
        param = ["Parameter 1", "Parameter 2"] 
        
        # Then we need to ask the user which inverted parameter he/she wants to plot
        dataplot = [(1,2)] # by default, the two first
        # Path of files
        self.fileResults = fileResults
        self.fileNAB = fileNAB
        self.FolderNA = FolderNA
        self.GraphWin = GraphWin
        self.inputParam = self.GraphWin.ParametersInput.DParameters
        # Double check if files exists
        if not os.path.exists(self.fileResults) and not os.path.exists(self.fileNAB):
            QErrorMessage(self).showMessage("File not found: NA_results.csv. Cannot plot inversion results.")
            return
        # if not os.path.exists(self.fileNAB):
        #     QErrorMessage(self).showMessage("File not found: nab.out. Cannot plot NAB results.")
       
        # Choose if you want the PDFs (Probability Density Function) 
        # plots 1D or 2D contour with the misfit plots (True of False)
        self.PDF_1D = False
        self.PDF_2D = False
        # Parameters to find peaks
        peak_thres = 0.05
        peak_min_dist = 30
        self.logValue = 0
        
        ###### Set parameters and boxes for the user ########
        self.ChoosePlotLabel = QLabel('Choose plot: ')
        self.ChoosePlotLabel.setFont(conf.fontBold12)
        self.ChoosePlotCombo = QComboBox()
        self.ChoosePlotCombo.addItems(['Misfit evolution','2D parameter space','2D parameter space + 1D-pdf','1D-pdf single parameter'])
        self.ChoosePlotCombo.setMinimumContentsLength(15)
        self.VBox = QVBoxLayout()
        HBox = QHBoxLayout()
        HBox.addWidget(self.ChoosePlotLabel)
        HBox.addWidget(self.ChoosePlotCombo)
        HBox.addStretch(1)
        self.VBox.addLayout(HBox)
        
        # For 2D parameter space plot
        self.widget1 = QWidget()
        self.Param1Label = QLabel('Parameter ID (x):')
        self.Param1Label.setFont(conf.font10)
        self.Param1Label.setToolTip('The parameter ID number you want along X axis')
        self.Param1Label.setAlignment(Qt.AlignCenter)
        self.Param1Edit = QLineEdit("1")
        self.Param1Edit.setValidator(QIntValidator())
        self.Param1Edit.setAlignment(Qt.AlignCenter)
        self.Param1Edit.setEnabled(False)
        self.Param2Label = QLabel('Parameter ID (y):')
        self.Param2Label.setFont(conf.font10)
        self.Param2Label.setToolTip('The parameter ID number you want along Y axis')
        self.Param2Label.setAlignment(Qt.AlignCenter)
        self.Param2Edit = QLineEdit("2")
        self.Param2Edit.setValidator(QIntValidator())
        self.Param2Edit.setAlignment(Qt.AlignCenter)
        self.Param2Edit.setEnabled(False)
        self.Name1Label = QLabel('Parameter Label (x):')
        self.Name1Label.setFont(conf.font10)
        self.Name1Label.setToolTip('How do you want to name the parameter along X axis')
        self.Name1Label.setAlignment(Qt.AlignCenter)
        self.Name1Edit = QLineEdit("Parameter 1")
        self.Name1Edit.setAlignment(Qt.AlignCenter)
        self.Name1Edit.setEnabled(False)
        self.Name2Label = QLabel('Parameter Label (y):')
        self.Name2Label.setFont(conf.font10)
        self.Name2Label.setToolTip('How do you want to name the parameter along Y axis')
        self.Name2Label.setAlignment(Qt.AlignCenter)
        self.Name2Edit = QLineEdit("Parameter 2")
        self.Name2Edit.setAlignment(Qt.AlignCenter)
        self.Name2Edit.setEnabled(False)
        self.tickSpaceXLabel = QLabel('Tick space (x):')
        self.tickSpaceXLabel.setFont(conf.font10)
        self.tickSpaceXLabel.setAlignment(Qt.AlignCenter)
        self.tickSpaceXLabel.setToolTip("Tick space along X axis.")
        self.tickSpaceXEdit = QLineEdit("0.1")
        self.tickSpaceXEdit.setValidator(conf.DoubleValidator)
        self.tickSpaceXEdit.setAlignment(Qt.AlignCenter)
        self.tickSpaceYLabel = QLabel('Tick space (y):')
        self.tickSpaceYLabel.setFont(conf.font10)
        self.tickSpaceYLabel.setAlignment(Qt.AlignCenter)
        self.tickSpaceYLabel.setToolTip("Tick space along Y axis.")
        self.tickSpaceYEdit = QLineEdit("0.1")
        self.tickSpaceYEdit.setValidator(conf.DoubleValidator)
        self.tickSpaceYEdit.setAlignment(Qt.AlignCenter)
        self.PeakThresholdLabel = QLabel("Peak threshold:")
        self.PeakThresholdLabel.setFont(conf.font10)
        self.PeakThresholdLabel.setAlignment(Qt.AlignCenter)
        self.PeakThresholdLabel.setToolTip("Normalized threshold to find peaks. Only the peaks with amplitude higher than the threshold will be detected (value: 0-1)")
        self.PeakThreshold = QLineEdit(conf.Variable_names["Peak threshold"])
        self.PeakThreshold.setAlignment(Qt.AlignCenter)
        self.PeakThreshold.setEnabled(False)
        self.PeakMinDistLabel = QLabel("Peak min distance:")
        self.PeakMinDistLabel.setFont(conf.font10)
        self.PeakMinDistLabel.setAlignment(Qt.AlignCenter)
        self.PeakMinDistLabel.setToolTip("Minimum distance between each detected peak. The peak with the highest amplitude is preferred to satisfy this constraint.")
        self.PeakMinDist = QLineEdit(conf.Variable_names["Peak min distance"])
        self.PeakMinDist.setAlignment(Qt.AlignCenter)
        self.PeakMinDist.setEnabled(False)
        self.NormalizationLabel = QLabel("Normalized to: ")
        self.NormalizationLabel.setFont(conf.font10)
        self.NormalizationLabel.setAlignment(Qt.AlignRight)
        self.NormalizationLabel.setToolTip("Normalize the misfit values to (e.g., number of data)")
        self.NormalizationValue = QLineEdit("1")
        self.NormalizationValue.setAlignment(Qt.AlignCenter)
        self.logBox = QCheckBox('log')
        self.logBox.setEnabled(True)
        self.SummaryInversionTitle = QLabel('Summary of inversion:')
        self.SummaryInversionTitle.setFont(conf.fontBold12)
        self.SummaryInversion = QLabel()
        self.SummaryInversion.setFont(conf.font8)
        self.SummaryInversion.setEnabled(False)
        self.read_summary_inversion()
        
        # Set layout
        Grid1 = QGridLayout()
        Grid1.setSpacing(10)
        Grid1.addWidget(self.Param1Label,0,1)
        Grid1.addWidget(self.Param2Label,0,2)
        Grid1.addWidget(self.Param1Edit,1,1)
        Grid1.addWidget(self.Param2Edit,1,2)
        Grid1.addWidget(self.Name1Label,2,1)
        Grid1.addWidget(self.Name2Label,2,2)
        Grid1.addWidget(self.Name1Edit,3,1)
        Grid1.addWidget(self.Name2Edit,3,2)
        Grid1.addWidget(self.NormalizationLabel,4,0)
        Grid1.addWidget(self.NormalizationValue,4,1)
        Grid1.addWidget(self.logBox,4,2)
        Grid1.addWidget(self.PeakThresholdLabel,5,1)
        Grid1.addWidget(self.PeakMinDistLabel,5,2)
        Grid1.addWidget(self.PeakThreshold,6,1)
        Grid1.addWidget(self.PeakMinDist,6,2)
        Vbox_SumInv = QVBoxLayout()
        Vbox_SumInv.addWidget(self.SummaryInversionTitle)
        Vbox_SumInv.addWidget(self.SummaryInversion)
        # Grid1.addWidget(self.tickSpaceXLabel,4,1)
        # Grid1.addWidget(self.tickSpaceYLabel,4,2)
        # Grid1.addWidget(self.tickSpaceXEdit,5,1)
        # Grid1.addWidget(self.tickSpaceYEdit,5,2)
        Grid1.setColumnStretch(3,1)
        Grid1.setRowStretch(7,1)
        self.widget1.setLayout(Grid1)
        self.VBox.addWidget(self.widget1)
        self.VBox.addLayout(Vbox_SumInv)
        
        #Signals
        self.ChoosePlotCombo.currentIndexChanged.connect(lambda: self.update_plotOption())
        self.logBox.stateChanged.connect(lambda: self.update_logBox())
        
    
    #----------------------------------------------------------
    def update_logBox(self):
        """ plot inversion result in log scale."""
        if self.logBox.isChecked():
            self.logValue = 1
        else:
            self.logValue = 0
        
    #----------------------------------------------------------
    def update_plotOption(self):
        """ Update the user interface according to the option chosen to plot inversion results."""
        if self.ChoosePlotCombo.currentIndex() == 3: # 1D PDF only
            # disable parameter 2 options
            self.Name2Edit.setEnabled(False)
            self.Param2Edit.setEnabled(False)
            self.Name1Edit.setEnabled(True)
            self.Param1Edit.setEnabled(True)
            self.PeakThreshold.setEnabled(True)
            self.PeakMinDist.setEnabled(True)
            
        elif self.ChoosePlotCombo.currentIndex() == 2: # 2D Parameter space + 1D pdf
            self.Name2Edit.setEnabled(True)
            self.Param2Edit.setEnabled(True)
            self.Name1Edit.setEnabled(True)
            self.Param1Edit.setEnabled(True)
            self.PeakThreshold.setEnabled(False)
            self.PeakMinDist.setEnabled(False)
        
        elif self.ChoosePlotCombo.currentIndex() == 1: # 2D Parameter space
            self.Name2Edit.setEnabled(True)
            self.Param2Edit.setEnabled(True)
            self.Name1Edit.setEnabled(True)
            self.Param1Edit.setEnabled(True)
            self.PeakThreshold.setEnabled(False)
            self.PeakMinDist.setEnabled(False)
            
        elif self.ChoosePlotCombo.currentIndex() == 0: # Misfit evolution
            self.Name2Edit.setEnabled(False)
            self.Param2Edit.setEnabled(False)
            self.Name1Edit.setEnabled(False)
            self.Param1Edit.setEnabled(False)
            self.PeakThreshold.setEnabled(False)
            self.PeakMinDist.setEnabled(False)
        
    #----------------------------------------------------------
    def read_summary_inversion(self,inv_results = 'NA/NA_int_results.csv',):
        """ Read NA_results.csv file and show inversion summary to the user. """
        
        # Check that NA_results exists
        if not os.path.isfile(self.fileResults):
            print('\033[91mERROR:\033[00m File %s does not exist' %(str(inv_results)))
            # update text
            self.SummaryInversion.setText("Inversion file NA_result.csv not found...")
        
        # Loading the data, calling as many variables as given in param + 1 for the misfit
        if inv_results[-3:] == 'txt':
            nb_var = np.loadtxt(self.fileResults, unpack = True)
        elif inv_results[-3:] == 'csv':
            nb_var = np.genfromtxt(self.fileResults, delimiter = ',', skip_header = 1, unpack = True)
        else:
            # raise ImportError('NA_results format not supported...\n\n')
            QErrorMessage(self).showMessage('NA_results format not supported...\n\n')
            return
        
        # Check if misfits are not with NaN or infinite values...
        # Raise error if this is the case, and recheck your inversion parameters
        if np.any(np.isnan(nb_var)) or not np.all(np.isfinite(nb_var)):
            print("\nInfinity: ", np.inf in nb_var)
            QErrorMessage(self).showMessage('NA_results.csv contains NaN or infinite numbers...\n\n')
            return
            # raise ValueError('NA_results.csv contains NaN or infinite numbers...\n\n')

        # Sort nb_var in function of misfit to get the best misfits over the lower ones 
        #          --> Better scatter plot visualisation
        nb_var_nosort = nb_var.T
        nb_var = nb_var.T[nb_var.T[:,0].argsort()[::-1]].T
        
        # Number of inverted parameters
        NbParam = len(nb_var[:,0]) -1 # Datascore is last column
        
        # update text
        index_min_misfit = np.where(nb_var[0,:] == min(nb_var[0,:]))[0][0]
        text = []
        for i in range(NbParam):
            if i < 9:
                text.append('   - Param00'+str(i+1)+' = ' + str(round(nb_var[i+1,index_min_misfit],3))+'\n')
            else:
                text.append('   - Param0'+str(i+1)+' = ' + str(round(nb_var[i+1,index_min_misfit],3))+'\n')
        text = np.asarray(text)
        
        self.SummaryInversion.setText(" Lowest misfit value: " + str(nb_var[0,index_min_misfit]) + '\n' +
                                      " Best fit parameters: " + '\n' + ' ' +
                                      ' '.join(text))
                                      
        
        
    #----------------------------------------------------------
    def plot_results(self,PDF_1D=False,PDF1D_only=False,PDF_2D=False,Misfit_evol=False):
        """ To read and plot data according to user's choice of parameters """
        
        # if self.ChoosePlotCombo.currentIndex() == 0:
        # Plot 2D parameter space
        param = [self.Name1Edit.text(),self.Name2Edit.text()]
        dataplot = [(int(self.Param1Edit.text()),int(self.Param2Edit.text()))]
        if dataplot[0][0] == dataplot[0][1]:
            dataplot = [(int(self.Param1Edit.text()),int(self.Param2Edit.text())-1)]
 
        tick_space = [float(self.tickSpaceXEdit.text()),float(self.tickSpaceYEdit.text())]
        # Set the size of the font for the x and y axes label
        size_x = 15
        size_y = 15
        size_m = 15

        # Set the size of the markers of the scatter plot and of the misfit
        size_mis = 50
        
        self.MisfitEvol = Misfit_evol

        # Parameters to find peaks
        peak_thres = float(self.PeakThreshold.text())
        peak_min_dist = int(self.PeakMinDist.text())
        self.NAplot_Pecube(param = param,
                      dataplot = dataplot,
                      inv_results = self.fileResults,
                      data_nab = self.fileNAB,
                      graph_path  = self.FolderNA,
                      PDF_1D = PDF_1D, PDF_2D = PDF_2D,
                      PDF1D_only = PDF1D_only,
                      tick_space = tick_space,
                      size_x = size_x, size_y = size_y, size_m = size_m,
                      size_mis = size_mis,
                      peak_thres = peak_thres, peak_min_dist = peak_min_dist,
                      misfitEvol = self.MisfitEvol)
        
        return
    
    #----------------------------------------------------------
    def read_NABout(self,data_nab, dataplot, PDF1D_only=False, pdf1d_results = None, pdf2d_results = None):
        """
        Read the nad.out file and write the PDF_DATA.txt and PDF_2D_DATA.txt
        
        INPUTS:
        	data_nab       : path and name of the nab.out file
            dataplot       : list of couple of variables we want to plot against each other
            pdf1d_results  : path and name of the file with the PDF 1D data if file saved ;
                             None by default
            pdf2d_results  : path and name of the file with the PDF 2D data if file saved ;
                             None by default
        
        OUTPUTS:
        	data1D : 1D marginals
            data2D : 2D marginals
        		
        USAGE:
            data1D, data2D = read_NABout(data_nab, dataplot, pdf1d_results, pdf2d_results)
        
        Author: Xavier Robert, Grenoble 2021/02/10
        
        """

        # Open the nab.out file and store the content in lines
        with open(data_nab, 'r') as f1r:
            lines = f1r.readlines()
            
        print("after read file")
        test_couple = []
        # # Get number of parameters
        # Nparam = next((int(l[-13:len(l)]) for l in lines if u'Number of dimensions' in l),None)
        # Nbins1D = next((int(l[-13:len(l)]) for l in lines if u'Number of bins per axis for 1D marginals' in l),None)
        # # Get number of parameters 2D
        # Nparam2D = next((int(l[-13:len(l)]) for l in lines if u'Number of 2D marginal pdfs to be calculated' in l),None)
        # # Check if the number of plots asked is compatible with the number of 2D marginals computed
        # if len(dataplot) > Nparam2D:
        #     raise ValueError('\033[91mERROR:\033[00m Number of plots (%s) > Number of 2D marginals computed (%s)' %(len(dataplot), Nparam2D))
        # # Record the 2D marginals couple present in the nab.out file
        # test_couple = [(int(line.split()[-2]), int(line.split()[-1])) for line in lines if u'2D Marginal for parameters :' in line]
        # # Find the number of bins for 2D marginals
        # Nbins2D = next((int(line[-13:len(line)]) for line in lines if u'Number of 2D marginal pdfs to be calculated' in line),None)
        for line in lines:
            # Get number of parameters
            if u'Number of dimensions' in line: 
                Nparam = int(line[-13:len(line)])
            if u'Number of bins per axis for 1D marginals' in line: 
                Nbins1D = int(line[-13:len(line)])

            # Get number of parameters 2D
            if u'Number of 2D marginal pdfs to be calculated' in line: 
                Nparam2D = int(line[-13:len(line)])
                # Check if the number of plots asked is compatible with the number of 2D marginals computed
                if len(dataplot) > Nparam2D:
                    QErrorMessage(self).showMessage('\033[91mERROR:\033[00m Number of plots (%s) > Number of 2D marginals computed (%s)' %(len(dataplot), Nparam2D))
                    return
                    # raise ValueError('\033[91mERROR:\033[00m Number of plots (%s) > Number of 2D marginals computed (%s)' %(len(dataplot), Nparam2D))
            # Record the 2D marginals couple present in the nab.out file
            if u'2D Marginal for parameters :' in line:
                test_couple.append((int(line.split()[-2]), int(line.split()[-1])))
            # Find the number of bins for 2D marginals
            if u'Number of bins per axis for 2D marginals' in line: 
                Nbins2D = int(line[-13:len(line)])
                
        print("after read through lines")
        # check if the couples asked to plot are OK with the marginals computed
        if not PDF1D_only:
            for elem in dataplot:
                if elem not in test_couple:
                    print(elem, test_couple)
                    QErrorMessage(self).showMessage('\033[91mERROR:\033[00m You did nto read all the comments in the beginning !!\n2D plot %s not in 2D marginals computed...\n\t\tPlease, rerun NAD with correct couples\n' %(str(elem)))
                    return
                    # raise ValueError('\033[91mERROR:\033[00m You did nto read all the comments in the beginning !!\n2D plot %s not in 2D marginals computed...\n\t\tPlease, rerun NAD with correct couples\n' %(str(elem)))
        
        print("after test couples")
        # intitiate variable
        # data1D = np.zeros((Nbins1D, 2*Nparam))
        # for i in range(0, Nparam, 1):
        #     # Find index of the beginning of the last iteration for the parameter i+1
        #     param_str = f'Marginalforparameter:{i + 1}\n'
        
        #     # Use a list comprehension to find the index of the parameter line
        #     zzz = 0
        #     indexparam = 0
        #     param_indices = [zzz for zzz, line in enumerate(lines) if line.replace(' ', '') == param_str]
        #     indexparam = param_indices[-1]
        #     # Copy the lines until the marginal i+1 in the variable data1D
        #     startline = indexparam + 2
        #     endline = startline + Nbins1D
        #     for k in range (startline, endline):
        #             data1D[k - indexparam -2, 2*i:2*i+2] = lines[k].split()[0:2]
        
        data1D = np.zeros((Nbins1D, 2*Nparam))
        for i in range(0, Nparam, 1):
            # Find index of the beginning of the last iteration for the parameter i+1
            zzz = 0
            indexparam = 0
            for zzz, z_str in enumerate([item.replace(u' ', u'') for item in lines]):
                if z_str == 'Marginalforparameter:' + str(i+1) +'\n':
                    indexparam = zzz
                    
            # Copy the lines until the marginal i+1 in the variable data1D
            for k in range (indexparam + 2, indexparam + 2 + Nbins1D):
                data1D[k - indexparam -2, 2*i:2*i+2] = lines[k].split()[0:2]
                
        print("after find indexes and copy lines")
        # If asked, save the data1D variable in a text file
        if pdf1d_results:
            if os.path.isfile(pdf1d_results):
                print('\t\t\033[91mWARNING:\033[00m File %s already exists \n \t\tI do not erase it\n' %(pdf1d_results))
            else:
                np.savetxt(pdf1d_results, data1D)
                print('\tFile %s written\n' %(pdf1d_results))

        # intitiate variable
        if not PDF1D_only:
            data2D = np.zeros((Nbins2D, Nbins2D*len(dataplot)))
            for i in range (0, len(dataplot)):
                # Find index of the beginning of the last iteration for the couple of parameters in dataplot
                try:
                    # We kept this line because it permits to raise an error we do not find the couple do plot.
                    # Normaly, this is already done before the call to thsi function...
                    #indexparam = [item.replace(u' ', u'') for item in lines].index('2DMarginalforparameters:' + str(dataplot[i][0]) + str(dataplot[i][1]) + '\n')
                    zzz = 0
                    indexparam = 0
                    # Find index of the beginning of the last iteration for the couple of parameters
                    for zzz, z_str in enumerate([item.replace(u' ', u'') for item in lines]):
                        if z_str == '2DMarginalforparameters:' + str(dataplot[i][0]) + str(dataplot[i][1]) + '\n':
                            indexparam = zzz
                            
                except ValueError:
                    QErrorMessage(self).showMessage('2D pdf couple asked do not agree with 2D pdf computed (%s, %s) ' %(str(dataplot[i][0]), str(dataplot[i][1])))
                    return
                    # raise ValueError('2D pdf couple asked do not agree with 2D pdf computed (%s, %s) ' %(str(dataplot[i][0]), str(dataplot[i][1])))
                # Copy the lines until the marginal i+1 in the variable ?
                for k in range (indexparam + 4, indexparam + 4 + Nbins2D):
                    data2D[k - indexparam - 4, (Nbins2D * i):(Nbins2D * i + Nbins2D)] = lines[k].split()
            
            # If asked, save the data2D variable in a text file
            if pdf2d_results:
                if os.path.isfile(pdf2d_results):
                    print('\t\t\033[91mWARNING:\033[00m File %s already exists \n \t\tI do not erase it\n' %(pdf2d_results))
                else:
                    np.savetxt(pdf2d_results, data2D)
                    print('\tFile %s written \n' %(pdf2d_results))
        else:
            data2D = data1D

        # Find parameter ranges (Uncomment the block if needed)
        #param_ranges = np.zeros((Nparam, 3))
        #indexparam = lines.index('  Parameter ranges\n')
        #for i in range (0, Nparam):
        #    param_ranges[i, 0:3] = lines[indexparam+2+i].split()[0:3]
        #
        #return data1D, data2D, param_ranges

        return data1D, data2D
    
    #----------------------------------------------------------
    def NAplot_Pecube(self,param, dataplot, tick_space,
                      inv_results = 'NA/NA_int_results.csv', data_nab = 'NA/NAB/nab.out',
                      graph_path  = 'NA/Graphs/',
                      PDF_1D = True, PDF_2D = False, PDF1D_only=False,pdf1d_results = None, pdf2d_results = None,
                      size_x = 15, size_y = 15, size_m = 15,
                      size_mis = 50,
                      peak_thres = 0.05, peak_min_dist = 30,misfitEvol = False):
        
        """
        Main code of NA Plot for Pecube inversions with NA

        INPUTS:
            param (list of str): Define as many variable as you have, with their unit
                                 Check the order in na.sum (open with text editor) or in NA_Results
                                 If you use the later, first column always is the misfit
                                 The others paramaters are from topo_parameters.txt then from fault_parameters.txt
                                 Exemple : param=['Offset (km)','Basal Temperature (°C)','Slip rate (km∕Ma)']
                                 /!\ IF YOU WANT TO USE A SLASH USE THIS ONE --> '∕' <-- , IT'S A UNICODE DIVISION SYMBOL
                                WINDOWS AND OSX DON'T ALLOW THE USE OF THE REGULAR SLASH 

            dataplot (list of couple of integers): Set the couple of variables to plot against each other
                                                   Exemple : Offset vs Slip rates
                                                   dataplot = [(1,3)] -->  plot=(1,3)
                                                   If you plot 2D pdfs (contours), please, CHECK that the couple of parameters to plot
                                                   are the same and in the same order than in the nab.in file.
                                                   This python script checks it and will insult you if this is not compatible !!!

            tick_space (array of floats): Set the space between ticks for x and y axes for each parameters
                                          (same order than the list param).
                                          If the tick format does not fit your variables, 
                                          you may need to modify the dictionnary tick_order in the function multiplot

            inv_results (str, optional): Name of the NA file with the inversion results; Usually NA_Results.
                                         Defaults to 'NA/NA_results.csv'.

            data_nab (str, optional): Name of the NAB file with the inversion results. Usually nab.out.
                                      Defaults to 'NA/NAB/nab.out'.

            graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                         Do not forget the '/' at the end.
                                         Defaults to 'NA/Graphs/'.

            PDF_1D (bool, optional): Choose if you want the 1-PDFs (Probability Density Function) 
                                     Defaults to True.

            PDF_2D (bool, optional): Choose if you want the 2-PDFs (Probability Density Function) 
                                     Defaults to False.

            pdf1d_results (string, optional): Print the 1-pdfs in a text file.
                                              Defaults to None.

            pdf2d_results (string, optional): Print the 2-pdfs in a text file.
                                              Defaults to None.

            size_x (int, optional): Size of the font for the x axes label. 
                                    Defaults to 15.

            size_y (int, optional): Size of the font for the y axes label. 
                                    Defaults to 15.

            size_m (int, optional): Size of the font for the markers label. 
                                    Defaults to 15.

            size_mis (int, optional): Size of the markers of the misfits. 
                                        Defaults to 50.

            peak_thres (float, optional): Threshold to find peaks; between 0. and 1.
                                          Default = 0.05.

            peak_min_dist (interger, optional): Minimum distance between the peaks
                                                Default = 30.
        
        Raises:
            NameError  : Input file names does not exists
            ValueError : Problems with number of values or with NaN
            ImportError: File format not supported, please, change the format
        """

        # If you want to print the pdfs in a text file, just modify the 2 next lines
        pdf1d_results = None #'NA/NAB/PDF_DATA.txt'
        pdf2d_results = None #'NA/NAB/PDF_2D_DATA.txt'

        print('###########################################################################\n')
        print('\t\tPlot results from NA inversions\n')
        print('\t\t  \xa9Xavier Robert - IRD-ISTerre\n')
        print('###########################################################################\n')
        # Check that NA_results exists
        if not os.path.isfile(inv_results):
            QErrorMessage(self).showMessage('\033[91mERROR:\033[00m File %s does not exist' %(str(inv_results)))
            return
            # raise NameError('\033[91mERROR:\033[00m File %s does not exist' %(str(inv_results)))
        # Check if the number of parameters is compatible with the number ticks_space
        if len(param) != len(tick_space):
            QErrorMessage(self).showMessage('\033[91mERROR:\033[00m the number of parameters is different than the number of ticks_space asked !\n Check it !')
            return
            # raise ValueError('\033[91mERROR:\033[00m the number of parameters is different than the number of ticks_space asked !\n Check it !')
        
        # Loading the data, calling as many variables as given in param + 1 for the misfit
        if inv_results[-3:] == 'txt':
            nb_var = np.loadtxt(inv_results, unpack = True)
        elif inv_results[-3:] == 'csv':
            nb_var = np.genfromtxt(inv_results, delimiter = ',', skip_header = 1, unpack = True)
        else:
            QErrorMessage(self).showMessage('NA_results format not supported...\n\n')
            return
            # raise ImportError('NA_results format not supported...\n\n')
        
        # Check if misfits are not with NaN or infinite values...
        # Raise error if this is the case, and recheck your inversion parameters
        if np.any(np.isnan(nb_var)) or not np.all(np.isfinite(nb_var)):
            print("\nInfinity: ", np.inf in nb_var)
            QErrorMessage(self).showMessage('NA_results.csv contains NaN or infinite numbers...\n\n')
            return
            # raise ValueError('NA_results.csv contains NaN or infinite numbers...\n\n')

        # Sort nb_var in function of misfit to get the best misfits over the lower ones 
        #          --> Better scatter plot visualisation
        nb_var_nosort = nb_var.T
        nb_var = nb_var.T[nb_var.T[:,0].argsort()[::-1]].T
        #nb_var.T[:] = nb_var.T[::-1]

        # Check that nab.out exists
        if not os.path.isfile(data_nab): 
            print ('\033[91mWARNING:\033[00m File %s does not exist\n\n' %(str(data_nab)))
            # PDF_1D = False
            # PDF_2D = False
            
        if not PDF_1D: 
            pdf1d_results = None
            data1D = None
        if not PDF_2D: 
            pdf2d_results = None
            data2D = None

        print('\x1b[32;1mBuilding NAB files for plotting...\x1b[0m')
        # Build Nab files for plotting
        self.isNAB = 0
        if PDF_1D or PDF_2D:
            # Check if nab.out exits
            if os.path.exists(data_nab):
                data1D, data2D = self.read_NABout(data_nab, dataplot,PDF1D_only, pdf1d_results, pdf2d_results)
                self.isNAB = 1
            else: # Read data from NA_results.csv or NA_int_results.csv
                # Read data, nb_var = data, already sorted
                shape = nb_var.shape
                nparam = shape[0]-1
                nbin = 30
                # Format of data1d is 2 columns per parameter (xpdf, ypdf)
                data1D = np.zeros((nbin,nparam*2))
                for p in range(nparam):
                    data = nb_var[p+1,:]
                    # nbin = self.get_bin_number(data)
                    heights,bins = np.histogram(data,bins=nbin)
                    indices_bin = np.digitize(data,bins)
                    # take misfit values for each bin
                    misfit_values = nb_var[0,:]
                    if int(self.inputParam['misfit_function']) == 3: # Means misfit is L2-norm
                        misfit_values = misfit_values*misfit_values 
                    binMids=bins[:-1]+np.diff(bins)/2.
                    misfit_bins = np.zeros((len(binMids),))
                    for k in range(nbin):
                        values = misfit_values[indices_bin==k+1]
                        # misfit_bins[k] = sum(values)/np.count_nonzero(values) #Mean misfit value for bin k
                        PDF_prior = 1/nbin
                        PDF_m = PDF_prior * np.exp(-0.5 * values)
                        PDF_m = PDF_m / sum(PDF_m + 1e-10)
                        if gmean(values) == gmean(values):
                            misfit_bins[k] = gmean(values)
                        else:
                            misfit_bins[k] = 1000
                        # misfit_bins[k] = sum(PDF_m) / (len(PDF_m)+1e-10) # Mean Posterior probability
                        
                    data1D[:nbin,2*p] = binMids
                    misfits =1/misfit_bins
                    data1D[:nbin,2*p+1] = (misfits / sum(misfits))
                    # data1D[:nbin,2*p+1] = misfit_bins/sum(misfit_bins)
                
                data2D = data1D
        # Plot data
        # Check if there is a folder to store the outputs
        print('\x1b[32;1m\nPlotting results...\x1b[0m')
        # if not os.path.exists(graph_path): os.mkdir(graph_path)
        if misfitEvol: # Plot Misfit evolution
            nb_models = len(nb_var_nosort[:,0])
            models = np.linspace(1,nb_models,nb_models)
            # Create the canva
            self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
            # Create a customed toolbar to show/hide ages from specific thermochronometer
            self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
            if self.logValue:
                self.plotSpace.axes.scatter(x=models,y=np.log10(nb_var_nosort[:,0]/int(self.NormalizationValue.text())),s=5,c=[0,0,0])
            else:
                self.plotSpace.axes.scatter(x=models,y=nb_var_nosort[:,0]/int(self.NormalizationValue.text()),s=5,c=[0,0,0])
            self.plotSpace.axes.set_xlabel('Model ID', fontsize=size_x)
            self.plotSpace.axes.set_ylabel('Misfit', fontsize=size_y)

            # close the output text file
            # fstats_w.close()
            # Add plot to interface
            VBox = QVBoxLayout()
            VBox.addWidget(self.toolbar)
            VBox.addWidget(self.plotSpace)
            Q = QWidget()
            Q.setLayout(VBox)
            self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'Misfit_evolution_')
            
        elif len(dataplot) > 1:
            for i in range(len(dataplot)):
                # Call the plot function for each parameter value to plot
                self.multiplot(param, nb_var, dataplot[0], 
                          tick_space, 
                          data1D, data2D,
                          size_x, size_y, size_m,
                          size_mis, 
                          PDF_1D, PDF_2D, graph_path, i)
        elif not PDF1D_only:
            self.multiplot(param, nb_var, dataplot[0],
                      tick_space, 
                      data1D, data2D,
                      size_x, size_y, size_m,
                      size_mis, 
                      PDF_1D, PDF_2D, graph_path)
        
        # Compute the stats, print them
        print('\x1b[32;1mComputing statistics...\x1b[0m')
        if not PDF_1D:
            print('\033[91mWARNING:\033[00m No 1D-pdf computed --> No stats computed !')
        elif PDF1D_only:
            self.statsparam(data1D, dataplot,param, graph_path, peak_thres, peak_min_dist, size_x, size_y)

        print('###########################################################################\n')
        print('\t\tEnd...\n')
        print('###########################################################################\n')

        return
    
    #----------------------------------------------------------
    def get_bin_number(self,data):
        """ Calculate the optimal number of bin for the dataset. """
        
        # calculate number of bins
        q3, q1 = np.percentile(data,[75, 25])
        IQR = q3 - q1
        hbin = 2* (IQR/np.cbrt(len(data)))
        nbin = round((max(data) - min(data)) / hbin)
        if nbin > 50:
            nbin = 50
        
        return nbin
    
    #----------------------------------------------------------
    def statsparam(self,data1D, dataplot,param, graph_path  = 'NA/Graphs/', peak_thres = 0.05, peak_min_dist = 30, size_x = 15, size_y = 15):
        """
        Function to find the number of gaussian functions, to compute their mean and sigma,
        and to plot the gaussians fit

        INPUTS:
            data1D ([np.array])   : Array of 1D-pdfs for each parameter

            param (list of string): List the name of the parameters 
                                    in the same order than define in Pecube inversion

            graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                         Do not forget the '/' at the end.
                                         Defaults to 'NA/Graphs/'.

            peak_thres (float, optional): Threshold to find peaks; between 0. and 1.
                                          Default = 0.05.

            peak_min_dist (interger, optional): Minimum distance between the peaks
                                                Default = 30.

            size_x (int, optional): Size of the font for the x axes label. 
                                    Defaults to 15.

            size_y (int, optional): Size of the font for the y axes label. 
                                    Defaults to 15.

        """
        
        # # Open the text file where to put the results
        # fstats_w = open(graph_path + 'Nab-stats.txt', 'w')
        # # Write the Header
        # fstats_w.write('Param \t Mean \t Mean_stdev \t Std_Err \t Std_Err_Stdev \n')
        
        # Create the canva
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        # Create a customed toolbar to show/hide ages from specific thermochronometer
        self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
        ax = self.plotSpace.axes

        # Faire la boucle sur i in range (0, len(param)):
        count = 0
        i = dataplot[0][0]-1
        print('\033[96m\t%s\033[00m' %(param[count]))
        # initiate fig
        # self.plotSpace.axes.clf()
        # print the 1D-pdf
        self.plotSpace.axes.plot(data1D.T[2*i], data1D.T[2*i+1], "mo", label = 'NA results')

        # find number of peaks/gaussian
        print("Peak threshold = ", peak_thres, peak_min_dist)
        indexes = peakutils.indexes(data1D.T[2*i+1], 
                                    thres = peak_thres, # Normalized threshold. Only the peaks with amplitude higher than the threshold will be detected.
                                    min_dist = peak_min_dist)
                                    #min_dist = (max(data1D.T[2*i])-min(data1D.T[2*i]))/2) # Minimum distance between each detected peak.
                                    #min_dist = max(data1D.T[2*i])) # Minimum distance between each detected peak.
        #pplot(data1D.T[2*i], data1D.T[2*i+1], indexes)
        #plt.show()
        print('\t\tNumber of peak(s) : %s' %(indexes.shape[0]))

        # Fit the Gaussian data
        # Compute the mean and 1 sigma erro used for the first guess of the fit
        mean = sum(data1D.T[2*i] * data1D.T[2*i+1]) / sum(data1D.T[2*i+1])
        sigma = np.sqrt(sum(data1D.T[2*i+1] * (data1D.T[2*i] - mean) ** 2) / sum(data1D.T[2*i+1]))
        # REM: to avoid RuntimeError from curve_fit, it is easier to divide sigma by ~10;
        #      this is done inside the curve_fit parameters

        # Do the fitting
        try:
            pars, cov = curve_fit(f = self.gaussians,
                                xdata = data1D.T[2*i], ydata = data1D.T[2*i+1], 
                                p0 = np.array([np.stack([data1D.T[2*i+1][indexes[k]], data1D.T[2*i][indexes[k]], sigma/10]) 
                                              for k in range(0, indexes.shape[0])]).reshape(-1),
                                bounds=(-np.inf, np.inf),
                                check_finite = True)
                                #maxfev = 5000)
            i_opt = True
        except RuntimeError:
            print(u'\t\t\033[91mWarning:\033[00m No optimization found with the least-square method; I am trying the dogbox method')
            try:
                pars, cov = curve_fit(f = self.gaussians,
                                    xdata = data1D.T[2*i], ydata = data1D.T[2*i+1], 
                                    p0 = np.array([np.stack([data1D.T[2*i+1][indexes[k]], data1D.T[2*i][indexes[k]], sigma/10]) 
                                                  for k in range(0, indexes.shape[0])]).reshape(-1),
                                    bounds=(-np.inf, np.inf),
                                    check_finite = True,
                                    method = 'dogbox',
                                    maxfev = 5000)
                i_opt = True
            except RuntimeError:
                print(u'\t\t\033[91mWarning:\033[00m No optimization found also with the dogbox method; No graph will be produced')
                i_opt = False
                pass
        except TypeError:
            print(u'\t\t\033[91mWarning:\033[00m No optimization found')
            i_opt = False

        if i_opt:
            # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
            stdevs = np.sqrt(np.diag(cov))
            # Calculate the residuals
            res = data1D.T[2*i+1] - self.gaussians(data1D.T[2*i], *pars)

            for k in range(0,indexes.shape[0]):
                # Compute the result for each peak
                pars_1 = pars[(3*k) : (3*k+3)]
                xerror = (max(data1D.T[2*i]) - min(data1D.T[2*i])) / 2
                if (pars_1[1] > (min(data1D.T[2*i]) - xerror)) & (pars_1[1] < (max(data1D.T[2*i] + xerror))):
                    # print on screen results
                    print('\t\tMean & sigma (Gaussian fit %s, %0.2f %s) : %0.2f +/- %0.2f' %(str(k+1), pars_1[0]*100, chr(37), pars_1[1], abs(pars_1[2])))
                    # Save the results in text file for each peak
                    line = (str(param[count]) + 'peak ' + str(k+1) + '\t' + str(pars[1]) + '\t' + str(stdevs[1]) +
                                           '\t' + str(pars[2]) + '\t' + str(stdevs[2]) +
                                           '\n')
                    # fstats_w.write(line)
                    # Plot the gaussian fitted
                    if pars_1[0] > 0:
                        self.plotgfit(data1D = data1D, i = i, pars = pars_1, ipeak = k+1)
                    else:
                        print(u'\t\t\t\033[91mWarning:\033[00m Negative amplitude for the fit %s...\n\t\t\t\033[91m-->\033[00m Not plotted' %(str(i)))
                else:
                    print(u'\t\t\t\033[91mWarning:\033[00m peak value outside of acceptable bounds ...\n\t\t\t\033[91m-->\033[00m Not plotted' )


            # Saving the plots as a pdf file
            # plt.savefig(graph_path + 'PDF-1D_param' + str(i+1) + '.pdf')
            print('\t\tPlotting PDF: '+ graph_path + 'PDF-1D_param' + str(count) + '.pdf\n')
            
        # Set axis names
        self.plotSpace.axes.set_xlabel(param[count], fontsize=size_x)
        self.plotSpace.axes.set_ylabel('Probability', fontsize=size_y)
        self.plotSpace.axes.set_xlim(min(data1D.T[2*i]), max(data1D.T[2*i]))
        self.plotSpace.axes.set_ylim(0, None)
        self.plotSpace.axes.legend(loc = 'best')

        # close the output text file
        # fstats_w.close()
        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'PDF_1D_')

        return
    
    #--------------------------------------------------------------------------
    def plotgfit(self,data1D, i, pars, ipeak = None):
        """
        Plot the gaussian fitted with the peak and the 1 sigma error

        Args:
            data1D (array)           : 1D-pdf data
            i (integer)              : Number of the parameter to plot
            pars (array)             : Gaussian(s) fit results for the parameter i
            ipeak (integer, optional): Number of peaks/gaussians. 
                                        Defaults to None.
        """
        
        if ipeak == 1 : labelname = '1-peak Gaussian fit'
        else: labelname = None
        self.plotSpace.axes.plot(data1D.T[2*i], 
                self.gaussian(data1D.T[2*i], pars[0], pars[1], pars[2]),
                "-b",
                label = labelname)

        # Colorize the gaussian between the 1 sigmas
        if ipeak == 1 : labelname = 'Acceptable values\n' + '(1\u03C3)'
        else: labelname = None
        self.plotSpace.axes.fill_between(data1D.T[2*i], 
                         self.gaussian(data1D.T[2*i], *pars),
                         0,
                         where = ((data1D.T[2*i] >= pars[1] - abs(pars[2])) & (data1D.T[2*i] <= pars[1] + abs(pars[2]))),
                         alpha=0.30, 
                         color='green', 
                         interpolate=True,
                         label = labelname)
                
        # Plot the mean
        self.plotSpace.axes.vlines(x = pars[1], 
                   ymin = 0, 
                   ymax = max(self.gaussian(data1D.T[2*i], *pars)),
                   colors = 'red',
                   label = 'Mean value %s (%0.1f %s)\n%0.2f +/- %0.2f' %(ipeak, pars[0]*100, chr(37), pars[1], pars[2]))
        #           label = 'Mean value %s (%0.2f %s)\n%0.2f +/- %0.2f' %(ipeak, pars[0]*100, chr(37), pars[1], pars[2]))

        # End of plot
        return
    
    #----------------------------------------------------------
    def gaussians(self,x, *gaussians_param):#amp1, cen1, sigma1, amp2, cen2, sigma2):
        """
        Compute the sum of several Gaussian functions

        INPUTS:
            x ([float])                               : data to desccribe
            gaussians_param (list of floats triplets) : list of tripplets that describes each gaussian pdf
                                                        with the amplitude, the mean and the 1 sigma for each pdf

        RETURNS:
            gaussians_results (np.array of floats): results of the sum of several Gaussian functions.
                                                    Th size is the same than the input vector x.
        """
        
        # Clear the variable
        gaussians_results = np.zeros(x.shape)
        # Do a loop on the number of peaks --> Determines the number of gaussians to stack
        for h in range(0, int(len(gaussians_param)/3)):
            gaussians_results = gaussians_results + self.gaussian(x, gaussians_param[3*h], 
                                                                gaussians_param[3*h+1], 
                                                                gaussians_param[3*h+2])

        return gaussians_results
    
    #----------------------------------------------------------
    def gaussian(self,x, ampg, meang, sigma1g):
        """
        Function to calculate the Gaussian with constants ampg, meang, and sigmag

        Args:
            x ([float])      : data to describe
            ampg ([float])   : amplitude of the Gaussian
            meang ([float])  : mean of the Gaussian
            sigma1g ([float]): 1-sigma1 of the Gaussian

        Returns:
            1 Gaussian pdf
        """
        
        return ampg * np.exp(-np.power(x - meang, 2)/(2 * np.power(sigma1g, 2)))
    
    #----------------------------------------------------------
    def multiplot(self,param, nb_var, plot,
                  tick_space, 
                  data1D = None, data2D = None,
                  size_x = 15, size_y = 15, size_m = 15,
                  size_mis = 50, 
                  PDF_1D = True, PDF_2D = True,
                  graph_path = 'NA/Graphs/',
                  i_param = None):
        """
        Function to plot the scatter plots with the 1D and 2D pdfs if asked.

        INPUTS:
            param       : name of the parameters
            nb_var      : NA_Results data
            plot        : parameter to plot
            tick_space  : Space between the ticks
            data1D      : 1D marginal. Defaults to None.
            data2D      : 2D marginal. Defaults to None.
            size_x      : Font size for x-axis. Defaults to 15.
            size_y      : Font size for y-axis. Defaults to 15.
            size_m      : Font size for ???. Defaults to 15.
            size_mis    : Size of the plot. Defaults to 50.
            PDF_1D      : Boolean to tell if we plot 1D-pdf (True) or not (False).
                          Defaults to True. Optional.
            PDF_2D      : Boolean to tell if we plot 2D-pdf (True) or not (False). 
                          Defaults to True. Optional.
            graph_path (str, optional) : Path where to save graphs and results Usually NA/Graphs.
                                         Do not forget the '/' at the end.
                                         Defaults to 'NA/Graphs/'.
            i_param     : Index of couple of parameters to plot. 
                          Defaults to None. Optional.

        OUTPUTS:
            Just pdf graphs !
        """
        
        # Define dictionnary to manage tick format
        # You may add you own settings if needed
        tick_order = {1000  : u'%4.0f',
                      100   : u'%4.0f',
                      10    : u'%4.0f',
                      1     : u'%4.1f',
    	              0.1   : u'%4.1f',
                      0.01  : u'%4.2f',
                      0.001 : u'%4.3f'}

        if i_param or i_param == 0:
            print('\033[96m\tPlotting :\033[00m ' + str(param[plot[0]-1]) + ' vs ' + str((param[plot[1]-1])) + ' (' + graph_path + 'PDF_' + str(i_param+1) + '.pdf)')
        else: 
            print('\033[96m\tPlotting :\033[00m ' + str(param[0]) + ' vs ' + str((param[1])) + ' (' + graph_path + 'PDF.pdf)')

        # Create the canva
        self.plotSpace = pgu.PlotCanvas(self, width=5, height=4)
        # Create a customed toolbar to show/hide ages from specific thermochronometer
        self.toolbar = NavigationToolbar(self.plotSpace,self.GraphWin)
        
        # Plot the results as a scatter plot

        # set the misfit value as a variable for the color
        ax = self.plotSpace.axes
        self.plotSpace.fig.clear()
        
        
        # Find the lowest misfit run and return the line/run number
        if np.isnan(min(nb_var[0])):
            QErrorMessage(self).showMessage('\t\t\033[91mERROR:\033[00m one of the misifit is NaN; \n\tthis is probably because you had an issue with the multithread inversion process...\n\tRe-run the inversion with 1 core only.')
            return
            # raise ValueError('\t\t\033[91mERROR:\033[00m one of the misifit is NaN; \n\tthis is probably because you had an issue with the multithread inversion process...\n\tRe-run the inversion with 1 core only.')
            # This error is known with clusters using OAR protocols.
            # It seems that there are errors when nodes discuss with themselves, 
            # Some misfit became NaN, and the misfits printed in the results files 
            # does not trully correspond to the parameters values associated in those files.
            # BE CAREFULL !
        else:
            issue_zero = 0
            if  min(nb_var[0]) != 0.0:
                bfit = np.where(nb_var[0] == min(nb_var[0]))[0][0]
            else: # find first value non zero
                issue_zero = 1
                bfit = np.where(nb_var[0] == 0.0)
                nb_var2 = nb_var[0]
                nb_var2[bfit]  = max(nb_var[0]) 
                nb_var[0][bfit] = max(nb_var[0])
                bfit = np.where(nb_var2 == min(nb_var2))[0][0]
                
            print("\t\tBest model is number %s with misfit of %s" %(str(bfit + 1), str(min(nb_var[0]))))
                  
        if self.logValue:
            nb_var[0] = np.log10(nb_var[0]/int(self.NormalizationValue.text()))
        else:
            nb_var[0] = nb_var[0]/int(self.NormalizationValue.text())
        
        if issue_zero:
            nb_var = nb_var.T[nb_var.T[:,0].argsort()[::-1]].T
            
        gs = GridSpec(4,5)
        ax_scatter = self.plotSpace.fig.add_subplot(gs[1:4, 0:3])
        ax_hist_x = self.plotSpace.fig.add_subplot(gs[0,0:3])
        ax_hist_y = self.plotSpace.fig.add_subplot(gs[1:4, 3])
        ax_cb = self.plotSpace.fig.add_subplot(gs[2:3,4])
        scatter_mis = ax_scatter.scatter(nb_var[plot[0]], 
                                  nb_var[plot[1]], 
                                  s = size_mis,            
                                  linewidth = 0.1,   
                                  edgecolor = 'k',   
                                  c = nb_var[0],  # data use for the color bar, here as a log10
                                  cmap = 'Spectral')
                               
        # Plot the best (lowest misfit) run as a yellow star
        ax_scatter.plot(nb_var[plot[0]][bfit],
                  nb_var[plot[1]][bfit],
                  '*',
                  markeredgewidth = 1,
                  markeredgecolor = 'k',
                  markersize = 25,
                  c = '#ffff00')

        # Set the misfit legend as a color bar
        cbar = self.plotSpace.fig.colorbar(scatter_mis,ax = ax_cb, fraction = 1,shrink = 2, aspect=10, label = 'Misfit')
        cbar.ax.set_ylabel('Misfit',
                           rotation = 270,
                           size = size_m,
                           labelpad = 20)
        # Find the right misfit digits to print
        for x_order in (1000, 100, 10, 1, 0.1, 0.01, 0.001):
            if int(max(nb_var[0])-min(nb_var[0])/x_order) in range (1,10):
                cbar.ax.yaxis.set_major_formatter(FormatStrFormatter(tick_order[x_order/10]))
        #cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
        
        ################
        # Axes settings
        ################

        # set the name of the axes
        ax_scatter.set_xlabel(param[0], fontsize=size_x)
        ax_scatter.set_ylabel(param[1], fontsize=size_y)
        
        ############################
        # PDFs Histogram plot code #
        ############################

        # Start 1D PDF or/and 2D PDF    
        if PDF_2D and self.isNAB:
            print('\t\tBuild 2D pdfs countours')  

            # Load the 2D pdf data and assing variables to it
            val_z2d = np.zeros((np.shape(data2D.T)[1], np.shape(data2D.T)[1]))  # empty array to stock the pdf results
            val_x2d = data1D.T[(plot[0]-1)*2]
            val_y2d = data1D.T[(plot[1]-1)*2]
            if i_param or i_param == 0:
                val_z2d[0:100, 0:100] = data2D.T[(i_param+1) * np.shape(data2D.T)[1] - np.shape(data2D.T)[1]:(i_param+2) * np.shape(data2D.T)[1] - np.shape(data2D.T)[1], :]
            else:
                val_z2d = data2D.T

            val_xx2d, val_yy2d = np.meshgrid(val_x2d, val_y2d)
            try:
                CS = ax_scatter.contour(val_xx2d, val_yy2d, val_z2d.T, 
                                 #levels = [0.4], cmap = 'CMRmap')
                                 #levels = 3, cmap = 'CMRmap')
                                 #levels = 3, cmap = 'gray_r')
                                 levels = 3, cmap = 'hot')
                ax_scatter.clabel(CS, inline = 2, fontsize = 12)
            except UserWarning:
                print('\t\t\033[91mWARNING:\033[00m No contour levels plotted')
        
        # Start 1D PDF or/and 2D PDF    
        if PDF_1D:
            print('\t\tBuild 1D pdfs')
            
            # Position the 2 new figures     
        
            # Disable the labels for the axes they share with the scatter plot    

            # Loop to select the right columns for the right parameters
            print(plot)
            if plot[0] == 1:
                ax_hist_x.plot(data1D.T[0], data1D.T[1],color='k')
            else:
                ax_hist_x.plot(data1D.T[plot[0] + (plot[0] - 2)], 
                             data1D.T[plot[0] + (plot[0] - 1)],color='k')
            if plot[1] == 1:
                ax_hist_y.plot(data1D.T[1], data1D.T[0],color='k')
            else:
                ax_hist_y.plot(data1D.T[plot[1] + (plot[1] - 1)], 
                             data1D.T[plot[1] + (plot[1] - 2)],color='k')

            ax_hist_x.set_ylim(ymin = 0)
            ax_hist_y.set_xlim(xmin = 0)
            ax_hist_y.set_yticks([])
            ax_hist_x.set_xticks([])
            ax_cb.axis('off')
        self.plotSpace.axes.set_box_aspect(1)
        
        # Add plot to interface
        VBox = QVBoxLayout()
        VBox.addWidget(self.toolbar)
        VBox.addWidget(self.plotSpace)
        Q = QWidget()
        Q.setLayout(VBox)
        self.GraphWin.add_plot2win(Q,self.plotSpace,self.toolbar,'Inversion2D')
        
        return