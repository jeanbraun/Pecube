""" 

 This module gathers all parameters, class and functions related to helium
thermochronometry.

@author: maxime Bernard

"""



from PyQt5.QtWidgets import (QPushButton,QWidget,QLabel,QCheckBox,QVBoxLayout,
                             QGridLayout,QComboBox,QGroupBox,QHBoxLayout,
                             QErrorMessage,QLineEdit,QTableWidgetItem,QRadioButton,
                             )
from PyQt5.QtGui import (QIntValidator)
import configs as conf
from PyQt5.QtCore import Qt
import main as pgui
import Utils.PGUI_utils as pgu
import Utils.interface as interface
import os
import numpy as np
import xarray as xr
import Thermochronology.helium_43 as He43


#############################################################################################
################################## apatite ##################################################
#############################################################################################
class apatite(QWidget):
    """
      This class set the parameters for apatite (U-Th)/He thermochronometry.
      
    """
    def __init__(self, AgeParam, d, Param, PFolder):
        super().__init__()

        # Project folder name
        self.PFolder = PFolder
        self.Param = Param
        self.d = d
        self.AgeParam = AgeParam
        self.parent = AgeParam
        self.Winparent = self.parent.parent
        self.grains = {}
        self.ID = 'AHE'
         
        # Check if old input
        if self.Param.DParameters['oldInput'] == '1':
            self.oldInput = 1
        else:
            self.oldInput = 0
        # To store samples
        self.LabelHe = QLabel('Set model parameters for He age computation:')
        self.LabelHe.setFont(conf.fontBold12)
        self.LabelHe.setAlignment(Qt.AlignLeft)
        # Choose Production-diffusion model
        self.ModelPDLabel = QLabel('Production-diffusion model:')
        self.ModelPDLabel.setFont(conf.font10)
        self.ModelPDLabel.setAlignment(Qt.AlignRight)
        self.ModelPD1 = QRadioButton("Finite difference")
        pgui.store_input(self.Winparent, {'PDModel_signal': str(1)})
        # self.ModelPD2 = QRadioButton("Monte Carlo")
        #self.ModelPD2.setToolTip("Production-diffusion model from Gautheron et al. (2010)")
        self.ModelDiff = QLabel('Diffusion model: ')
        self.ModelDiff.setFont(conf.font8)
        self.ModelDiff.setAlignment(Qt.AlignRight)
        self.MDCombo = QComboBox()
        self.MDCombo.addItem('Farley et al. (2000)')
        self.MDCombo.addItem('Shuster et al. (2006)')
        self.MDCombo.addItem('Gautheron et al. (2009)')
        self.MDCombo.addItem('Flowers et al. (2009) - RDAAM')
        # self.MDCombo.addItem('Willett et al. (2017) - ADAM')
        
        self.MDCombo.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Diffusion_model_AHe']]))
        self.D0 = QLabel('D0 (cm²/s):')
        self.D0.setFont(conf.font8)
        self.D0.setAlignment(Qt.AlignCenter)
        self.D0Value = QLineEdit(self)
        self.D0Value.setText(Param.DParameters[conf.Variable_names['D0_Apatite']])
        self.D0Value.setAlignment(Qt.AlignCenter)
        self.D0Value.setValidator(conf.DoubleValidator)
        self.ActEner = QLabel('Ea (kJ/mol):')
        self.ActEner.setFont(conf.font8)
        self.ActEner.setAlignment(Qt.AlignCenter)
        self.ActEnerValue = QLineEdit(self)
        self.ActEnerValue.setText(Param.DParameters[conf.Variable_names['Ea_Apatite']])
        self.ActEnerValue.setAlignment(Qt.AlignCenter)
        self.ActEnerValue.setValidator(conf.DoubleValidator)
        self.NbIter = QLabel('Number of iterations:')
        self.NbIter.setFont(conf.font8)
        self.NbIter.setToolTip("Number of iteration for the Monte Carlo simulation.")
        self.NbIter.setAlignment(Qt.AlignCenter)
        self.NbIterEdit = QLineEdit()
        self.NbIterEdit.setValidator(QIntValidator())
        self.NbIterEdit.setText(self.Param.DParameters['NMonteIter'])
        self.alphaEjec = QLabel('Alpha ejection:')
        self.alphaEjec.setFont(conf.font8)
        self.alphaEjec.setAlignment(Qt.AlignRight)
        self.alphaCombo = QComboBox()
        self.alphaCombo.addItem('no ejection')
        self.alphaCombo.addItem('ejection')
        self.alphaCombo.addItem('redistribution')
        # self.alphaCombo.setCurrentIndex(int(self.Param.DParameters['Alpha_Ejec']))
        self.alphaDist = QLabel('stopping distance:')
        self.alphaDist.setFont(conf.font8)
        self.alphaDist.setAlignment(Qt.AlignRight)
        self.alphaDistCombo = QComboBox()
        self.alphaDistCombo.addItem('none')
        self.alphaDistCombo.addItem('Farley et al. (1996)')
        self.alphaDistCombo.addItem('Ketcham et al. (2011)')
        self.alphaDistCombo.setCurrentIndex(int(self.Param.DParameters['Alpha_Ejec_Flag']))
        self.alphaDistCombo.setFont(conf.font8)
        self.SaveButton = QPushButton('save samples file...')
        self.SaveButton.setFont(conf.font8)
        self.CheckGrains = QCheckBox('Grains characteristics')
        self.Check43He = QCheckBox('4He/3He predictions')
        self.Check43He.setToolTip("Enabling 4He/3He predictions")
        data = {'Grain ID': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Radius (µm)": ["0"]*int(Param.DParameters['ngrains']),
                "U (ppm)": ["0"]*int(Param.DParameters['ngrains']),
                "Th (ppm)": ["0"]*int(Param.DParameters['ngrains']),
                "rmr0 ": ["0"]*int(Param.DParameters['ngrains'])}
        #How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['AHe_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        self.GrainTable = interface.WinTable(data, self.TotGrain, 7)
        self.updateTableGrain(self.GrainTable, 1, Param)
        #Grain characterisics
        # self.GrainParam = self.AgeParam.GrainParam
        self.He43Param = He43.apatite(self.AgeParam,self, self.Param,self.d.nbSampleValue,self.d.Samples,self.PFolder)

        #if old input load values
        if self.oldInput:
            self.MDCombo.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Diffusion_model_AHe']]))
            self.D0Value.setText(self.Param.DParameters[conf.Variable_names['D0_Apatite']])
            self.ActEnerValue.setText(self.Param.DParameters[conf.Variable_names['Ea_Apatite']])
            # self.alphaCombo.setCurrentIndex(int(self.Param.DParameters['Alpha_Ejec_Mod_Flag']))
            self.alphaDistCombo.setCurrentIndex(int(self.Param.DParameters['Alpha_Ejec_Flag']))
            # self.CheckGrains.setChecked(True)
            if self.Param.DParameters[conf.Variable_names['4He/3He']] == '1':
                self.Check43He.setChecked(True)
                self.He43Param = He43.apatite(self.AgeParam, self, self.Param,self.d.nbSampleValue,self.d.Samples,self.PFolder)
            # self.GrainParam = self.AgeParam.GrainParam
            
        # Check for edited text
        self.D0Value.editingFinished.connect(lambda: pgui.store_input(
            d, {conf.Variable_names['D0_Apatite']: str(self.D0Value.text())}))
        self.ActEnerValue.editingFinished.connect(lambda: pgui.store_input(
            d, {conf.Variable_names['Ea_Apatite']: str(self.ActEnerValue.text())}))
        self.NbIterEdit.editingFinished.connect(lambda: pgui.store_input(
            d, {'NMonteIter': str(self.NbIterEdit.text())}))
        
        #Define Splitters
        self.ParamSplit = QWidget()
        
        # Extra box
        self.VBox1 = QVBoxLayout()
        self.VBox2 = QVBoxLayout()
        self.GroupBox1 = QGroupBox("Specific kinetics")
        self.GridHe = QGridLayout()
        self.GridHe.setSpacing(10)
        self.Hbox1 = QHBoxLayout()
        self.Hbox1.addWidget(self.ModelPDLabel)
        self.Hbox1.addWidget(self.ModelPD1)
        # self.Hbox1.addWidget(self.ModelPD2)
        self.Hbox1.addStretch(1)
        self.GridHe.addWidget(self.ModelDiff, 1, 0)
        self.GridHe.addWidget(self.MDCombo, 1, 1)
        self.GridHe.addWidget(self.D0, 2, 1)
        self.GridHe.addWidget(self.D0Value, 3, 1)
        self.GridHe.addWidget(self.ActEner, 2, 2)
        self.GridHe.addWidget(self.ActEnerValue,3, 2)
        # self.GridHe.addWidget(self.rmr0, 2, 3)
        # self.GridHe.addWidget(self.rmr0Value, 3, 3)
        self.GridHe.addWidget(self.alphaDist, 4, 0)
        self.GridHe.addWidget(self.alphaDistCombo, 4, 1)
        # self.GridHe.addWidget(self.CheckGrains, 6, 1)
        self.GridHe.addWidget(self.Check43He, 6, 1)
        self.GridHe.setColumnStretch(4, 2)
        # self.VBox2.addLayout(self.Hbox1)
        self.VBox2.addLayout(self.GridHe)
        self.GroupBox1.setLayout(self.VBox2)
        self.VBox1.addWidget(self.GroupBox1)
        self.VBox1.addWidget(self.GrainTable)
        self.VBox1.setStretch(1,2)
        
        self. updateWidgets()
        self.Q = QWidget()
        self.Q.setLayout(self.VBox1)
        
        # Check for updates
        # When window quit on 'ok' open table for grains
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.CheckGrains.stateChanged.connect(lambda: self.CheckGrainsChanged())
        self.Check43He.stateChanged.connect(lambda: self.Check43HeChanged())
        # self.ModelPD1.toggled.connect(lambda: self.update_PDModel(self.ModelPD1))
        # self.ModelPD2.toggled.connect(lambda: self.update_PDModel(self.ModelPD2))
        self.MDCombo.currentIndexChanged.connect(
            lambda: self.store_input_ComboBoxDiffModel(self.MDCombo, conf.Variable_names['Diffusion_model_AHe'], d))
        # self.alphaCombo.currentIndexChanged.connect(
        #     lambda: self.store_input_ComboBoxAlpha(self.alphaCombo, 'alpha_ejection', d))
        self.alphaDistCombo.currentIndexChanged.connect(
            lambda: self.store_input_ComboBoxAlphaDist(self.alphaDistCombo, 'Alpha_Ejec_Flag', d))
        self.SaveButton.clicked.connect(lambda: self.savefile())
        if int(self.Param.DParameters['PDModel_signal']) == 1:
            self.ModelPD1.setChecked(True)
            # self.ModelPD2.setChecked(False)
        else:
            self.ModelPD1.setChecked(False)
            # self.ModelPD2.setChecked(True)
    
    #--------------------------------------------------------
    def updateWidgets(self):
        """To enabke/disable widgets according to predictions wanted ("for all nodes" or
        sample specific")"""
        if self.Winparent.AgeCombo.currentIndex() == 2:#Means we want grain-specific
            self.GroupBox1.setEnabled(True)
            self.MDCombo.clear()
            self.MDCombo.addItem('Farley et al. (2000)')
            self.MDCombo.addItem('Shuster et al. (2006)')
            self.MDCombo.addItem('Gautheron et al. (2009)')
            self.MDCombo.addItem('Flowers et al. (2009) - RDAAM')
            self.MDCombo.addItem('Willett et al. (2017) - ADAM')
            self.MDCombo.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Diffusion_model_AHe']]))
            # self.updateTableGrain(self.GrainTable, 1,self.Param)
            
        else:
            #We are not in inversion mode
            self.MDCombo.clear()
            self.MDCombo.addItem('Farley et al. (2000)')
            self.MDCombo.addItem('Shuster et al. (2006)')
            self.MDCombo.addItem('Gautheron et al. (2009)')
            self.MDCombo.addItem('Flowers et al. (2009) - RDAAM')
            self.MDCombo.addItem('Willett et al. (2017) - ADAM')
            self.MDCombo.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Diffusion_model_AHe']]))
            
    #--------------------------------------------------------
    def updateTableGrain(self, t, h, p):
        """
          To store input value for the Grain table
          t = table
          h : parameter
          p : parameters
        """
        self.grains = {}
        h_temp = 1
        nSample = 0
        if h == 2: # Update from Table Samples
            try:
                self.GrainTable.cellChanged.disconnect() #Shut down update from the user
            except TypeError:
                print('TypeError in disconnect cell changed in age table for AHe.')
                pass
            # if (self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 6 or
            # self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 0):
            #     h_temp = 2
            #     h = 1
            # Check if we are in "for all nodes"
            if self.parent.parent.AgeCombo.currentIndex() == 1:
                self.TotGrain = 1
                t.setRowCount(1)
                p.TableParameters["SampleName1"] = "Grain"
            elif self.parent.parent.AgeCombo.currentIndex() == 2:
                # Get total number of grains
                self.TotGrain = self.Winparent.Samples['AHe_TotGrain']
                t.setRowCount(int(self.TotGrain))
                # Get Sample ID for each grain, read from sample table
                for j in range(int(self.Winparent.nbSampleValue.text())):
                        nb_grains = int(self.Winparent.Samples['AHe_ngrains'+str(j+1)])
                        if nb_grains > 0:
                            sampleName = self.Winparent.Samples['SampleName'+str(j+1)]
                            p.TableParameters["SampleName"+str(j+1)] = sampleName
        
        # Update Table
        if h >= 1:
          # dParam = DefaultParameterValues()
          # Check if we are in "for all nodes"
          nb_grains = []
          if self.parent.parent.AgeCombo.currentIndex() == 1:
              nSample = 1
              nb_grains.append(1)
              # Grain Id
              t.setItem(0, 0, QTableWidgetItem(p.DParameters['SampleName']+'1_1'))
              sampleName = p.DParameters['SampleName']+'1_1'
              #Obs Age, here always zero
              t.setItem(0, 1, QTableWidgetItem(p.DParameters['Age_grain1']))
              #Obs error, here always zero
              t.setItem(0, 2, QTableWidgetItem(p.DParameters['Error_grain1']))
              #Radius
              t.setItem(0, 3, QTableWidgetItem(p.DParameters['ASIZE']))
              # Uranium
              t.setItem(0, 4, QTableWidgetItem(p.DParameters['AUPPM']))
              # Thorium
              t.setItem(0, 5, QTableWidgetItem(p.DParameters['AThPPM']))
              # Rmr0
              t.setItem(0, 6, QTableWidgetItem(p.DParameters[conf.Variable_names['FTL_kinetic_parameter_value_AHe']]))

          elif self.parent.parent.AgeCombo.currentIndex() == 2:
              nSample = int(self.Winparent.nbSampleValue.text())
              for j in range(nSample):
                  nb_grains.append(int(self.Winparent.Samples['AHe_ngrains'+str(j+1)]))
              
              count = 0
              for j in range(nSample):
                  if nb_grains[j] > 0:
                      for i in range(nb_grains[j]):
                          if "SampleName"+str(j+1) in p.TableParameters:
                              t.setItem(count, 0, QTableWidgetItem(p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)))
                              sampleName = p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)
                          else:
                              t.setItem(count, 0, QTableWidgetItem(p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)))
                              sampleName = p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)
                          if "AHE_"+sampleName in p.TableParameters:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.TableParameters["AHE_"+sampleName]))
                          else:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.DParameters['Age_grain1']))
                          if "DAHE_"+sampleName in p.TableParameters:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.TableParameters['DAHE_'+sampleName]))
                          else:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.DParameters['Error_grain1']))
                          if "ASIZE_"+sampleName in p.TableParameters:
                              t.setItem(count, 3, QTableWidgetItem(
                              str(p.TableParameters['ASIZE_'+sampleName])))
                          else:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.DParameters['ASIZE']))
                          if "AUPPM_"+sampleName in p.TableParameters:
                              t.setItem(count, 4, QTableWidgetItem(
                                  str(p.TableParameters['AUPPM_'+sampleName])))
                          else:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.DParameters['AUPPM']))
                          if "ATHPPM_"+sampleName in p.TableParameters:
                              t.setItem(count, 5, QTableWidgetItem(
                                  str(p.TableParameters['ATHPPM_'+sampleName])))
                          else:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.DParameters['AThPPM']))
                          if conf.Variable_names['FTL_kinetic_parameter_value_AHe']+"_"+sampleName in p.TableParameters:
                              t.setItem(count, 6, QTableWidgetItem(
                                  str(p.TableParameters[conf.Variable_names['FTL_kinetic_parameter_value_AHe']+"_"+sampleName])))
                          else:
                              t.setItem(count,6, QTableWidgetItem(
                                  p.DParameters[conf.Variable_names['FTL_kinetic_parameter_value_AHe']]))
                          count += 1

        # Read the table
        num_rows, num_cols = t.rowCount(), t.columnCount()
        temp_dict = {}
        temp_GrainID = []
        try:
            for row in range(num_rows):
                for col in range(num_cols):
                    item = t.item(row, col)
                    if item is None:
                        break
                    if col == 0:
                        sampleName = item.text()
                        temp_GrainID.append(sampleName)
                    if col == 1:
                        value = float(item.text())
                        temp_dict['AHE_'+sampleName] = str(item.text())
                        p.TableParameters['AHE_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DAHE_'+sampleName] = str(item.text())
                        p.TableParameters['DAHE_'+sampleName] = str(item.text())
                    if col == 3:
                        value = float(item.text())
                        temp_dict['ASIZE_'+sampleName] = str(item.text())
                        p.TableParameters['ASIZE_'+sampleName] = str(item.text())
                    elif col == 4:
                        value = float(item.text())
                        temp_dict['AUPPM_'+sampleName] = str(item.text())
                        p.TableParameters['AUPPM_'+sampleName] = str(item.text())
                    elif col == 5:
                        value = float(item.text())
                        temp_dict['ATHPPM_'+sampleName] = str(item.text())
                        p.TableParameters['ATHPPM_'+sampleName] = str(item.text())
                    elif col == 6:
                        value = float(item.text())
                        temp_dict[conf.Variable_names['FTL_kinetic_parameter_value_AHe']+"_"+sampleName] = str(item.text())
                        p.TableParameters[conf.Variable_names['FTL_kinetic_parameter_value_AHe']+"_"+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            print("Issue in reading AHe table")
            QErrorMessage(self).showMessage("In AHe observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        temp_dict['GrainIDs'] = temp_GrainID
        self.grains.update(temp_dict)
        
        # # Update AFT rmr0
        # if self.parent.AFTParam:
        #     try:
        #         self.parent.AFTParam.GrainTable.cellChanged.disconnect()
        #         self.parent.AFTParam.updateTableGrain(self.parent.AFTParam.GrainTable, 1, p)
        #         self.parent.AFTParam.GrainTable.cellChanged.connect(lambda:self.parent.AFTParam.updateTableGrain(self.parent.AFTParam.GrainTable, 0, p))
        #     except:
        #         print('Line 3092 in Thermochronometers.py: disconnect() failed')
                
        # We write U, Th for all nodes
        if self.parent.parent.AgeCombo.currentIndex() == 1:
            for k, v in self.grains.items():
                if k.startswith('AUPPM'):
                    pgui.store_input(self.parent.parent, {'AUPPM': str(v)})
                elif k.startswith('ATHPPM'):
                    pgui.store_input(self.parent.parent, {'AThPPM': str(v)})
                elif k.startswith('ASIZE'):
                    pgui.store_input(self.parent.parent, {'ASIZE': str(v)})
                elif k.startswith(conf.Variable_names['FTL_kinetic_parameter_value_AHe']):
                    pgui.store_input(self.parent.parent, {conf.Variable_names['FTL_kinetic_parameter_value_AHe']: str(v)})

        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))
        
    #----------------------------------------------------------
    def update_PDModel(self,model):
        """
          To update the parameters for the production-diffusion model.
        
        """
          
        if model.text() == "Finite difference" and model.isChecked():
            # Extra box
            self.alphaDist = QLabel('stopping distance:')
            self.alphaDist.setFont(conf.font8)
            self.alphaDist.setAlignment(Qt.AlignRight)
            self.alphaDistCombo = QComboBox()
            self.alphaDistCombo.addItem('none')
            self.alphaDistCombo.addItem('Farley et al. (1996)')
            self.alphaDistCombo.addItem('Ketcham et al. (2011)')
            self.alphaDistCombo.setCurrentIndex(int(self.Param.DParameters['Alpha_Ejec_Flag']))
            self.alphaDistCombo.setFont(conf.font8)
            self.GridHe.replaceWidget(self.NbIter,self.alphaDist)
            self.GridHe.replaceWidget(self.NbIterEdit,self.alphaDistCombo)
            self.NbIter.deleteLater()
            self.NbIterEdit.deleteLater()
            pgui.store_input(self.Winparent, {'PDModel_signal': str(1)})
        elif model.text() == "Monte Carlo" and model.isChecked():
            # Extra box
            self.NbIter = QLabel('Number of iterations:')
            self.NbIter.setFont(conf.font8)
            self.NbIter.setToolTip("Number of iteration for the Monte Carlo simulation.")
            self.NbIter.setAlignment(Qt.AlignCenter)
            self.NbIterEdit = QLineEdit()
            self.NbIterEdit.setValidator(QIntValidator())
            self.NbIterEdit.setText(self.Param.DParameters['NMonteIter'])
            self.GridHe.replaceWidget(self.alphaDist,self.NbIter)
            self.GridHe.replaceWidget(self.alphaDistCombo,self.NbIterEdit)
            self.alphaDist.deleteLater()
            self.alphaDistCombo.deleteLater()
            pgui.store_input(self.Winparent, {'PDModel_signal': str(0)})
            
    #----------------------------------------------------------
    # def CheckGrainsChanged(self):
    #     """ To show the window with the characteristics of each grain """
    #     if self.CheckGrains.isChecked() == True:
    #         self.wgrains = pgu.WinExtraParameters(800, 1000)
    #         self.wgrains.OkButton.clicked.connect(lambda: self.wgrains.close())
    #         self.AgeParam.GrainParam = GrainsParameters(self.AgeParam, self.AgeParam.Param,self.AgeParam.parent.nbSampleValue,self.AgeParam.parent.Samples,'AHE')
    #         # We are not in inversion mode
    #         if self.Winparent.AgeCombo.currentIndex() == 1:# Means we want for all nodes
    #             self.GrainParam.TotGrain = 1 #One line to set U, Th, grainsize
    #             # Updtae Gain characteristic table, one row with default values
    #             self.GrainParam.GrainTable.setRowCount(1)
    #             self.GrainParam.GrainTable.setItem(0, 0, QTableWidgetItem('Apatite1_1'))
    #             self.GrainParam.GrainTable.setItem(0, 1, QTableWidgetItem(self.Param.DParameters['grainradius']))
    #             self.GrainParam.GrainTable.setItem(0, 2, QTableWidgetItem(self.Param.DParameters['Uppm']))
    #             self.GrainParam.GrainTable.setItem(0, 3, QTableWidgetItem(self.Param.DParameters['Thppm']))
    #             self.wgrains.layout.addLayout(self.GrainParam.Layout, 0, 0, 1, 3)
    #         else:
    #             self.wgrains.layout.addLayout(self.AgeParam.GrainParam.Layout, 0, 0, 1, 3)
    #         Q2 = QWidget()
    #         Q2.setLayout(self.wgrains.layout)
    #         self.wgrains.setCentralWidget(Q2)
    #         conf.WindowsOpen.append(self.wgrains)
    #         self.wgrains.show()
    
    #----------------------------------------------------------
    def Check43HeChanged(self):
        """ To show the window to set the parameters for the prediction of the
        4He/3He release spectrum. """
        
        if self.Check43He.isChecked() == True:
            # Store the flag
            self.d.Samples.update({conf.Variable_names['4He/3He']:'1'})
            pgui.store_input(self.Winparent, {conf.Variable_names['4He/3He']: '1'})
            if self.Winparent.AgeCombo.currentIndex() == 2: #Means compute for sample-specific
                # Open a window to provide Heating schedule
                self.w43He = interface.WinExtraParameters(1000, 800,text="4He/3He parameters")
                self.w43He.OkButton.clicked.connect(lambda: self.Close_43He())
                self.He43Param = He43.apatite(self.AgeParam, self, self.Param,self.d.nbSampleValue,self.d.Samples,self.PFolder)
                self.w43He.layout.addLayout(self.He43Param.VBox, 0, 0)
                Q = QWidget()
                Q.setLayout(self.w43He.layout)
                self.w43He.setCentralWidget(Q)
                conf.WindowsOpen.append(self.w43He)
                self.w43He.show()
        else:
            # Check if any former 43 file and remove it
            path = os.path.join(conf.PecubeFolderPath,self.PFolder,"data",self.Winparent.dataFolderEdit1.text(),conf.DefaultParameterValues().DParameters['fileName43'])
            if os.path.exists(path):
                os.remove(path)
            self.d.Samples.update({conf.Variable_names['4He/3He']:'0'}) 
            pgui.store_input(self.Winparent, {conf.Variable_names['4He/3He']: '0'})
    
    def Close_43He(self):
        try:
            self.He43Param.write_in_file()
        except Exception as E:
            QErrorMessage(self).showMessage("Problem in writing file for 4He/3He dataset: "+str(E)+" Please, check your input data.")
            return
        self.w43He.close()
        
    #----------------------------------------------------------
    def store_input_ComboBoxDiffModel(self, b, s, d):
        """ Update the value of the Combo box each it is changed. """
        flag=0
        ValueD0 = 50
        ValueEa = 137.6536
        if b.currentIndex() == 0:  # Farley
            flag = 0
            ValueD0 = 50
            ValueEa = 137.6536
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        elif b.currentIndex() == 1:  # Shuster
            flag = 1
            ValueD0 = 0.569
            ValueEa = 119.9971
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        elif b.currentIndex() == 2:  # Gautheron
            flag = 2
            ValueD0 = 0.002
            ValueEa = 109.250
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        elif b.currentIndex() == 3:  # Flowers
            flag = 3
            ValueD0 = 0.6071
            ValueEa = 122.2983
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        elif b.currentIndex() == 4:  # Willett
            flag = 4
            ValueD0 = 0.6071
            ValueEa = 122.2983
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        pgui.store_input(d, {s: str(flag)})
        pgui.store_input(d, {conf.Variable_names['D0_Apatite']: str(ValueD0)})
        pgui.store_input(d, {conf.Variable_names['Ea_Apatite']: str(ValueEa)})

    #----------------------------------------------------------
    def store_input_ComboBoxAlpha(self, b, s, d):
        """ Update the value of the Combo box each it is changed. """
        if b.currentIndex() == 0:  # no ejection
            flag = 0
        elif b.currentIndex() == 1:  # ejection
            flag = 1
        elif b.currentIndex() == 2:  # redistribution
            flag = 2
        pgui.store_input(d, {s: str(flag)})

    #----------------------------------------------------------
    def store_input_ComboBoxAlphaDist(self, b, s, d):
        """ Update the value of the Combo box each it is changed. """
        if b.currentIndex() == 0:  # no stopping distance
            flag = 0
        elif b.currentIndex() == 1:  # Farley et al. (1996)
            flag = 1
        elif b.currentIndex() == 2:  # Ketcham et al. (2011)
            flag = 2
        pgui.store_input(d, {s: str(flag)})



#############################################################################################
################################## zircon ###################################################
#############################################################################################
class zircon(QWidget):
    """
      This class set the parameters for zircon (U-Th)/He thermochronometry.
      
    """
    def __init__(self, parent, Winparent, Param, PFolder):
        super().__init__()
        
        # Parameters
        if Param.DParameters['oldInput'] == '1':
            self.oldInput = 1
        else:
            self.oldInput = 0
        self.parent = parent
        self.Winparent = Winparent
        self.grains = {}
        self.Param = Param
        self.DiffusionLabel = QLabel('Diffusion model:')
        self.DiffusionLabel.setFont(conf.font10)
        self.DiffusionModel = QComboBox()
        items = ['Reiners et al., 2004','Guenthner et al (2013)']
        self.DiffusionModel.addItems(items)
        self.DiffusionModel.setMinimumContentsLength(15)
        self.DiffusionModel.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Diffusion_model_ZHe']])-5)
        self.D0 = QLabel('D0 (cm²/s):')
        self.D0.setFont(conf.font10)
        self.D0.setAlignment(Qt.AlignCenter)
        self.D0Value = QLineEdit(self)
        self.D0Value.setText(Param.DParameters[conf.Variable_names['D0_Zircon']])
        self.ActEner = QLabel('Ea (kJ/mol):')
        self.ActEner.setFont(conf.font10)
        self.ActEner.setAlignment(Qt.AlignCenter)
        self.ActEnerValue = QLineEdit(self)
        self.ActEnerValue.setText(Param.DParameters[conf.Variable_names['Ea_Zircon']])
        self.alphaDist = QLabel('stopping distance:')
        self.alphaDist.setFont(conf.font8)
        self.alphaDist.setAlignment(Qt.AlignRight)
        self.alphaDistCombo = QComboBox()
        self.alphaDistCombo.addItem('none')
        self.alphaDistCombo.addItem('Farley et al. (1996)')
        self.alphaDistCombo.addItem('Ketcham et al. (2011)')
        self.alphaDistCombo.setCurrentIndex(int(self.Param.DParameters['Alpha_Ejec_Flag']))
        data = {'Grain ID': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Radius (µm)": ["0"]*int(Param.DParameters['ngrains']),
                "U (ppm)": ["0"]*int(Param.DParameters['ngrains']),
                "Th (ppm)": ["0"]*int(Param.DParameters['ngrains'])}
        # self.GrainParam = self.parent.ZirconParam
        # How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['ZHe_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        self.GrainTable = interface.WinTable(data, self.TotGrain, 6)
        self.updateTableGrain(self.GrainTable, 1, Param)
        # self.CheckGrains = QCheckBox('Grains characteristics')
        
        # if self.oldInput:
        #     self.GrainParam = self.parent.ZirconParam
        
        # Add in Layout
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.VBox = QVBoxLayout()
        self.GroupBox = QGroupBox("Specific kinetics")
        self.Grid.addWidget(self.DiffusionLabel, 0, 0)
        self.Grid.addWidget(self.DiffusionModel, 0, 1)
        self.Grid.addWidget(self.D0,1,0)
        self.Grid.addWidget(self.D0Value,1,1)
        self.Grid.addWidget(self.ActEner,1,2)
        self.Grid.addWidget(self.ActEnerValue,1,3)
        self.Grid.addWidget(self.alphaDist,2,0)
        self.Grid.addWidget(self.alphaDistCombo,2,1)
        # self.Grid.addWidget(self.CheckGrains,2,1)
        self.Grid.setRowStretch(5,1)
        self.Grid.setColumnStretch(4,1)
        self.VBox.addLayout(self.Grid)
        self.GroupBox.setLayout(self.VBox)
        self.VBox2 = QVBoxLayout()
        self.VBox2.addWidget(self.GroupBox)
        self.VBox2.addWidget(self.GrainTable)
        self.VBox2.setStretch(1,2)
        
        # Signals
        self.D0Value.editingFinished.connect(lambda: pgui.store_input(
            Winparent, {conf.Variable_names['D0_Zircon']: str(self.D0Value.text())}))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        self.ActEnerValue.editingFinished.connect(lambda: pgui.store_input(
            Winparent, {conf.Variable_names['Ea_Zircon']: str(self.ActEnerValue.text())}))
        self.DiffusionModel.currentIndexChanged.connect(
            lambda: self.store_input_ComboBoxDiffModel(self.DiffusionModel, conf.Variable_names['Diffusion_model_ZHe'], self.Winparent))
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.CheckGrains.stateChanged.connect(lambda: self.CheckGrainsChanged())
        
        self.Q = QWidget()
        self.Q.setLayout(self.VBox2)
            
    #----------------------------------------------------------
    def store_input_ComboBoxDiffModel(self, diff, text, Winparent):
        # Update the value of the Combo box each it is changed
        nbAHeRDmodel = 5
        if diff.currentIndex() == 0:  # Reiners et al. (2004)
            flag = 0+nbAHeRDmodel
            ValueD0 = 0.46
            ValueEa = 169.0336
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        elif diff.currentIndex() == 1:  # Guenthner et al. (2013)
            flag = 1+nbAHeRDmodel
            ValueD0 = 0.006367 # cm2/s
            ValueEa = 71 # kJ/mol
            self.D0Value.setText(str(ValueD0))
            self.ActEnerValue.setText(str(ValueEa))
        pgui.store_input(Winparent, {text: str(flag)})
        pgui.store_input(Winparent, {conf.Variable_names['D0_Zircon']: str(ValueD0)})
        pgui.store_input(Winparent, {conf.Variable_names['Ea_Zircon']: str(ValueEa)})
    
    # #----------------------------------------------------------
    # def CheckGrainsChanged(self):
    #     """ To show the window with the characteristics of each grain. """
    #     if self.CheckGrains.isChecked() == True:
    #         self.wgrains = pgu.WinExtraParameters(800, 1000)
    #         self.wgrains.OkButton.clicked.connect(lambda: self.wgrains.close())
    #         self.parent.ZirconParam = ZirconGrainsParameters(self.parent, self.parent.Param,self.parent.parent.nbSampleValue,self.parent.parent.Samples,'ZHE')
    #         self.wgrains.layout.addLayout(self.parent.ZirconParam.Layout, 0, 0, 1, 3)
    #         Q2 = QWidget()
    #         Q2.setLayout(self.wgrains.layout)
    #         self.wgrains.setCentralWidget(Q2)
    #         conf.WindowsOpen.append(self.wgrains)
    #         self.wgrains.show()
    
    #--------------------------------------------------------
    def updateTableGrain(self, t, h, p):
        """ 
          To store input value for the Grain table
          t = table
          h : parameter
          p : parameters
          d : pecube parameters (window)
        """
        # Get total number of grains
        self.grains = {}
        h_temp = 1
        nSample = 0
        if h == 2: # Update from Table Samples
            try:
                self.GrainTable.cellChanged.disconnect() # Shut down update from the user
            except TypeError:
                pass
            # if (self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 4 or
            # self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 0):
            #     h_temp = 2
            #     h = 1
            # Check if we are in "for all nodes"
            if self.parent.parent.AgeCombo.currentIndex() == 1:
                self.TotGrain = 1
                t.setRowCount(1)
                p.TableParameters["SampleName1"] = "Grain"
            elif self.parent.parent.AgeCombo.currentIndex() == 2:
                # Get total number of grains
                self.TotGrain = self.Winparent.Samples['ZHe_TotGrain']
                t.setRowCount(int(self.TotGrain))
                # Get Sample ID for each grain, read from sample table
                for j in range(int(self.Winparent.nbSampleValue.text())):
                        nb_grains = int(self.Winparent.Samples['ZHe_ngrains'+str(j+1)])
                        if nb_grains > 0:
                            sampleName = self.Winparent.Samples['SampleName'+str(j+1)]
                            p.TableParameters["SampleName"+str(j+1)] = sampleName
        
        # Update Table
        if h >= 1:
          nb_grains = []
          if self.parent.parent.AgeCombo.currentIndex() == 1:
              nSample = 1
              nb_grains.append(1)
          elif self.parent.parent.AgeCombo.currentIndex() == 2:
              nSample = int(self.Winparent.nbSampleValue.text())
              for j in range(nSample):
                  nb_grains.append(int(self.Winparent.Samples['ZHe_ngrains'+str(j+1)]))
                  
          count = 0
          for j in range(nSample):
              if nb_grains[j] > 0:
                  for i in range(nb_grains[j]):
                      if "SampleName"+str(j+1) in p.TableParameters:
                          t.setItem(count, 0, QTableWidgetItem(p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)))
                          sampleName = p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)
                      else:
                          t.setItem(count, 0, QTableWidgetItem(p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)))
                          sampleName = p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)
                      if "ZHE_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.TableParameters['ZHE_'+sampleName]))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DZHE_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.TableParameters['DZHE_'+sampleName]))
                      else:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.DParameters['Error_grain1']))
                      if "ZSIZE_"+sampleName in p.TableParameters:
                          t.setItem(count, 3, QTableWidgetItem(
                          str(p.TableParameters['ZSIZE_'+sampleName])))
                      else:
                          t.setItem(count, 3, QTableWidgetItem(
                              p.DParameters['ZSIZE']))
                      if "ZUPPM_"+sampleName in p.TableParameters:
                          t.setItem(count, 4, QTableWidgetItem(
                              str(p.TableParameters['ZUPPM_'+sampleName])))
                      else:
                          t.setItem(count, 4, QTableWidgetItem(
                              p.DParameters['ZUPPM']))
                      if "ZTHPPM_"+sampleName in p.TableParameters:
                          t.setItem(count, 5, QTableWidgetItem(
                              str(p.TableParameters['ZTHPPM_'+sampleName])))
                      else:
                          t.setItem(count, 5, QTableWidgetItem(
                              p.DParameters['ZThPPM']))
                      count += 1
        # Read the table
        num_rows, num_cols = t.rowCount(), t.columnCount()
        temp_dict = {}
        try:
            for row in range(num_rows):
                for col in range(num_cols):
                    item = t.item(row, col)
                    if item is None:
                        break
                    if col == 0:
                        sampleName = item.text()
                    if col == 1:
                        value = float(item.text())
                        temp_dict['ZHE_'+sampleName] = str(item.text())
                        p.TableParameters['ZHE_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DZHE_'+sampleName] = str(item.text())
                        p.TableParameters['DZHE_'+sampleName] = str(item.text())
                    elif col == 3:
                        value = float(item.text())
                        temp_dict['ZSIZE_'+sampleName] = str(item.text())
                        p.TableParameters['ZSIZE_'+sampleName] = str(item.text())
                    elif col == 4:
                        value = float(item.text())
                        temp_dict['ZUPPM_'+sampleName] = str(item.text())
                        p.TableParameters['ZUPPM_'+sampleName] = str(item.text())
                    elif col == 5:
                        value = float(item.text())
                        temp_dict['ZTHPPM_'+sampleName] = str(item.text())
                        p.TableParameters['ZTHPPM_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In ZHe observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        # We write U, Th for all nodes
        if self.parent.parent.AgeCombo.currentIndex() == 1:
            for k, v in self.grains.items():
                if k.startswith('ZUPPM'):
                    pgui.store_input(self.parent.parent, {'ZUPPM': str(v)})
                elif k.startswith('ZTHPPM'):
                    pgui.store_input(self.parent.parent, {'ZThPPM': str(v)})
                elif k.startswith('ZSIZE'):
                    pgui.store_input(self.parent.parent, {'ZSIZE': str(v)})
                    
                    
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))


