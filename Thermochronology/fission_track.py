""" 

 This module gathers all parameters, class and functions related to fission tracks
thermochronometry.

@author: maxime Bernard

"""

from PyQt5.QtWidgets import (QWidget,QLabel,QCheckBox,QVBoxLayout,
                             QGridLayout,QComboBox,QGroupBox,QHBoxLayout,
                             QErrorMessage,QLineEdit,QTableWidgetItem,QRadioButton,
                             )
import configs as conf
from PyQt5.QtCore import Qt
import Utils.PGUI_utils as pgu
import Utils.interface as interface
import main as pgui


#############################################################################################
################################## apatite ##################################################
#############################################################################################
class apatite(QWidget):
    """
      This class set the parameters for apatite FT thermochronometer
      Winparent = parameters in WindowParameters class 
    
    """
    def __init__(self, parent, Winparent, Param, PFolder):
        super().__init__()
        
        #Parameters
        self.parent = parent
        self.Winparent = Winparent
        self.grains = {}
        self.Param = Param
        self.rhoSTLabel = QLabel('Density reduction in standard:')
        self.rhoSTLabel.setFont(conf.font10)
        self.rhoSTLabel.setAlignment(Qt.AlignLeft)
        self.rhoST = QLineEdit()
        self.rhoST.setText(Param.DParameters['rhoST'])
        # self.rmr0 = QLabel('rmr0:')
        # self.rmr0.setFont(conf.font10)
        # self.rmr0.setAlignment(Qt.AlignCenter)
        # self.rmr0Value = QLineEdit(self)
        # self.rmr0Value.setText(Param.DParameters['rmr0_spec'])
        self.AnnealingLabel = QLabel('Annealing model:')
        self.AnnealingLabel.setFont(conf.font10)
        self.AnnealingLabel.setAlignment(Qt.AlignRight)
        self.AnnealingModel = QComboBox()
        items = ['van der Beek (1995)','Ketcham et al. (1999)','Ketcham et al. (2007)']
        items = ['Ketcham et al. (2007)']
        self.AnnealingModel.addItems(items)
        self.AnnealingModel.setMinimumContentsLength(15)
        # self.AnnealingModel.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['FissionTrackModel']]))
        self.AnnealingModel.setCurrentIndex(0) # Temporary access to Ketcham et al. (2007) only
        self.ageFTL = QLabel('Use FTL: ')
        self.ageFTL.setFont(conf.font10)
        self.ageFTL.setAlignment(Qt.AlignRight)
        self.ageFTL.setToolTip("Use fission track length in the misfit calculation.")
        self.ageFTLb = QComboBox()
        self.ageFTLb.setMinimumContentsLength(15)
        items = ['no','Mean FTL']
        self.ageFTLb.addItems(items)
        self.ageFTLb.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Use Fission Track']]))
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "rmr0": ["0"]*int(Param.DParameters['ngrains']),
                "MFTL (µm)":["0"]*int(Param.DParameters['ngrains']),
                "Error (µm)":["0"]*int(Param.DParameters['ngrains'])}
        # How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['AFT_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        self.GrainTable = interface.WinTable(data, self.TotGrain,6)
        self.updateTableGrain(self.GrainTable, 1, Param)
        self.Kinetic_parameter_label = QLabel('Kinetic parameter: ')
        self.Kinetic_parameter_label.setToolTip('Kinetic parameter for fission track annealing')
        self.Kinetic_parameter_label.setFont(conf.font10)
        self.Kinetic_parameter_label.setAlignment(Qt.AlignRight)
        self.Kinetic_parameter_list = QComboBox()
        items = ['Dpar (µm)',' Cl (apfu)', 'OH (apfu)','Cl (wt%)','rmr0']
        self.Kinetic_parameter_list.addItems(items)
        self.Kinetic_parameter_list.setMinimumContentsLength(15)
        self.Kinetic_parameter_list.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['FTL_kinetic_parameter']]))
        self.Initial_track_length = QLabel('Initial track length: ')
        self.Initial_track_length.setToolTip('Choose from which to calculate the initial non-annealed track length.')
        self.Initial_track_length.setFont(conf.font10)
        self.Initial_track_length.setAlignment(Qt.AlignRight)
        self.Initial_track_length_list = QComboBox()
        items = ['Single value (µm)','Dpar (µm)']
        self.Initial_track_length_list.addItems(items)
        self.Initial_track_length_list.setMinimumContentsLength(15)
        self.Initial_track_length_list.setCurrentIndex(int(self.Param.DParameters[conf.Variable_names['Initial_FTL_model']]))
        self.Initial_track_length_value = QLineEdit()
        self.Initial_track_length_value.setText(Param.DParameters[conf.Variable_names['Unannealed_FTL_value']])
        self.MFTL_error_label = QLabel('MFTL error (µm):')
        self.MFTL_error_label.setFont(conf.font10)
        self.MFTL_error_label.setAlignment(Qt.AlignRight)
        self.MFTL_error_label.setToolTip('Set an uncertainty on the MFTL. This will be use in the misfit calculation (for inversion).')
        self.MFTL_error = QLineEdit()
        self.MFTL_error.setText(Param.DParameters[conf.Variable_names['Mean_fission_track_length_error_inversion']])
        self.MFTL_error.setEnabled(False)
        self.MFTL_std_error_label = QLabel('MFTL std error (µm):')
        self.MFTL_std_error_label.setFont(conf.font10)
        self.MFTL_std_error_label.setAlignment(Qt.AlignRight)
        self.MFTL_std_error_label.setToolTip('Set an uncertainty on the std of MFTL. This will be use in the misfit calculation (for inversion).')
        self.MFTL_std_error = QLineEdit()
        self.MFTL_std_error.setText(Param.DParameters[conf.Variable_names['Mean_fission_track_length_std_error_inversion']])
        self.MFTL_std_error.setEnabled(False)
        if self.ageFTLb.currentIndex() == 1:
            self.MFTL_error.setEnabled(True)
            self.MFTL_std_error.setEnabled(True)
        
        # Add in Layout
        self.GroupBox = QGroupBox("Specific kinetics")
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.AnnealingLabel, 0, 0)
        self.Grid.addWidget(self.AnnealingModel, 0, 1)
        self.Grid.addWidget(self.ageFTL, 0,3)
        self.Grid.addWidget(self.ageFTLb, 0,4)
        self.Grid.addWidget(self.Initial_track_length, 1, 0)
        self.Grid.addWidget(self.Initial_track_length_list, 1, 1)
        self.Grid.addWidget(self.Initial_track_length_value, 1, 2)
        self.Grid.addWidget(self.MFTL_error_label, 1, 3)
        self.Grid.addWidget(self.MFTL_error, 1, 4)
        self.Grid.addWidget(self.Kinetic_parameter_label, 2, 0)
        self.Grid.addWidget(self.Kinetic_parameter_list, 2, 1)
        self.Grid.addWidget(self.MFTL_std_error_label, 2, 3)
        self.Grid.addWidget(self.MFTL_std_error, 2, 4)
        self.Grid.addWidget(self.rhoSTLabel, 3, 0)
        self.Grid.addWidget(self.rhoST, 3, 1)
        self.Grid.setColumnStretch(5,2)
        self.GroupBox.setLayout(self.Grid)
        self.vBox = QVBoxLayout()
        self.vBox.addWidget(self.GroupBox)
        self.vBox.addWidget(self.GrainTable)
        
        # Signals
        self.AnnealingModel.currentIndexChanged.connect(
            lambda: self.store_input_ComboBoxFT())
        self.rhoST.editingFinished.connect(lambda: pgui.store_input(self.Winparent,{'rhoST': str(self.rhoST.text())}))
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        self.ageFTLb.currentIndexChanged.connect(lambda: self.update_FTLcomputation())
        self.Initial_track_length_value.textChanged.connect(lambda: self.update_FTLInitValue())
        self.Initial_track_length_list.currentIndexChanged.connect(lambda: self.update_FTLInit())
        self.Kinetic_parameter_list.currentIndexChanged.connect(lambda: self.update_FTLkinetic())
        self.MFTL_error.editingFinished.connect(lambda:pgui.store_input(self.Winparent,{conf.Variable_names['Mean_fission_track_length_error_inversion']: str(self.MFTL_error.text())}))
        self.MFTL_std_error.editingFinished.connect(lambda:pgui.store_input(self.Winparent,{conf.Variable_names['Mean_fission_track_length_std_error_inversion']: str(self.MFTL_std_error.text())}))
        
        self.Q = QWidget()
        self.Q.setLayout(self.vBox)
    
    
    #--------------------------------------------------------
    def update_FTLInitValue(self):
        """ Update the Initial track length."""
        pgui.store_input(self.Winparent, {conf.Variable_names['Unannealed_FTL_value']: str(self.Initial_track_length_value.text())})
                                          
    #--------------------------------------------------------
    def update_FTLInit(self):
        """ Update the Initial track length calculation."""
        pgui.store_input(self.Winparent, {conf.Variable_names['Initial_FTL_model']: str(self.Initial_track_length_list.currentIndex())})
        if self.Initial_track_length_list.currentIndex() == 1:
            self.Initial_track_length_value.setEnabled(False)
        else:
            self.Initial_track_length_value.setEnabled(True)
                                          
                                          
    #--------------------------------------------------------
    def update_FTLkinetic(self):
        """ Update the kinetic parameter."""
        pgui.store_input(self.Winparent, {conf.Variable_names['FTL_kinetic_parameter']: str(self.Kinetic_parameter_list.currentIndex())})
        # Update initial track list
        if self.Kinetic_parameter_list.currentIndex() == 0: #Dpar
            self.Initial_track_length_list.setItemText(1,'Dpar (µm)')
            Headers = ['Grain','Obs. age(Ma)','Error (Ma)','Dpar (µm)']
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.Kinetic_parameter_list.currentIndex() == 1: #Cl (apfu)
            self.Initial_track_length_list.setItemText(1,'Cl (apfu)')
            Headers = ['Grain','Obs. age(Ma)','Error (Ma)','Cl (apfu)']
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.Kinetic_parameter_list.currentIndex() == 2: #OH (apfu)
            self.Initial_track_length_list.setItemText(1,'OH (apfu)')
            Headers = ['Grain','Obs. age(Ma)','Error (Ma)','OH (apfu)']
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.Kinetic_parameter_list.currentIndex() == 3: #Cl (wt %)
            self.Initial_track_length_list.setItemText(1,'Cl (wt %)')
            Headers = ['Grain','Obs. age(Ma)','Error (Ma)','Cl (wt %)']
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.Kinetic_parameter_list.currentIndex() == 4: #rmr0
            self.Initial_track_length_list.setItemText(1,'none')
            Headers = ['Grain','Obs. age(Ma)','Error (Ma)','rmr0']
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        
    #--------------------------------------------------------
    def update_FTLcomputation(self):
        """ Flag for computation of Fission track length distribution. """
        pgui.store_input(self.Winparent, {conf.Variable_names["Use Fission Track"]: str(self.ageFTLb.currentIndex())})
        if self.ageFTLb.currentIndex() == 1:
            self.MFTL_error.setEnabled(True)
            self.MFTL_std_error.setEnabled(True)
        else:
            self.MFTL_error.setEnabled(False)
            self.MFTL_std_error.setEnabled(False)
        
    #--------------------------------------------------------
    def store_input_ComboBoxFT(self):
        """ Store input for fission track annealing. """
        if self.AnnealingModel.currentIndex() == 0:
            flag = 0
        elif self.AnnealingModel.currentIndex() == 1:
            flag = 1
        elif self.AnnealingModel.currentIndex() == 2:
            flag = 2
        pgui.store_input(self.Winparent, {conf.Variable_names['FissionTrackModel']: str(flag+2)}) #☻+2 temporary limit to Ketcham et al. 2007!
    
    #--------------------------------------------------------
    def updateTableGrain(self, t, h, p):
        """ 
          To store input value for the Grain table
          t = table
          h : parameter
          p : parameters
          d : pecube parameters (window)
        """
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
               self.TotGrain = self.Winparent.Samples['AFT_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['AFT_ngrains'+str(j+1)])
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
                    nb_grains.append(int(self.Winparent.Samples['AFT_ngrains'+str(j+1)]))

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
                      if "AFT_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.TableParameters['AFT_'+sampleName]))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DAFT_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.TableParameters['DAFT_'+sampleName]))
                      else:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.DParameters['Error_grain1']))
                      if conf.Variable_names['FTL_kinetic_parameter_value_AFT']+"_"+sampleName in p.TableParameters:
                          t.setItem(count, 3, QTableWidgetItem(
                              p.TableParameters[conf.Variable_names['FTL_kinetic_parameter_value_AFT']+"_"+sampleName]))
                      else:
                          t.setItem(count, 3, QTableWidgetItem(
                              p.DParameters[conf.Variable_names['FTL_kinetic_parameter_value_AFT']]))
                      # MFTL
                      if conf.Variable_names['Mean Fission Track Length']+"_"+sampleName in p.TableParameters:
                          t.setItem(count,4, QTableWidgetItem(
                              p.TableParameters[conf.Variable_names['Mean Fission Track Length']+"_"+sampleName]))
                      else:
                          t.setItem(count, 4, QTableWidgetItem(
                              p.DParameters[conf.Variable_names['Mean Fission Track Length']]))
                      # MFTL error
                      if conf.Variable_names['Mean Fission Track Length error']+"_"+sampleName in p.TableParameters:
                          t.setItem(count,5, QTableWidgetItem(
                              p.TableParameters[conf.Variable_names['Mean Fission Track Length error']+"_"+sampleName]))
                      else:
                          t.setItem(count, 5, QTableWidgetItem(
                              p.DParameters[conf.Variable_names['Mean Fission Track Length error']]))
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
                        temp_dict['AFT_'+sampleName] = str(item.text())
                        p.TableParameters['AFT_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DAFT_'+sampleName] = str(item.text())
                        p.TableParameters['DAFT_'+sampleName] = str(item.text())
                    elif col == 3:
                        value = float(item.text())
                        temp_dict[conf.Variable_names['FTL_kinetic_parameter_value_AFT']+"_"+sampleName] = str(item.text())
                        p.TableParameters[conf.Variable_names['FTL_kinetic_parameter_value_AFT']+"_"+sampleName] = str(item.text())
                    # MFTL
                    elif col == 4:
                        value = float(item.text())
                        temp_dict[conf.Variable_names['Mean Fission Track Length']+"_"+sampleName] = str(item.text())
                        p.TableParameters[conf.Variable_names['Mean Fission Track Length']+"_"+sampleName] = str(item.text())
                    # MFTL error
                    elif col == 5:
                        value = float(item.text())
                        temp_dict[conf.Variable_names['Mean Fission Track Length error']+"_"+sampleName] = str(item.text())
                        p.TableParameters[conf.Variable_names['Mean Fission Track Length error']+"_"+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In AFT observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)

         # We write U, Th for all nodes
        if self.parent.parent.AgeCombo.currentIndex() == 1:
            for k, v in self.grains.items():
                if k.startswith(conf.Variable_names['FTL_kinetic_parameter_value_AFT']):
                    pgui.store_input(self.parent.parent, {conf.Variable_names['FTL_kinetic_parameter_value_AFT']: str(v)})
       
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))



#############################################################################################
################################## zircon ###################################################
#############################################################################################
class zircon(QWidget):
    """
      This class set the parameters for zircon FT thermochronometer.
      Winparent = parameters in WindowParameters class.
    
    """
    def __init__(self, parent, Winparent, Param, PFolder):
        super().__init__()
        
        # Parameters
        self.parent = parent
        self.Winparent = Winparent
        self.grains = {}
        self.Param = Param
        # How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['ZFT_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains'])}
        self.GrainTable = interface.WinTable(data, self.TotGrain, 3)
        self.updateTableGrain(self.GrainTable, 1, Param)
        self.CheckGrains = QCheckBox('Grains characteristics')
        
        # Add in Layout
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        # if self.parent.AgeCombo.currentIndex() == 2:#Means we want grain-specific
        #     self.Grid.addWidget(self.AnnealingLabel, 0, 0)
        #     self.Grid.addWidget(self.AnnealingModel, 0, 1)
        #     self.Grid.addWidget(self.rmr0, 1, 0)
        #     self.Grid.addWidget(self.rmr0Value, 1, 1)
        #     self.Grid.addWidget(self.CheckGrains,2,1)
        self.Grid.addWidget(self.GrainTable, 3, 0, 2, 3)
        self.Grid.setRowStretch(3,1)
        self.Grid.setColumnStretch(2,2)
        
        # Signals
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        
        self.Q = QWidget()
        self.Q.setLayout(self.Grid)
    
    #--------------------------------------------------------
    def updateTableGrain(self, t, h, p):
        """ 
          To store input value for the Grain table
          t = table
          h : parameter
          p : parameters
          d : pecube parameters (window).
        """
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
               self.TotGrain = self.Winparent.Samples['ZFT_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['ZFT_ngrains'+str(j+1)])
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
                    nb_grains.append(int(self.Winparent.Samples['ZFT_ngrains'+str(j+1)]))
          # dParam = DefaultParameterValues()
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
                      if "ZFT_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.TableParameters['ZFT_'+sampleName]))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DZFT_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.TableParameters['DZFT_'+sampleName]))
                      else:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.DParameters['Error_grain1']))
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
                        temp_dict['ZFT_'+sampleName] = str(item.text())
                        p.TableParameters['ZFT_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DZFT_'+sampleName] = str(item.text())
                        p.TableParameters['DZFT_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In ZFT observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))