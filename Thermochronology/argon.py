""" 

 This module gathers all parameters, class and functions related to argon
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
################################## K-Feldspar ###############################################
#############################################################################################
class K_Feldspar(QWidget):
    """
      This class set the parameters for K-Feldspar Ar thermochronometer
      Winparent = parameters in WindowParameters class 
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
                nb_grains = self.Winparent.Samples['KAr_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains'])}
        self.GrainTable = interface.WinTable(data, self.TotGrain, 3)
        self.updateTableGrain(self.GrainTable, 1, Param)
        self.CheckGrains = QCheckBox('Grains characteristics')
        self.D0 = QLabel('D0 (cm²/s):')
        self.D0.setFont(conf.font10)
        self.D0.setAlignment(Qt.AlignCenter)
        self.D0Value = QLineEdit(self)
        self.D0Value.setText(Param.DParameters[conf.Variable_names['D0_Feldspar']])
        self.ActEner = QLabel('Ea (kJ/mol):')
        self.ActEner.setFont(conf.font10)
        self.ActEner.setAlignment(Qt.AlignCenter)
        self.ActEnerValue = QLineEdit(self)
        self.ActEnerValue.setText(Param.DParameters[conf.Variable_names['Ea_Feldspar']])
        
        # Check for edited text
        self.D0Value.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['D0_Feldspar']: str(self.D0Value.text())}))
        self.ActEnerValue.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['Ea_Feldspar']: str(self.ActEnerValue.text())}))
        
        # Add in Layout
        self.MainLayout = QVBoxLayout()
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.D0,0,0)
        self.Grid.addWidget(self.D0Value,0,1)
        self.Grid.addWidget(self.ActEner,0,2)
        self.Grid.addWidget(self.ActEnerValue,0,3)
        self.Grid.setRowStretch(5,1)
        self.Grid.setColumnStretch(4,1)
        self.MainLayout.addLayout(self.Grid)
        self.MainLayout.addWidget(self.GrainTable)
        self.Grid.setRowStretch(3,1)
        self.Grid.setColumnStretch(2,2)
        self.MainLayout.setStretch(1,2)
        
        # Signals
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        
        self.Q = QWidget()
        self.Q.setLayout(self.MainLayout)
    
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
               self.TotGrain = self.Winparent.Samples['KAr_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['KAr_ngrains'+str(j+1)])
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
                    nb_grains.append(int(self.Winparent.Samples['KAr_ngrains'+str(j+1)]))
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
                      if "KAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.TableParameters['KAR_'+sampleName]))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DKAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.TableParameters['DKAR_'+sampleName]))
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
                        temp_dict['KAR_'+sampleName] = str(item.text())
                        p.TableParameters['KAR_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DKAR_'+sampleName] = str(item.text())
                        p.TableParameters['DKAR_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In KAr observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        if h == 2: # Update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))
            
            
            
#############################################################################################
################################## Biotite ##################################################
#############################################################################################
class Biotite(QWidget):
    """
       This class set the parameters for Biotite Ar thermochronometer
       Winparent = parameters in WindowParameters class 
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
                nb_grains = self.Winparent.Samples['BAr_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains'])}
        self.GrainTable = interface.WinTable(data, self.TotGrain, 3)
        self.updateTableGrain(self.GrainTable, 1, Param)
        self.CheckGrains = QCheckBox('Grains characteristics')
        self.D0 = QLabel('D0 (cm²/s):')
        self.D0.setFont(conf.font10)
        self.D0.setAlignment(Qt.AlignCenter)
        self.D0Value = QLineEdit(self)
        self.D0Value.setText(Param.DParameters[conf.Variable_names['D0_Biotite']])
        self.ActEner = QLabel('Ea (kJ/mol):')
        self.ActEner.setFont(conf.font10)
        self.ActEner.setAlignment(Qt.AlignCenter)
        self.ActEnerValue = QLineEdit(self)
        self.ActEnerValue.setText(Param.DParameters[conf.Variable_names['Ea_Biotite']])
        
        # Check for edited text
        self.D0Value.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['D0_Biotite']: str(self.D0Value.text())}))
        self.ActEnerValue.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['Ea_Biotite']: str(self.ActEnerValue.text())}))
        
        # Add in Layout
        self.MainLayout = QVBoxLayout()
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.D0,0,0)
        self.Grid.addWidget(self.D0Value,0,1)
        self.Grid.addWidget(self.ActEner,0,2)
        self.Grid.addWidget(self.ActEnerValue,0,3)
        self.Grid.setRowStretch(5,1)
        self.Grid.setColumnStretch(4,1)
        self.MainLayout.addLayout(self.Grid)
        self.MainLayout.addWidget(self.GrainTable)
        self.MainLayout.setStretch(1,2)
        
        # Signals
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        
        self.Q = QWidget()
        self.Q.setLayout(self.MainLayout)
    
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
               self.TotGrain = self.Winparent.Samples['BAr_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['BAr_ngrains'+str(j+1)])
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
                    nb_grains.append(int(self.Winparent.Samples['BAr_ngrains'+str(j+1)]))
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
                      if "BAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.TableParameters['BAR_'+sampleName]))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DBAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.TableParameters['DBAR_'+sampleName]))
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
                        temp_dict['BAR_'+sampleName] = str(item.text())
                        p.TableParameters['BAR_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DBAR_'+sampleName] = str(item.text())
                        p.TableParameters['DBAR_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In BAr observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
      
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))
            
            
#############################################################################################
################################## Muscovite ################################################
#############################################################################################
class Muscovite(QWidget):
    """
      This class set the parameters for Muscovite Ar 
      Winparent = parameters in WindowParameters class 
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
                nb_grains = self.Winparent.Samples['MAr_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains'])}
        self.GrainTable = interface.WinTable(data, self.TotGrain, 3)
        self.updateTableGrain(self.GrainTable, 1, Param)
        self.CheckGrains = QCheckBox('Grains characteristics')
        self.D0 = QLabel('D0 (cm²/s):')
        self.D0.setFont(conf.font10)
        self.D0.setAlignment(Qt.AlignCenter)
        self.D0Value = QLineEdit(self)
        self.D0Value.setText(Param.DParameters[conf.Variable_names['D0_Muscovite']])
        self.ActEner = QLabel('Ea (kJ/mol):')
        self.ActEner.setFont(conf.font10)
        self.ActEner.setAlignment(Qt.AlignCenter)
        self.ActEnerValue = QLineEdit(self)
        self.ActEnerValue.setText(Param.DParameters[conf.Variable_names['Ea_Muscovite']])
        
        # Check for edited text
        self.D0Value.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['D0_Muscovite']: str(self.D0Value.text())}))
        self.ActEnerValue.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['Ea_Muscovite']: str(self.ActEnerValue.text())}))
        
        # Add in Layout
        self.MainLayout = QVBoxLayout()
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.D0,0,0)
        self.Grid.addWidget(self.D0Value,0,1)
        self.Grid.addWidget(self.ActEner,0,2)
        self.Grid.addWidget(self.ActEnerValue,0,3)
        self.Grid.setRowStretch(5,1)
        self.Grid.setColumnStretch(4,1)
        self.MainLayout.addLayout(self.Grid)
        self.MainLayout.addWidget(self.GrainTable)
        self.MainLayout.setStretch(1,2)
        
        # Signals
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        
        self.Q = QWidget()
        self.Q.setLayout(self.MainLayout)
    
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
               self.TotGrain = self.Winparent.Samples['MAr_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['MAr_ngrains'+str(j+1)])
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
                    nb_grains.append(int(self.Winparent.Samples['MAr_ngrains'+str(j+1)]))
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
                      if "MAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.TableParameters['MAR_'+sampleName]))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DMAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.TableParameters['DMAR_'+sampleName]))
                      else:
                          t.setItem(count, 2, QTableWidgetItem(
                              p.DParameters['Error_grain1']))
                      count += 1
        #Read the table
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
                        temp_dict['MAR_'+sampleName] = str(item.text())
                        p.TableParameters['MAR_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DMAR_'+sampleName] = str(item.text())
                        p.TableParameters['DMAR_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In MAr observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))
            
        
#############################################################################################
################################## Hornblende ###############################################
#############################################################################################
class Hornblende(QWidget):
    """
      This class set the parameters for Hornblende Ar thermochronometer
      Winparent = parameters in WindowParameters class 
    
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
                nb_grains = self.Winparent.Samples['HAr_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Obs. age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error (Ma)": ["0"]*int(Param.DParameters['ngrains'])}
        self.GrainTable = interface.WinTable(data, self.TotGrain, 3)
        self.updateTableGrain(self.GrainTable, 1, Param)
        self.CheckGrains = QCheckBox('Grains characteristics')
        self.D0 = QLabel('D0 (cm²/s):')
        self.D0.setFont(conf.font10)
        self.D0.setAlignment(Qt.AlignCenter)
        self.D0Value = QLineEdit(self)
        self.D0Value.setText(Param.DParameters[conf.Variable_names['D0_Hornblende']])
        self.ActEner = QLabel('Ea (kJ/mol):')
        self.ActEner.setFont(conf.font10)
        self.ActEner.setAlignment(Qt.AlignCenter)
        self.ActEnerValue = QLineEdit(self)
        self.ActEnerValue.setText(Param.DParameters[conf.Variable_names['Ea_Hornblende']])
        
        # Check for edited text
        self.D0Value.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['D0_Hornblende']: str(self.D0Value.text())}))
        self.ActEnerValue.editingFinished.connect(lambda: pgui.store_input(
            self.Winparent, {conf.Variable_names['Ea_Hornblende']: str(self.ActEnerValue.text())}))
        
        # Add in Layout
        self.MainLayout = QVBoxLayout()
        self.Grid = QGridLayout()
        self.Grid.setSpacing(10)
        self.Grid.addWidget(self.D0,0,0)
        self.Grid.addWidget(self.D0Value,0,1)
        self.Grid.addWidget(self.ActEner,0,2)
        self.Grid.addWidget(self.ActEnerValue,0,3)
        self.Grid.setRowStretch(5,1)
        self.Grid.setColumnStretch(4,1)
        self.MainLayout.addLayout(self.Grid)
        self.MainLayout.addWidget(self.GrainTable)
        self.MainLayout.setStretch(1,2)
        
        # Signals
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        
        self.Q = QWidget()
        self.Q.setLayout(self.MainLayout)
    
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
               self.TotGrain = self.Winparent.Samples['HAr_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['HAr_ngrains'+str(j+1)])
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
                    nb_grains.append(int(self.Winparent.Samples['HAr_ngrains'+str(j+1)]))
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
                      if "HAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 1, QTableWidgetItem(
                              str(p.TableParameters['HAR_'+sampleName])))
                      else:
                          t.setItem(count, 1, QTableWidgetItem(
                              p.DParameters['Age_grain1']))
                      if "DHAR_"+sampleName in p.TableParameters:
                          t.setItem(count, 2, QTableWidgetItem(
                              str(p.TableParameters['DHAR_'+sampleName])))
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
                        temp_dict['HAR_'+sampleName] = str(item.text())
                        p.TableParameters['HAR_'+sampleName] = str(item.text())
                    elif col == 2:
                        value = float(item.text())
                        temp_dict['DHAR_'+sampleName] = str(item.text())
                        p.TableParameters['DHAR_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In HAr observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))