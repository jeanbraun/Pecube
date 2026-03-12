""" 

 This module gathers all parameters, class and functions related to trapped-charge
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
import main as pgui
import Utils.interface as interface

#############################################################################################
################################## Thermoluminescence #######################################
#############################################################################################
class ThermoLuminescence(QWidget):
    """
      This class set the parameters for Thermo-Luminescence
      Winparent = parameters in WindowParameters class.
    
    """
    def __init__(self, parent, Winparent, Param, PFolder):
        super().__init__()
        
        # Parameters
        self.parent = parent
        self.Winparent = Winparent
        self.grains = {}
        self.Param = Param
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "N/N": ["0"]*int(Param.DParameters['ngrains']),
                "Error": ["0"]*int(Param.DParameters['ngrains']),
                "Dose rate (Gy/kyr)": ["5.0"]*int(Param.DParameters['ngrains']),
                "D0 (Gy)": ["800"]*int(Param.DParameters['ngrains']),
                "a": ["1.8"]*int(Param.DParameters['ngrains']),
                "b": ["1.8"]*int(Param.DParameters['ngrains']),
                "Et (eV)": ["1.4"]*int(Param.DParameters['ngrains']),
                "log_s (1/s)": ["0"]*int(Param.DParameters['ngrains']),
                "log_r": ["-5.5"]*int(Param.DParameters['ngrains'])}
        # How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['TL_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        self.GrainTable = interface.WinTable(data, self.TotGrain, 10)
        self.updateTableGrain(self.GrainTable, 1, Param)
        
        self.TLDoser = QLabel('Dose rate (Gy/kyr): ')
        self.TLDoser.setFont(conf.font8)
        self.TLDoser.setAlignment(Qt.AlignCenter)
        self.TLDoser.setToolTip(
            "Dose rate (Doser) used for TL data computation.")
        self.TLD0 = QLabel('D0 (Gy): ')
        self.TLD0.setAlignment(Qt.AlignCenter)
        self.TLD0.setFont(conf.font8)
        self.TLD0.setToolTip(
            "Onset of dose saturation (D0) for TL data computation.")
        self.TLa = QLabel('a: ')
        self.TLa.setFont(conf.font8)
        self.TLa.setAlignment(Qt.AlignCenter)
        self.TLa.setToolTip(
            "Kinetic orders of trapping for TL data computation.")
        self.TL_b = QLabel('b: ')
        self.TL_b.setFont(conf.font8)
        self.TL_b.setAlignment(Qt.AlignCenter)
        self.TL_b.setToolTip(
            "Kinetic orders of detrapping for TL data computation.")
        self.TL_Et = QLabel('Activation energy (eV): ')
        self.TL_Et.setFont(conf.font8)
        self.TL_Et.setAlignment(Qt.AlignCenter)
        self.TL_Et.setToolTip(
            "Activation energy (Et) used for TL data computation.")
        self.TL_logs = QLabel('log_s (1/s):')
        self.TL_logs.setFont(conf.font8)
        self.TL_logs.setAlignment(Qt.AlignCenter)
        self.TL_logs.setToolTip(
            "Log of thermal frequancy factor for TL data computation")
        self.TL_logrho = QLabel('log(\u03C1):')
        self.TL_logrho.setFont(conf.font8)
        self.TL_logrho.setAlignment(Qt.AlignCenter)
        self.TL_logrho.setToolTip(
            "log of dimensionless recombination center density for TL data computation")

        # Create text boxes
        self.TLDoserEdit1 = QLineEdit(self)
        self.TLDoserEdit1.setText(Param.DParameters['TL_doser'])
        self.TLDoserEdit1.setAlignment(Qt.AlignCenter)
        self.TLD0Edit2 = QLineEdit(self)
        self.TLD0Edit2.setText(Param.DParameters['TL_D0'])
        self.TLD0Edit2.setAlignment(Qt.AlignCenter)
        self.TLaEdit3 = QLineEdit(self)
        self.TLaEdit3.setText(Param.DParameters['TL_a'])
        self.TLaEdit3.setAlignment(Qt.AlignCenter)
        self.TL_bEdit4 = QLineEdit(self)
        self.TL_bEdit4.setText(Param.DParameters['TL_b'])
        self.TL_bEdit4.setAlignment(Qt.AlignCenter)
        self.TL_EtEdit5 = QLineEdit(self)
        self.TL_EtEdit5.setText(Param.DParameters['TL_Et'])
        self.TL_EtEdit5.setAlignment(Qt.AlignCenter)
        self.TL_logsEdit6 = QLineEdit(self)
        self.TL_logsEdit6.setText(Param.DParameters['TL_logs'])
        self.TL_logsEdit6.setAlignment(Qt.AlignCenter)
        self.TL_logrhoEdit7 = QLineEdit(self)
        self.TL_logrhoEdit7.setText(Param.DParameters['TL_logrho'])
        self.TL_logrhoEdit7.setAlignment(Qt.AlignCenter)
        self.TL_Model_Label = QLabel("Kinetic model:")
        self.TL_Model_Label.setAlignment(Qt.AlignLeft)
        self.TL_Model = QComboBox()
        items = ['GOK_FAD']
        self.TL_Model.addItems(items)
        

        # Extra box
        fakeWidget = QLabel()
        self.GroupBox1 = QGroupBox()
        self.VBox = QVBoxLayout()
        self.GridTL = QGridLayout()
        self.GridTL.setColumnStretch(3, 2)
        # self.GridTL.setSpacing(5)
        self.GridTL.addWidget(fakeWidget, 0, 0)
        self.GridTL.addWidget(self.TL_Model_Label,1,1)
        self.GridTL.addWidget(self.TL_Model,1,2)
        self.GridTL.addWidget(fakeWidget, 3, 3)
        self.GroupBox1.setLayout(self.GridTL)
        self.VBox.addWidget(self.GroupBox1)
        self.VBox.addWidget(self.GrainTable)
        
        # Signals
        self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
            self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        self.TL_Model.currentIndexChanged.connect(lambda: self.updateTLModel())

        
        self.Q = QWidget()
        self.Q.setLayout(self.VBox)
    
    #--------------------------------------------------------
    def updateTLModel(self):
        """ Update kinetic model for TL."""
        
        pgui.store_input(self.parent.parent,{conf.Variable_names['ThL Model Label']: str(self.TL_Model.currentIndex())})
        
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
        if h == 2: #update from Table Samples
            try:
                self.GrainTable.cellChanged.disconnect() #Shut down update from the user
            except TypeError:
                pass
            # if (self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 4 or
            # self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 0):
            #     h_temp = 2
            #     h = 1
                
            #Check if we are in "for all nodes"
            if self.parent.parent.AgeCombo.currentIndex() == 1:
               self.TotGrain = 1
               t.setRowCount(1)
               p.TableParameters["SampleName1"] = "Grain"
            elif self.parent.parent.AgeCombo.currentIndex() == 2:
               # Get total number of grains
               self.TotGrain = self.Winparent.Samples['TL_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['TL_ngrains'+str(j+1)])
                       if nb_grains > 0:
                           sampleName = self.Winparent.Samples['SampleName'+str(j+1)]
                           p.TableParameters["SampleName"+str(j+1)] = sampleName
        
        # Update Table
        if h >= 1:
          nb_grains = []
          if self.parent.parent.AgeCombo.currentIndex() == 1:
              nSample = 1
              nb_grains.append(1)
              # Grain Id
              t.setItem(0, 0, QTableWidgetItem(p.DParameters['SampleName']+'1_1'))
              sampleName = p.DParameters['SampleName']+'1_1'
              #Obs Age, here always zero
              t.setItem(0, 1, QTableWidgetItem(p.DParameters['TL_NN']))
              #Obs error, here always zero
              t.setItem(0, 2, QTableWidgetItem(p.DParameters['TL_DNN']))
              #Radius
              t.setItem(0, 3, QTableWidgetItem(p.DParameters['TL_doser']))
              t.setItem(0, 4, QTableWidgetItem(p.DParameters['TL_D0']))
              # Uranium
              t.setItem(0, 5, QTableWidgetItem(p.DParameters['TL_a']))
              # Thorium
              t.setItem(0, 6, QTableWidgetItem(p.DParameters['TL_b']))
              # Rmr0
              t.setItem(0, 7, QTableWidgetItem(p.DParameters['TL_Et']))
              t.setItem(0, 8, QTableWidgetItem(p.DParameters['TL_logs']))
              t.setItem(0, 9, QTableWidgetItem(p.DParameters['TL_logrho']))

          elif self.parent.parent.AgeCombo.currentIndex() == 2:
              nSample = int(self.Winparent.nbSampleValue.text())
              for j in range(nSample):
                  nb_grains.append(int(self.Winparent.Samples['TL_ngrains'+str(j+1)]))
              
              count = 0
              for j in range(nSample):
                  if nb_grains[j] > 0:
                      for i in range(nb_grains[j]):
                          if "SampleName"+str(j+1) in p.TableParameters:
                              t.setItem(count, 0, QTableWidgetItem(p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)))
                              sampleName = p.TableParameters['SampleName'+str(j+1)]
                          else:
                              t.setItem(count, 0, QTableWidgetItem(p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)))
                              sampleName = p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)
                          if sampleName+"_"+str(i+1)+"_N/N_TL" in p.TableParameters:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_N/N_TL"]))
                          else:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.DParameters['TL_NN']))
                          if sampleName+"_"+str(i+1)+"_DNN_TL" in p.TableParameters:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_DNN_TL"]))
                          else:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.DParameters['TL_DNN']))
                          if sampleName+"_"+str(i+1)+"_DOSER_TL" in p.TableParameters:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_DOSER_TL"]))
                          else:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.DParameters['TL_doser']))
                          if sampleName+"_"+str(i+1)+"_D0_TL" in p.TableParameters:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_D0_TL"]))
                          else:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.DParameters['TL_D0']))
                          if sampleName+"_"+str(i+1)+"_ATL_TL" in p.TableParameters:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_ATL_TL"]))
                          else:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.DParameters['TL_a']))
                          if sampleName+"_"+str(i+1)+"_BTL_TL" in p.TableParameters:
                              t.setItem(count, 6, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_BTL_TL"]))
                          else:
                              t.setItem(count, 6, QTableWidgetItem(
                                  p.DParameters['TL_b']))
                          if sampleName+"_"+str(i+1)+"_ET_TL" in p.TableParameters:
                              t.setItem(count, 7, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_ET_TL"]))
                          else:
                              t.setItem(count, 7, QTableWidgetItem(
                                  p.DParameters['TL_Et']))
                          if sampleName+"_"+str(i+1)+"_LOGS_TL" in p.TableParameters:
                              t.setItem(count, 8, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_LOGS_TL"]))
                          else:
                              t.setItem(count, 8, QTableWidgetItem(
                                  p.DParameters['TL_logs']))
                          if sampleName+"_"+str(i+1)+"_LOGRHO_TL" in p.TableParameters:
                              t.setItem(count, 9, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_LOGRHO_TL"]))
                          else:
                              t.setItem(count, 9, QTableWidgetItem(
                                  p.DParameters['TL_logrho']))
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
                        float(item.text())
                        temp_dict['TL_NN_'+sampleName] = str(item.text())
                        p.TableParameters['TL_NN_'+sampleName] = str(item.text())
                    elif col == 2:
                        float(item.text())
                        temp_dict['DTL_NN_'+sampleName] = str(item.text())
                        p.TableParameters['TL_DNN_'+sampleName] = str(item.text())
                    elif col == 3:
                        float(item.text())
                        temp_dict['TL_Doser_'+sampleName] = str(item.text())
                        p.TableParameters['TL_Doser_'+sampleName] = str(item.text())
                    elif col == 4:
                        float(item.text())
                        temp_dict['TL_D0_'+sampleName] = str(item.text())
                        p.TableParameters['TL_D0_'+sampleName] = str(item.text())
                    elif col == 5:
                        float(item.text())
                        temp_dict['TL_a_'+sampleName] = str(item.text())
                        p.TableParameters['TL_a_'+sampleName] = str(item.text())
                    elif col == 6:
                        float(item.text())
                        temp_dict['TL_b_'+sampleName] = str(item.text())
                        p.TableParameters['TL_b_'+sampleName] = str(item.text())
                    elif col == 7:
                        float(item.text())
                        temp_dict['TL_Et_'+sampleName] = str(item.text())
                        p.TableParameters['TL_Et_'+sampleName] = str(item.text())
                    elif col == 8:
                        float(item.text())
                        temp_dict['TL_logs_'+sampleName] = str(item.text())
                        p.TableParameters['TL_logs_'+sampleName] = str(item.text())
                    elif col == 9:
                        float(item.text())
                        temp_dict['TL_logrho_'+sampleName] = str(item.text())
                        p.TableParameters['TL_logrho_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In Th.L observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        # We write for all nodes
        if self.parent.parent.AgeCombo.currentIndex() == 1:
            for k, v in self.grains.items():
                if k.startswith('TL_Doser'):
                    pgui.store_input(self.parent.parent, {'TL_doser': str(v)})
                elif k.startswith('TL_D0'):
                    pgui.store_input(self.parent.parent, {'TL_D0': str(v)})
                elif k.startswith('TL_Et'):
                    pgui.store_input(self.parent.parent, {'TL_Et': str(v)})
                elif k.startswith('TL_a'):
                    pgui.store_input(self.parent.parent, {'TL_a': str(v)})
                elif k.startswith('TL_b'):
                    pgui.store_input(self.parent.parent, {'TL_b': str(v)})
                elif k.startswith('TL_logs'):
                    pgui.store_input(self.parent.parent, {'TL_logs': str(v)})
                elif k.startswith('TL_lorho'):
                    pgui.store_input(self.parent.parent, {'TL_logrho': str(v)})
       
        if h == 2: # Update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))
            

#############################################################################################
########################################## OSL ##############################################
#############################################################################################
class OSL(QWidget):
    """
      This class set the extra parameters for OSL/TL sample specific age computation
      Winparent = parameters in WindowParameters class
    
    """
    def __init__(self, parent, Winparent, Param, PFolder):
        super().__init__()
        
        # Parameters
        self.parent = parent
        self.Winparent = Winparent
        self.grains = {}
        self.Param = Param
        data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                "Age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                "Error": ["0"]*int(Param.DParameters['ngrains']),
                "Dose rate (Gy/kyr)": [Param.DParameters['OSL_doser']]*int(Param.DParameters['ngrains']),
                "D0 (Gy)": [Param.DParameters['OSL_D0']]*int(Param.DParameters['ngrains']),
                "Eu (eV)": [Param.DParameters['OSL_Eu']]*int(Param.DParameters['ngrains']),
                "Et (eV)": [Param.DParameters['OSL_Et']]*int(Param.DParameters['ngrains']),
                "log_s (1/s)": [Param.DParameters['OSL_logs']]*int(Param.DParameters['ngrains']),
                "log_r": [Param.DParameters['OSL_logrho']]*int(Param.DParameters['ngrains']),
                "a":  [Param.DParameters['OSL_a']]*int(Param.DParameters['ngrains']),
                "Lmax":  [Param.DParameters['OSL_Lmax']]*int(Param.DParameters['ngrains'])}
        # How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['OSL_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        self.GrainTable = interface.WinTable(data, self.TotGrain, 11)
        self.updateTableGrain(self.GrainTable, 1, Param)
        
        self.OSLDoser = QLabel('Dose rate (Gy/kyr): ')
        self.OSLDoser.setFont(conf.font8)
        self.OSLDoser.setAlignment(Qt.AlignCenter)
        self.OSLDoser.setToolTip(
            "Dose rate (Doser) used for OSL data computation.")
        self.OSLD0 = QLabel('D0 (Gy): ')
        self.OSLD0.setAlignment(Qt.AlignCenter)
        self.OSLD0.setFont(conf.font8)
        self.OSLD0.setToolTip(
            "Onset of dose saturation (D0) for OSL data computation.")
        self.OSL_Et = QLabel('Activation energy Et (eV): ')
        self.OSL_Et.setFont(conf.font8)
        self.OSL_Et.setAlignment(Qt.AlignCenter)
        self.OSL_Et.setToolTip(
            "Activation energy (Et) used for OSL data computation.")
        self.OSL_Eu = QLabel('Activation energy Eu (eV): ')
        self.OSL_Eu.setFont(conf.font8)
        self.OSL_Eu.setAlignment(Qt.AlignCenter)
        self.OSL_Eu.setToolTip(
            "Activation energy (Eu) used for OSL data computation.")
        self.OSL_logs = QLabel('log_s (1/s):')
        self.OSL_logs.setFont(conf.font8)
        self.OSL_logs.setAlignment(Qt.AlignCenter)
        self.OSL_logs.setToolTip(
            "Log of thermal frequancy factor for OSL data computation")
        self.OSL_logrho = QLabel('log(\u03C1):')
        self.OSL_logrho.setFont(conf.font8)
        self.OSL_logrho.setAlignment(Qt.AlignCenter)
        self.OSL_logrho.setToolTip(
            "log of dimensionless recombination center density for OSL data computation")
        self.OSLmisfitTarget_Label = QLabel("Calculate misfit on: ")
        self.OSLmisfitTarget_Label.setFont(conf.font8)
        self.OSLmisfitTarget_Label.setAlignment(Qt.AlignCenter)
        self.OSLmisfitTarget_Label.setToolTip(
            "Choose if you want to calculate misfit on OSL/IRSL ages or n/N values.")

        # Create text boxes
        self.OSLDoserValue = QLineEdit(self)
        self.OSLDoserValue.setText(Param.DParameters['OSL_doser'])
        self.OSLDoserValue.setAlignment(Qt.AlignCenter)
        self.OSLD0Value = QLineEdit(self)
        self.OSLD0Value.setText(Param.DParameters['OSL_D0'])
        self.OSLD0Value.setAlignment(Qt.AlignCenter)
        self.OSL_EtValue = QLineEdit(self)
        self.OSL_EtValue.setText(Param.DParameters['OSL_Et'])
        self.OSL_EtValue.setAlignment(Qt.AlignCenter)
        self.OSL_EuValue = QLineEdit(self)
        self.OSL_EuValue.setText(Param.DParameters['OSL_Eu'])
        self.OSL_EuValue.setAlignment(Qt.AlignCenter)
        self.OSL_logsValue = QLineEdit(self)
        self.OSL_logsValue.setText(Param.DParameters['OSL_logs'])
        self.OSL_logsValue.setAlignment(Qt.AlignCenter)
        self.OSL_logrhoValue = QLineEdit(self)
        self.OSL_logrhoValue.setText(Param.DParameters['OSL_logrho'])
        self.OSL_logrhoValue.setAlignment(Qt.AlignCenter)
        self.OSL_Model_Label = QLabel("Kinetic model:")
        self.OSL_Model_Label.setAlignment(Qt.AlignLeft)
        self.OSL_Model = QComboBox()
        items = ['BandTail_FAD','GOK_GAUSS_FAD','GOK_FAD']
        self.OSL_Model.addItems(items)
        self.OSL_Model.setCurrentIndex(int(Param.DParameters[conf.Variable_names['OSL Model Label']]))
        self.OSLMisfitCombo = QComboBox()
        items = ['Age','n/N']
        self.OSLMisfitCombo.addItems(items)
        self.OSLMisfitCombo.setCurrentIndex(int(Param.DParameters[conf.Variable_names['OSL_misfit_target']]))

        # Extra box
        fakeWidget = QLabel()
        self.GroupBox1 = QGroupBox()
        self.VBox = QVBoxLayout()
        self.Grid = QGridLayout()
        self.Grid.setColumnStretch(3,2)
        self.Grid.setSpacing(5)
        self.Grid.addWidget(fakeWidget, 0, 0)
        self.Grid.addWidget(self.OSL_Model_Label,1,1)
        self.Grid.addWidget(self.OSL_Model,1,2)
        self.Grid.addWidget(self.OSLmisfitTarget_Label,2,1)
        self.Grid.addWidget(self.OSLMisfitCombo,2,2)
        self.Grid.addWidget(fakeWidget, 5, 3)
        self.GroupBox1.setLayout(self.Grid)
        self.VBox.addWidget(self.GroupBox1)
        self.VBox.addWidget(self.GrainTable)
        
        # Signals
        # self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 0, Param))
        # self.Winparent.SampleTable.cellChanged.connect(lambda: self.updateTableGrain(
        #     self.GrainTable, 2, Param))
        self.OSL_Model.currentIndexChanged.connect(lambda: self.updateOSLModel())
        self.OSLMisfitCombo.currentIndexChanged.connect(lambda: self.updateOSLmisfit())

        
        self.Q = QWidget()
        self.Q.setLayout(self.VBox)
    
    #--------------------------------------------------------
    def updateOSLmisfit(self):
        """Update the misfit calculation: on ages or on n/N values ? """
        pgui.store_input(self.parent.parent,{conf.Variable_names['OSL_misfit_target']: str(self.OSLMisfitCombo.currentIndex())})  
        if self.OSLMisfitCombo.currentIndex() == 0 and self.OSL_Model.currentIndex() == 1: # Age + GOK_Gauss_fad model
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.OSLMisfitCombo.currentIndex() == 1 and self.OSL_Model.currentIndex() == 1: # n/N + GOK_Gauss_fad model
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        
        elif self.OSLMisfitCombo.currentIndex() == 0 and self.OSL_Model.currentIndex() == 2: # Age + GOK_fad model
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK_b","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.OSLMisfitCombo.currentIndex() == 1 and self.OSL_Model.currentIndex() == 2: # n/N + GOK_fad model
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK_b","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
            
        elif self.OSLMisfitCombo.currentIndex() == 1 and self.OSL_Model.currentIndex() == 0: # n/N + band-tail model
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","Eu (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.OSLMisfitCombo.currentIndex() == 0 and self.OSL_Model.currentIndex() == 0: # Age + band-tail model
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","Eu (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
            
    #--------------------------------------------------------
    def updateOSLModel(self):
        """ Update kinetic model for OSL."""
        
        # save parameter for Pecube
        pgui.store_input(self.parent.parent,{conf.Variable_names['OSL Model Label']: str(self.OSL_Model.currentIndex())})
        
        # Change Column Labels according to the model 
        if self.OSL_Model.currentIndex() == 1 and self.OSLMisfitCombo.currentIndex() == 0: # GOK_Gauss_FAD model + age
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.OSL_Model.currentIndex() == 1 and self.OSLMisfitCombo.currentIndex() == 1: # GOK_Gauss_FAD model + n/N
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)

        elif self.OSL_Model.currentIndex() == 2 and self.OSLMisfitCombo.currentIndex() == 0: # GOK_FAD model + age
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK_b","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.OSL_Model.currentIndex() == 2 and self.OSLMisfitCombo.currentIndex() == 1: # GOK_FAD model + n/N
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK_b","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)

        elif self.OSL_Model.currentIndex() == 0  and self.OSLMisfitCombo.currentIndex() == 0: # Band-tail mode + Age
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","Eu (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.OSL_Model.currentIndex() == 0  and self.OSLMisfitCombo.currentIndex() == 1: # Band-tail mode + n/N
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","Eu (eV)","Et (eV)","log_s (1/s)","log_r","GOK_a","Lmax"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
            
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
        try: 
            OSL_Model = self.OSL_Model.currentIndex()
        except:
            OSL_Model = 0

        if h == 2: # update from Table Samples
            try:
                self.GrainTable.cellChanged.disconnect() # Shut down update from the user
            except TypeError:
                pass
            # if (self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 4 or
            # self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 0):
            #     h_temp = 2
            #     h = 1
                
            #Check if we are in "for all nodes"
            if self.parent.parent.AgeCombo.currentIndex() == 1:
               self.TotGrain = 1
               t.setRowCount(1)
               p.TableParameters["SampleName1"] = "Grain"
            elif self.parent.parent.AgeCombo.currentIndex() == 2:
               # Get total number of grains
               self.TotGrain = self.Winparent.Samples['OSL_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['OSL_ngrains'+str(j+1)])
                       if nb_grains > 0:
                           sampleName = self.Winparent.Samples['SampleName'+str(j+1)]
                           p.TableParameters["SampleName"+str(j+1)] = sampleName
        
        
        # Update Table
        if h >= 1:
          # dParam = DefaultParameterValues()
          nb_grains = []
          if self.parent.parent.AgeCombo.currentIndex() == 1:
              nSample = 1
              nb_grains.append(1)
              # Grain Id
              t.setItem(0, 0, QTableWidgetItem(p.DParameters['SampleName']+'1_1'))
              sampleName = p.DParameters['SampleName']+'1_1'
              #Obs Age or nN
              t.setItem(0, 1, QTableWidgetItem(p.DParameters['OSL_NN']))
              #Obs error
              t.setItem(0, 2, QTableWidgetItem(p.DParameters['OSL_DNN']))
              # Dose rate
              t.setItem(0, 3, QTableWidgetItem(p.DParameters['OSL_doser']))
              # Characteristic dose
              t.setItem(0, 4, QTableWidgetItem(p.DParameters['OSL_D0']))
              # Model-specific parameter
              if OSL_Model == 0: # SSE-BTS
                t.setItem(0, 5, QTableWidgetItem(p.DParameters['OSL_Eu']))
              elif OSL_Model == 1: # GOK-GAUSS
                t.setItem(0, 5, QTableWidgetItem(p.DParameters['OSL_Eu']))
              elif OSL_Model == 2: # GOK
                  t.setItem(0, 5, QTableWidgetItem(p.DParameters['OSL_b']))
              # Activation energy
              t.setItem(0, 6, QTableWidgetItem(p.DParameters['OSL_Et']))
              # Escape Frequency factor
              t.setItem(0, 7, QTableWidgetItem(p.DParameters['OSL_logs']))
              # Density of recombination centre
              t.setItem(0, 8, QTableWidgetItem(p.DParameters['OSL_logrho']))
              # GOk a
              t.setItem(0, 9, QTableWidgetItem(p.DParameters['OSL_a']))
              # Maximum intensity
              t.setItem(0, 10, QTableWidgetItem(p.DParameters['OSL_Lmax']))
              

          elif self.parent.parent.AgeCombo.currentIndex() == 2:
              nSample = int(self.Winparent.nbSampleValue.text())
              for j in range(nSample):
                  nb_grains.append(int(self.Winparent.Samples['OSL_ngrains'+str(j+1)]))
              
              count = 0
              for j in range(nSample):
                  if nb_grains[j] > 0:
                      for i in range(nb_grains[j]):
                          if "SampleName"+str(j+1) in p.TableParameters:
                              t.setItem(count, 0, QTableWidgetItem(p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)))
                              sampleName = p.TableParameters['SampleName'+str(j+1)]
                          else:
                              t.setItem(count, 0, QTableWidgetItem(p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)))
                              sampleName = p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)
                          if sampleName+"_"+str(i+1)+"_N/N_OSL" in p.TableParameters:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_N/N_OSL"]))
                          else:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.DParameters['OSL_NN']))
                          if sampleName+"_"+str(i+1)+"_DNN_OSL" in p.TableParameters:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_DNN_OSL" ]))
                          else:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.DParameters['OSL_DNN']))
                          if sampleName+"_"+str(i+1)+"_DOSER_OSL" in p.TableParameters:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_DOSER_OSL"]))
                          else:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.DParameters['OSL_doser']))
                          if sampleName+"_"+str(i+1)+"_D0_OSL" in p.TableParameters:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_D0_OSL"]))
                          else:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.DParameters['OSL_D0']))
                          if sampleName+"_"+str(i+1)+"_EU_OSL" in p.TableParameters:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_EU_OSL"]))
                          else:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.DParameters['OSL_Eu']))
                          if sampleName+"_"+str(i+1)+"_ET_OSL" in p.TableParameters:
                              t.setItem(count, 6, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_ET_OSL"]))
                          else:
                              t.setItem(count, 6, QTableWidgetItem(
                                  p.DParameters['OSL_Et']))
                          if sampleName+"_"+str(i+1)+"_LOGS_OSL" in p.TableParameters:
                              t.setItem(count, 7, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_LOGS_OSL"]))
                          else:
                              t.setItem(count, 7, QTableWidgetItem(
                                  p.DParameters['OSL_logs']))
                          if sampleName+"_"+str(i+1)+"_LOGRHO_OSL" in p.TableParameters:
                              t.setItem(count, 8, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_LOGRHO_OSL"]))
                          else:
                              t.setItem(count, 8, QTableWidgetItem(
                                  p.DParameters['OSL_logrho']))
                          if sampleName+"_"+str(i+1)+"_a_OSL" in p.TableParameters:
                              t.setItem(count, 9, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_a_OSL"]))
                          else:
                              t.setItem(count, 9, QTableWidgetItem(
                                  p.DParameters['OSL_a']))
                          if sampleName+"_"+str(i+1)+"_Lmax_OSL" in p.TableParameters:
                              t.setItem(count, 10, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_Lmax_OSL"]))
                          else:
                              t.setItem(count, 10, QTableWidgetItem(
                                  p.DParameters['OSL_Lmax']))
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
                        float(item.text())
                        temp_dict['OSL_NN_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_NN_'+sampleName] = str(item.text())
                    elif col == 2:
                        float(item.text())
                        temp_dict['OSL_DNN_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_DNN_'+sampleName] = str(item.text())
                    elif col == 3:
                        float(item.text())
                        temp_dict['OSL_Doser_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_Doser_'+sampleName] = str(item.text())
                    elif col == 4:
                        float(item.text())
                        temp_dict['OSL_D0_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_D0_'+sampleName] = str(item.text())
                    elif col == 5:
                        float(item.text())
                        temp_dict['OSL_Eu_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_Eu_'+sampleName] = str(item.text())
                    elif col == 6:
                        float(item.text())
                        temp_dict['OSL_Et_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_Et_'+sampleName] = str(item.text())
                    elif col == 7:
                        float(item.text())
                        temp_dict['OSL_logs_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_logs_'+sampleName] = str(item.text())
                    elif col == 8:
                        float(item.text())
                        temp_dict['OSL_logrho_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_logrho_'+sampleName] = str(item.text())
                    elif col == 9:
                        float(item.text())
                        temp_dict['OSL_a_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_a_'+sampleName] = str(item.text())
                    elif col == 10:
                        float(item.text())
                        temp_dict['OSL_Lmax_'+sampleName] = str(item.text())
                        p.TableParameters['OSL_Lmax_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In OSL observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        # We write for all nodes
        if self.parent.parent.AgeCombo.currentIndex() == 1:
            for k, v in self.grains.items():
                if k.startswith('OSL_Doser'):
                    pgui.store_input(self.parent.parent, {'OSL_doser': str(v)})
                elif k.startswith('OSL_D0'):
                    pgui.store_input(self.parent.parent, {'OSL_D0': str(v)})
                elif k.startswith('OSL_Et'):
                    pgui.store_input(self.parent.parent, {'OSL_Et': str(v)})
                elif k.startswith('OSL_Eu'):
                    pgui.store_input(self.parent.parent, {'OSL_Eu': str(v)})
                elif k.startswith('OSL_b'):
                    pgui.store_input(self.parent.parent, {'OSL_b': str(v)})
                elif k.startswith('OSL_logs'):
                    pgui.store_input(self.parent.parent, {'OSL_logs': str(v)})
                elif k.startswith('OSL_lorho'):
                    pgui.store_input(self.parent.parent, {'OSL_logrho': str(v)})
                elif k.startswith('OSL_a'):
                    pgui.store_input(self.parent.parent, {'OSL_a': str(v)})
                elif k.startswith('OSL_Lmax'):
                    pgui.store_input(self.parent.parent, {'OSL_Lmax': str(v)})
       
       
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))
            


#############################################################################################
########################################## ESR ##############################################
#############################################################################################
class ESR(QWidget):
    """
      This class set the parameters for ESR thermochronometer
      Winparent = parameters in WindowParameters class
    
    """
    def __init__(self, parent, Winparent, Param, PFolder):
        super().__init__()
        
        # Parameters
        self.parent = parent
        self.Winparent = Winparent
        self.grains = {}
        self.Param = Param
        if int(Param.DParameters[conf.Variable_names['ESR Model Label']]) == 0: # SSE-Gauss
            data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                    "Age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                    "Error": ["0"]*int(Param.DParameters['ngrains']),
                    "Dose rate (Gy/kyr)": [Param.DParameters['ESR_doser']]*int(Param.DParameters['ngrains']),
                    "D0 (Gy)": [Param.DParameters['ESR_D0']]*int(Param.DParameters['ngrains']),
                    "sigmaEt": [Param.DParameters['ESR_sigmaEt']]*int(Param.DParameters['ngrains']),
                    "Et (eV)": [Param.DParameters['ESR_Et']]*int(Param.DParameters['ngrains']),
                    "log_s (1/s)": [Param.DParameters['ESR_logs']]*int(Param.DParameters['ngrains']),
                    "Lmax": [Param.DParameters['ESR_Lmax']]*int(Param.DParameters['ngrains']),
                    "Nil": [Param.DParameters['ESR_b']]*int(Param.DParameters['ngrains'])
                    }
        elif int(Param.DParameters[conf.Variable_names['ESR Model Label']]) == 1: # GOK-GOK
            data = {'Grain': ['1']*int(Param.DParameters['ngrains']),
                    "Age (Ma)": ["0"]*int(Param.DParameters['ngrains']),
                    "Error": ["0"]*int(Param.DParameters['ngrains']),
                    "Dose rate (Gy/kyr)": [Param.DParameters['ESR_doser']]*int(Param.DParameters['ngrains']),
                    "D0 (Gy)": [Param.DParameters['ESR_D0']]*int(Param.DParameters['ngrains']),
                    "GOK a": [Param.DParameters['ESR_a']]*int(Param.DParameters['ngrains']),
                    "Et (eV)": [Param.DParameters['ESR_Et']]*int(Param.DParameters['ngrains']),
                    "log_s (1/s)": [Param.DParameters['ESR_logs']]*int(Param.DParameters['ngrains']),
                    "Lmax": [Param.DParameters['ESR_Lmax']]*int(Param.DParameters['ngrains']),
                    "GOK b": [Param.DParameters['ESR_b']]*int(Param.DParameters['ngrains'])
                    }
        # How many grains in total?
        self.TotGrain = 0
        for j in range(int(self.Winparent.nbSampleValue.text())):
                nb_grains = self.Winparent.Samples['ESR_ngrains'+str(j+1)]
                self.TotGrain += int(nb_grains)   
        self.GrainTable = interface.WinTable(data, self.TotGrain, 10)
        self.updateTableGrain(self.GrainTable, 1, Param)
        
        self.ESRDoser = QLabel('Dose rate (Gy/kyr): ')
        self.ESRDoser.setFont(conf.font8)
        self.ESRDoser.setAlignment(Qt.AlignCenter)
        self.ESRDoser.setToolTip(
            "Dose rate (Doser) used for ESR data computation.")
        self.ESRD0 = QLabel('D0 (Gy): ')
        self.ESRD0.setAlignment(Qt.AlignCenter)
        self.ESRD0.setFont(conf.font8)
        self.ESRD0.setToolTip(
            "Onset of dose saturation (D0) for ESR data computation.")
        self.ESR_Et = QLabel('Activation energy Et (eV): ')
        self.ESR_Et.setFont(conf.font8)
        self.ESR_Et.setAlignment(Qt.AlignCenter)
        self.ESR_Et.setToolTip(
            "Activation energy (Et) used for ESR data computation.")
        self.OSL_Eu = QLabel('Activation energy Eu (eV): ')
        self.OSL_Eu.setFont(conf.font8)
        self.OSL_Eu.setAlignment(Qt.AlignCenter)
        self.OSL_Eu.setToolTip(
            "Activation energy (Eu) used for OSL data computation.")
        self.ESR_logs = QLabel('log_s (1/s):')
        self.ESR_logs.setFont(conf.font8)
        self.ESR_logs.setAlignment(Qt.AlignCenter)
        self.ESR_logs.setToolTip(
            "Log of thermal frequancy factor for OSL data computation")
        self.ESR_sigmaEt = QLabel('Sigma_Et (eV):')
        self.ESR_sigmaEt.setFont(conf.font8)
        self.ESR_sigmaEt.setAlignment(Qt.AlignCenter)
        self.ESR_sigmaEt.setToolTip(
            "log of dimensionless recombination center density for OSL data computation")
        self.ESRmisfitTarget_Label = QLabel("Calculate misfit on: ")
        self.ESRmisfitTarget_Label.setFont(conf.font8)
        self.ESRmisfitTarget_Label.setAlignment(Qt.AlignCenter)
        self.ESRmisfitTarget_Label.setToolTip(
            "Choose if you want to calculate misfit on ESr ages or n/N values.")

        # Create text boxes
        self.ESRDoserValue = QLineEdit(self)
        self.ESRDoserValue.setText(Param.DParameters['ESR_doser'])
        self.ESRDoserValue.setAlignment(Qt.AlignCenter)
        self.ESRD0Value = QLineEdit(self)
        self.ESRD0Value.setText(Param.DParameters['ESR_D0'])
        self.ESRD0Value.setAlignment(Qt.AlignCenter)
        self.ESR_EtValue = QLineEdit(self)
        self.ESR_EtValue.setText(Param.DParameters['ESR_Et'])
        self.ESR_EtValue.setAlignment(Qt.AlignCenter)
        self.ESR_sigmaEtValue = QLineEdit(self)
        self.ESR_sigmaEtValue.setText(Param.DParameters['ESR_sigmaEt'])
        self.ESR_sigmaEtValue.setAlignment(Qt.AlignCenter)
        self.ESR_logsValue = QLineEdit(self)
        self.ESR_logsValue.setText(Param.DParameters['ESR_logs'])
        self.ESR_logsValue.setAlignment(Qt.AlignCenter)
        self.ESR_Model_Label = QLabel("Kinetic model:")
        self.ESR_Model_Label.setAlignment(Qt.AlignLeft)
        self.ESR_Model = QComboBox()
        items = ['SSE-Gauss','GOK-GOK']
        self.ESR_Model.addItems(items)
        self.ESR_Model.setCurrentIndex(int(Param.DParameters[conf.Variable_names['ESR Model Label']]))
        self.ESRMisfitCombo = QComboBox()
        items = ['Ages','n/N']
        self.ESRMisfitCombo.addItems(items)
        self.ESRMisfitCombo.setCurrentIndex(int(Param.DParameters[conf.Variable_names['ESR_misfit_target']]))

        # Extra box
        fakeWidget = QLabel()
        self.GroupBox1 = QGroupBox()
        self.VBox = QVBoxLayout()
        self.Grid = QGridLayout()
        self.Grid.setColumnStretch(3,2)
        self.Grid.setSpacing(5)
        self.Grid.addWidget(fakeWidget, 0, 0)
        self.Grid.addWidget(self.ESR_Model_Label,1,1)
        self.Grid.addWidget(self.ESR_Model,1,2)
        self.Grid.addWidget(self.ESRmisfitTarget_Label,2,1)
        self.Grid.addWidget(self.ESRMisfitCombo,2,2)
        self.Grid.addWidget(fakeWidget, 5, 3)
        self.GroupBox1.setLayout(self.Grid)
        self.VBox.addWidget(self.GroupBox1)
        self.VBox.addWidget(self.GrainTable)
        
        # Signals
        self.ESR_Model.currentIndexChanged.connect(lambda: self.updateESRModel())
        self.ESRMisfitCombo.currentIndexChanged.connect(lambda: self.updateESRmisfit())
        
        self.Q = QWidget()
        self.Q.setLayout(self.VBox)
    
    #--------------------------------------------------------
    def updateESRmisfit(self):
        """Update the misfit calculation: on ages or on n/N values ? """

        if self.ESRMisfitCombo.currentIndex() == 0 and self.ESR_Model.currentIndex() == 0: # Age + GOK-GAUSS model
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","Lmax","Nil"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)

        elif self.ESRMisfitCombo.currentIndex() == 1 and self.ESR_Model.currentIndex() == 0: # n/N + GOK-GAUSS model
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","Lmax","Nil"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        
        elif self.ESRMisfitCombo.currentIndex() == 0 and self.ESR_Model.currentIndex() == 1: # Age + GOK-GOK model
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","Lmax","Nil"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)

        elif self.ESRMisfitCombo.currentIndex() == 1 and self.ESR_Model.currentIndex() == 1: # n/N + GOK-GOK model
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK a","Et (eV)","log_s (1/s)","Lmax","GOK b"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        
        pgui.store_input(self.parent.parent,{conf.Variable_names['ESR_misfit_target']: str(self.ESRMisfitCombo.currentIndex())})
        
    #--------------------------------------------------------
    def updateESRModel(self):
        """ Update kinetic model for ESR."""
        
        pgui.store_input(self.parent.parent,{conf.Variable_names['ESR Model Label']: str(self.ESR_Model.currentIndex())})

         # Change Column Labels according to the model 
        if self.ESR_Model.currentIndex() == 0 and self.ESRMisfitCombo.currentIndex() == 0: # SSE_Gauss model + age
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","Lmax","Nil"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.ESR_Model.currentIndex() == 0 and self.ESRMisfitCombo.currentIndex() == 1: # SSE_Gauss model + n/N
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","SigmaEt (eV)","Et (eV)","log_s (1/s)","Lmax","Nil"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.ESR_Model.currentIndex() == 1 and self.ESRMisfitCombo.currentIndex() == 0: # GOK_GOK model + age
            Headers = ["Grain", "Age (Ma)", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK a","Et (eV)","log_s (1/s)","Lmax","GOK b"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        elif self.ESR_Model.currentIndex() == 1 and self.ESRMisfitCombo.currentIndex() == 1: # GOK_Gauss model + n/N
            Headers = ["Grain", "n/N", "Error","Dose rate (Gy/kyr)","D0 (Gy)","GOK a","Et (eV)","log_s (1/s)","Lmax","GOK b"]
            self.GrainTable.setHorizontalHeaderLabels(Headers)
        
    
    #--------------------------------------------------------
    def updateTableGrain(self, t, h, p):
        """
          To store input value for the Grain table
          t = table
          h : parameter
          p : parameters
          d : pecube parameters (window)
        
        """
        try: 
            ESR_Model = self.ESR_Model.currentIndex()
        except:
            ESR_Model = 0

        self.grains = {}
        h_temp = 1
        nSample = 0
        if h == 2: # update from Table Samples
            try:
                self.GrainTable.cellChanged.disconnect() # Shut down update from the user
            except TypeError:
                pass
            # if (self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 4 or
            # self.Winparent.SampleTable.column(self.Winparent.SampleTable.currentItem()) == 0):
            #     h_temp = 2
            #     h = 1
            #Check if we are in "for all nodes"
            if self.parent.parent.AgeCombo.currentIndex() == 1:
               self.TotGrain = 1
               t.setRowCount(1)
               p.TableParameters["SampleName1"] = "Grain"
            elif self.parent.parent.AgeCombo.currentIndex() == 2:
               # Get total number of grains
               self.TotGrain = self.Winparent.Samples['ESR_TotGrain']
               t.setRowCount(int(self.TotGrain))
               # Get Sample ID for each grain, read from sample table
               for j in range(int(self.Winparent.nbSampleValue.text())):
                       nb_grains = int(self.Winparent.Samples['ESR_ngrains'+str(j+1)])
                       if nb_grains > 0:
                           sampleName = self.Winparent.Samples['SampleName'+str(j+1)]
                           p.TableParameters["SampleName"+str(j+1)] = sampleName 
        
        # Update Table
        if h >= 1:
          # dParam = DefaultParameterValues()
          nb_grains = []
          if self.parent.parent.AgeCombo.currentIndex() == 1:
              nSample = 1
              nb_grains.append(1)
              # Grain Id
              t.setItem(0, 0, QTableWidgetItem(p.DParameters['SampleName']+'1_1'))
              sampleName = p.DParameters['SampleName']+'1_1'
              #Obs Age, here always zero
              t.setItem(0, 1, QTableWidgetItem(p.DParameters['ESR_NN']))
              #Obs error, here always zero
              t.setItem(0, 2, QTableWidgetItem(p.DParameters['ESR_DNN']))
              #Radius
              t.setItem(0, 3, QTableWidgetItem(p.DParameters['ESR_doser']))
              t.setItem(0, 4, QTableWidgetItem(p.DParameters['ESR_D0']))
              if ESR_Model == 0:
                t.setItem(0, 5, QTableWidgetItem(p.DParameters['ESR_sigmaEt']))
              elif ESR_Model == 1:
                t.setItem(0, 5, QTableWidgetItem(p.DParameters['ESR_a']))
              t.setItem(0, 6, QTableWidgetItem(p.DParameters['ESR_Et']))
              t.setItem(0, 7, QTableWidgetItem(p.DParameters['ESR_logs']))
              t.setItem(0, 8, QTableWidgetItem(p.DParameters['ESR_Lmax']))
              t.setItem(0, 9, QTableWidgetItem(p.DParameters['ESR_b']))

          elif self.parent.parent.AgeCombo.currentIndex() == 2:
              nSample = int(self.Winparent.nbSampleValue.text())
              for j in range(nSample):
                  nb_grains.append(int(self.Winparent.Samples['ESR_ngrains'+str(j+1)]))
              
              count = 0
              for j in range(nSample):
                  if nb_grains[j] > 0:
                      for i in range(nb_grains[j]):
                          if "SampleName"+str(j+1) in p.TableParameters:
                              t.setItem(count, 0, QTableWidgetItem(p.TableParameters['SampleName'+str(j+1)]+'_'+str(i+1)))
                              sampleName = p.TableParameters['SampleName'+str(j+1)]
                          else:
                              t.setItem(count, 0, QTableWidgetItem(p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)))
                              sampleName = p.DParameters['SampleName']+str(j+1)+'_'+str(i+1)
                          if sampleName+"_"+str(i+1)+"_N/N_ESR" in p.TableParameters:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_N/N_ESR" ]))
                          else:
                              t.setItem(count, 1, QTableWidgetItem(
                                  p.DParameters['Age_grain1']))
                          if sampleName+"_"+str(i+1)+"_DNN_ESR" in p.TableParameters:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_DNN_ESR"]))
                          else:
                              t.setItem(count, 2, QTableWidgetItem(
                                  p.DParameters['Error_grain1']))
                          if sampleName+"_"+str(i+1)+"_DOSER_ESR" in p.TableParameters:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_DOSER_ESR"]))
                          else:
                              t.setItem(count, 3, QTableWidgetItem(
                                  p.DParameters['ESR_doser']))
                          if sampleName+"_"+str(i+1)+"_D0_ESR" in p.TableParameters:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_D0_ESR"]))
                          else:
                              t.setItem(count, 4, QTableWidgetItem(
                                  p.DParameters['ESR_D0']))
                              
                          if sampleName+"_"+str(i+1)+"_SIGMAET_ESR" in p.TableParameters and ESR_Model == 0:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_SIGMAET_ESR"]))
                          elif sampleName+"_"+str(i+1)+"_a_ESR" in p.TableParameters and ESR_Model == 1:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_a_ESR"]))
                          else:
                              t.setItem(count, 5, QTableWidgetItem(
                                  p.DParameters['ESR_sigmaEt']))
                              
                          if sampleName+"_"+str(i+1)+"_ET_ESR" in p.TableParameters:
                              t.setItem(count, 6, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_ET_ESR"]))
                          else:
                              t.setItem(count, 6, QTableWidgetItem(
                                  p.DParameters['ESR_Et']))
                          if sampleName+"_"+str(i+1)+"_LOGS_ESR" in p.TableParameters:
                              t.setItem(count, 7, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_LOGS_ESR"]))
                          else:
                              t.setItem(count, 7, QTableWidgetItem(
                                  p.DParameters['ESR_logs']))
                          if sampleName+"_"+str(i+1)+"_Lmax_ESR" in p.TableParameters:
                              t.setItem(count, 8, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_Lmax_ESR"]))
                          else:
                              t.setItem(count, 8, QTableWidgetItem(
                                  p.DParameters['ESR_Lmax']))
                          if sampleName+"_"+str(i+1)+"_b_ESR" in p.TableParameters:
                              t.setItem(count, 9, QTableWidgetItem(
                                  p.TableParameters[sampleName+"_"+str(i+1)+"_b_ESR"]))
                          else:
                              t.setItem(count, 9, QTableWidgetItem(
                                  p.DParameters['ESR_b']))
                    
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
                        float(item.text())
                        temp_dict['ESR_NN_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_NN_'+sampleName] = str(item.text())

                    elif col == 2:
                        float(item.text())
                        temp_dict['ESR_DNN_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_DNN_'+sampleName] = str(item.text())

                    elif col == 3:
                        float(item.text())
                        temp_dict['ESR_Doser_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_Doser_'+sampleName] = str(item.text())

                    elif col == 4:
                        float(item.text())
                        temp_dict['ESR_D0_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_D0_'+sampleName] = str(item.text())

                    elif col == 5:
                        float(item.text())
                        temp_dict['ESR_sigmaEt_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_sigmaEt_'+sampleName] = str(item.text())

                        temp_dict['ESR_a_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_a_'+sampleName] = str(item.text())

                    elif col == 6:
                        float(item.text())
                        temp_dict['ESR_Et_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_Et_'+sampleName] = str(item.text())

                    elif col == 7:
                        float(item.text())
                        temp_dict['ESR_logs_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_logs_'+sampleName] = str(item.text())

                    elif col == 8:
                        float(item.text())
                        temp_dict['ESR_Lmax_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_Lmax_'+sampleName] = str(item.text())
                        
                    elif col == 9:
                        float(item.text())
                        temp_dict['ESR_b_'+sampleName] = str(item.text())
                        p.TableParameters['ESR_b_'+sampleName] = str(item.text())
                if item is None:
                    break
        except ValueError:
            QErrorMessage(self).showMessage("In ESR observed ages table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        self.grains.update(temp_dict)
        
        # We write for all nodes
        if self.parent.parent.AgeCombo.currentIndex() == 1:
            for k, v in self.grains.items():
                if k.startswith('ESR_Doser'):
                    pgui.store_input(self.parent.parent, {'ESR_doser': str(v)})
                elif k.startswith('ESR_D0'):
                    pgui.store_input(self.parent.parent, {'ESR_D0': str(v)})
                elif k.startswith('ESR_Et'):
                    pgui.store_input(self.parent.parent, {'ESR_Et': str(v)})
                elif k.startswith('ESR_sigmaEt'):
                    pgui.store_input(self.parent.parent, {'ESR_sigmaEt': str(v)})
                elif k.startswith('ESR_logs'):
                    pgui.store_input(self.parent.parent, {'ESR_logs': str(v)})
                elif k.startswith('ESR_Lmax'):
                    pgui.store_input(self.parent.parent, {'ESR_Lmax': str(v)})
                elif k.startswith('ESR_a'):
                    pgui.store_input(self.parent.parent, {'ESR_a': str(v)})
                elif k.startswith('ESR_b'):
                    pgui.store_input(self.parent.parent, {'ESR_b': str(v)})
       
        if h == 2: #update from Table Samples
            self.GrainTable.cellChanged.connect(lambda: self.updateTableGrain(
                self.GrainTable, 0, self.Param))