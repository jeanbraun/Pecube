""" 

 This module gathers all parameters, class and functions related to 4He/3He
thermochronometry.

@author: maxime Bernard

"""

from PyQt5.QtWidgets import (QPushButton,QWidget,QLabel,QCheckBox,QVBoxLayout,
                             QGridLayout,QComboBox,QHBoxLayout,QTreeView,QListView,
                             QErrorMessage,QLineEdit,QTableWidgetItem,
                             )
from PyQt5.QtGui import (QIntValidator)
import configs as conf
from PyQt5.QtCore import Qt
import Utils.PGUI_utils as pgu
import Utils.interface as interface
import os
import numpy as np
import xarray as xr
from PyQt5.Qt import QStandardItemModel, QStandardItem


################################################################################################
################################## apatite ###################################################
################################################################################################
class apatite(QWidget):
    """
      This class set the extra parameters for 4He/3He profiles computation. Provide
      Observed profiles with following information:
          Sample,LAT,LON,HEIGHT,GRAINSIZE,AGE43,DAGE43,DUR43,TEMP43,F3He,DF3He,RSRB,DRSRB
    
      parent = AgeParameters class
      system = parent (either AHe or ZHe)
    
    """
    def __init__(self, parent, system, Param, nSamples,Samples,PFolder):
        super().__init__()
        
        # to store parameters for 43He
        self.He43Storage = {}
        
        ########## Parameters ###########
        # self.nbSampleValue = nSamples
        self.Samples = Samples
        self.Param = Param
        self.PFolder = PFolder
        self.parent = parent
        self.system = system
        # self.HSParameters = QLabel("Parameters:")
        # self.HSParameters.setFont(conf.fontBold12)
        self.SampleListLabel = QLabel("Choose grain:")
        self.SampleListLabel.setFont(conf.font10)
        self.SampleListLabel.setAlignment(Qt.AlignCenter)
        self.SampleList = QComboBox()
        self.SampleList.setMinimumContentsLength(15)
        self.HeatingSchedule = QLabel("Observations/Schedule:")
        self.HeatingSchedule.setFont(conf.fontBold12)
        self.NbstepsLabel = QLabel('Number of steps:')
        self.NbstepsLabel.setFont(conf.font10)
        self.NbstepsLabel.setAlignment(Qt.AlignCenter)
        self.NbstepsLabel.setToolTip("Number of steps for the heating schedule.")
        self.NbstepsEdit = QLineEdit()
        self.NbstepsEdit.setValidator(QIntValidator(0,100))
        self.NbstepsEdit.setText(str(self.Param.DParameters["NstepsHe43"]))
        self.CheckAllSampleLabel = QLabel("Set sample specific:")
        self.CheckAllSampleLabel.setAlignment(Qt.AlignRight)
        self.CheckAllSampleLabel.setFont(conf.font10)
        self.CheckAllSampleLabel.setToolTip("Set the heating schedule for each sample.")
        self.CheckAllSample = QCheckBox()
        self.ExplainingText = QLabel('Before providing 4He/3He data, make sure '+
                                     'to first provide the grain characteristics (radius, Uppm etc...).'+
                                     '\nSelect the grain you want to provide 4He/3He data, fill the table in, '+
                                     'then click on "Add grain" to add the data on the list for predictions.'+
                                     '\nTo delete data, check the grain label on the list on the left and click "Remove grain".')
        self.ExplainingText.setFont(conf.font11)
        self.ExplainingText.setAlignment(Qt.AlignLeft)
        self.AddButton = QPushButton("Add grain")
        self.RemoveButton = QPushButton("Remove grain")
        # List of samples, the sample IDs are taken from the GrainsParameters class
        self.SListLabel = QLabel('Select grain(s): ')
        self.SListLabel.setFont(conf.font10)
        self.SListLabel.setAlignment(Qt.AlignRight)
        try:
            self.SampleList_temp = self.parent.apatiteHelium.grains['GrainIDs']
        except AttributeError:
            self.SampleList_temp = []
        items = [name for name in self.SampleList_temp]
        self.SampleList.addItems(items)
        self.SampleList.setCurrentIndex(0)
        
        # self.SampleListCombo.addItems(SampleList)
        # self.SampleListCombo.setMinimumContentsLength(10)
        # Initiate the dataset
        self.dicoData = {}
        
        # How many samples in total?
        data = {
                'Duration (hours)': ['0']*int(Param.DParameters['NstepsHe43']),
                "Heating (°C)": ['0']*int(Param.DParameters['NstepsHe43']),
                "F3He":['0']*int(Param.DParameters['NstepsHe43']),
                "Error F":['0']*int(Param.DParameters['NstepsHe43']),
                "Rs/Rb":['0']*int(Param.DParameters['NstepsHe43']),
                "Error R":['0']*int(Param.DParameters['NstepsHe43'])}
        self.HeatingTable = interface.WinTable(data, int(self.NbstepsEdit.text()), 6)
        self.updateTableHeating(1)
        
        #self.HeatingTable.cellChanged.connect(lambda: self.updateTableHeating(0))
        self.NbstepsEdit.editingFinished.connect(lambda: self.updateTableHeating(1))
        self.SampleList.currentIndexChanged.connect(lambda: self.updateSampleList())
        self.AddButton.clicked.connect(lambda: self.add_data())
        self.RemoveButton.clicked.connect(lambda: self.remove_data())
    
        # Tree
        self.treeView = QTreeView()
        self.treeView.setHeaderHidden(True)
        self.rootNode = QListView()
        self.treeModel = QStandardItemModel(self.rootNode)
        self.treeView.setModel(self.treeModel)
        self.treeView.expandAll()
        self.itemsTree = []
       
        #Grid
        Grid1 = QGridLayout()
        Grid1.setSpacing(10)
        Grid1.addWidget(self.HeatingTable,0,3,1,4)
        Grid1.addWidget(self.treeView,0,0,1,2)

        hbox = QHBoxLayout()
        hbox.addWidget(self.SampleListLabel)
        hbox.addWidget(self.SampleList)
        hbox.addWidget(self.NbstepsLabel)
        hbox.addWidget(self.NbstepsEdit)
        hbox.addWidget(self.AddButton)
        hbox.addWidget(self.RemoveButton)
        hbox.addStretch(2)
        
        self.VBox = QVBoxLayout()
        self.VBox.setSpacing(20)
        self.VBox.addWidget(self.HeatingSchedule)
        self.VBox.addWidget(self.ExplainingText)
        self.VBox.addLayout(hbox)
        self.VBox.addLayout(Grid1,5)
    
    
    #----------------------------------------------
    def add_data(self):
        
        """ To add data to the list of predictions..."""
        # First check if item already exists
        GrainsID = []
        n = 0
        while self.treeModel.item(n):
            if self.treeModel.item(n).isCheckable():
                GrainsID.append(self.treeModel.item(n).text())
                n+=1
        if self.SampleList.currentText() in GrainsID:
            QErrorMessage(self).showMessage("This grain has already been provided. Please, remove it before updating.")
        else:
            item = pgu.StandardItem(self.SampleList.currentText(), 10,set_checked=False)
            self.treeModel.appendRow(item)
            self.treeView.setModel(self.treeModel)
            #Read data in the table
            self.read_table(self.HeatingTable)
    
    #----------------------------------------------
    def remove_data(self):
        
        """ To remove data from the list of predictions..."""
        #Delete from tree view
        i = 0
        try:
            while self.treeModel.item(i):
                if self.treeModel.item(i).isCheckable():
                    if self.treeModel.item(i).checkState() == Qt.Checked:
                        self.treeModel.removeRow(i)
                        i -= 1
                i+=1
        except IndexError:
            pass
        
    #----------------------------------------------
    def get_data(self): 
        """
            To count the number of items that are checked and get Grain ID.
            It means the number of grain the user wants to predict 4He/3He data

        """
        
        # Get nstep
        nstep = int(self.NbstepsEdit.text())
        grainIDList = []
        LatList = []
        LonList = []
        ElevList = []
        RadiusList = []
        AgeList = []
        DAgeList = []
        UList = []
        ThList = []
        rmr0List = []
        dictionary = {}
        SampleName_temp = self.SampleList.currentText()
        lengthSampleName = len(SampleName_temp)
        index = pgu.findOccurrences(SampleName_temp,"_")
        SampleName = SampleName_temp[:index[-1]] # last occurence is grain ID
        print('Sample_ID : ', SampleName)
        for i in range(nstep):
            index = self.SampleList.currentIndex() # Sample ID 
            grainIDList.append(self.SampleList.currentText())
            LatList.append(self.parent.parent.Samples[SampleName+'_lat'])
            LonList.append(self.parent.parent.Samples[SampleName+'_lon'])
            ElevList.append(self.parent.parent.Samples[SampleName+'_Elev'])

            RadiusList.append( self.parent.apatiteHelium.grains['ASIZE_'+self.SampleList_temp[index]])
            UList.append(self.parent.apatiteHelium.grains['AUPPM_'+self.SampleList_temp[index]])
            ThList.append(self.parent.apatiteHelium.grains['ATHPPM_'+self.SampleList_temp[index]])
            rmr0List.append(self.parent.apatiteHelium.grains[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]+"_"+self.SampleList_temp[index]])
            AgeList.append(self.system.grains[self.system.ID+'_'+self.SampleList_temp[index]])
            DAgeList.append(self.system.grains['D'+self.system.ID+'_'+self.SampleList_temp[index]])
                
        dictionary['GrainID'] = grainIDList
        dictionary['Lat'] = LatList
        dictionary['Lon'] = LonList
        dictionary['Elev'] = ElevList
        dictionary['Radius'] = RadiusList
        dictionary['Age'] = AgeList
        dictionary['DAge'] = DAgeList
        dictionary['U'] = UList
        dictionary['Th'] = ThList
        dictionary[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]] = rmr0List
        
        return dictionary
    
    #----------------------------------------------
    def updateSampleList(self):
        # Update nstep
        try:
            self.dicoData["SAMPLE_"+self.SampleList.currentText()]
            self.NbstepsEdit.setText(str(self.dicoData["Nstep_"+self.SampleList.currentText()]))
        except KeyError:
            self.NbstepsEdit.setText('20')
        # update the table
        self.updateTableHeating(1)
    
    #----------------------------------------------
    def read_table(self,table):
        """ Read the table's values when a cell is changed"""
        
        num_rows, num_cols = self.HeatingTable.rowCount(), self.HeatingTable.columnCount()
        grainId = self.SampleList.currentText()
        DUR43 = []
        TEMP43 = []
        rel = []
        drel = []
        arel = []
        darel = []
        try:
            for j in range(num_rows):
                DUR43.append(table.item(j, 0).text())
                TEMP43.append(table.item(j, 1).text())
                rel.append(table.item(j, 2).text())
                drel.append(table.item(j, 3).text())
                arel.append(table.item(j, 4).text())
                darel.append(table.item(j, 5).text())
            self.dicoData["DUR43_"+grainId] = np.asarray(DUR43)
            self.dicoData["TEMP43_"+grainId] = np.asarray(TEMP43)
            self.dicoData["rel_"+grainId] = np.asarray(rel)
            self.dicoData["drel_"+grainId] = np.asarray(drel)
            self.dicoData["arel_"+grainId] = np.asarray(arel)
            self.dicoData["darel_"+grainId] = np.asarray(darel)
        except:
            print('Error in read_table 43HEParameters')
            
    #----------------------------------------------
    def updateTableHeating(self, h):
        """
          To store input value for the Grain table
          h : if 1, update table according to default parameters
          The data are stored in a xarray dataset with attribute '43He data'
        """

        if self.SampleList.currentIndex() > -1:
            dico = self.get_data()
        else:
            dico = {}
        self.nGrains = 1
        nSteps = int(self.NbstepsEdit.text())
        #Update number of lines
        nrows = self.nGrains*nSteps
        # Build dataset
        try:
            dico['GrainID']
            GrainNames = np.transpose(np.asarray(dico['GrainID']))
            Latitude = np.asarray(dico['Lat'])
            Longitude = np.asarray(dico['Lon'])
            Elevation = np.asarray(dico['Elev']) 
            Radius = np.asarray(dico['Radius'])
            Age = np.asarray(dico['Age'])
            DAge = np.asarray(dico['DAge'])
            UPPM = np.asarray(dico['U'])
            THPPM = np.asarray(dico['Th'])
            RMR0 = np.asarray(dico[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]])
            
        except KeyError:
            GrainNames = np.transpose(np.asarray(nSteps*[self.Param.DParameters['Grain43']]))
            Latitude = np.asarray([self.Param.DParameters['Lat43'] for i in range(nSteps)])
            Longitude = np.asarray([self.Param.DParameters['Lon43'] for i in range(nSteps)])
            Elevation = np.asarray([self.Param.DParameters['Elev43'] for i in range(nSteps)]) 
            Radius = np.asarray([self.Param.DParameters['Radius43'] for i in range(nSteps)])
            Age = np.asarray([self.Param.DParameters['Age43'] for i in range(nSteps)])
            DAge = np.asarray([self.Param.DParameters['DAge43'] for i in range(nSteps)])
            UPPM = np.asarray([self.Param.DParameters['AUPPM'] for i in range(nSteps)])
            THPPM = np.asarray([self.Param.DParameters['AThPPM'] for i in range(nSteps)])
            RMR0 = np.asarray([self.Param.DParameters[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]] for i in range(nSteps)])
        
        #Check if data already have been provided
        try:
            self.dicoData["SAMPLE_"+str(GrainNames[0])]
            Duration = self.dicoData["DUR43_"+str(GrainNames[0])]
            Heating = self.dicoData["TEMP43_"+str(GrainNames[0])]
            F3He = self.dicoData["rel_"+str(GrainNames[0])]
            DF3He = self.dicoData["drel_"+str(GrainNames[0])]
            RsRb = self.dicoData["arel_"+str(GrainNames[0])]
            DRsRb = self.dicoData["darel_"+str(GrainNames[0])]
        except KeyError: 
            Array = np.zeros((nSteps,1))
            Duration = np.asarray(self.Param.DParameters['Duration1'])
            Heating = np.asarray(self.Param.DParameters['Heating1'])
            F3He = np.transpose(Array*float(self.Param.DParameters['F3He']))[0] 
            DF3He = np.transpose(Array*float(self.Param.DParameters['DF3He']))[0] 
            RsRb = np.transpose(Array*float(self.Param.DParameters['RsRb']))[0]
            DRsRb = np.transpose(Array*float(self.Param.DParameters['DRsRb']))[0]

        self.dicoData["SAMPLE_"+str(GrainNames[0])] = GrainNames
        self.dicoData["LON_"+str(GrainNames[0])]= Longitude
        self.dicoData["LAT_"+str(GrainNames[0])] = Latitude
        self.dicoData["HEIGHT_"+str(GrainNames[0])] = Elevation
        self.dicoData["SIZE_"+str(GrainNames[0])] = Radius
        self.dicoData["AGE43_"+str(GrainNames[0])] = Age
        self.dicoData["DAGE43_"+str(GrainNames[0])] = DAge
        self.dicoData["UPPM_"+str(GrainNames[0])] = UPPM
        self.dicoData["THPPM_"+str(GrainNames[0])] = THPPM
        self.dicoData["DUR43_"+str(GrainNames[0])] = Duration
        self.dicoData["TEMP43_"+str(GrainNames[0])] = Heating
        self.dicoData["rel_"+str(GrainNames[0])] = F3He
        self.dicoData["drel_"+str(GrainNames[0])] = DF3He
        self.dicoData["arel_"+str(GrainNames[0])] = RsRb
        self.dicoData["darel_"+str(GrainNames[0])] = DRsRb
        self.dicoData[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]+"_"+str(GrainNames[0])] = RMR0
        self.dicoData["Nstep_"+str(GrainNames[0])] = nSteps
        
        # Feed Table with data values
        if h:
            self.HeatingTable.setRowCount(nrows)
            num_rows, num_cols = self.HeatingTable.rowCount(), self.HeatingTable.columnCount()
            for j in range(self.nGrains):
                for n, value in enumerate(self.dicoData["DUR43_"+str(GrainNames[0])]):
                    self.HeatingTable.setItem(n,0, QTableWidgetItem(str(value)))
                    
                for n, value in enumerate(self.dicoData["TEMP43_"+str(GrainNames[0])]):
                    self.HeatingTable.setItem(n,1, QTableWidgetItem(str(value)))
                    
                for n, value in enumerate(self.dicoData["rel_"+str(GrainNames[0])]):
                    self.HeatingTable.setItem(n,2, QTableWidgetItem(str(value)))
                    
                for n, value in enumerate(self.dicoData["drel_"+str(GrainNames[0])]):
                    self.HeatingTable.setItem(n,3, QTableWidgetItem(str(value)))
                    
                for n, value in enumerate(self.dicoData["arel_"+str(GrainNames[0])]):
                    self.HeatingTable.setItem(n,4, QTableWidgetItem(str(value)))
                    
                for n, value in enumerate(self.dicoData["darel_"+str(GrainNames[0])]):
                    self.HeatingTable.setItem(n,5, QTableWidgetItem(str(value)))
                    
        
    #----------------------------------------------
    def read_in_file(self):
        """
            Read values in the table and update dataset.
            
        """
        Duration = np.zeros((int(self.NbstepsEdit.text()),self.nGrains))
        Heating = np.zeros((int(self.NbstepsEdit.text()),self.nGrains))
        F3He = np.zeros((int(self.NbstepsEdit.text()),self.nGrains))
        DF3He = np.zeros((int(self.NbstepsEdit.text()),self.nGrains)) 
        RsRb = np.zeros((int(self.NbstepsEdit.text()),self.nGrains)) 
        DRsRb = np.zeros((int(self.NbstepsEdit.text()),self.nGrains)) 
        num_rows, num_cols = self.HeatingTable.rowCount(), self.HeatingTable.columnCount()
        try:
            count_steps = 0
            count_grain = 0
            for row in range(num_rows):
                for col in range(num_cols):
                    item = self.HeatingTable.item(row, col)
                    if col == 0:
                        Duration[count_steps][count_grain] = float(item.text())
                    if col == 1:
                        Heating[count_steps][count_grain] = float(item.text())
                    elif col == 2:
                        F3He[count_steps][count_grain] = float(item.text())
                    elif col == 3:
                        DF3He[count_steps][count_grain] = float(item.text())
                    elif col == 4:
                        RsRb[count_steps][count_grain] = float(item.text())
                    elif col == 5:
                        DRsRb[count_steps][count_grain] = float(item.text())
                count_steps += 1
                if count_steps == int(self.NbstepsEdit.text()):
                    count_steps = 0
                    count_grain += 1

        except ValueError:
            QErrorMessage(self).showMessage("In 4He/3He table: \n'"+item.text()+"' is not a float number. Please, provide a float number.")
            return
        # Update Dataset
        self.dicoData["DUR43_"+str(GrainNames[0])] = Duration
        self.ds.TEMP43.values = Heating
        self.ds.rel.values = F3He
        self.ds.drel.values = DF3He
        self.ds.arel.values = RsRb
        self.ds.darel.values = DRsRb
        
    #----------------------------------------------
    def write_in_file(self):
        """
          To store the 4He/3He data in a xarray dataset
          attribute is "43He data"
        """
        # Make sure to read the latest values in the table
        # self.read_in_file()
        # Get grain Id of selected grains
        n = 0
        GrainsID = []
        while self.treeModel.item(n):
            if self.treeModel.item(n).isCheckable():
                GrainsID.append(self.treeModel.item(n).text())
                n += 1
        # Combine data from all selected grains
        dico = {
            "SAMPLE":[],
            "LON":[],
            "LAT":[],
            "HEIGHT":[],
            "SIZE":[],
            "AGE43":[],
            "DAGE43":[],
            "UPPM":[],
            "THPPM":[],
            conf.Variable_names["FTL_kinetic_parameter_value_AHe"]:[],
            "DUR43":[],
            "TEMP43":[],
            "rel":[],
            "drel":[],
            "arel":[],
            "darel":[],
            }
        for name in GrainsID:
            dico["SAMPLE"].append(self.dicoData["SAMPLE_"+name])
            dico["LON"].append(self.dicoData["LON_"+name])
            dico["LAT"].append(self.dicoData["LAT_"+name])
            dico["HEIGHT"].append(self.dicoData["HEIGHT_"+name])
            dico["SIZE"].append(self.dicoData["SIZE_"+name])
            dico["AGE43"].append(self.dicoData["AGE43_"+name])
            dico["DAGE43"].append(self.dicoData["DAGE43_"+name])
            dico["UPPM"].append(self.dicoData["UPPM_"+name])
            dico["THPPM"].append(self.dicoData["THPPM_"+name])
            dico[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]].append(self.dicoData[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]+"_"+name])
            dico["DUR43"].append(self.dicoData["DUR43_"+name])
            dico["TEMP43"].append(self.dicoData["TEMP43_"+name])
            dico["rel"].append(self.dicoData["rel_"+name])
            dico["drel"].append(self.dicoData["drel_"+name])
            dico["arel"].append(self.dicoData["arel_"+name])
            dico["darel"].append(self.dicoData["darel_"+name])
        
        # indexes = [j+1+(nsteps+1)*l for j in range(nsteps) for l in range(self.nGrains)]
        name = os.path.join(conf.PecubeFolderPath,self.PFolder,"data",self.parent.parent.dataFolderEdit1.text(),conf.DefaultParameterValues().DParameters['fileName43'])
        self.ds = xr.Dataset(
            data_vars=dict(
                SAMPLE = (["x"], np.concatenate(dico["SAMPLE"])),
                LON = (["x"], np.concatenate(dico["LON"])),
                LAT = (["x"], np.concatenate(dico["LAT"])),
                HEIGHT = (["x"], np.concatenate(dico["HEIGHT"])),
                SIZE = (["x"], np.concatenate(dico["SIZE"])),
                AGE43 = (["x"], np.concatenate(dico["AGE43"])),
                DAGE43 = (["x"], np.concatenate(dico["DAGE43"])),
                UPPM = (["x"], np.concatenate(dico["UPPM"])),
                THPPM = (["x"], np.concatenate(dico["THPPM"])),
                KINAHE = (["x"], np.concatenate(dico[conf.Variable_names["FTL_kinetic_parameter_value_AHe"]])),
                DUR43 = (["x"], np.concatenate(dico["DUR43"])),
                TEMP43 = (["x"], np.concatenate(dico["TEMP43"])),
                rel = (["x"], np.concatenate(dico["rel"])),
                drel = (["x"], np.concatenate(dico["drel"])),
                arel = (["x"], np.concatenate(dico["arel"])),
                darel = (["x"], np.concatenate(dico["darel"]))
                ),
            attrs = dict(description="43He data")
            )
        try:
            df = self.ds.to_dataframe()
            # Rework the dataframe for Pecube, some parts have to be empty
            # df.loc[df.index[indexes],"SAMPLE"] = None
            # df.loc[df.index[indexes],"LON"] = None
            # df.loc[df.index[indexes],"LAT"] = None
            # df.loc[df.index[indexes],"HEIGHT"] = None
            # df.loc[df.index[indexes],"SIZE"] = None
            # df.loc[df.index[indexes],"AGE43"] = None
            # df.loc[df.index[indexes],"DAGE43"] = None
            # df.loc[df.index[indexes],"UPPM"] = None
            # df.loc[df.index[indexes],"THPPM"] = None
            # df.loc[df.index[indexes],"RMR0"] = None
            # Write csv file    
            dummy = df.to_csv(name,index=False)
        except PermissionError:
            if os.path.exists(name):
                # os.chmod(os.path.join(self.PFolder,"data",Fname),0o777)
                os.chmod(name,0o666)
                file = open(name, 'w+')
        except OSError:
            os.mkdir(os.path.join(conf.PecubeFolderPath,self.PFolder,"data",self.parent.parent.dataFolderEdit1.text()))
            df = self.ds.to_dataframe().to_csv(name,index=False)