#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to plot output results in the interface. 
It contains classes to read and load model output from Pecube.

@author: maxime Bernard

"""


import os
import xarray as xr
import pandas as pd
import configs as conf

#-----------------------------------------------------------------------------
#-------------------------------- Classes ------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

class dataSets:
    """ 
    This class stores the data read in the output files of Pecube.
    
    """
    def __init__(self,fileList,folderName,SampleFile):
        super().__init__()
        
        self.fileList = fileList
        self.folderName = folderName
        self.SampleFile = SampleFile+".csv"
        self.Data_path = os.path.join(conf.PecubeFolderPath,self.folderName,"Data",SampleFile)
        # self.GrainID = GrainID
        # Datasets array
        self.MainDataset = xr.Dataset()
        
        self.loadData()
        
    #----------------------------------------------------------
    def loadData(self):
        """ 
        To load available data from files provided by the interface.
        The data are read once and stored in a xarray Dataset.
        
        """
        
        # Path to output directory
        self.Output_path = os.path.join(conf.PecubeFolderPath,self.folderName,"Output")
        # print(self.Output_path)
        #### Search file name and load data in xarray dataset ####
        if os.path.exists(os.path.join(self.Output_path,"CompareAGE.csv")):
            # Start with CompareAGE.csv file
            if "CompareAGE.csv" in self.fileList:
                # Read file 
                dataFrame = pd.read_csv(os.path.join(self.Output_path,"CompareAGE.csv"))
                # Transform to xarray dataset
                self.MainDataset = dataFrame.to_xarray()
                # Rename variables
                self.MainDataset = self.MainDataset.rename({"SAMPLE": 'SID',
                                                            "LON":"LONOBS",
                                                            "LAT":"LATOBS",})
                # self.MainDataset.coords["latitude"] = (("latitude"),self.MainDataset.Latitude.data)
                # self.MainDataset.coords["longitude"] = (("longitude"),self.MainDataset.Longitude.data)
                # Get variable names
                varNames = list(self.MainDataset.keys())
                # Assign variable attributes to filter later by name
                for i in varNames:
                    self.MainDataset[i].attrs = {"description": 'CompareAges'}
            
            # Read sample data if any
            if self.SampleFile in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Data_path,self.SampleFile))
                # Transform to xarray dataset
                dataset = dataFrame.to_xarray()
                # Get variable names
                varNames = list(dataset.keys())
                # Assign variable attributes to filter later by name
                for i in varNames:
                    dataset[i].attrs = {"description": 'InputSample'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset])
        
            # Time-temperature paths
            if "TimeTemperaturePaths.csv" in self.fileList:
                dataFrame = pd.read_csv(os.path.join(self.Output_path,"TimeTemperaturePaths.csv"))
                dataset2 = dataFrame.set_index("Time").to_xarray()
                dataset2 = dataset2.rename({"Time": "Time_Tt"})
                # Get variable names
                varNames = list(dataset2.keys())
                # Change variable names according to grain ID
                res = dict(zip(varNames, self.MainDataset.SID.values))
                dataset2 = dataset2.rename_vars(res)
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'Age evolution'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            if "TimeAgeAHe.csv" in self.fileList:
                dataFrame = pd.read_csv(os.path.join(self.Output_path,"TimeAgeAHe.csv"),index_col=False)
                dataset2 = dataFrame.set_index("Time").to_xarray()
                dataset2 = dataset2.rename({"Time": "Time_AgeAHe"})
                # Add Thermochronometer ID in front of variable (Sample) name
                varNames = list(dataset2.keys())
                dataset2 = dataset2.rename_vars(dict(zip(varNames,["AHE"+n for n in varNames])))
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'AHe age-specific'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            if "TimeAgeAFT.csv" in self.fileList:
                dataFrame = pd.read_csv(os.path.join(self.Output_path,"TimeAgeAFT.csv"),index_col=False)
                dataset2 = dataFrame.set_index("Time").to_xarray()
                dataset2 = dataset2.rename({"Time": "Time_AgeAFT"})
                # Add Thermochronometer ID in front of variable (Sample) name
                varNames = list(dataset2.keys())
                dataset2 = dataset2.rename_vars(dict(zip(varNames,["AFT"+n for n in varNames])))
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'AFT age-specific'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            if "TimeAgeZHe.csv" in self.fileList:
                dataFrame = pd.read_csv(os.path.join(self.Output_path,"TimeAgeZHe.csv"),index_col=False)
                dataset2 = dataFrame.set_index("Time").to_xarray()
                dataset2 = dataset2.rename({"Time": "Time_AgeZHe"})
                # Add Thermochronometer ID in front of variable (Sample) name
                varNames = list(dataset2.keys())
                dataset2 = dataset2.rename_vars(dict(zip(varNames,["ZHE"+n for n in varNames])))
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'ZHe age-specific'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
                
            if "Compare43He.csv" in self.fileList:
                dataFrame = pd.read_csv(os.path.join(self.Output_path,"Compare43He.csv"),index_col=False)
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"SAMPLE": 'Sample43Pred'})
                # Add Thermochronometer ID in front of variable (Sample) name
                # varNames = list(dataset2.keys())
                # dataset2 = dataset2.rename_vars(dict(zip(varNames,["ZHE"+n for n in varNames])))
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'Compare43He'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
             
             # Read input 43 profiles
            if conf.Variable_names["43_inputFile"] in self.fileList:
                dataFrame = pd.read_csv(os.path.join(self.Data_path,conf.Variable_names["43_inputFile"]),index_col=False)
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"SAMPLE": 'Sample43',
                                            "LON":"LON43",
                                            "LAT":"LAT43",
                                            "HEIGHT":"HEIGHT43",
                                            "SIZE":"SIZE43",
                                            "UPPM":"UPPM43",
                                            "THPPM":"THPPM43",
                                            conf.Variable_names['FTL_kinetic_parameter_value_AHe']:"FTLKIN43",
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'input43He'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
                
            # Read sample ThL data if any
            if conf.Variable_names['ThL_file'] in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Output_path,conf.Variable_names["ThL_file"]))
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"SAMPLE":"SAMPLE_THL",
                                            "LON":"LON_TL",
                                            "LAT":"LAT_TL",
                                            "HEIGHTOBS":"HEIGHTOBS_THL",
                                            "HEIGHTPRED":"HEIGHTPRED_THL",
                                            "NNOBS":"ThL_NN_OBS",
                                            "NNPRED":"ThL_NN_PRED"
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'ThL_data'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            # Read sample OSL data if any
            if conf.Variable_names['OSL_file'] in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Output_path,conf.Variable_names["OSL_file"]))
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"SAMPLE":"SAMPLE_OSL",
                                            "LON":"LON_OSL",
                                            "LAT":"LAT_OSL",
                                            "HEIGHTOBS":"HEIGHTOBS_OSL",
                                            "HEIGHTPRED":"HEIGHTPRED_OSL",
                                            "OSLOBS":"OSLOBS",
                                            "OSLPRED":"OSLPRED"
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'OSL_data'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            # Read sample ESR data if any
            if conf.Variable_names['ESR_file'] in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Output_path,conf.Variable_names["ESR_file"]))
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"SAMPLE":"SAMPLE_ESR",
                                            "LON":"LON_ESR",
                                            "LAT":"LAT_ESR",
                                            "HEIGHTOBS":"HEIGHTOBS_ESR",
                                            "HEIGHTPRED":"HEIGHTPRED_ESR",
                                            "ESROBS":"ESROBS",
                                            "ESRPRED":"ESRPRED"
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'ESR_data'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            # Read sample ThL data if any
            if conf.Variable_names['ThL_file_input'] in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Data_path,conf.Variable_names["ThL_file_input"]))
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"SAMPLETL":"SAMPLETHL",
                                            "LON":"LON_TL2",
                                            "LAT":"LAT_TL2",
                                            "HEIGHT":"HEIGHT_THL",
                                            "DOSER":"DOSER_THL",
                                            "D0":"D0_THL",
                                            "ET":"ET_THL",
                                            "ATL":"ATL_THL",
                                            "BTL":"BTL_THL",
                                            "LOGS":"LOGS_THL",
                                            "LOGRHO":"LOGRHO_THL",
                                            "N/N" : "ThL_NN",
                                            "DNN" : "ThL_DNN",
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'ThL_data_input'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            # Read sample OSL data if any
            if conf.Variable_names['OSL_file_input'] in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Data_path,conf.Variable_names["OSL_file_input"]))
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"LON":"LON_OSL2",
                                            "LAT":"LAT_OSL2",
                                            "HEIGHT":"HEIGHT_OSL",
                                            "DOSER":"DOSER_OSL",
                                            "D0":"D0_OSL",
                                            "ET":"ET_OSL",
                                            "EU":"EU_OSL",
                                            "LOGS":"LOGS_OSL",
                                            "LOGRHO":"LOGRHO_OSL",
                                            "N/N" : "OSLOBS_input",
                                            "DNN" : "DOSL",
                                            "Lmax": "Lmax_OSL",
                                            "a":"a_OSL"
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'OSL_data_input'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])
            
            # Read sample ESR data if any
            if conf.Variable_names['ESR_file_input'] in self.fileList:
                # Read file
                dataFrame = pd.read_csv(os.path.join(self.Data_path,conf.Variable_names["ESR_file_input"]))
                dataset2 = dataFrame.to_xarray()
                dataset2 = dataset2.rename({"LON":"LON_ESR2",
                                            "LAT":"LAT_ESR2",
                                            "HEIGHT":"HEIGHT_ESR",
                                            "DOSER":"DOSER_ESR",
                                            "D0":"D0_ESR",
                                            "ET":"ET_ESR",
                                            "SIGMAET":"SIGMAET_ESR",
                                            "LOGS":"LOGS_ESR_OBS",
                                            "N/N" : "ESROBS_input",
                                            "DNN" : "DESR",
                                            "Lmax": "Lmax_ESR",
                                            "a":"a_ESR",
                                            "b":"b_ESR"
                                            })
                # Assign variable attributes to filter later by name
                for i in list(dataset2.keys()):
                    dataset2[i].attrs = {"description": 'ESR_data_input'}
                # Merge datasets
                self.MainDataset = xr.merge([self.MainDataset,dataset2])