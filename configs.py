#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This file contains global variables to be used across python modules of the package,
as well as the default values for the input parameter for Pecube.

@author: Maxime Bernard
"""

from PyQt5.QtGui import (QFont, QDoubleValidator,QBrush,QColor)
from PyQt5.QtCore import (QLocale, QCoreApplication,QFile,QIODevice,QTextStream)
from PyQt5.QtWidgets import QMessageBox,QFileDialog
import os
import sys

# Set a list to store open windows
WindowsOpen = []

##############################################################################
############################ Global parameters ##############################
##############################################################################

########### Path to directories #############
if sys.platform == 'darwin' or sys.platform=='linux':
    #Handle whether application is open from command line or executable
    if getattr(sys, 'frozen', False):
        FolderPath = sys._MEIPASS
    else:
        FolderPath = os.path.abspath(sys.executable)[:-10]
    # FolderPath = "/Users/maxime/Documents/Post-doc/"
    encoding_label = 'cp858'
    NEWLINE_SIZE_IN_BYTES = 1
elif sys.platform == 'win32' or sys.platform =='cygwin':
    FolderPath = os.path.realpath(os.getcwd())
    encoding_label = 'ANSI'
    NEWLINE_SIZE_IN_BYTES = 2
    import win32con, win32api
    
# Locate Icones directory
IconPath = os.path.join(FolderPath, "Icones")


# For the format of doubles
# Force dot notation for the decimals
DoubleValidator = QDoubleValidator()
locale = QLocale(QLocale.English, QLocale.UnitedStates)
DoubleValidator.setLocale(locale)
DoubleValidator.setNotation(QDoubleValidator.StandardNotation)

# All the fonts to used in the interface
font6 = QFont("Helvetica", 6)
font8 = QFont("Helvetica", 8)
fontBold8 =  QFont("Helvetica", 8, QFont.Bold)
font10 = QFont("Helvetica", 10)
fontBold10 = QFont("Helvetica", 10, QFont.Bold)
fontItalic10 = QFont("Helvetica", 10, italic=True)
font11 =  QFont("Helvetica", 11)
fontBold11 =  QFont("Helvetica", 11, QFont.Bold)
font12 = QFont("Helvetica", 12)
fontLine12 = QFont("Helvetica", 12)
fontLine12.setUnderline(True)
fontBold12 =  QFont("Helvetica", 12, QFont.Bold)
fontBold12.setUnderline(True)

# Colors
ColorObservations = QBrush(QColor(183,134,32))
ColorTextDefault = QBrush(QColor(255,255,255))

# The path where the Pecube directory is located
# PecubePath = os.path.join(QCoreApplication.applicationDirPath(),"Pecube")
PecubeFolderPath = os.path.join(FolderPath,"Pecube")
ProjectName = ''

# Get the style of the interface from Combinear.qss file
stream = QFile(os.path.join(FolderPath,"Combinear.qss"))
stream.open(QIODevice.ReadOnly)
styleSheet = QTextStream(stream).readAll()

# Read the preferences file to get the path of the Pecube folder
PreferencesPath = os.path.join(FolderPath,"preferences.txt")
file = open(PreferencesPath, 'r')
PrefList = {}
for line in file:
    if "=" in line:
        (key, _, par) = line.split()
        PrefList[str(key)] = par
if PrefList['Theme'] == 'Dark':
    styleSheet = styleSheet
    appStyle = "Fusion"
else:
    styleSheet = ""
    appStyle = "Fusion"
file.close()
    
PecubeFolderPath = PrefList['PecubePath']
PecubeFolderPath = os.path.abspath(PecubeFolderPath)
                
#############################################################
# Parameters for 3D plots
if PrefList['Theme'] == 'Dark':
    Colors3Dplot = {'axes' : 'w',
                'background':[0.3,0.3,0.3],
                'colorbarLabels':'w'
    	  }
else:
    Colors3Dplot = {'axes' : 'k',
                    'background':[0.9,0.9,0.9],
                    'colorbarLabels':'k'
        	  }    
        
#############################################################
# Default parameters
# Variable_names: where to change/add new parameter labels for Pecube
# Do not change the labels on the left side, this is the labels that are
# read everywhere in PecubeGUI.
Variable_names = {
    'Default_age':'default_age',
    '4He/3He': 'He43_flag',
    'Diffusion_model_AHe':'RDmodel',
    'Mode_age_computation':'age_computation',
    'Edge Age for 4He/3He':'EdgeAge',
    "43_inputFile":"He43_data.csv",
    "Inversion_43":"He43_inv_flag",
    "D0_Biotite" : "D0_BAr",
    "Ea_Biotite" : "Ea_BAr",
    "D0_Feldspar" : "D0_KAr",
    "Ea_Feldspar" : "Ea_KAr",
    "D0_Muscovite" : "D0_MAr",
    "Ea_Muscovite" : "Ea_MAr",
    "D0_Hornblende" : "D0_HAr",
    "Ea_Hornblende" : "Ea_HAr",
    "D0_Apatite" : "D0_AHe",
    "Ea_Apatite" : "Ea_AHe",
    "D0_Zircon" : "D0_ZHe",
    "Ea_Zircon" : "Ea_ZHe",
    'Diffusion_model_ZHe':'RDmodelz',
    "FissionTrackModel":"FT_code_flag",
    "FTL_kinetic_parameter":"Kinetic_FTL_Parameter",
    "FTL_kinetic_parameter_value_AFT":"KINAFT",
    "FTL_kinetic_parameter_value_AHe":"KINAHE",
    "Initial_FTL_model":'Init_FTL_Model',
    "Unannealed_FTL_value": 'Init_FTL_value',
    "Use Fission Track":'age_FTL_flag',
    "Mean Fission Track Length":"MFTL",
    "Mean Fission Track Length error":"DMFTL",
    "Mean_fission_track_length_error_inversion":"MFTL_error",
    "Mean_fission_track_length_std_error_inversion":"MFTL_std_error",
    "Rmr0": "rmr0",
    "eFoldingDepth":"D_heat",
    "Peak threshold": "0.05",
    "Peak min distance":"30",
    "ThL_file":"CompareTL.csv",
    "OSL_file":"CompareOSL.csv",
    "ESR_file":"CompareESR.csv",
    "ThL_file_input":"ThL_data.csv",
    "OSL_file_input":"OSL_data.csv",
    "ESR_file_input":"ESR_data.csv",
    "4He/3He observation file":"43_data.csv",
    "Starting incision": "time_incision_start",
    "Ending incision": 'time_incision_stop',
    'Incision depth':'depth_incision',
    'Scaling_propagation_rate':"tauH",
    'Headward erosion flag': 'do_headward',
    'Scale erosion time':'erosional_time_scale',
    'ThL vs time file':'TLTime.csv',
    'OSL vs time file':'OSLTime.csv',
    'ESR vs time file':'ESRTime.csv',
    'ThL Model Label':'TL_Model',
    'OSL Model Label':'OSL_Model',
    'ESR Model Label':'ESR_Model',
    'Custom ref topo':'topo_ref_custom',
    'ESR_misfit_target': 'ESR_misfit_target',
    'OSL_misfit_target': 'OSL_misfit_target',
    'Save age for batch': 'save_ages_inversion',
    'use inversion or batch': 'inversion_mode',
    'Number intervals for range values': 'ninter_batch',
    'File Age ESR batch': 'AgeESR.csv',
    'File Age AHe batch': 'AgeApatiteHelium.csv',
    'File Age AFT batch': 'AgeApatiteFT.csv',
    'File Age ZHe batch': 'AgeZirconHelium.csv',
    'File Age ZFT batch': 'AgeZirconFT.csv',
    'File Age KAr batch': 'AgeKSparArgon.csv',
    'File Age BAr batch': 'AgeBiotiteArgon.csv',
    'File Age MAr batch': 'AgeMuscoviteArgon.csv',
    'File Age HAr batch': 'AgeHornblendeArgon.csv',
    'File int inversion': 'NA_int_res.csv',
    'Synthetic topo file name': 'Synth',
    'minimum time step':'min_dt',
    'Heat_production_flag':'heatprod_flag'
    }


class DefaultParameterValues:
    """
    set the default parameters here
    
    """
    
    def __init__(self):
        super().__init__()
        # Topographic Parameters
        self.DParameters = {
            'topo_file_name': 'Nil',
            'nx': '31',
            'ny': '31',
            'lon0': '0',
            'lat0': '0',
            'dlon': '0.0083333',
            'dlat': '0.0083333',
            'nskip': '1',
            'Raster_flag':'0',
            # SPM files default
            'initTime_spm':'0.0',
            'endTime_spm':'0.0',
            'nstep_spm':'0',
            'npreStep':'0',
            'uRateSpm':'0.0',
            'Tsl_spm':'0.0',
            'lrate_spm':'0.0',
            # Time evolution parameters
            'ntime': '1',
            Variable_names['Scale erosion time']: '0',
            'time_topo1': '0',
            'topo_wavelength':'10.0',
            'topo_offset':'0.0',
            'topo_amp':'1.0',
            'topo_phase':'0.0',
            'amplification1': '1',
            'offset1': '0',
            'output1': '0',
            'phase1':'0',
            'topo_ref': '0',
            Variable_names['Headward erosion flag']: '0.0',
            Variable_names['Starting incision']: '0.0',
            Variable_names['Ending incision']: '0.0',
            Variable_names['Incision depth']: '0.0',
            Variable_names['Scaling_propagation_rate']: '0.0',
            Variable_names['Custom ref topo']: '0.0',
            # Thermal parameters
            'thickness': '35',
            'basal_temperature': '700',
            Variable_names['Heat_production_flag'] : '0',
            'nz': '21',
            'sea_level_temperature': '0',
            'lapse_rate': '0',
            'thermal_diffusivity': '25',
            'heat_production': '0',
            'Radioactive_heat':'0',
            'Specific_heat':'880',
            'surfHeatFlux':'0',
            Variable_names['minimum time step']: '1',
            Variable_names['eFoldingDepth']: '20',
            # Data parameters
            'data_folder': 'Nil',
            Variable_names['Default_age']: '-1',
            # Tectonic parameters
            'nfault': '0',
            'lon1': '0',
            'lat1': '0',
            'lon2': '0',
            'lat2': '0',
            'npoint1': '0',
            'bottom_left': '1',
            'bottom_right': '1',
            'top_right': '1',
            'top_left': '1',
            'nstep1': '0',
            'uplift':'0',
            'fault_advect_flag':'0',
            'logarithmic_velocity':'0',
            'uplift':'-1',
            # Output parameters
            'PDModel_signal':'1',
            'age_AHe_flag':'0',
            'age_AFT_flag':'0',
            'age_ZHe_flag':'0',
            'age_ZFT_flag':'0',
            'age_KAr_flag':'0',
            'age_BAr_flag':'0',
            'age_MAr_flag':'0',
            'age_HAr_flag':'0',
            Variable_names['Use Fission Track']:'0',
            Variable_names["Mean Fission Track Length"]:'0',
            Variable_names["Mean Fission Track Length error"]:'0',
            Variable_names["Mean_fission_track_length_error_inversion"]:'0.4',
            Variable_names["Mean_fission_track_length_std_error_inversion"]:'0.3',
            'age_OSL_flag':'0',
            'age_TL_flag':'0',
            'age_ESR_flag':'0',
            'TL_NN':'0.0',
            'TL_DNN':'0.0',
            'OSL_NN':'0.0',
            'OSL_DNN':'0.0',
            'TL_doser': '5.0',
            'TL_D0':'800.0',
            'TL_a':'1',
            'TL_b':'1',
            'TL_Et':'1.4',
            'TL_logs':'12.0',
            'TL_logrho':'-5.5',
            'OSL_doser':'5.0',
            'OSL_D0':'800.0',
            'OSL_Et':'1.4',
            'OSL_Eu':'0.12',
            'OSL_logs':'12.0',
            'OSL_logrho':'-5.5',
            'OSL_SigmaEt':'0.0',
            'OSL_a':'1.0',
            'OSL_b':'1.0',
            'OSL_Lmax': '1.0',
            Variable_names['OSL_misfit_target']:'0',
            'ESR_doser':'4.265',
            'ESR_NN':'0.0',
            'ESR_DNN':'0.0',
            'ESR_D0':'5186.8',
            'ESR_Et': '1.763',
            'ESR_sigmaEt': '0.096',
            'ESR_logs': '14.58',
            'ESR_a': '1',
            'ESR_b': '1',
            'ESR_Lmax': '1',
            'echo_input_file': '0',
            Variable_names['ESR_misfit_target']:'0',
            Variable_names['ThL Model Label']: '0',
            Variable_names['OSL Model Label']: '0',
            Variable_names['ESR Model Label']: '0',
            Variable_names['Mode_age_computation']: '0',
            Variable_names['FissionTrackModel']:'2',
            Variable_names['Initial_FTL_model']:'0',
            Variable_names['Unannealed_FTL_value']:'16.3',
            Variable_names['FTL_kinetic_parameter']:'4',
            Variable_names['FTL_kinetic_parameter_value_AFT']:'0.79',
            Variable_names['FTL_kinetic_parameter_value_AHe']:'0.79',
            Variable_names['Diffusion_model_AHe']: '0',
            Variable_names['D0_Apatite']: '50',
            Variable_names['Ea_Apatite']: '137.6536',
            Variable_names['Diffusion_model_ZHe']:'5',
            Variable_names['D0_Zircon']: '0.46',
            Variable_names['Ea_Zircon']: '169.0336',
            'rhoST':'0.893',
            Variable_names['Rmr0']: '0.79',
            'ASIZE':'60',
            'AUPPM':'10.0',
            'AThPPM':'40.0',
            'ZSIZE':'60',
            'ZSize':'60',
            'ZUPPM':'100.0',
            'ZUppm':'100.0',
            'ZThPPM':'400.0',
            'ZThppm':'400.0',
            'grainradius':'0.0',
            'Alpha_Ejec_Flag': '2',
            'Alpha_Ejec_Mod_Flag': '1',
            Variable_names['D0_Feldspar'] : '0.014',
            Variable_names['Ea_Feldspar'] : '120.0',
            Variable_names['D0_Biotite'] : '0.4',
            Variable_names['Ea_Biotite'] : '211.09',
            Variable_names['D0_Muscovite'] : '0.04',
            Variable_names['Ea_Muscovite'] : '217.36',
            Variable_names['D0_Hornblende'] : '0.06',
            Variable_names['Ea_Hornblende'] : '276.0',
            'ngrains': '0',
            'ngrains1':'0',
            'lat_grain1': '0',
            'lon_grain1': '0',
            'Elev_grain1': '0',
            'Age_grain1': '0',
            'Error_grain1': '0',
            'SFe3':'0',
            'SampleName':'sample',
            Variable_names['4He/3He']:'0',
            'save_PTT_paths':'0',
            'save_cooling_rates':'0',
            'save_eroded_volume':'0',
            'Nb_Samples': '0',
            'NMonteIter': '10000',
            # 43He parameters
            'Duration1':[0.38,0.38,0.51,0.66,0.66,0.46,0.45,
                        0.48,0.66,0.53,0.48,0.50,0.56,0.63,
                        0.50,0.50,0.50,0.50,0.50,0.50],
            'Heating1':[200.0,270.0,290.0,300.0,310.0,330.0,
                          340.0,350.0,350.0,370.0,400.0,410.0,
                          420.0,440.0,460.0,475.0,500.0,600.0,
                          700.0,900.0],
            'Grain43':'Grain',
            'Lat43':0.0,
            'Lon43':0.0,
            'Radius43':60.0,
            'Age43':0.0,
            'DAge43':0.0,
            'F3He':0.0,
            'DF3He':0.0,
            'RsRb':0.0,
            'DRsRb':0.0,
            'NstepsHe43':20,
            'Elev43':0.0,
            'fileName43':Variable_names["43_inputFile"],
            Variable_names['Inversion_43']: '0',
            # Isostasy parameters
            'EET': '20',
            'isostasy':'0',
            'rho_crust': '2700.0',
            'rho_astenosphere': '3150.0',
            'young_modulus': '1e11',
            'poisson_ratio': '0.25',
            'nx_isostasy': '1024',
            'ny_isostasy': '1024',
            # Inversion parameters
            Variable_names['Save age for batch']: '0',
            Variable_names['use inversion or batch']: '0',
            Variable_names['Number intervals for range values']: '1',
            'error_predictions':'0.0',
            'NbCores':'0',
            'maxNIt': '4',
            'sampleSizeFirst': '8',
            'sampleSizeOther': '8',
            'nbCellResample': '4',
            'misfit_weight_AHE': '1.0',
            'misfit_weight_FTLD': '1.0',
            'misfit_weight_TH': '1.0',
            'misfit_weight_43He': '1.0',
            'misfit_weight_TL': '1.0',
            'misfit_function':'1',
            'misfit_slope':'1.0',
            'misfit_weight_ESR':'1.0',
            'misfit_weight_OSL':'1.0',
            'misfit_function':'1',
            'maximum_number_of_iterations': '4',
            'sample_size_for_first_iteration': '8',
            'sample_size_for_all_other_iterations':'8',
            'number_of_cells_to_resample':'4',
            # Miscelllaneous flags
            'shear_heating': '0', 
            'debug':'0',
            #Extra GUI parameters
            'oldInput':'0',
            'TurnOffPecube':'0',
            #Zonation
            'zonation_flag':'0',
            'nLayers':'0'}
        self.TableParameters = {'nb_grains': '1',}

