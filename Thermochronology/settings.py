"""This file gather all global variables for use with thermochronometers

This is where you can change for instance:
- the marker style for each thermochronometers for plotting
- their maker face color
- etc...

"""

import configs as conf

###################### Parameters for plotting #########################

# Assign each thermochronometer a marker shape
Marker_dict = {
    'AHE':'o',
    'AFT':'d',
    'ZHE':'s',
    'ZFT':'o',
    'KAR':'o',
    'BAR':'v',
    'MAR':'^',
    'HAR':'o',
    'ThL':'o',
    'OSL':'o',
    'ESR':'o'
    
    }

#################################################
def dict_pecube():
    """
	Definition of dictionnaries

	You may have to change it depending on your Pecube version, and on your settings...
    
    Author: from Xavier Robert, University of Grenoble
	"""

    agecol = {'alt' : 'HEIGHT',
    	      'AHE' : 'AHE',
        	  'AFT' : 'AFT',
			  'ZHE' : 'ZHE',
        	  'ZFT' : 'ZFT',
	          'KAR' : 'KAR',
    	      'BAR' : 'BAR',
        	  'MAR' : 'MAR',
	          'HAR' : 'HAR',
    	      'FTL' : 'FT',
              'THL' : 'ThL_NN_',
              'ESR' : 'ESR',
              'OSL' : 'OSL',
              'MFTL': conf.Variable_names["Mean Fission Track Length"],
              'DMFTL': conf.Variable_names["Mean Fission Track Length error"]
        	  }

    # errname: respective column number of the error on data in the data input file 
    errname = {'AHE' : 'DAHE',
			   'AFT' : 'DAFT',
			   'ZHE' : 'DZHE',
			   'ZFT' : 'DZFT',
			   'KAR' : 'DKAR',
        	   'BAR' : 'DBAR',
	           'MAR' : 'DMAR',
    	       'HAR' : 'DHAR',
               'HE43':'darel',
               'S3HE':'drel',
               'THL' : 'ThL_DNN',
               'ESR' : 'DESR',
               'OSL' : 'DOSL',
               'MFTL': conf.Variable_names["Mean Fission Track Length error"]
			   }

	# agename: legend of each data system         
    agename = {'AHE' : 'AHe (Ma)',
    	       'AFT' : 'AFT (Ma)',
        	   'ZHE' : 'ZHe (Ma)',
               'ZFT' : 'ZFT (Ma)',
    	       'KAR' : 'KAr (Ma)',
        	   'BAR' : 'Biot. Ar (Ma)',
               'MAR' : 'Musc. Ar (Ma)',
    	       'HAR' : 'Hb Ar (Ma)',
        	   'FTL' : 'FT length (µm)',
               'alt': 'elevation (m)',
               'AHE43':'AGE43',
               'THL' : 'ThL',
               'OSL' : 'OSL',
               'ESR' : 'ESR',
               'MFTL': conf.Variable_names["Mean Fission Track Length"] + '(µm)',
               'DMFTL': conf.Variable_names["Mean Fission Track Length error"]
               }

	# predname: legend of each predicted system         
    predname = {'AHE' : 'Predicted AHe',
    	        'AFT' : 'Predicted AFT',
        	    'ZHE' : 'Predicted ZHe',
	            'ZFT' : 'Predicted ZFT',
    	        'KAR' : 'Predicted KAr',
        	    'BAR' : 'Predicted Biot. Ar',
	            'MAR' : 'Predicted Musc. Ar',
    	        'HAR' : 'Predicted Hb Ar',
        	    'FTL' : 'Predicted FT length',
				'alt' : 'Predicted elevation',
                'THL' : 'Predicted ThL',
                'OSL' : 'Predicted OSL',
                'ESR' : 'Predicted ESR',
                'MFTL': 'Predicted '+conf.Variable_names["Mean Fission Track Length"] + '(µm)',
                'DMFTL': 'Predicted '+ conf.Variable_names["Mean Fission Track Length error"]+ '(µm)'
	           }

    # colores: Colors used for the different age system
    colores = {'AHE' : 'y',
    	       'AFT' : 'r',
        	   'ZHE' : 'g',
	           'ZFT' : 'b',
    	       'KAR' : 'k',
        	   'BAR' : 'c',
	           'MAR' : 'm',
    	       'HAR' : '0.75',
        	   'FTL' : 'y',
			   'alt' : 'y',
               'He43obs':'k',
               'He43pred':'b',
               'THL' : '#cca300',
               'OSL' : '#cc2900',
               'ESR' : '#0033cc',
               'MFTL': 'm',
               'DMFTL': 'orange'
	           }

    profdict = {'Latitude'  : 'LAT', 
                'Longitude' : 'LON', 
                'Altitude'  : 'HEIGHT', 
                'Projected' : 'PROJ'}
    
    # For 4He/3He thermochronometry
    He43dict = {
              '4HE' : 'arel', 
              'S3He' : 'rel',
              'S3HePred': 'RELEASEDPRED',
              '4HePred': 'AGERELEASEDPRED',
              'Duration'  : 'DUR43',
              'Temperature' : 'TEMP43',
              'Ejection':'ARELEJEC'}
    
    # coloramp: for gradient colors
    coloramp = {'AHE' : 'y',
    	       'AFT' : 'r',
        	   'ZHE' : 'g',
	           'ZFT' : 'b',
    	       'KAR' : 'k',
        	   'BAR' : 'c',
	           'MAR' : 'm',
    	       'HAR' : '0.75',
        	   'FTL' : 'y',
			   'alt' : 'y',
               'He43obs':'k',
               'He43pred':'b',
               'THL' : ['#cca300','#ffeb99'],
               'OSL' : ['#cc2900','#ffc2b3'],
               'ESR' : ['#0033cc','#ccd9ff'],
               'MFTL' : 'm',
	           }
    
    ageMarker = {
        'AHE':'o',
        'AFT':'d',
        'ZHE':'s',
        'ZFT':'s',
        'KAR':'o',
        'BAR':'v',
        'MAR':'^',
        'HAR':'o',
        'THL':'o',
        'OSL':'o',
        'ESR':'o',
        'MFTL':'d',
        'DMFTL':'d'
        
        }

    return agecol, errname, agename, predname, colores, profdict,He43dict,coloramp,ageMarker