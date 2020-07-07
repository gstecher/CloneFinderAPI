from config.CloneFinderParams import CloneFinderParams
from parsers.ParamsParser import ParamsParser
import os.path
import sys
from config.FormatInput import FormatInput
clone_finder_params = None # used as a global instance of the params object

class ParamsLoader(object):
    """
        Loads command-line parameters and parses the options.ini file
    """
    
    def __init__(self):
        global clone_finder_params
        self._params_file = 'options.ini'
        
    @property 
    def params_file(self):
        return self._params_file
    
    @params_file.setter
    def params_file(self, value):
        self._params_file = value
    
    def load_params(self):    
        parser = ParamsParser()
        adjust_format = FormatInput()		
        if os.path.isfile(self._params_file) == False:
            Exception("The required config.ini file is missing.")        
        result = parser.parse_config_file(self._params_file)        
        print("parsing command-line parameters...")
        Data = sys.argv[1]
        result.data_format = Data
        if Data=='snv': 
            result.snv_data_file = sys.argv[2]  
            result.cnv_data_file = sys.argv[2][:-4]+'snv-CNV.txt'			
            adjust_format.snv2snv(result.snv_data_file,'withCNVfile')			
            result.input_data_file = sys.argv[2][:-4]+'snv.txt'			
        else: 
            print('the command should be python CloneFinder.py snv [input]\npython CloneFinder.py ')   	

        result.input_id = sys.argv[2][:-4] +Data		
        clone_finder_params = result
        return clone_finder_params
    