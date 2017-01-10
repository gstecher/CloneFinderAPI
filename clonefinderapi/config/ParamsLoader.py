from config.CloneFinderParams import CloneFinderParams
from parsers.ParamsParser import ParamsParser
import os.path
import sys

class ParamsLoader(object):
    """
        Loads command-line parameters and parses the options.ini file
    """
    
    def __init__(self):
        self._params_file = 'options.ini'
        
    @property 
    def params_file(self):
        return self._params_file
    
    @params_file.setter
    def params_file(self, value):
        self._params_file = value
    
    def load_params(self):    
        parser = ParamsParser()
        if os.path.isfile(self._params_file) == False:
            Exception("The required config.ini file is missing.")        
        result = parser.parse_config_file(self._params_file)        
        print "parsing command-line parameters..."
        Data = sys.argv[1]
        result.data_format = Data
        if Data=='snv': 
            result.input_data_file = sys.argv[2]            
        elif Data=='ccf': 
            result.ccf_data_file = sys.argv[2]
            result.read_coverage= sys.argv[3]	
            result.snv_data_file = sys.argv[2][:-4]+'snv.txt'
        elif Data=='cnv':
            result.input_data_file0 = sys.argv[2]
            result.input_data_file = sys.argv[2][:-4]+'Adj.txt'
        elif Data=='cnv-post':
            result.input_data_file = sys.argv[2]    	
            result.cnv_data_file = sys.argv[2][:-4]+'-CNV.txt'
        else: 
            print 'the command should be python CloneFinder.py snv [input]\npython CloneFinder.py ccf [input]\n or\npython CloneFinder.py cnv [input]\n or \npython CloneFinder.py cnv-post [input]'   	
        result.input_id = sys.argv[2].split('/')[-1][:-4]        
        return result
    