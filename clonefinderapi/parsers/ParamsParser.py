from config.CloneFinderParams import CloneFinderParams
import ConfigParser
import os.path

class ParamsParser(object):
    
    def __init__(self):        
        self._messages = []
        
    @property 
    def messages(self):
        return self._messages

    def parse_config_file(self, config_file):        
        if os.path.isfile(config_file) == False:
            IOError('input config file not found')        
        parser = ConfigParser.RawConfigParser()
        try:
            print "loading config file: " + config_file
            params = CloneFinderParams()
            parser.read(config_file)            
            params.bino_num = parser.getint('parameters', 'BinoNum')
            params.bino_cutoff = parser.getfloat('parameters', 'BinoCutoff')
            params.error_rate = parser.getfloat('parameters', 'ErrorRate')
            params.freq_cutoff = parser.getfloat('parameters', 'FreqCutoff')
            params.option_a = parser.get('options', 'OptionA')
            params.option_b = parser.get('options', 'OptionB')
            params.option_c = parser.get('options', 'OptionC')
            params.option_d = parser.get('options', 'OptionD')
            params.option_e = parser.get('options', 'OptionE')
            params.option_f = parser.get('options', 'OptionF')
            params.option_g = parser.get('options', 'OptionG')
            params.option_h = parser.get('options', 'OptionH')
            print 'config parameters loaded successfully'
            return params
        except Exception as e:
            print(str(e))
            self._messages.append(str(e))
            return False