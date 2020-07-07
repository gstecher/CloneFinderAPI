from abc import ABCMeta, abstractmethod
from tsp_profiles.TumorSampleProfileList import TumorSampleProfileList

class AbstractTSPParser(object, metaclass=ABCMeta):
    
    def __init__(self):
        self._messages = []
        self._input_data_file = ''
        self.tumor_sample_profiles = TumorSampleProfileList()
        
    @abstractmethod
    def parse(self): pass
        
    @property 
    def messages(self):
        return self._messages
    
    @property
    def input_data_file(self):
        return self._input_data_file
    
    @input_data_file.setter
    def input_data_file(self, value):
        self._input_data_file = value
        
