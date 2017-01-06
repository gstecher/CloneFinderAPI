from abc import ABCMeta, abstractmethod
from cancer_cell_fractions.CancerCellFractionList import CancerCellFractionList

class AbstractCCFParser(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        self._messages = []
        self._input_data_file = ''
        self._fraction_list = CancerCellFractionList()
        
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
    