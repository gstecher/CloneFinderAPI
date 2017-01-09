from parsers.AbstractCCFParser import AbstractCCFParser
from cancer_cell_fractions.CancerCellFraction import CancerCellFraction
import os.path

class DefaultCCFParser(AbstractCCFParser):
    """
        Parser for ccf text files wich generates a CancerCellFractionList
    """
    
    def _parse_header(self, header):
        header = header.strip().split("\t")
        
    def _parse_row(self, row):
        row = row.strip().split("\t")
        if len(row) != 3:
            Exception('invalid format. Expected 3 columns, found ' + str(len(row)))
        fraction = CancerCellFraction()
        fraction.ecm1 = float(row[0])
        fraction.ecm2 = float(row[1])
        fraction.ecpt = float(row[2])
        self._fraction_list.add(fraction)
                
    def parse(self):
        try:
            if os.path.isfile(self._input_data_file) == False:
                IOError('input data file not found')
            data = open(self._input_data_file, 'r').readlines()
            self._parse_header(data[0])
            index = 1
            while index < len(data):
                self._parse_row(data[index])
                index += 1
            return True
        except Exception as e:
            self.messages.append(str(e))
            return False
        
    def get_cancer_cell_fraction_list(self):
        return self._fraction_list

