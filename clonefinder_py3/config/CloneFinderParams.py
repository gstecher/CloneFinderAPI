import os.path

class CloneFinderParams(object):
    """
        Encapsulation of input paramters
    """
    
    def __init__(self):
        self._data_format = ''
        self._freq_cutoff = 0.02
        self._total_read_cut=50
        self._mutant_read_cut=2	
        self._input_id = ''             
        self._input_data_file = ''
        self._input_data_file0 = ''
        self._ccf_data_file = ''
        self._cnv_data_file = ''
        self._snv_data_file = ''
        
        
    @property 
    def input_data_file(self):
        return self._input_data_file
    @input_data_file.setter
    def input_data_file(self, value):
        if os.path.isfile(value) == False:
            IOError("input data file not found")
        self._input_data_file = value;

    @property 
    def input_data_file0(self):
        return self._input_data_file0
    @input_data_file0.setter
    def input_data_file0(self, value):
        if os.path.isfile(value) == False:
            IOError("input data file not found")
        self._input_data_file0 = value;
        
    @property 
    def ccf_data_file(self):
        return self._ccf_data_file
    @ccf_data_file.setter
    def ccf_data_file(self, value):
        if os.path.isfile(value) == False:
            IOError("ccf data file not found")
        self._ccf_data_file = value;

    @property 
    def cnv_data_file(self):
        return self._cnv_data_file
    @cnv_data_file.setter
    def cnv_data_file(self, value):
        if os.path.isfile(value) == False:
            IOError("cnv data file not found")
        self._cnv_data_file = value;
        
    @property 
    def snv_data_file(self):
        return self._snv_data_file
    @snv_data_file.setter
    def snv_data_file(self, value):
        if os.path.isfile(value) == False:
            IOError("snv data file not found")
        self._snv_data_file = value;
        
    @property
    def input_id(self):
        return self._input_id
    @input_id.setter
    def input_id(self, value):
        self._input_id = value
        
    @property 
    def read_coverage(self):
        return self._read_coverage
    @read_coverage.setter
    def read_coverage(self, value):
        self._read_coverage = value
        
    @property
    def data_format(self):
        return self._data_format    
    @data_format.setter
    def data_format(self, value):
        if (value.lower() != 'cnv') and (value.lower() != 'cnv-post') and (value.lower() != 'snv') and (value.lower() != 'ccf'):
            Exception('invalid data format. Must be one of the following: cnv cnv-post snv ccf')
        self._data_format = value
        
    @property
    def bino_num(self):
        return self._bino_num    
    @bino_num.setter
    def bino_num(self, value):
        self._bino_num = value
	
    @property
    def total_read_cut(self):
        return self._total_read_cut    
    @total_read_cut.setter
    def total_read_cut(self, value):
        self._total_read_cut = value

    @property
    def mutant_read_cut(self):
        return self._mutant_read_cut    
    @mutant_read_cut.setter
    def mutant_read_cut(self, value):
        self._mutant_read_cut = value		
        
    @property
    def bino_cutoff(self):
        return self._bino_cutoff    
    @bino_cutoff.setter
    def bino_cutoff(self, value):
        self._bino_cutoff = value
        
    @property
    def error_rate(self):
        return self._error_rate    
    @error_rate.setter
    def error_rate(self, value):
        self._error_rate = value
        
    @property 
    def freq_cutoff(self):
        return self._freq_cutoff    
    @freq_cutoff.setter
    def freq_cutoff(self, value):
        self._freq_cutoff = value
 
    def _is_integer(self, value):
        try:
            int(value)
            return True
        except ValueError:
            return False
        
    def to_string(self):
        result = ''
        my_vars = vars(self)
        for key, value in my_vars.items():
            result += key[1:] + '\t\t' + str(value) + '\n'
        return result
