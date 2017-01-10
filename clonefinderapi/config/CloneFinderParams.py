import os.path

class CloneFinderParams(object):
    """
        Encapsulation of input paramters
    """
    
    def __init__(self):
        self._bino_num = 50
        self._data_format = ''
        self._bino_cutoff = 0.95
        self._error_rate = 0.01
        self._freq_cutoff = 0.01
        self._read_coverage = 100.0
        self._input_id = ''
        self._option_a='On'     #On or Off (Cut extra mutations from a clone) 
        self._option_b='5,3,On' #maximum_number_of_clusters_to_be_produced,Minumum_number_of_SNVs_per_cluster,Decomposed clone has trunc mutaion and is not an ancestral clone)
        self._option_c='On,5'   #On or Off,Number_of_SNVs_per_cluster' (decompose clones based on clone phylogeny) 
        self._option_d='On,3'   #On or Off,Number_of_SNVs_per_cluster (make ancestral clones at middle of a branch)
        self._option_e='On'     #On or Off (Filter backward/paralle mutations)
        self._option_f='On'     #On or Off (Add unassigned SNVs to clones)
        self._option_g='On'     #On or Off (Flag bad data)
        self._option_h='On'     #On or Off (Combine similar clones)                
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

    #On or Off (Cut extra mutations from a clone)
    @property
    def option_a(self):
        return self._option_a    
    @option_a.setter
    def option_a(self, value):
        if (value.lower() != 'on') and (value.lower() != 'off'):
            Exception("invalid option given for OptionA. Must be 'On' or 'Off'")
        self._option_a = value;
        
    @property
    def option_b(self):
        return self._option_b    
    @option_b.setter
    def option_b(self, value):
        tokens = value.split(',')
        if len(tokens) != 3:
            ValueError("invalid value given for OptionB. Should be something like: '5,3,On'")
        if self._is_integer(tokens[0]) == False:
            ValueError("invalid value given for OptionB. Should be something like: '5,3,On'")
        if self._is_integer(tokens[1]) == False:
            ValueError("invalid value given for OptionB. Should be something like: '5,3,On'")  
        if (tokens[2].lower() != 'on') and (tokens[2].lower() != 'off'):
            ValueError("invalid value given for OptionB. Should be something like: '5,3,On'")            
        self._option_b = value;
        
    @property
    def option_c(self):
        return self._option_c
    @option_c.setter
    def option_c(self, value):
        tokens = value.split(',')
        if len(tokens) != 2:
            ValueError("invalid value given for OptionC. Should be something like: 'On,5'")        
        if (tokens[0].lower() != 'on') and (tokens[0].lower() != 'off'):
            ValueError("invalid value given for OptionC. Should be something like: 'On,5'")            
        if self._is_integer(tokens[1]) == False:
            ValueError("invalid value given for OptionC. Should be something like: 'On,5'")          
        self._option_c = value;
        
    @property
    def option_d(self):
        return self._option_d
    @option_d.setter
    def option_d(self, value):
        tokens = value.split(',')
        if len(tokens) != 2:
            ValueError("invalid value given for OptionD. Should be something like: 'On,3'")        
        if (tokens[0].lower() != 'on') and (tokens[0].lower() != 'off'):
            ValueError("invalid value given for OptionD. Should be something like: 'On,3'")            
        if self._is_integer(tokens[1]) == False:
            ValueError("invalid value given for OptionD. Should be something like: 'On,3'")                  
        self._option_d = value;
        
    @property
    def option_e(self):
        return self._option_e    
    @option_e.setter
    def option_e(self, value):
        if (value.lower() != 'on') and (value.lower() != 'off'):
            ValueError("invalid value given for OptionE. Should be something like: '5,3,On'")      
        self._option_e = value;
        
    @property
    def option_f(self):
        return self._option_f    
    @option_f.setter
    def option_f(self, value):
        if (value.lower() != 'on') and (value.lower() != 'off'):
            ValueError("invalid value given for OptionF. Should be something like: '5,3.On'")       
        self._option_f = value;
        
    @property
    def option_g(self):
        return self._option_g    
    @option_g.setter
    def option_g(self, value):
        if (value.lower() != 'on') and (value.lower() != 'off'):
            ValueError("invalid value given for OptionG. Should be something like: '5,3,On'")        
        self._option_g = value;     
        
    @property
    def option_h(self):
        return self._option_h    
    @option_h.setter
    def option_h(self, value):
        if (value.lower() != 'on') and (value.lower() != 'off'):
            ValueError("invalid value given for OptionH. Should be something like: '5,3,On'")        
        self._option_h = value;     
        
    def _is_integer(self, value):
        try:
            int(value)
            return True
        except ValueError:
            return False
        
    def to_string(self):
        result = ''
        my_vars = vars(self)
        for key, value in my_vars.iteritems():
            result += key[1:] + '\t\t' + str(value) + '\n'
        return result
