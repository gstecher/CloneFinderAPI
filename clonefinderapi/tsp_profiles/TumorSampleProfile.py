class TumorSampleProfile(object):
    
    """
        encapsulation of a tumor sample profile
    """
            
    def __init__(self, name):
        if name != '':
            self.name = name
        else:
            self.name = 'SNVs'
        self._read_counts = []
        
    def __iter__(self):
        return iter(self._read_counts)
        
    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value
        
    @property
    def count(self):
        return len(self._read_counts)
    
    def ref_count(self, index):
        if (index >= 0) and (index < len(self._read_counts)):
            return self._read_counts[index].num_ref
        IndexError('index out of bounds')        

    def alt_count(self, index):
        if (index >= 0) and (index < len(self._read_counts)):
            return self._read_counts[index].num_alt
        IndexError('index out of bounds') 
        
    def count_id(self, index):
        if (index >= 0) and (index < len(self._read_counts)):
            return self._read_counts[index].id
        IndexError('index out of bounds')           
        
    def ref_frequency(self, index):
        if (index >= 0) and (index < len(self._read_counts)):
            return self._read_counts[index].ref_frequency()
        IndexError('index out of bounds')
        
    def alt_frequency(self, index):
        if (index >= 0) and (index < len(self._read_counts)):
            return self._read_counts[index].alt_frequency()
        IndexError('index out of bounds')
        
    def get_read_count(self, index):
        if (index >= 0) and (index < len(self._read_counts)):
            return self._read_counts[index]
        IndexError('index out of bounds')
        
    def add(self, read_count):
        self._read_counts.append(read_count)
        
    def get_alignment_string(self):
        result = ''
        if len(self._read_counts) > 0:
            for read_count in self._read_counts:
                if read_count.num_alt > 0:
                    result += 'T'
                else:
                    result += 'A'
        return result
            
       