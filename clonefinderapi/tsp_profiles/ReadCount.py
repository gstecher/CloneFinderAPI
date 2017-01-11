class ReadCount(object):
    
    """
        A data structure encapsulating the number of reference and alternate
        alleles for a given locus for a single tumor sample profile
    """
    
    def __init__(self):
        self._id = -1
        self._num_alt = 0
        self._num_ref = 0
        
    def total(self): 
        result = self._num_alt + self._num_ref
        return result

    def ref_frequency(self):
        if self.total() == 0:
            return 0.0
        return 1.0 * self._num_ref / self.total()

    def alt_frequency(self):
        if self.total() == 0:
            return 0.0		
        return 1.0 * self._num_alt / self.total()
    
    @property
    def id(self):
        return self._id
    
    @id.setter
    def id(self, value):
        self._id = value
        
    @property
    def num_alt(self):
        return self._num_alt
    @num_alt.setter
    def num_alt(self, value):
        self._num_alt = value
        
    @property
    def num_ref(self):
        return self._num_ref
    
    @num_ref.setter
    def num_ref(self, value):
        self._num_ref = value        
              
      