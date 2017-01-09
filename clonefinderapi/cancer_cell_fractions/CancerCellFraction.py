class CancerCellFraction(object):
    
    """
        Encapsulation of data contained in ccf text files
    """
    
    def __init__(self):
        self._ecm1 = 0.0
        self._ecm2 = 0.0
        self._ecpt = 0.0

        
    @property
    def ecm1(self):
        return self._ecm1
    
    @ecm1.setter
    def ecm1(self, value):
        self._ecm1 = value
        
    @property
    def ecm2(self):
        return self._ecm2
    
    @ecm2.setter
    def ecm2(self, value):
        self._ecm2 = value   
        
    @property
    def ecpt(self):
        return self._ecpt
    
    @ecpt.setter
    def ecpt(self, value):
        self._ecpt = value        