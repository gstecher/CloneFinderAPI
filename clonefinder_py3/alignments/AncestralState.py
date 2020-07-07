class AncestralState(object):
    
    def __init__(self):
        self._index = -1
        self._label = ''
        self._des1 = -1
        self._des2 = -1
        self._seq_data = ''
        
    def __str__(self):
        result = ''
        my_vars = vars(self)
        for key, value in my_vars.items():
            result += key[1:] + '\t\t' + str(value) + '\n'
        return result        
        
    @property 
    def index(self):
        return self._index
    @index.setter
    def index(self, value):
        self._index = value
        
    @property
    def label(self):
        return self._label
    @label.setter
    def label(self, value):
        self._label = value
        
    @property 
    def des1(self):
        return self._des1
    @des1.setter
    def des1(self, value):
        self._des1 = value
        
    @property 
    def des2(self):
        return self._des2
    @des2.setter
    def des2(self, value):
        self._des2 = value        
                
    @property 
    def num_sites(self):
        return len(self._seq_data)
        
    @property 
    def seq_data(self):
        return self._seq_data
    @seq_data.setter
    def seq_data(self, value):
        self._seq_data = value
        
    @property
    def is_leaf_node(self):
        return ((self._des1 > 0) and (self._des2 > 0))