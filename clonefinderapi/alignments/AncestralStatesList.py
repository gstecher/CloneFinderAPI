class AncestralStatesList(object):
    
    def __init__(self):
        self._items = []
        self._filename = ''
        self._newick_string = ''
        
    def __iter__(self):
        return iter(self._items)
    
    def __str__(self):
        result = ''
        my_vars = vars(self)
        for key, value in my_vars.iteritems():
            result += key[1:] + '\t\t' + str(value) + '\n'
        return result        
        
    def add_state(self, state):
        self._items.append(state)
        
    def item(self, index):
        return self._items[index]
    
    @property 
    def filename(self):
        return self._filename
    @filename.setter
    def filename(self, value):
        self._filename = value
        
    @property
    def num_sites(self):
        if len(self._items) > 0:
            return self._items[0].num_sites
        return 0
    
    @property
    def num_taxa(self):
        return len(self._items)
    
    @property 
    def newick_string(self):
        return self._newick_string
    @newick_string.setter
    def newick_string(self, value):
        self._newick_string = value
    
    