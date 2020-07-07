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
        for key, value in my_vars.items():
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
    def num_nodes(self):
        return len(self._items)
    
    @property 
    def newick_string(self):
        return self._newick_string
    @newick_string.setter
    def newick_string(self, value):
        self._newick_string = value
        
    def get_alignment(self, remove_duplicates = True):
        temp = []
        result = ''
        temp.append('#MEGA')
        temp.append('!Title ' + self._filename)
        temp.append('!Format datatype=dna;')
        temp.append(' ')
        for item in self._items:
            if item.is_leaf_node:
                name = '#' + item.label
            else:
                name = '#Node' + str(item.des1) + 'S0T' + str(item.des2)
            if remove_duplicates == True:
                if item.seq_data in temp:
                        continue
            temp.append(name)
            temp.append(item.seq_data)
            
        for item in temp:
            result += item + "\n"
        return result            
            
        
    
    