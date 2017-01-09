class MegaSimpleMP(object):
    
    def __init__(self):
        self.mao_file = 'infer_MP_nucleotide.mao'
        self._alignment_file = ''
        
        self._summary = ''
        self._ancestral_states_files = []        
        self._newick_file = ''
        self._num_trees = 0
        
        
    def do_mega_mp(self, alignment_file_name):
        result = False
        
        
    @property
    def mao_file(self):
        return self._mao_file
    
    @mao_file.setter
    def mao_file(self, value):
        self._mao_file = value
        
    @property 
    def alignment_file(self):
        return self._alignment_file
    
    @alignment_file.setter
    def alignment_file(self, value):
        self._alignment_file = value
        
    @property
    def newick_file(self):
        return self._newick_file
    
    @property
    def ancestral_states_files(self):
        return self._ancestral_states_files
    
    @property
    def num_trees(self):
        return len(self._ancestral_states_files)
        
