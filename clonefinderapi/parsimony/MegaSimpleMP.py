from alignments.FreqToMegaSeq import FreqToMegaSeq
from parsers.AncestralStatesParser import AncestralStatesParser
import os
import tempfile

"""
    A wrapper around MEGA-CC for Maximum Parsimony tree construction
    
    Given a MEGA sequence alignment, an instance of this class can generate
    an MP tree (or multiple equally parsimonious trees) using MEGA-CC and 
    construct ancestral sequences for all trees generated. The ancestral sequences
    are stored in the ancestral_states_list property, where each element in the list
    also stores the newick tree for that set of ancestral states
"""

class MegaSimpleMP(object):
    
    def __init__(self):
        self.mao_file = 'infer_MP_nucleotide.mao'
        self._alignment_file = ''        
        self._summary = ''
        self._newick_file = ''
        self._num_trees = 0
        self._mega_id = ''
        self._ancestral_states_list = []
        self._newick_trees = []
        self._temp_dir = tempfile.gettempdir() + os.sep
        
    def do_mega_mp(self, alignment_builder, mega_id):
        print 'constructing MP tree'
        result = False
        self._update_file_names(mega_id)        
        alignment_builder.save_mega_alignment_to_file(self._alignment_file)
        cl = self._command_line_string()
        os.system(cl)
        if os.path.isfile(self._newick_file) == True:
            result = True
            nf = open(self._newick_file, 'r')
            ns = nf.readlines()
            print 'MP tree(s):'
            for line in ns:
                print line
                self._newick_trees.append(line)
            nf.close()
            self._retrieve_ancestral_states()
            self._cleanup_temp_files()
        return result
        
    def _update_file_names(self, mega_id):        
        
        print 'executing megacc parsimony tree construction in ' + self._temp_dir
        self._mega_id = mega_id
        self._alignment_file = self._temp_dir + mega_id + '.meg'
        self._newick_file = self._temp_dir + mega_id + '.nwk'
        
    def _command_line_string(self):
        return 'megacc -a ' + self._mao_file + ' -d ' + self._alignment_file + ' -o ' + self._newick_file
    
    def _retrieve_ancestral_states(self):        
        files = self._get_ancestral_states_files()
        for file in files:
            parser = AncestralStatesParser()
            parser.input_file_name = file
            if not parser.parse() == True:
                IOError('failed to parse ancestral states file')
            states_list = parser.get_ancestral_states()
            self._ancestral_states_list.append(states_list)
            os.remove(file)
        
    def _cleanup_temp_files(self):    
        os.remove(self._alignment_file)
        os.remove(self._newick_file)
        summary_file = self._temp_dir + self._mega_id + '_summary.txt'
        os.remove(summary_file)
        
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
    
    
    def _get_ancestral_states_files(self):
        result = []
        filename = self._temp_dir + self._mega_id + '_ancestral_states.txt'
        if os.path.isfile(filename) == True:
            result.append(filename)
        index = 0
        filename = filename = self._temp_dir + self._mega_id + '_ancestral_states_' + str(index) + '.txt'
        while os.path.isfile(filename) == True:
            result.append(filename)
            index += 1
            filename = filename = self._temp_dir + self._mega_id + '_ancestral_states_' + str(index) + '.txt'
        return result
    
    @property
    def num_trees(self):
        return len(self._ancestral_states_files)
    
    @property 
    def ancestral_states_list(self):
        return self._ancestral_states_list
    
    @property 
    def newick_trees(self):
        return self._newick_trees
        
