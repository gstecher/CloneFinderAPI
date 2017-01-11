from parsers.AncestralStatesParser import AncestralStatesParser
from alignments.AncestralStatesList import AncestralStatesList
from parsers.DefaultTSPParser import DefaultTSPParser
from alignments.FreqToMegaSeq import FreqToMegaSeq
from parsimony.MegaSimpleMP import MegaSimpleMP

import os.path
import sys


if __name__ == "__main__":
#    print sys.path
    filename = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/copyNumberVariationAdj.txt'
#    filename = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/ancestral_states.txt'
#    parser = AncestralStatesParser()
#    parser.input_file_name = filename
#    if parser.parse() == False:
#        print parser.messages
#    else:
#        states = parser.get_ancestral_states()
#        for state in states:
#            print state
    parser = DefaultTSPParser()
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
        tsp_list = parser.get_tumor_sample_profile_list()
        align_builder = FreqToMegaSeq()
        align_builder.initialize(tsp_list, True) # pass in True to remove duplicates, False to keep duplicates
        mega_seqs = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/mega_alignment.meg'
        tree_builder = MegaSimpleMP()
        tree_builder.mao_file = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/infer_MP_nucleotide.mao'
        id = 'mega_alignment' # id will be used internally for file names
        status = tree_builder.do_mega_mp(align_builder, id)
        if status == True:
            alignment = tree_builder.alignment_least_back_parallel_muts()
            print 'best alignment'
            print alignment
        else:
            print 'failed to run megaMP'
        
