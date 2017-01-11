from parsers.DefaultTSPParser import DefaultTSPParser
from parsers.DefaultCCFParser import DefaultCCFParser
from alignments.FreqToMegaSeq import FreqToMegaSeq
from parsimony.MegaSimpleMP import MegaSimpleMP
#from binomial.BinomialReplicater import BinomialReplicater
from binomial.SignificanceTester import SignificanceTester
from regression.CloneFrequencyComputer import CloneFrequencyComputer
BinoRepNum=50
significance_cutoff=0.95
ErroRate=0.01
if __name__ == "__main__":
   # filename = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\ccf.txt'
   # parser = DefaultCCFParser()
    filename = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\copyNumberVariationAdj.txt'
    parser = DefaultTSPParser()
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
        tsp_list = parser.get_tumor_sample_profile_list()
        align_builder = FreqToMegaSeq()
        align_builder.initialize(tsp_list, True) # pass in True to remove duplicates, False to keep duplicates
        mega_seqs = align_builder.get_mega_alignment()
      #  MPtree = MegaSimpleMP()
       # MPtree.do_mega_mp(align_builder, 'test')		
        clone_frequency = CloneFrequencyComputer(mega_seqs, tsp_list)	
        significant_clones = SignificanceTester(BinoRepNum, significance_cutoff, ErroRate, tsp_list, mega_seqs, clone_frequency)
        significant_clone_in_tumor, significant_clone_list = significant_clones.count_detection()
        print 	significant_clone_in_tumor, significant_clone_list	
		
   #     clusters = SNPClusterGenerator(mega_seqs, tsp_list, clone_frequencies, OptionB[0].strip())		
        		