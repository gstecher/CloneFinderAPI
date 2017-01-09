from parsers.DefaultTSPParser import DefaultTSPParser
from parsers.DefaultCCFParser import DefaultCCFParser
#from binomial.MakeBinomialReplicates import MakeBinomialReplicates
from binomial.BinomialReplicater import BinomialReplicater
BinoRepNum=2
if __name__ == "__main__":
   # filename = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\ccf.txt'
   # parser = DefaultCCFParser()
    filename = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\copyNumberVariationAdj.txt'
    parser = DefaultTSPParser()	
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
      #  print parser._parse_row()	
      #  ccf_list = parser.get_cancer_cell_fraction_list()
#        data = ccf_list.to_string()
       # print ccf_list
       # ccf_list.save_to_file('C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\ccf2.tsv')
      
        tsp_list = parser.get_tumor_sample_profile_list()
        input_header, A = parser._parse_header(open(filename,'r').readlines()[0].strip())
	
			 
        A = BinomialReplicater(BinoRepNum, tsp_list, input_header)
        A.make_rep(BinoRepNum, tsp_list)	

    
