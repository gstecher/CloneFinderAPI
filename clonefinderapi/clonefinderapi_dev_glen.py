from parsers.DefaultCCFParser import DefaultCCFParser
from alignments.FreqToMegaSeq import FreqToMegaSeq

if __name__ == "__main__":
    filename = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/ccf.txt'
    parser = DefaultCCFParser()
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
        ccf_list = parser.get_cancer_cell_fraction_list()
#        data = ccf_list.to_string()
        print ccf_list
        ccf_list.save_to_file('/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/ccf2.tsv')
    
