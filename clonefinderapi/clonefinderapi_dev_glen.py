from parsers.AncestreeTSPParser import AncestreeTSPParser
from alignments.FreqToMegaSeq import FreqToMegaSeq

if __name__ == "__main__":
    filename = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/ancestree_input.txt'
    parser = AncestreeTSPParser()
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
        profile_list = parser.get_tumor_sample_profile_list()
        alignment_builder = FreqToMegaSeq()
        alignment_builder.initialize(profile_list)
        mega_alignment = alignment_builder.get_mega_alignment_string()
        print mega_alignment
        alignment_builder.save_mega_alignment_to_file('/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/ancestree_input.meg')
    
