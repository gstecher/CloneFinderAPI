from tsp_profiles.TumorSampleProfileList import TumorSampleProfileList

class FreqToMegaSeq:
    """
        Generate a MEGA sequence alignment file from a TumorSampleProfile
        
        A MEGA dna sequence alignment will be generated based on presence/absence
        of SNVs in the tumor sample profile. A sequence is generated for each
        tumor sample where a an 'A' represents absence of SNV at a given site
        and a 'T' represents presence of SNV at a given site
    """
    
    def __init__(self):
        self._mega_seqs = []
        
    def initialize(self, tumor_sample_profiles, remove_duplicates = True):
        self._mega_seqs.append('#MEGA')
        self._mega_seqs.append('!Title ' + tumor_sample_profiles.name + ';')
        self._mega_seqs.append('!Format datatype=dna' + ';')
        self._mega_seqs.append(' ')
        self._mega_allseqs = []		
        num_sites = tumor_sample_profiles.num_read_counts()
        hg19 = self._get_hg19_sequence(num_sites)
        self._mega_seqs.append('#hg19')
        self._mega_seqs.append(hg19)
        
        for profile in tumor_sample_profiles:
            name = profile.name
            seqdata = profile.get_alignment_string()
            self._mega_allseqs.append('#' + name)
            self._mega_allseqs.append(seqdata) 			
            if remove_duplicates == True:
                if seqdata in self._mega_seqs:
                    continue
            self._mega_seqs.append('#' + name)
            self._mega_seqs.append(seqdata)                            
 
			
    def _get_hg19_sequence(self, num_sites):
        result = ''
        if num_sites > 0:
            index = 1
            while index <= num_sites:
                result = result + 'A'
                index += 1
        return result
        
    def get_mega_alignment(self):
        return self._mega_seqs
    
    def get_mega_alignment_string(self):
        result = ''
        for item in self._mega_seqs:
            result += item + "\n"
        return result
    
    def save_mega_alignment_to_file(self, filename):
        destination = open(filename,'w')
        destination.write(self.get_mega_alignment_string())
        destination.close()        
             
    def get_mega_allalignment(self):
        return self._mega_allseqs        
        

