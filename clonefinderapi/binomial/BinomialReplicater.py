from tsp_profiles.TumorSampleProfileList import TumorSampleProfileList
from tsp_profiles.TumorSampleProfile import TumorSampleProfile
from parsers.DefaultTSPParser import DefaultTSPParser
from tsp_profiles.ReadCount import ReadCount
import numpy

class BinomialReplicater(object):
    """
        Generate TumorSampleProfile replicates
        
        Total read count and SNV frequency is obtained from Original TumorSampleProfile.
        Binomial distribution is used.		
    """
    def __init__(self, tsp_list):
        self.tsp_list = tsp_list	
		

    def make_rep(self):
            self.newtumor_sample_profiles = TumorSampleProfileList()        		 
            for profile in self.tsp_list: 
                tumor = profile._name
                if self.newtumor_sample_profiles.profile_exists(tumor) == False:
                    newprofile = TumorSampleProfile(tumor)
                    self.newtumor_sample_profiles.add(newprofile)
                read_count = ReadCount()								
                for read_count in profile:
                     total_count=read_count.total()
                     alt_freq=read_count.alt_frequency()
                     newalt_count =numpy.random.binomial(total_count, alt_freq, size=None)
                     newref_count = total_count - newalt_count
                     read_count.num_ref = newref_count
                     read_count.num_alt = newalt_count               
                     newprofile.add(read_count)
                     					 
            return self.newtumor_sample_profiles 					 
   