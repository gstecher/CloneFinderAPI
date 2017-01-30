from tsp_profiles.TumorSampleProfileList import TumorSampleProfileList
from tsp_profiles.TumorSampleProfile import TumorSampleProfile



class tsp_information():
    """
        Generate TumorSampleProfile replicates
        
        Total read count and SNV frequency is obtained from Original TumorSampleProfile.
        Binomial distribution is used.		
    """
    def __init__(self, tsp_list):
        self.tsp_list = tsp_list

    def make_single_tsp_list(self, target_tsp):
            single_sample_profiles = TumorSampleProfileList()        		 
            for profile in self.tsp_list: 
                tumor = profile._name
                if tumor == target_tsp:				
                    if single_sample_profiles.profile_exists(tumor) == False:
                        newprofile = TumorSampleProfile(tumor)
                        single_sample_profiles.add(newprofile)
							
                    for read_count in profile:             
                         newprofile.add(read_count)
                     					 
            return single_sample_profiles

    def tumor2alt_frequency(self):
        self.v_obs = {}
        for profile in self.tsp_list: 
            tumor = profile.name			
            v_list=[]				
            for read_count in profile:
                 v_list.append(read_count.alt_frequency())
            self.v_obs[tumor]=v_list
        return self.v_obs			