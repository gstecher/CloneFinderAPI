from tsp_profiles.TumorSampleProfileList import TumorSampleProfileList
from tsp_profiles.TumorSampleProfile import TumorSampleProfile
from parsers.DefaultTSPParser import DefaultTSPParser
import numpy

class BinomialReplicater(object):
    """
        Generate TumorSampleProfile replicates
        
        Total read count and SNV frequency is obtained from Original TumorSampleProfile.
        Binomial distribution is used.		
    """
    def __init__(self, rep_number, tsp_list, input_order):
        self.RepNum = rep_number 
        self.tsp_list = tsp_list
        self.input_order = input_order		
		

    def make_rep(self, tsp_list, input_order):
        RepID=1
        input_header = ''
       # print self.input_order		
        for input in self.input_order:		
            input_header += input+'\t'
        input_header = input_header[:-1]+'\n'			
		
        while RepID<=self.RepNum:
            newinput={}
            out=input_header			
         #   print out
            for profile in self.tsp_list: 
                tumor = vars(profile)['_name']	
             #   print tumor				
                ref_list=[]
                alt_list=[]				
                for read_count in profile:
                     alt_count = vars(read_count)['_num_alt']
                     ref_count = vars(read_count)['_num_ref']
                     total_count=alt_count + ref_count
                     alt_freq=1.0*alt_count/total_count	
                     newalt_count =numpy.random.binomial(total_count, alt_freq, size=None)
                     newref_count = total_count - newalt_count
                     alt_list.append(str(newalt_count))
                     ref_list.append(str(newref_count))
                newinput[tumor+':ref']=ref_list	
                newinput[tumor+':alt']=alt_list	
            snv_num = len(ref_list)
            snv_id = 0
            while snv_id <snv_num:
                for input in self.input_order:
                    out += newinput[input][snv_id]+'\t'
                out = out[:-1]+'\n'
                snv_id += 1	

            destination = open('BinoRep'+str(RepID)+'.txt','w')
            destination.write(out)
            destination.close()   			
            RepID+=1

     	
