from tsp_profiles.TumorSampleProfileList import TumorSampleProfileList
from tsp_profiles.TumorSampleProfile import TumorSampleProfile
from parsers.DefaultTSPParser import DefaultTSPParser
from binomial.BinomialReplicater import BinomialReplicater
from regression.CloneFrequencyComputer import CloneFrequencyComputer

class SignificanceTester():
    """
        Compute clone frequency (F) for each replicates and find significant clones.
        
        List clones that are significant within a tumor. Detect clones for each tumor.	
    """
	
    def __init__(self, rep_number, significance_cutoff, clone_frequency_cutoff, tsp_list, sequences, original_clone_frequency):
        self.RepNum = rep_number
        self.Sig = significance_cutoff		
        self.CutOff = clone_frequency_cutoff		
        self.tsp_list = tsp_list
        self.sequences = sequences		
        self.OriCloFre = original_clone_frequency
	
    def count_detection(self):
        RepID=1
        self.hit_count={}		
        while RepID<=self.RepNum: 
            replicate = BinomialReplicater(self.tsp_list)
            tsp_rep_list = replicate.make_rep()
            clone_frequency = CloneFrequencyComputer(self.sequences, tsp_rep_list, self.CutOff)
            hitseq, clone_frequency_table = clone_frequency.regress()			
            for tumor in clone_frequency_table:
                if RepID==1: self.hit_count[tumor] = {}
                clone_frequencies = clone_frequency_table[tumor]
                for clone in clone_frequencies:
                    if self.hit_count[tumor].has_key(clone)!=True:	self.hit_count[tumor][clone] = 0				
                    if clone_frequencies[clone]>0: 
                          self.hit_count[tumor][clone]+=1.0/self.RepNum								
            RepID += 1			
        self.hit_count1={}
        self.hit_clone_list=[]		
        for tumor in self.hit_count:
            self.hit_count1[tumor]={}
            CloneHit=0
            hit_count=self.hit_count[tumor]
            for clone in hit_count:
                if hit_count[clone]>self.Sig: 
                     CloneHit+=1
                     self.hit_count1[tumor][clone]=1
                     self.hit_clone_list.append(clone)					 
                else: 
                    self.hit_count1[tumor][clone]=0
                    print 'binomial is insignificant',tumor,clone,hit_count[clone]					
            if CloneHit==0:
                LargeCloneFre = self.find_large_clone_frequency(tumor)
                for clone in LargeCloneFre:	
                    self.hit_count1[tumor][clone]=1				
                    self.hit_clone_list.append(clone)
        self.hit_clone_list=list(set(self.hit_clone_list))					
        return self.hit_count1, self.hit_clone_list	
		
    def find_large_clone_frequency(self, tumor):
        	
        if self.OriCloFre.has_key(tumor)==True: OriClo2Fre=self.OriCloFre[tumor]
        else: 	OriClo2Fre=self.OriCloFre[tumor[2:]]	
        Lar=0
        LarClo=[]
        for Clo in OriClo2Fre:
            Fre=OriClo2Fre[Clo]
            if Fre>Lar:
                Lar=Fre
                LarClo=[Clo]
            elif Fre==Lar: LarClo.append(Clo)
        return LarClo	
		