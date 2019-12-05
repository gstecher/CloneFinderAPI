from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
#from output.CloneFrequencyWriter import CloneFrequencyWriter
from alignments.MegaAlignment import MegaAlignment
#from parsimony.MegaMP import MegaMP
#from parsimony.MegaMPtiming import MegaMPtiming
#from parsimony.TreeAnalizer import TreeAnalizer
#from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
#import scipy.optimize
#import numpy
#import scipy
#import os
#import random
#import itertools

class CloneFrequencyComputer_EachTumor(object):
    """
        Compute clone frequency (F) from mutant allele fraction matrix (C) and observed variant frequencies (v_obs).
        
        C x F = v_obs	
    """
    def __init__(self, seqs_with_ancestor, tsp_list, CNV_info, freq_cutoff, ReadCountTable, IniCloFre):
      self.CutOff = freq_cutoff	
      if seqs_with_ancestor!={}:	
        self.IniCloFre=IniCloFre	  
        Align=MegaAlignment()	
        self.ini_seq_builder = seqs_with_ancestor 
        self.ini_clone_order, self.ini_clone_seq = Align.name2seq(self.ini_seq_builder)
        self.tsp_list = tsp_list        
        
        self._CNV_file = CNV_info
        self.ReadCountTable = ReadCountTable
        self.SNVnum = len(ReadCountTable[ReadCountTable.keys()[0]])
        self.CNVnum = len(CNV_info[CNV_info.keys()[0]])		
		
    def	ComputeHitCloneOnly(self):
        Align=MegaAlignment()	
        NewCloFre={}
        NewCloSeq={}
        print 	'initial clone frequency',self.IniCloFre	
        for Tu in self.IniCloFre:	
             sub_seq_builder = ['#MEGA','!Title SNVs;','!Format datatype=dna;']	
             CloFre=self.IniCloFre[Tu]
             for Clo in CloFre:
                 if CloFre[Clo]>self.CutOff: sub_seq_builder+=['#'+Clo,self.ini_clone_seq['#'+Clo]]
             print Tu,sub_seq_builder				 
             clone_frequency_cnv = CloneFrequencyComputer_cnv1(sub_seq_builder, self.tsp_list, self._CNV_file, self.CutOff, self.ReadCountTable)
             clone_frequency_cnv.regress_cnv()
             NewCloFre[Tu]=	 clone_frequency_cnv.Tumor2Clone_frequency[Tu]	
             Hit_clone_order, Hit_clone_seq = Align.name2seq(clone_frequency_cnv.hitclone_seq_builder)
             print 'new clofre',clone_frequency_cnv.Tumor2Clone_frequency[Tu]			 
             for HitClo in Hit_clone_seq:			 
                 NewCloSeq[HitClo]=Hit_clone_seq[HitClo]			 
        				
        print NewCloFre
        NewSeq_buil = Align.UpMeg(NewCloSeq, NewCloSeq.keys())		
        return NewSeq_buil, NewCloFre
				  
 
	  
 