from alignments.MegaAlignment import MegaAlignment
import scipy.optimize
import numpy

class CloneFrequencyComputer_cnv1(object):
    """
        Compute clone frequency (F) from mutant allele fraction matrix (C) and observed variant frequencies (v_obs).
        
        C x F = v_obs	
    """
    def __init__(self, seqs_with_ancestor, v_obs, CNV_info, freq_cutoff):
        self.CutOff = freq_cutoff	
        Align=MegaAlignment()	
        self.ini_clone_order, self.ini_clone_seq = Align.name2seq(seqs_with_ancestor)
        self._CNV_file = CNV_info
        self.v_obs = v_obs		
							  
    def make_Cmatrix(self, SNV_seq):
        Min=''	
        snv_num = len(SNV_seq[self.ini_clone_order[0]])
        snv=0	
        while snv < snv_num:
            for Name in self.ini_clone_order:					
                if Name!='#hg19' and Name!='#Normal':
                   Nuc=SNV_seq[Name][snv]
                   if Nuc=='T': 
                          Min+='0.5 '						  
                   else: 
                         Min+='0 '				   		
            Min=Min[:-1]+'; '
            snv+=1	     		
        Min=Min[:-2]  
        Min=numpy.matrix(Min)	
        return Min	
  
    def regress_cnv(self): 
        Align=MegaAlignment()	
        self.Tumor2Clone_frequency = {}
        HitCloSeq_dic={}
        for tumor in self.v_obs:
             v_obs_single = self.v_obs[tumor]
             v_obs_single_sub = []	
             Seq_dic_sub={}			 
             RmSNVPosi=[]
             CNVls=  self._CNV_file[tumor]
             Len=len(CNVls)
             c=0	
             while c<Len:
                 if CNVls[c]=='normal':
                      v_obs_single_sub.append(v_obs_single[c])
				  
                 else: RmSNVPosi.append(c)
                 c+=1
             for Clo in self.ini_clone_order:
                  NewSeq=''
                  OldSeq=self.ini_clone_seq[Clo]				  
                  c=0	
                  while c<Len:
                      if RmSNVPosi.count(c)==0: NewSeq+=OldSeq[c]
                      c+=1					  
                  Seq_dic_sub[Clo]=NewSeq	
             print(tumor,'exclude bad SNVs for clone frequency computation', RmSNVPosi)				  		 
             Cmatrix_noCNV = self.make_Cmatrix(Seq_dic_sub)				
             Clone2Freq = self.do_nnls0(Cmatrix_noCNV, v_obs_single_sub)
             self.Tumor2Clone_frequency['T-'+tumor]=Clone2Freq
             for Clo in Clone2Freq:
                  if Clone2Freq[Clo]>0:
                            HitCloSeq_dic['#'+Clo]= self.ini_clone_seq['#'+Clo]					   
        self.hitclone_seq_builder = Align.UpMeg(HitCloSeq_dic, [])       			

    def do_nnls0(self,Cmatrix, v_obs_ls):	
            Clone2Freq={}	   	
            clone_frequency = scipy.optimize.nnls(Cmatrix,v_obs_ls)[0] 
            clone_id = 0
            for clone in self.ini_clone_order:		
                if 	clone_frequency[clone_id] > self.CutOff	: Fre = clone_frequency[clone_id]
                else: Fre = 0				
                Clone2Freq[clone[1:]] = Fre
                clone_id += 1
            return Clone2Freq
