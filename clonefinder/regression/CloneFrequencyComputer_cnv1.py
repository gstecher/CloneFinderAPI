from output.CloneFrequencyWriter import CloneFrequencyWriter
from alignments.MegaAlignment import MegaAlignment
from parsimony.MegaMP import MegaMP
from parsimony.MegaMPtiming import MegaMPtiming
from parsimony.TreeAnalizer import TreeAnalizer
from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
import scipy.optimize
import numpy
import scipy
import os
import random
import itertools

class CloneFrequencyComputer_cnv1(object):
    """
        Compute clone frequency (F) from mutant allele fraction matrix (C) and observed variant frequencies (v_obs).
        
        C x F = v_obs	
    """
    def __init__(self, seqs_with_ancestor, tsp_list, CNV_info, freq_cutoff, ReadCountTable):
      self.CutOff = freq_cutoff	
      if seqs_with_ancestor!={}:	
        Align=MegaAlignment()	
        self.ini_seq_builder = seqs_with_ancestor 
        self.ini_clone_order, self.ini_clone_seq = Align.name2seq(self.ini_seq_builder)
        self.tsp_list = tsp_list        
        self.make_readcount()
        self._CNV_file = CNV_info
        self.ReadCountTable = ReadCountTable
        self.SNVnum = len(ReadCountTable[ReadCountTable.keys()[0]])
        self.CNVnum = len(CNV_info[CNV_info.keys()[0]])		
		
    def	make_readcount(self):	
        self.v_obs = {}
        self.readCount = {}		
	
        for profile in self.tsp_list: 
            tumor = profile.name			
            v_list=[]	
            ReadCount_list=[]			
            for read_count in profile:
                 v_list.append(read_count.alt_frequency())
                 ReadCount_list.append([read_count.num_ref,read_count.num_alt])				 
            self.v_obs[tumor]=v_list
            self.readCount[tumor]=ReadCount_list				


				  
    def make_mut_wild_allele_count_noCNV(self, PreAbsCNV, clone_order, SNV_seq):
        CNV_file={}	
        SNVnum=len(SNV_seq[clone_order[0]])		
        Align=MegaAlignment()
        TreeAna=TreeAnalizer()		
        for Clo in 	clone_order:
              CNV_file[Clo.replace('#','')+':m']=[]	
              CNV_file[Clo.replace('#','')+':w']=[]				  
        if 	PreAbsCNV == {}:
          c=0
          while c<SNVnum:		  
            for Clo in 	clone_order:
              Nuc=SNV_seq[Clo][c]
              if Nuc=='T':			  
                  CNV_file[Clo.replace('#','')+':m'].append(1)	
                  CNV_file[Clo.replace('#','')+':w'].append(1)
              else:				  
                  CNV_file[Clo.replace('#','')+':m'].append(0)	
                  CNV_file[Clo.replace('#','')+':w'].append(2)
            c+=1
        return CNV_file			
			
  
   
  
		
    def make_Min(self, clone_order, clone_seq, MutWildAlleleCount):	
        Min=''
        Cmat_dic={}	
        for Name in clone_order:
               # if snv==0:print Name, clone_seq[Name]			
                if Name!='#hg19' and Name!='#Normal': Cmat_dic[Name[1:]]=[]		
        snv_num = len(clone_seq[clone_order[0]])
        snv=0
      #  print clone_order		
        while snv < snv_num:
            for Name in clone_order:
               # if snv==0:print Name, clone_seq[Name]			
                if Name!='#hg19' and Name!='#Normal' and MutWildAlleleCount.has_key(Name[1:]+':m')==True:
                       # print Name				
                        CNVm=float(MutWildAlleleCount[Name[1:]+':m'][snv])		
                        CNVw=float(MutWildAlleleCount[Name[1:]+':w'][snv])
                        if CNVm+CNVw==0: value=0	
                        else: value=CNVm/(CNVm+CNVw)
                        Min+=str(value)+' '
                        Cmat_dic[Name[1:]].append(value)						
                elif Name!='#hg19' and Name!='#Normal':
                   Nuc=clone_seq[Name][snv]
                   if Nuc=='T': 
                          Min+='0.5 '
                          Cmat_dic[Name[1:]].append(0.5)							  
                   else: 
                         Min+='0 '				   
                         Cmat_dic[Name[1:]].append(0)				
            Min=Min[:-1]+'; '
            snv+=1	     		
        Min=Min[:-2]
       # print Min
       # print Min	   
        Min=numpy.matrix(Min)
        #print Min		
        return Min, Cmat_dic		
  
 
    def regress_cnv(self):
        Align=MegaAlignment()	
        CloFreAna = CloneFrequencyAnalizer()
        self.Tumor2Clone_frequency = {}
        HitCloSeq_dic={}
        self.tumor2CNVSNVposi={}		
      #  print 'nnls removing SNV-CNVs'
        for tumor in self.v_obs:
           #  print tumor
             			 
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
				  
             self.tumor2CNVSNVposi[tumor]=RmSNVPosi			 
             MutWildAlleleCount_noCNV = self.make_mut_wild_allele_count_noCNV({}, self.ini_clone_order, Seq_dic_sub)#PreAbsCNV, clone_order, SNV_seq, Tu2CloFre		
             Cmatrix_noCNV, Cmatrix_noCNV_dic = self.make_Min(self.ini_clone_order, Seq_dic_sub, MutWildAlleleCount_noCNV)	
             self.Cmatrix_noCNV_mat = Cmatrix_noCNV
             self.Cmatrix_noCNV_dic = Cmatrix_noCNV_dic		
             Clone2Freq = self.do_nnls0(Cmatrix_noCNV, self.ini_clone_order, v_obs_single_sub)
             self.Tumor2Clone_frequency['T-'+tumor]=Clone2Freq
             for Clo in Clone2Freq:
                  if Clone2Freq[Clo]>0:
                       if HitCloSeq_dic.has_key('#'+Clo)!=True:
                            HitCloSeq_dic['#'+Clo]= self.ini_clone_seq['#'+Clo]					   
        self.hitclone_seq_builder = Align.UpMeg(HitCloSeq_dic, [])       	
        CloFreAna.save_frequency_table_to_file('Ini_freq.txt', self.Tumor2Clone_frequency, [])	
        Align.save_mega_alignment_to_file('Ini.meg', self.hitclone_seq_builder)
		
     #   self.Tu2EstSNVFre_noCNV, self.Tu2Est_minus_ObsLs_noCNV = self.EstimateSNVfre(Clone_frequency_file_noCNV, M_noCNV, self.ReadCountTable)	
      #  self.getout('Est.txt',self.Tu2EstSNVFre_noCNV)	
       # self.getout('delta.txt',self.Tu2Est_minus_ObsLs_noCNV)		
 
    def getout(self,Out,out):
        In=''
        Ls=out.keys()
        Len=len(out[Ls[0]])
        for i in Ls:
           In+=i+'\t'
        In=In[:-1]+'\n'
        c=0
        while c<Len:
           for i in Ls:
               In+=str(out[i][c])+'\t'
           In=In[:-1]+'\n'
           c+=1	
        OutF=open(Out,'w')
        OutF.write(In)
        OutF.close()	

  
		

    def do_nnls0(self,Cmatrix, clone_order, v_obs_ls):	#use this	
            Clone2Freq={}	   

	
            clone_frequency = scipy.optimize.nnls(Cmatrix,v_obs_ls)[0] #self.do_ScipyNnls(tumor)#float list	

            clone_id = 0
            for clone in clone_order:
		
                if 	clone_frequency[clone_id] > self.CutOff	: Fre = clone_frequency[clone_id]
                else: Fre = 0				
                Clone2Freq[clone[1:]] = Fre#*2
                clone_id += 1
         #   print 	'clone frequencies\n',clone_order,'\n',clone_frequency	
            return Clone2Freq

    def EstimateSNVfre(self, Tu2CloFre,clone_seq0,ReadCount):
      Align=MegaAlignment()
      cloorder,clone_seq=Align.name2seq(clone_seq0)	

      tumor2estSNV={}
      tumor2diff={}	
    
      for tumor in Tu2CloFre:
        clone2frequency=Tu2CloFre[tumor]	
        tumor=tumor.split('-')[-1]		
        estSNVfreLs=[]
        DiffLs=[]		
        snv_num=len(ReadCount[tumor+':ref'])		


        c=0
        while c<snv_num:
          estSNVfre=0		
          for Clo in clone2frequency:
          	  
            S=clone_seq['#'+Clo]
           		
            if str(clone2frequency[Clo]).find('e')!=-1: F=0		
            else: F=clone2frequency[Clo]/2		
            if S[c]=='T': estSNVfre+=F

          estSNVfreLs.append(estSNVfre)
          Obs=1.0*float(ReadCount[tumor+':alt'][c])/(float(ReadCount[tumor+':alt'][c])+float(ReadCount[tumor+':ref'][c]))	
          Dif=estSNVfre-Obs		
          DiffLs.append(Dif)		  
          c+=1		  
        tumor2estSNV[tumor]=estSNVfreLs
        tumor2diff[tumor]=DiffLs
      return tumor2estSNV,tumor2diff	
	  
 
	  
 