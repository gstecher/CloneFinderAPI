from alignments.MegaAlignment import MegaAlignment
from tsp_profiles.tsp_information import tsp_information
#from regression.CloneFrequencyComputer import CloneFrequencyComputer
#from regression.CloneFrequencyComputer_cnv import CloneFrequencyComputer_cnv
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
from output.CloneFrequencyWriter import CloneFrequencyWriter
from scipy import stats
import numpy as np
import scipy
import numpy
from sklearn.cluster import KMeans

class Calculation:
    def compute_estimated_SNVfrequency1(self, CloFre, CloSeq):
        TMP={}
        c=0
     	
        snv_num=len(CloSeq[CloSeq.keys()[0]])		
        while c<snv_num:
            TMP[c]=0
            c+=1
        for C in CloFre:
          if CloFre[C]>0:		
            S=self.clone_seq['#'+C]
	
      	
            F=CloFre[C]
            c=0
            while c<snv_num:
              				
                if S[c]=='T': TMP[c]+=F*0.5
                c+=1
        return TMP	
	
    def GetClu2Seq_cnv(self,Tu,NormalPosi): 
        Clu2Seq = {}
        Clu2EstSNVfre = {}	

        ID=1
        for Exp in self.Exp2CluSNV:
        #   print Exp		
           CluSNVs=self.Exp2CluSNV[Exp]  
           NegPos='Neg'
           Tot=0		   
           for CluSNV in CluSNVs:		   
               Tot+=len(CluSNV)		   
           for CluSNV in CluSNVs:
            TuClu=Tu+'Clu'+str(ID)
            Seq=''
            c=0
            CluMutC=0		
            while c<self.snv_num:
              if NormalPosi.count(c)!=0:			
                SNV='S'+str(c)		
                Code=SNV in CluSNV
                if Code==True: 
                    Seq+='T'
                    CluMutC+=1					
                else: Seq+='A'
              c+=1
         #   print CluMutC,self.MinSNVnum				
            if CluMutC>=(0.05*self.snv_num):			
              Clu2Seq[TuClu]=Seq      
           #   print 'In'
              Clu2EstSNVfre['Clu'+str(ID)] = str(Exp)+'-'+NegPos	  
              NegPos='Pos'		
              ID+=1
	
        return Clu2Seq, Clu2EstSNVfre		

    def compute_diff1(self, ObsSNV, EstSNV):		
        TMP={}
        c=0
        snv_num=len(ObsSNV)		
        while c<snv_num:	
            Ob=ObsSNV[c]
            Ex=EstSNV[c]		  
            diff=(Ob-Ex)
            Code=Ex in TMP
            if Code!=True: TMP[Ex]={}
            Code=diff in TMP[Ex]
            if Code!=True: TMP[Ex][diff]={}
            TMP[Ex][diff]['S'+str(c)]=Ob
            c+=1			
        return TMP

	
class SNPClusterGenerator_cnv1(Calculation):
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, ini_seq_builder, tsp_list, clone_frequencies, CNV, freq_cutoff ):
        self.tsp_list = tsp_list
        self.freq_cutoff = freq_cutoff			
        self.Tu2CloFre = clone_frequencies
        self.all_tsp = tsp_information(tsp_list)
        self.CloFreCutOff = self.freq_cutoff        				
        self.v_obs = self.all_tsp.tumor2alt_frequency()
        Align = MegaAlignment()
        self.clone_order, self.clone_seq = Align.name2seq(ini_seq_builder)
        self._CNV_file =CNV		

    def GetExp2CluSNV(self): 
        Exp2CluSNV={}
        for Exp in self.Exp2Diff2ObsSNVID:
            Diff2ObsSNVID=self.Exp2Diff2ObsSNVID[Exp]
            if Exp>0:     
                Neg=-99999 #largest negative difference
                Pos=99999 #smallest positive difference
                NegSNV=[]
                PosSNV=[]
                ZeroSNV=[]
                ObsNegLs=[]
                ObsPosLs=[]	
                ObsZeroLs=[]			
                for Diff in Diff2ObsSNVID:
                    ObsSNVID=Diff2ObsSNVID[Diff]
                    for SNV in ObsSNVID: 
                      if Diff<0:
                         ObsNegLs.append(Diff+Exp)				  
                         if Diff>Neg: Neg=Diff
                         NegSNV.append(SNV)
                      elif Diff>0:
                         ObsPosLs.append(Diff+Exp)				  
                         if Diff<Pos: Pos=Diff
                         PosSNV.append(SNV)
                      elif Diff==0: 
                         ZeroSNV.append(SNV)
                         ObsZeroLs.append(Diff+Exp)					 
                if Neg!=-99999 and Pos!=99999:
                    if ObsZeroLs!=[]:
                         PosSq=Pos*Pos
                         NegSq=Neg*Neg
                         if PosSq<NegSq: 
                             ObsPosLs+=	ObsZeroLs
                             PosSNV+=ZeroSNV						 
                         else:
                             ObsNegLs+=ObsZeroLs
                             NegSNV+=ZeroSNV						 
                    if len(ObsNegLs)>=1: NegMed=numpy.mean(ObsNegLs)
                    else: NegMed=0.0
                    if len(ObsPosLs)>=1: PosMed=numpy.mean(ObsPosLs)
                    else: PosMed=0.0
                    Dif=abs(NegMed-PosMed)			
                    pval1=1
                    pval2=1				
                    if len(ObsNegLs)>=3:
                          NegSE=NegMed*(1-NegMed)/len(ObsNegLs)
                          SE=(NegSE)**0.5
                          if SE==0: pval1=0	
                          else:					  
                            T=Dif/SE
                            pval1 = stats.t.sf(abs(T), len(ObsNegLs)-1)  					  
                    if len(ObsPosLs)>=3:
                          PosSE=PosMed*(1-PosMed)/len(ObsPosLs)
                          SE=(PosSE)**0.5
                          if SE==0: pval1=0	
                          else:							  
                            T=Dif/SE
                            pval2 = stats.t.sf(T, len(ObsPosLs)-1)  
                    if pval1<0.05 or pval2<0.05: Clusters=[NegSNV,PosSNV]
                    else: Clusters=[NegSNV+PosSNV,[]]
                else: Clusters=[NegSNV+PosSNV,[]]				
                Exp2CluSNV[Exp]=Clusters	
        return Exp2CluSNV
		
    def add_cluster_Cmat(self, seq_list):
        Align = MegaAlignment()
        inCloLs, inClo2Seq = Align.name2seq(seq_list)
        Cmat_dic = {}
        Cmat_mat = ''
        for inClo1 in inCloLs:
                 inClo = inClo1[1:]		
                 Seq = inClo2Seq['#'+inClo]
                 Len = len(Seq)
                 c=0
                 C_val = []				 
                 while c<Len:
                     if Seq[c]=='T': C_val.append(0.5)
                     else: C_val.append(0)
                     c+=1	
                 Cmat_dic[inClo]=C_val				 
        Cmat_mat, Cmat_clone_order = self.convert_Cmatdic_to_mat(Cmat_dic)
        return Cmat_mat, Cmat_dic, Cmat_clone_order	
		
    def convert_Cmatdic_to_mat(self, Cmat_dic):
        Cin=''    
        snv_num = len(Cmat_dic[Cmat_dic.keys()[0]])
        clone_order0 = Cmat_dic.keys()
        clone_order=[]		
        for clone in clone_order0:
              if clone!='hg19': clone_order.append('#'+clone)			
        snv=0	
        while snv < snv_num:
            for Name in clone_order0:
                        Cval = Cmat_dic[Name][snv]					   
                        Cin+=str(Cval)+' '								
            Cin=Cin[:-1]+'; '
            snv+=1	     		
        Cin=Cin[:-2]   
        Cin=numpy.matrix(Cin)	
        return Cin, clone_order	
		
    def cluster_cnv(self):
        self.tumor2clusters = {}
        self.snv_num = len(self.clone_seq[self.clone_order[0]])	
        self.tumor2estSNV={}
        for tumor in self.v_obs:
           # print tumor		
            self.ObsSNVLs_all = self.v_obs[tumor]
            CloFre=	self.Tu2CloFre['T-'+tumor]
            if 	self.Tu2CloFre.has_key('T-'+tumor)==True:		  
                self.clone2frequency = self.Tu2CloFre['T-'+tumor]
            else:
                self.clone2frequency = self.Tu2CloFre[tumor]
            clone_seq_dic_sub={}

            ObsSNV_sub=[]			
            CNVlist = self._CNV_file[tumor]
            NoCNVPosi=[]
            c=0
            Len=len(CNVlist)			
            while c<Len:
                   if CNVlist[c]=='normal' : NoCNVPosi.append(c) #
                   c+=1            
            for Posi in NoCNVPosi:
                ObsSNV_sub.append(self.ObsSNVLs_all[Posi])	
                for Clone in self.clone_seq :
                     if clone_seq_dic_sub.has_key(Clone)!=True: clone_seq_dic_sub[Clone]=''
                     clone_seq_dic_sub[Clone]+=self.clone_seq[Clone][Posi]
         #   print clone_seq_dic_sub.keys()					 
            Site2EstSNV_sub = Calculation.compute_estimated_SNVfrequency1(self,CloFre, clone_seq_dic_sub)
					
            self.Exp2Diff2ObsSNVID = Calculation.compute_diff1(self, ObsSNV_sub, Site2EstSNV_sub)		

            self.Exp2CluSNV = self.GetExp2CluSNV() 	##get_candidate_decomposed_clones change
          #  print self.Exp2CluSNV			
            Clu2Seq, Clu2EstSNVfre=Calculation.GetClu2Seq_cnv(self, tumor, NoCNVPosi) #
          #  print Clu2EstSNVfre,#	Clu2Seq 
            if Clu2Seq!={}: 
                seq_list_sub=[]
                for HitClo in self.clone2frequency:
                    if self.clone2frequency[HitClo]	> 0:			
                        seq_list_sub+=['#'+HitClo, clone_seq_dic_sub['#'+HitClo]]
                for Clu in Clu2Seq:
                    seq_list_sub+=['#'+Clu, Clu2Seq[Clu]]
                new_Cmat_mat, new_Cmat_dic, new_clone_list = 	self.add_cluster_Cmat(seq_list_sub)	
                clone_frequency_cnv = CloneFrequencyComputer_cnv1({}, {}, {}, self.freq_cutoff, {})
	
              #  print new_clone_list			
                hitclone_frequency_dic = clone_frequency_cnv.do_nnls0(new_Cmat_mat, new_clone_list, ObsSNV_sub)				
                self.tumor2clusters[tumor] = [Clu2EstSNVfre, hitclone_frequency_dic, seq_list_sub]
             #   print hitclone_frequency_dic
            else: self.tumor2clusters[tumor]=[] 
      #  print 	self.tumor2clusters		
        return self.tumor2clusters				
		
			

