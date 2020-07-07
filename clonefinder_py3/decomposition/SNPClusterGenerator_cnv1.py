from alignments.MegaAlignment import MegaAlignment
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
from scipy import stats
import numpy

class Calculation:
    def compute_estimated_SNVfrequency1(self, CloFre, CloSeq): 
        eSNVfre={}
        c=0     	
        snv_num=len(CloSeq[list(CloSeq.keys())[0]])		
        while c<snv_num:
            eSNVfre[c]=0
            c+=1
        for C in CloFre:
          if CloFre[C]>0:		
            S=self.clone_seq['#'+C]     	
            F=CloFre[C]
            c=0
            while c<snv_num:              				
                if S[c]=='T': eSNVfre[c]+=F*0.5
                c+=1
        return eSNVfre	

    def GetClu2Seq_cnv1(self,CluSNV): 
            Seq=''
            c=0	
            while c<self.snv_num:
                if CluSNV.count('S'+str(c))!=0: 
                    Seq+='T'				
                else: Seq+='A'
                c+=1				
            return Seq			

    def compute_diff1(self, ObsSNV, EstSNV):	
        Dif={}
        c=0
        snv_num=len(ObsSNV)		
        while c<snv_num:	
            Ob=ObsSNV[c]
            Ex=EstSNV[c]		  
            diff=(Ob-Ex)
          
            if (Ex in Dif)!=True: Dif[Ex]={}
          
            if (diff in Dif[Ex])!=True: Dif[Ex][diff]={}
            Dif[Ex][diff]['S'+str(c)]=Ob
            c+=1			
        return Dif
	
class SNPClusterGenerator_cnv1(Calculation):

    def __init__(self, ini_seq_builder, v_obs, clone_frequencies, CNV, freq_cutoff ):
        self.freq_cutoff = freq_cutoff			
        self.Tu2CloFre = clone_frequencies
        self.CloFreCutOff = self.freq_cutoff        				
        self.v_obs = v_obs
        Align = MegaAlignment()
        self.clone_order, self.clone_seq = Align.name2seq(ini_seq_builder)
        self._CNV_file =CNV		
        self.snv_num = len(self.clone_seq[self.clone_order[0]])	

    def testClu(self,NegMed,ObsNegLs,PosMed,ObsPosLs):
                          Dif=abs(NegMed-PosMed)			
                          pval1=1
                          pval2=1				
 
                          NegSE=NegMed*(1-NegMed)/len(ObsNegLs)
                          SE=(NegSE)**0.5
                          if SE==0: pval1=0	
                          else:					  
                            T=Dif/SE
                            pval1 = stats.t.sf(abs(T), len(ObsNegLs)-1)  					  
                    
                          PosSE=PosMed*(1-PosMed)/len(ObsPosLs)
                          SE=(PosSE)**0.5
                          if SE==0: pval1=0	
                          else:							  
                            T=Dif/SE
                            pval2 = stats.t.sf(T, len(ObsPosLs)-1)  
                          return pval1,pval2	
						  
    def GetExp2CluSNV(self,Exp2Diff2ObsSNVID,Tu,NormalPosi): 
        CluID2CluSNV={}		
        ExpOrder=list(Exp2Diff2ObsSNVID.keys())
        ExpOrder.sort(reverse=True)		    	
        CluID=0 #start from the root		
        for Exp in ExpOrder:
            Diff2ObsSNVID=Exp2Diff2ObsSNVID[Exp]
            if Exp>0:     
            
                NegSNV=[]
                PosSNV=[]
                ZeroSNV=[]
                ObsNegLs=[]
                ObsPosLs=[]	
                ObsZeroLs=[]			
                for Diff in Diff2ObsSNVID:
                    ObsSNVID=Diff2ObsSNVID[Diff]
                    for SNV in ObsSNVID: 
                      if NormalPosi.count(int(SNV.replace('S','')))!=0 and Diff<0:
                         ObsNegLs.append(Diff+Exp)				  
                      
                         NegSNV.append(SNV)
                      elif NormalPosi.count(int(SNV.replace('S','')))!=0 and Diff>0:
                         ObsPosLs.append(Diff+Exp)				  
                        
                         PosSNV.append(SNV)
                      elif NormalPosi.count(int(SNV.replace('S','')))!=0 and Diff==0: 
                         ZeroSNV.append(SNV)
                         ObsZeroLs.append(Diff+Exp)					 
                if len(NegSNV)>=3 and len(PosSNV)>=3:
                    NegMed=numpy.mean(ObsNegLs)
                  
                    PosMed=numpy.mean(ObsPosLs)
                    			
                    if ObsZeroLs!=[]:
                         PosSq=PosMed*PosMed
                         NegSq=NegMed*NegMed
                         if PosSq<NegSq: 
                             ObsPosLs+=	ObsZeroLs
                             PosSNV+=ZeroSNV						 
                         else:
                             ObsNegLs+=ObsZeroLs
                             NegSNV+=ZeroSNV						 

                    pval1,pval2=self.testClu(NegMed,ObsNegLs,PosMed,ObsPosLs)							
                    if pval1<0.05 or pval2<0.05: 
                          CluID2CluSNV['#Clu'+str(CluID)]=Calculation.GetClu2Seq_cnv1(self,PosSNV) 
                          CluID2CluSNV['#Clu'+str(CluID+0.5)]=Calculation.GetClu2Seq_cnv1(self,NegSNV)						  
      
                    else: 
                          CluID2CluSNV['#Clu'+str(CluID)]=Calculation.GetClu2Seq_cnv1(self,PosSNV+NegSNV)					
               
                elif len(PosSNV+NegSNV)>=3: CluID2CluSNV['#Clu'+str(CluID)]=Calculation.GetClu2Seq_cnv1(self,PosSNV+NegSNV)	#Clusters=[NegSNV+PosSNV,[]]				
           
                CluID+=1				
	
        return CluID2CluSNV				
    def cluster_cnv(self):
        self.tumor2clusters = {}
        Align=MegaAlignment()
        for tumor in self.v_obs:		
            self.ObsSNVLs_all = self.v_obs[tumor]	  
            clone2frequency = self.Tu2CloFre['T-'+tumor]
            clone_seq_dic_sub={}
            ObsSNV_sub=[]			
            CNVlist = self._CNV_file[tumor]
            NoCNVPosi=[]
            c=0			
            while c<self.snv_num:
                   if CNVlist[c]=='normal' : 
                       NoCNVPosi.append(c) #
                       ObsSNV_sub.append(self.ObsSNVLs_all[c])	
                       for Clone in self.clone_seq :
                             if (Clone in clone_seq_dic_sub)!=True: clone_seq_dic_sub[Clone]=''
                             clone_seq_dic_sub[Clone]+=self.clone_seq[Clone][c]					   
                   c+=1            				 
            Site2EstSNV_sub = Calculation.compute_estimated_SNVfrequency1(self,clone2frequency, clone_seq_dic_sub)					
            Exp2Diff2ObsSNVID = Calculation.compute_diff1(self, ObsSNV_sub, Site2EstSNV_sub)		
            Clu2Seq = self.GetExp2CluSNV(Exp2Diff2ObsSNVID,tumor, NoCNVPosi) 	##get_candidate_decomposed_clones change		
        
            if Clu2Seq!={}: 
                seq_list_sub=[]			
                for HitClo in clone2frequency:
                    if clone2frequency[HitClo]	> 0:			
                        seq_list_sub+=['#'+HitClo, clone_seq_dic_sub['#'+HitClo]]					
                for Clu in Clu2Seq:
                    seq_list_sub+=[Clu, Clu2Seq[Clu]]		
                clone_frequency_clu = CloneFrequencyComputer_cnv1(seq_list_sub, {tumor:ObsSNV_sub}, {tumor:CNVlist}, self.freq_cutoff)
                clone_frequency_clu.regress_cnv()					
				
                self.tumor2clusters[tumor] = [clone_frequency_clu.hitclone_seq_builder, clone_frequency_clu.Tumor2Clone_frequency]
            else: self.tumor2clusters[tumor]=[] 
	
        return self.tumor2clusters				
		
			

