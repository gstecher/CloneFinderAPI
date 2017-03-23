from alignments.MegaAlignment import MegaAlignment
from tsp_profiles.tsp_information import tsp_information
from regression.CloneFrequencyComputer import CloneFrequencyComputer
from regression.CloneFrequencyComputer_cnv import CloneFrequencyComputer_cnv
from scipy import stats
import numpy

class Calculation:

    def compute_estimated_SNVfrequency(self):
        TMP={}
        c=0
        while c<self.snv_num:
            TMP[c]=0
            c+=1
        for C in self.clone2frequency:
            S=self.clone_seq['#'+C]
            if str(self.clone2frequency[C]).find('e')!=-1: F=0		
            else: F=self.clone2frequency[C]/2
            c=0
            while c<self.snv_num:
                if S[c]=='T': TMP[c]+=F
                c+=1
        return TMP	
		
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
	
    def GetClu2Seq(self,Tu): 
        Clu2Seq = {}
        Clu2EstSNVfre = {}	

        ID=1
        for Exp in self.Exp2CluSNV:
           CluSNVs=self.Exp2CluSNV[Exp]  
           NegPos='Neg'
           for CluSNV in CluSNVs:
            TuClu=Tu+'Clu'+str(ID)
            Seq=''
            c=0
            CluMutC=0		
            while c<self.snv_num:
                SNV='S'+str(c)		
                Code=SNV in CluSNV
                if Code==True: 
                    Seq+='T'
                    CluMutC+=1					
                else: Seq+='A'
                c+=1
            if CluMutC>=self.MinSNVnum:			
              Clu2Seq[TuClu]=Seq      

              Clu2EstSNVfre['Clu'+str(ID)] = str(Exp)+'-'+NegPos	  
              NegPos='Pos'		
              ID+=1
	
        return Clu2Seq, Clu2EstSNVfre

    def compute_diff(self):		
        TMP={}
        c=0
        while c<self.snv_num:	
            Ob=self.ObsSNVLs[c]
            Ex=self.Site2EstSNV[c]		  
            diff=(Ob-Ex)
            Code=Ex in TMP
            if Code!=True: TMP[Ex]={}
            Code=diff in TMP[Ex]
            if Code!=True: TMP[Ex][diff]={}
            TMP[Ex][diff]['S'+str(c)]=Ob
            c+=1			
        return TMP	
		
class SNPClusterGenerator_cnv(Calculation):
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, MEGAalignment, tsp_list, clone_frequencies, CNV_info, option_b, freq_cutoff):
        self.tsp_list = tsp_list
        self.CNV_info = CNV_info
        self.option_b = option_b
        self.freq_cutoff = freq_cutoff			
        self.Tu2CloFre = clone_frequencies
        self.MinSNVnum = int(option_b)
        self.all_tsp = tsp_information(tsp_list)
        self.CloFreCutOff = self.freq_cutoff        				
        self.v_obs = self.all_tsp.tumor2alt_frequency()
        Align = MegaAlignment()
        self.clone_order, self.clone_seq = Align.name2seq(MEGAalignment)
        #self.CNV_information= CNV_information
    #    print 'tumor sncestor clones wihtout CNV adjustment clone frequencies',self.Tu2CloFre		
		
    def cluster(self):
        self.tumor2clusters = {}
        self.snv_num = len(self.clone_seq[self.clone_order[0]])	
        self.tumor2estSNV={}
        for tumor in self.v_obs:
            self.ObsSNVLs = self.v_obs[tumor]
	
            if 	self.Tu2CloFre.has_key('T-'+tumor)==True:		  
                self.clone2frequency = self.Tu2CloFre['T-'+tumor]
            else:
                self.clone2frequency = self.Tu2CloFre[tumor]			
            self.Site2EstSNV = Calculation.compute_estimated_SNVfrequency(self)
            self.tumor2estSNV[tumor]=self.Site2EstSNV			
            self.Exp2Diff2ObsSNVID = Calculation.compute_diff(self)			
            self.Exp2CluSNV = Calculation.GetExp2CluSNV(self) 
            Clu2Seq, Clu2EstSNVfre=Calculation.GetClu2Seq(self, tumor) #
	
            if Clu2Seq!={}: 
                seq_list=[]
                for HitClo in self.clone2frequency:
                    if self.clone2frequency[HitClo]	> 0:			
                        seq_list+=['#'+HitClo, self.clone_seq['#'+HitClo]]
                for Clu in Clu2Seq:
                    seq_list+=['#'+Clu, Clu2Seq[Clu]]
                tsp_single = self.all_tsp.make_single_tsp_list(tumor)					
                clone_frequency = CloneFrequencyComputer_cnv(seq_list, tsp_single, self.CNV_info, self.option_b, self.freq_cutoff)
                hitclone_sequences, hitclone_frequency_dic = clone_frequency.snvGenotype_after_regress()
                self.tumor2clusters[tumor] = [Clu2EstSNVfre, hitclone_frequency_dic, hitclone_sequences]

            else: self.tumor2clusters[tumor]=[] 
       # print 	self.tumor2clusters		
        return self.tumor2clusters				

    def estimatedSNV_save_to_file(self):
        self.cluster()	
        out=''
        c=0
        TuOrder=[]
        out=''		
        for tumor in self.tumor2estSNV:
            TuOrder.append(tumor)
            out+=tumor+':obs\t'+tumor+':est\t'
        out=out[:-1]+'\n'			
        while c<self.snv_num:
           for tumor in TuOrder:		
                out+=str(self.v_obs[tumor][c])+'\t'+str(self.tumor2estSNV[tumor][c])+'\t'
           out=out[:-1]+'\n'
           c+=1
        return self.v_obs, self.tumor2estSNV 		   
      #  OutF=open('EstObsSNVfreqTuAnc.txt','w')
       # OutF.write(out)
        #OutF.close()		
           		   
