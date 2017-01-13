from alignments.MegaAlignment import MegaAlignment
from tsp_profiles.tsp_information import tsp_information
from regression.CloneFrequencyComputer import CloneFrequencyComputer
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
  #  outEst='CluID\tEst SNVfre\n'
  #  OutEst=Tu+'_EstSNVFreqClu.txt'	
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
         # outEst+='Clu'+str(ID)+'\t'+str(Exp)+'-'+NegPos+'\n'
              Clu2EstSNVfre['Clu'+str(ID)] = str(Exp)+'-'+NegPos	  
              NegPos='Pos'		
              ID+=1
  #  Functions.GetOut(OutEst,outEst)		
        return Clu2Seq, Clu2EstSNVfre
##
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
   #     print TMP,c			
        return TMP			
class SNPClusterGenerator(Calculation):
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, MEGAalignment, tsp_list, clone_frequencies, minimum_num_SNV_per_cluster, clone_frequency_cutoff):
        #self.M = MEGAalignment 
        self.tsp_list = tsp_list
        self.Tu2CloFre = clone_frequencies
        self.MinSNVnum = int(minimum_num_SNV_per_cluster)
        self.all_tsp = tsp_information(tsp_list)
        self.CloFreCutOff = clone_frequency_cutoff        		
############		
        self.v_obs = self.all_tsp.tumor2alt_frequency()
  #      for profile in self.tsp_list: 
   #         tumor = profile.name			
    #        v_list=[]				
     #       for read_count in profile:
      #           v_list.append(read_count.alt_frequency())
       #     self.v_obs[tumor]=v_list
#############
        Align = MegaAlignment()
        self.clone_order, self.clone_seq = Align.name2seq(MEGAalignment)		
   #     self.clone_order=[]
    #    self.clone_seq = {}
     #   Name=''		
      #  for Seq in self.M:
       #     if Seq[0]=='#' and Seq!='#MEGA':
        #        if Seq!='#hg19' and Seq!='#Normal':
         #           self.clone_order.append(Seq)
          #          self.out += '\t' + Seq[1:]					
           #         self.clone_seq[Seq]=''
            #        Name=Seq					
      #      elif Name!='':
       #         self.clone_seq[Name] += Seq	
#############
  #  def make_single_tsp_list(self, target_tsp):
   #         single_sample_profiles = TumorSampleProfileList()        		 
    #        for profile in self.tsp_list: 
     #           tumor = profile._name
      #          if tumor == target_tsp:				
       #             if single_sample_profiles.profile_exists(tumor) == False:
        #                newprofile = TumorSampleProfile(tumor)
         #               single_sample_profiles.add(newprofile)
          #          #read_count = ReadCount()								
           #         for read_count in profile:             
            #             newprofile.add(read_count)
                     					 
            #return single_sample_profiles 
#######			



		
    def cluster(self):
        self.tumor2clusters = {}
        self.snv_num = len(self.clone_seq[self.clone_order[0]])	

        for tumor in self.v_obs:
            self.ObsSNVLs = self.v_obs[tumor]
          #  print 'h', self.Tu2CloFre, tumor			
            self.clone2frequency = self.Tu2CloFre['T-'+tumor]
            self.Site2EstSNV = Calculation.compute_estimated_SNVfrequency(self)
            self.Exp2Diff2ObsSNVID = Calculation.compute_diff(self)			
            self.Exp2CluSNV = Calculation.GetExp2CluSNV(self) 
            Clu2Seq, Clu2EstSNVfre=Calculation.GetClu2Seq(self, tumor) #
          #  print tumor, Clu2Seq,Clu2EstSNVfre	,self.Exp2CluSNV, self.Exp2Diff2ObsSNVID		
            if Clu2Seq!={}: 
                seq_list=[]
                for HitClo in self.clone2frequency:
                    if self.clone2frequency[HitClo]	> 0:			
                        seq_list+=['#'+HitClo, self.clone_seq['#'+HitClo]]
                for Clu in Clu2Seq:
                    seq_list+=['#'+Clu, Clu2Seq[Clu]]
                tsp_single = self.all_tsp.make_single_tsp_list(tumor)					
                clone_frequency = CloneFrequencyComputer(seq_list, tsp_single, self.CloFreCutOff)
                hitclone_sequences, hitclone_frequency_dic = clone_frequency.regress()
                self.tumor2clusters[tumor] = [Clu2EstSNVfre, hitclone_frequency_dic, hitclone_sequences]
               # print self.tumor2clusters
            else: self.tumor2clusters[tumor]=[] 			   
        return self.tumor2clusters				
        #        Functions.GetOut(Tu+'AddClu.meg',outMeg)
         #       Functions.ExtObs(Turefmut2ReadC,Tu,Len)       		
          #      os.system('python '+'DoSLEPNoBoo.py '+Tu+'.txt '+Tu+'AddClu.meg '+SLEP_dir+' MakeClusterPy')

