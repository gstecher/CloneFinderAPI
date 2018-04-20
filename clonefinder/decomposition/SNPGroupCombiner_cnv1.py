from alignments.MegaAlignment import MegaAlignment
#from alignments.AlignmentAnnalizer import AlignmentAnnalizer
from tsp_profiles.tsp_information import tsp_information
#from regression.CloneFrequencyComputer import CloneFrequencyComputer
#from regression.CloneFrequencyComputer_cnv import CloneFrequencyComputer_cnv
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
import copy
class SNPGroupCombiner_cnv1():
    """
        Generate a MEGA sequence alignment file from a TumorSampleProfile
        
        A MEGA dna sequence alignment will be generated based on presence/absence
        of SNVs in the tumor sample profile. A sequence is generated for each
        tumor sample where a an 'A' represents absence of SNV at a given site
        and a 'T' represents presence of SNV at a given site
    """
    
    def __init__(self, cluster_information, original_seq, tumor_seq, tsp_list, clone_frequency_cutoff, CNV_info, ReadCountTable):
        self.Tu2Cluster = cluster_information
        self.CNV_info = CNV_info
        self.ReadCountTable=ReadCountTable
        Align = MegaAlignment()	
        self.OriAncOrder, self.OriAnc2Seq0 = Align.name2seq(original_seq)	
        self.TOrder, self.T2Seq = Align.name2seq(tumor_seq)	
        self.SharePosi= Align.GetSharePosi1(self.OriAnc2Seq0, 'T')
        self.all_tsp = tsp_information(tsp_list)
        self.CloFreCutOff = clone_frequency_cutoff        				
        self.v_obs = self.all_tsp.tumor2alt_frequency()		
        identical_seq_list = Align.identify_similar_seq(tumor_seq, 0)
        self.identical_seq = Align.make_similar_seq_dic(identical_seq_list)			
  
        self.freq_cutoff =self.CloFreCutOff

	
    def combinations(self,target,data,ID,TMP):
 
     for i in range(len(data)):
         new_target = copy.copy(target)
         new_data = copy.copy(data)
         new_target.append(data[i])
         new_data = data[i+1:]
         TMP[ID]=new_target
         ID+=1		 
       #  print 'add',new_target,ID		 
         TMP ,ID=self.combinations(new_target, new_data,ID,TMP)
     return TMP ,ID		 
    def get_decomposed_seq(self):
        print 'Decompose incorrect sample genotype clones'	
        RmCloLs=[]
        DecomRep='n'
        Tu2DecomSeq = {}
        Tu2NewCloFre = {}
        AllSeq=self.OriAnc2Seq0		
        for Tu in self.Tu2Cluster:
            self.ClusterInfo = self.Tu2Cluster[Tu]		
            if self.ClusterInfo != []:
                self.TuSeq = self.T2Seq['#'+Tu]
                self.LenSeq = len(self.TuSeq)				
                DecomSeq, CloneFre = self.get_candidate_decomposed_clones(Tu)
                Tu2DecomSeq[Tu]=DecomSeq			
                Tu2NewCloFre['T-'+Tu]=CloneFre							
                if DecomSeq!=[]:
                    RmCloLs.append(Tu)
                    DecomRep='y'
            else:  
                Tu2DecomSeq[Tu]=[]		                             				
        return 	Tu2DecomSeq, Tu2NewCloFre, RmCloLs				

    def get_candidate_decomposed_clones(self, target_tumor):
        Align = MegaAlignment()	
        CluInf_tu = self.ClusterInfo#[target_tumor]		
        NameOrder, Name2Seq = Align.name2seq(CluInf_tu[2])
      #  print target_tumor, CluInf_tu[0],CluInf_tu[1]	
        HitCloCluLs = CluInf_tu[1]#['T-'+target_tumor]	
        TuIdentical_seq = self.identical_seq['T-'+target_tumor]
        LenSeq = len(Name2Seq[NameOrder[0]])
        TuSeq=self.T2Seq['#'+target_tumor]		
        Clu2center=CluInf_tu[0]
        SigCluLs=[]
        HitCloLs=[]
        HitCloSeq_dic={}		
        RootClu=''		
        LarCen=0.0		
        for Hit in HitCloCluLs:		
          if HitCloCluLs[Hit]>0.02:          		  
            if Hit[:len(target_tumor+'Clu')]==target_tumor+'Clu' and Hit.find('REP')==-1: 
                SigCluLs.append(Hit)
                CluName='Clu'+Hit.split('Clu')[-1]
                Center=float(Clu2center[CluName].split('-')[0])
                				
                for CluN in Clu2center:
                    Center2=float(Clu2center[CluN].split('-')[0])
                    Sign=Clu2center[CluN].split('-')[1]					
                    if Center==Center2 and CluName!=CluN: SigCluLs.append(target_tumor+CluN)
                    if LarCen<Center2:# and Sign=='Pos':  #Pos for middle cut, Neg for K-means
                        LarCen=Center2
                        if Center==Center2: RootClu=target_tumor+CluN						
                    elif  LarCen<=Center2 and Sign=='Pos':  #Pos for middle cut, Neg for K-means
                        LarCen=Center2
                        if Center==Center2: RootClu=target_tumor+CluN                       					
                				
            else: 
               HitCloLs.append(Hit)
               HitCloSeq_dic['#'+Hit]=Name2Seq['#'+Hit]			   

      #  print 'cluls0',SigCluLs, HitCloLs, RootClu	
        if RootClu!='':
             SigCluLs.remove(RootClu)
             HitCloLs.append(RootClu)		
            # print 'cluls',SigCluLs, HitCloLs, RootClu			
        if SigCluLs!=[]: CluCombLs,IDend=self.combinations([],SigCluLs,0,{})
        else:CluCombLs={}
      #  print CluCombLs 			 
        if RootClu!='' or CluCombLs!={}:   
             print 'make cluster comb'		
             CloCan2Seq={}
             Got_Candidate='n'			 
             for Root in HitCloLs:
                RootSeq=Name2Seq['#'+Root]
                LenSeq=len(RootSeq)				
                RootMut=Align.GetMutPos(RootSeq)
                CloCan2Seq['#'+Root]=RootSeq
                Got_Candidate='y'
                if CluCombLs!={}:				
                 for ID in CluCombLs:
                    CluLs=CluCombLs[ID]	
                 #   print 'try make combo',Root,CluLs					
                    CluN=''
                    MutPosLs=[]						
                    for Clu in CluLs:  
                        Seq=Name2Seq['#'+Clu]
                        CluMut=Align.GetMutPos(Seq)
                        MutPosLs+=	CluMut							
                        CluN+=Clu.replace(target_tumor+'Clu','Clu')

                    MutPosLs=list(set(MutPosLs))
                    Go='y'					
                    for Mut in MutPosLs:
                         if RootMut.count(Mut)!=0: Go='n'
					 
                    if Go=='y':	
                         AllMutPosLs=MutPosLs+RootMut					
                         Seq=Align.ModSeq('A'*LenSeq,AllMutPosLs,'T',LenSeq)
                         Redun_ls=Align.find_redundant(Seq,HitCloSeq_dic)	
                       				
                         if Redun_ls==[]:			
                            CloCan2Seq['#'+target_tumor+Root.replace(target_tumor+'Clu','Clu')+CluN]=Seq
                            Got_Candidate='y'						
  
             if Got_Candidate=='y':	  
                      Can_list=CloCan2Seq.keys()
                   #   print 'find the good comb',Can_list		
         
                      new_seq = Align.UpMeg(CloCan2Seq,Can_list)			   
                      alt_frequency = []
                      CNVls=self.CNV_info[target_tumor]
                      Len=len(CNVls)
                      c=0	
                      TuMatPosi=[]
                      tumor_genotype=''					  
                      while c<Len:
                        if CNVls[c]=='normal':			 
                            alt_frequency.append(self.v_obs[target_tumor][c])
                            if 	self.v_obs[target_tumor][c]>0: 
                               TuMatPosi.append(c)
                               tumor_genotype+='T'
                            else: tumor_genotype+='A'							   
                        c+=1	
             		   
                      clone_frequency = CloneFrequencyComputer_cnv1({}, {}, {}, self.freq_cutoff, {})
         	
                      MutWildAlleleCount_noCNV = clone_frequency.make_mut_wild_allele_count_noCNV({}, Can_list, CloCan2Seq)#PreAbsCNV, clone_order, SNV_seq, Tu2CloFre		
                      Cmatrix_noCNV, Cmatrix_noCNV_dic = clone_frequency.make_Min(Can_list, CloCan2Seq, MutWildAlleleCount_noCNV)				 
                      Clone2Freq = clone_frequency.do_nnls0(Cmatrix_noCNV, Can_list, alt_frequency)			 
 
                      out2=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']	
                      AllMut=[]
                      NewClone2Freq={}					  
                      CluHit='n'					  
                      for Clo0 in Clone2Freq:
                            NewClone2Freq[Clo0]=Clone2Freq[Clo0]					  
                            if Clone2Freq[Clo0]>0.02:
                           
                                SeqMutPos=Align.GetMutPos(CloCan2Seq['#'+Clo0])
                                TuSeq='y'
                                for Mut in 	SeqMutPos:
                                    if TuMatPosi.count(Mut)!=0: AllMut.append(Mut)								
                                for Mut in TuMatPosi:
                                     if SeqMutPos.count(Mut)==0: TuSeq='n'
                                Iden='n'
                                for OriClo in self.OriAnc2Seq0:
                                      c=0	
                                      Dif='n'									  
                                      while c<Len:
                                           if self.OriAnc2Seq0[OriClo][c]!=	CloCan2Seq['#'+Clo0][c]: Dif='y'
                                           c+=1
                                      if Dif=='n': Iden=OriClo 
                                if Iden!='n': 
                                      out2+=[Iden, self.OriAnc2Seq0[Iden]] 
                                      NewClone2Freq[Iden[1:]]=Clone2Freq[Clo0]
                                      NewClone2Freq[Clo0]=0		                                      								
                                elif TuSeq=='n': 
                                			   
                                     out2+=['#'+Clo0.replace(target_tumor+target_tumor,target_tumor), CloCan2Seq['#'+Clo0]] 
                                     if Clo0.find('Clu')!=-1 and Clo0.find('REP')==-1: CluHit='y'									 
                                else:
                                      out2+=['#'+target_tumor, CloCan2Seq['#'+Clo0]] 
                                      NewClone2Freq[target_tumor]=Clone2Freq[Clo0]
                                      NewClone2Freq[Clo0]=0									  
                      AllMut=list(set(AllMut))
                      if len(AllMut)<len(TuMatPosi):out2+=['#'+target_tumor, tumor_genotype] 					  
                      if CluHit=='y':
                   		
                       #  print 'Decomposed!'	,target_tumor,NewClone2Freq,out2
                         return out2,NewClone2Freq			

        return [],{}		

	   


  