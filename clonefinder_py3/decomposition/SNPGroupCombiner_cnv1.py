from alignments.MegaAlignment import MegaAlignment
from decomposition.SNPClusterGenerator_cnv1 import SNPClusterGenerator_cnv1
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
from alignments.FreqToMegaSeq import FreqToMegaSeq
import copy
class SNPGroupCombiner_cnv1():
    
 
    def __init__(self, ini_seq_builder, v_obs, clone_frequencies, CNV, freq_cutoff, tumor_seqs):
        self.freq_cutoff = freq_cutoff			
        self.Tu2CloFre = clone_frequencies
        self.CloFreCutOff = self.freq_cutoff        				
        self.v_obs = v_obs
        Align = MegaAlignment()
        self.clone_order, self.clone_seq = Align.name2seq(ini_seq_builder)
        self.ini_seq_builder = ini_seq_builder		
        self._CNV_file =CNV		
        self.snv_num = len(self.clone_seq[self.clone_order[0]])	
        self.tumor_seqs = tumor_seqs	
    def combinations(self,target,data,ID,TMP):
 
     for i in range(len(data)):
         new_target = copy.copy(target)
         new_data = copy.copy(data)
         new_target.append(data[i])
         new_data = data[i+1:]
         TMP[ID]=new_target
         ID+=1		 
      	 
         TMP ,ID=self.combinations(new_target, new_data,ID,TMP)
     return TMP ,ID
    def extract_hitseq(self,seq_buil,CloFre,Cut):
        Align=MegaAlignment()
        CloLs,Clo2Seq=Align.name2seq(seq_buil)	
        Hit={}
        for Clo in CloFre:
             if CloFre[Clo]>Cut: Hit['#'+Clo]=Clo2Seq['#'+Clo]
        return Hit			 
    def get_decomposed_seq(self):
        Align=MegaAlignment()	
        TuLs, Tu2Seq = Align.name2seq(self.tumor_seqs)			
        print('make SNV clusters')
        clusters = SNPClusterGenerator_cnv1(self.ini_seq_builder, self.v_obs, self.Tu2CloFre, self._CNV_file, self.freq_cutoff)		
        Tumor_cluster_dic = clusters.cluster_cnv()	#Tu2Cluster={tumor:[[seq_builder,{tumor:{clone frequency}}]]}			
        print('Decompose incorrect sample genotype clones')	
 
        AllhitWithDecom={}	
        All_convol_tuseq=[]	
        DecomLs=[]		
        
        for Tu in Tumor_cluster_dic:
            ClusterInfo = Tumor_cluster_dic[Tu]		
            if ClusterInfo != []:
           
                			  
                HitWithDecomSeq_build,convol_tuseq = self.get_candidate_decomposed_clones(Tu,ClusterInfo,Tu2Seq['#'+Tu])
                if convol_tuseq!='': 				
                    A1,HitWithDecomSeq_dic=Align.name2seq(HitWithDecomSeq_build)
                    AllhitWithDecom.update(HitWithDecomSeq_dic)	
                    All_convol_tuseq.append(convol_tuseq)
                    DecomLs.append(Tu)					
                  
                else:
                    Original_hit_seq_dic = self.extract_hitseq(self.ini_seq_builder,self.Tu2CloFre['T-'+Tu],self.freq_cutoff)				
                    AllhitWithDecom.update(Original_hit_seq_dic)
            else:
                    Original_hit_seq_dic = self.extract_hitseq(self.ini_seq_builder,self.Tu2CloFre['T-'+Tu],self.freq_cutoff)				
                    AllhitWithDecom.update(Original_hit_seq_dic)			   
            
        if DecomLs==[]:
          	
             return self.clone_seq,'no decomposed clone was made'
        else:
         	
            for ConvTuSeq in All_convol_tuseq:
                Redun_ls=Align.find_redundant(ConvTuSeq,AllhitWithDecom) 
                if Redun_ls!=[]:
                    				
                     return self.clone_seq,'tumor genotype that was decomposed was hit in different tumor: failed decomposition'	
          
            return AllhitWithDecom,'decomposed'+str(DecomLs)			
     			
    def findcombohit(self,seq_builder):
        Align=MegaAlignment()
        SeqLs,SeqDic=Align.name2seq(seq_builder)
        Find='n'
        for i in SeqLs:
            if i.find('Clu')!=-1: Find='y'
        return Find			
    def get_candidate_decomposed_clones(self, target_tumor, CluInf_tu,Tuseq):
        Align = MegaAlignment()	
     	
        NameOrder, Name2Seq = Align.name2seq(CluInf_tu[0])
	
        LenSeq = len(Name2Seq[NameOrder[0]])
    
        SigCluLs=[]		
        for Name in NameOrder: #Root is the first cluster or initial candidate clone
               if Name!='#Clu0' and Name.find('Clu')!=-1: SigCluLs.append(Name)
        CluCombLs,IDend=self.combinations([],SigCluLs,0,{})   
        print(target_tumor,'make cluster comb',SigCluLs,CluCombLs,NameOrder)
   	
        if CluCombLs!={}:   
			 
             CloCan2Seq={}
             Got_Candidate='n'			 
             for Root in NameOrder: #Root is the first cluster or initial candidate clone
               if Root=='#Clu0' or Root.find('Clu')==-1:			 
                RootSeq=Name2Seq[Root]
                if Root=='#Clu0': CloCan2Seq['#'+target_tumor+'Clu0']=RootSeq  #Root is candidate clone              				
                RootMut=Align.GetMutPos(RootSeq)
                Got_Candidate='y'
                if CluCombLs!={}:				
                 for ID in CluCombLs:
                    CluLs=CluCombLs[ID]	
            			
                    CluN=''
                    MutPosLs=[]						
                    for Clu in CluLs:  
                        Seq=Name2Seq[Clu]
                        CluMut=Align.GetMutPos(Seq)
                        MutPosLs+=	CluMut							
                        CluN+=Clu.replace('#','')

                    Good='y'					
                    for Mut in MutPosLs:
                         if RootMut.count(Mut)!=0: Good='n'
					 
                    if Good=='y':	
                         AllMutPosLs=MutPosLs+RootMut					
                         Seq=Align.ModSeq('A'*LenSeq,AllMutPosLs,'T',LenSeq)
                         Redun_ls=Align.find_redundant(Seq,self.clone_seq) #all other clones ####	
                       				
                         if Redun_ls==[]:			
                            CloCan2Seq['#'+target_tumor+Root.replace('#','')+CluN]=Seq
                   					
  
             if CloCan2Seq!={}:	  
	
                      CloCan2Seq.update(self.clone_seq) 
                      Can_list=list(CloCan2Seq.keys())					  
                            					  
                      new_seq = Align.UpMeg(CloCan2Seq,Can_list)
                   						   
                      clone_frequency_combo = CloneFrequencyComputer_cnv1(new_seq, {target_tumor:self.v_obs[target_tumor]}, {target_tumor:self._CNV_file[target_tumor]}, self.freq_cutoff)
                      clone_frequency_combo.regress_cnv()					
                      CluComboHit=self.findcombohit(clone_frequency_combo.hitclone_seq_builder)
                      if CluComboHit=='y':
                            print('test the quality of clustercombo, by removing tumor seq (if any)')
                            hit_seq_ls,hit_seq_dic=Align.name2seq(clone_frequency_combo.hitclone_seq_builder) 							
                            Tuseq_ls=Align.find_redundant(Tuseq,hit_seq_dic)	
                            if Tuseq_ls==[]:
                                  print('tumor genotype did not hit, so clustercombo is good')							
                                  return clone_frequency_combo.hitclone_seq_builder,Tuseq
                            else:
                                  print('tumor genotype was hit, so test if clustercombo still hit without tumor genotype: testing if clustercombo genotypes fit well')
                                  Tuseq_ls=Align.find_redundant(Tuseq,CloCan2Seq)									  
                                  sub_hit_seq=[]
                                  for seqname in CloCan2Seq:
                                        if Tuseq_ls.count(seqname)==0:sub_hit_seq+=[seqname,CloCan2Seq[seqname]]
                            
                                  clone_frequency_combo_new = CloneFrequencyComputer_cnv1(sub_hit_seq, {target_tumor:self.v_obs[target_tumor]}, {target_tumor:self._CNV_file[target_tumor]}, self.freq_cutoff)
                                  clone_frequency_combo_new.regress_cnv()					
                                  CluComboHit=self.findcombohit(clone_frequency_combo_new.hitclone_seq_builder)
                                  if CluComboHit=='y': 
                           					  
                                     return clone_frequency_combo_new.hitclone_seq_builder,Tuseq 
                                  else: 
                                     return CluInf_tu[0],''								  
                      else: return CluInf_tu[0] ,''                                 								  
                                							
                      	
             else: return CluInf_tu[0],''
        return CluInf_tu[0],''		

	   


  