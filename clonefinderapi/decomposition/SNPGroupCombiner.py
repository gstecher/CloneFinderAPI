from alignments.MegaAlignment import MegaAlignment
#from alignments.AlignmentAnnalizer import AlignmentAnnalizer
from tsp_profiles.tsp_information import tsp_information
from regression.CloneFrequencyComputer import CloneFrequencyComputer

class SNPGroupCombiner():
    """
        Generate a MEGA sequence alignment file from a TumorSampleProfile
        
        A MEGA dna sequence alignment will be generated based on presence/absence
        of SNVs in the tumor sample profile. A sequence is generated for each
        tumor sample where a an 'A' represents absence of SNV at a given site
        and a 'T' represents presence of SNV at a given site
    """
    
    def __init__(self, cluster_information, original_seq, tumor_seq, tsp_list, clone_frequency_cutoff, optionOnOff, CluNumCutOff):
        self.Tu2Cluster = cluster_information
        Align = MegaAlignment()	
        self.OriAncOrder, self.OriAnc2Seq0 = Align.name2seq(original_seq)	
        self.TOrder, self.T2Seq = Align.name2seq(tumor_seq)	
        self.SharePosi= Align.GetSharePosi1(self.OriAnc2Seq0, 'T')
        self.all_tsp = tsp_information(tsp_list)
        self.CloFreCutOff = clone_frequency_cutoff        				
        self.v_obs = self.all_tsp.tumor2alt_frequency()		
        identical_seq_list = Align.identify_similar_seq(tumor_seq, 0)
        self.identical_seq = Align.make_similar_seq_dic(identical_seq_list)			
        self.TrunkAnc = optionOnOff	
        self.ClusterNum = CluNumCutOff
      	
    def get_decomposed_seq(self):
        print 'Decompose incorrect sample genotype clones'	
        RmCloLs=[]
        DecomRep='n'
        Tu2DecomSeq = {}		
        for Tu in self.Tu2Cluster:
            self.ClusterInfo = self.Tu2Cluster[Tu]		
            if self.ClusterInfo != []:
                self.TuSeq = self.T2Seq['#'+Tu]
                self.LenSeq = len(self.TuSeq)				
                DecomSeq = self.get_candidate_decomposed_clones(Tu)
                Tu2DecomSeq[Tu]=DecomSeq
                if DecomSeq!=[]:
                    RmCloLs.append(Tu)
                    DecomRep='y'
				
        return 	Tu2DecomSeq, RmCloLs				

    def get_candidate_decomposed_clones(self, target_tumor):
        Align = MegaAlignment()	
        NameOrder, Name2Seq = Align.name2seq(self.ClusterInfo[2])
        HitCloCluLs = self.ClusterInfo[1]['T-'+target_tumor]	
        TuIdentical_seq = self.identical_seq['T-'+target_tumor]
        LenSeq = len(Name2Seq[NameOrder[0]])
        TuSeq=self.T2Seq['#'+target_tumor]		

        SigCluLs=[]
        HitCloLs=[]
        for Hit in HitCloCluLs:
          if HitCloCluLs[Hit]>0:		
            if Hit[:len(target_tumor+'Clu')]==target_tumor+'Clu' and Hit.find('REP')==-1: SigCluLs.append(Hit)
            else: HitCloLs.append(Hit)
		
        CluTa = self.ClusterInfo[0]
        Freq2Clu={}
        FreqLs=[]
		
        for Clu in CluTa:
            		
            FreqNegPos=CluTa[Clu].split('-')		
            if FreqNegPos[0][-1]!='e':
		
              Freq=float(FreqNegPos[0])
              NegPos=FreqNegPos[1]
              Code=Freq in Freq2Clu
              if Code!=True: 
                Freq2Clu[Freq]=[]
                FreqLs.append(Freq)
              CloID=target_tumor+Clu
              Code=CloID in SigCluLs	
              if Code==True:Freq2Clu[Freq].append(Clu+NegPos)

        FreqLs.sort()
        FreqLs.reverse()
        CluOrder=[]
		
        for Freq in FreqLs:
                CluLs=Freq2Clu[Freq]
                NegPosDic={'Neg':'','Pos':''}
                for Clu in CluLs:
                    NegPos=Clu[-3:]
                    NegPosDic[NegPos]=Clu
                if NegPosDic['Pos']!='':CluOrder.append(NegPosDic['Pos'])
                if NegPosDic['Neg']!='': CluOrder.append(NegPosDic['Neg'])
        		
        TrunkClu=[]
	
        for Clu in CluOrder:  
            Seq=Name2Seq['#'+target_tumor+Clu[:-3]]
            CluMut=Align.GetMutPos(Seq)
            Trunk='n'	
            for i in CluMut:	
                Code=i in self.SharePosi
                if Code==True: Trunk='y' 
            if Trunk=='y':TrunkClu.append(Clu)
                		
        NewNameOrder=[]
        if TrunkClu!=[]:
	
            if self.TrunkAnc=='On':CluOrderNew=self.OrderClu(CluOrder,TrunkClu)
            else: CluOrderNew=CluOrder	
            if len(CluOrderNew)>self.ClusterNum: CluOrderNew=CluOrderNew[:self.ClusterNum]
            CloCan=self.GetCloCan(CluOrderNew)
            CloCan2Seq={}   
            for CluLs in CloCan:
                    CluN=''
                    MutPosLs=[]			
                    for Clu in CluLs:    		
                        CluN+=Clu[:-3]
                        Clu=Clu[:-3]
                        Seq=Name2Seq['#'+target_tumor+Clu]
                        CluMut=Align.GetMutPos(Seq)
                        MutPosLs+=	CluMut			       
                    MutPosLs=list(set(MutPosLs))
                    Seq=Align.ModSeq('A'*LenSeq,MutPosLs,'T',LenSeq)
                    DifNum=Align.CountDifNum(Seq,TuSeq)			        		 
                    if DifNum!=0:			
                        CloCan2Seq[CluN]=Seq
                        NewNameOrder.append(CluN)
        				
        NewCloLs=[]
        out2=[]
        AllSeq=self.OriAnc2Seq0
        if NewNameOrder!=[]:
		
          Added='n'	
          for Clu in NewNameOrder:
            Clu='#'+Clu
            Seq=CloCan2Seq[Clu[1:]]   				
            RmClu='n'	#Should not be ancestor of other clones
            if self.TrunkAnc=='On': 	
              for Ori in self.OriAnc2Seq0:
                
                if TuIdentical_seq.count(Ori[1:])==0:				

                  OriSeq=self.OriAnc2Seq0[Ori]
                  RmClu=self.DecideRmClu(Seq,OriSeq,RmClu)
         	 
            if RmClu=='n':
                Added='y'	
                AllSeq[Clu]=CloCan2Seq[Clu[1:]]		
                NewCloLs.append(Clu)		
        
          CluHit='n' 
          if Added=='y' :
	  
               NewCloLs+=HitCloLs 

               new_seq = Align.UpMeg(AllSeq,NewCloLs)			   

               alt_frequency = self.all_tsp.make_single_tsp_list(target_tumor)			   
               clone_frequency = CloneFrequencyComputer(new_seq, alt_frequency, self.CloFreCutOff)
               hitclone_sequences, hitclone_frequency_dic = clone_frequency.regress()			   
     	   
               CloLs0=hitclone_frequency_dic['T-'+target_tumor]
               		   
               for Clo0 in CloLs0:
                   if CloLs0[Clo0]>0:			   
                       if Clo0.find('Clu')!=-1 and Clo0.find('REP')==-1: CluHit='y'
                       if TuIdentical_seq.count(Clo0)==0: 
                            out2+=['#'+Clo0, AllSeq['#'+Clo0]] 			   
          if CluHit=='y':
            print 'Decomposed!'	,target_tumor
            return out2			

        return []		


    def GetCloCan(self, Ls):  
        ID2Name={}
        ID=0
        for Name in Ls:
            ID2Name[ID]=Name
            ID+=1
        IDMax=ID		
        CombLs=[[0]]
        NewCombLs=CombLs	
        Go='y'
        while Go=='y':
          NewCombLs0=[]
          for Comb in NewCombLs:
            ID=Comb[-1]+1
            if ID==IDMax: Go='n'		
            New=Comb		
            while ID<IDMax:
                New=Comb+[ID]
                NewCombLs0.append(New)		
                ID+=1 
          NewCombLs=NewCombLs0		
          CombLs+=NewCombLs   
        IDCombLs=CombLs[1:]
        TMP=[]
        for IDComb in IDCombLs:
            In=[]		
            for ID in IDComb:
                In.append(ID2Name[ID])         		
            TMP.append(In)
        TMP.append([Ls[0]])
        return TMP

    def  DecideRmClu(self, Seq,OthSeq,OriRmClu):	  
          c=0
          Bad='y'
          Len=len(Seq)  
          while c<Len:
             if Seq[c]=='T' and OthSeq[c]=='A': Bad='n'	 
             c+=1
          if Bad=='y': OriRmClu='y'
          return OriRmClu
		  
    def OrderClu(self, InOrder,Trunk):
        for i in InOrder:
           Code=i in Trunk
           if Code==True:
               First=i
               break
        New=[First]
        for i in InOrder:
            if i!=First: New.append(i)
        return New
##	

  