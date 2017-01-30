from alignments.MegaAlignment import MegaAlignment
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from tsp_profiles.tsp_information import tsp_information
from output.CloneFrequencyWriter import CloneFrequencyWriter
from decomposition.SNPClusterGenerator import SNPClusterGenerator

class MutationFilter:

    def __init__(self, OptionA, OptionC, OptionD, tumor_seq, tsp_list, ErroRate, CloFreCut, mao_file):
        self.OptionA = OptionA
        self.OptionC = OptionC
        self.OptionD = OptionD.split(',')
        self.ErroRate = ErroRate
        Align = MegaAlignment()
        self.tumor_list, self.tumor2seq = Align.name2seq(tumor_seq)
        self.Len= len(self.tumor2seq[self.tumor_list[0]])
        self.mao_file = mao_file
        self.tsp_list = tsp_list
        self.CloFreCut = CloFreCut
        TSPinfo=tsp_information(tsp_list)		
        self.Tu2SNV=TSPinfo.tumor2alt_frequency()		
       	
    def CutExtraMutClo(self, seq_list, clone_frequency):
        Align = MegaAlignment()
        CloOrder, Clo2Seq = Align.name2seq(seq_list)
        Trunk_mutation_position = Align.GetSharePosi1(Clo2Seq, 'T')		
        NewClo2Seq='n'
        newseq_list = ['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
        new_clone_frequency = {}		
        Done=[]	 	
        for Tu in clone_frequency:
	
            Cls=clone_frequency[Tu]			
            Tu=Tu[2:]		
            TuSeq=self.tumor2seq['#'+Tu]	
            Clo2Freq={}	
            for Clo in Cls:
               CloFre=Cls[Clo]			
               if CloFre>0:			
                 CloSeq=Clo2Seq['#'+Clo]
                 OthSeqDic = self.get_other_seq('#'+Clo,Clo2Seq)				 
                 ChangePos=[]
                 other_wild = Align.GetSharePosi1(OthSeqDic, 'A')
            #     print 'below other wild and trunk',other_wild,Trunk_mutation_position 				 
                 GoodPosi=0				 
                 c=0
                 while c<self.Len:			     
                    if TuSeq[c]=='A' and CloSeq[c]=='T' : 
                         ChangePos.append(c)
                    elif TuSeq[c]=='T' and CloSeq[c]=='T':
                         if other_wild.count(c)==0 and Trunk_mutation_position.count(c)==0: #not uniqu not trunl
                              GoodPosi+=1						 
                    c+=1
                 ChangePos_proportion = 1.0*len(ChangePos)/self.Len	
                 GoodPos_proportion = 1.0*GoodPosi/self.Len	
              #   print Tu, Clo, ChangePos, GoodPosi                				 
                 if ChangePos_proportion>self.ErroRate and GoodPos_proportion>self.ErroRate: 
                   # if GoodPos_proportion>self.ErroRate:
                     print 'h'					
                     NewSeq=Align.ModSeq(CloSeq,ChangePos,'A',self.Len)
                     Clo2Freq[Clo+'Mod'+Tu] = CloFre                    					 
                     newseq_list+=['#'+Clo+'Mod'+Tu, NewSeq]					 
                     NewClo2Seq='y'
                 else: 
                    Clo2Freq[Clo] = CloFre    				 
                    if Done.count('#'+Clo)==0:
                         newseq_list+=['#'+Clo, CloSeq]
                         Done.append('#'+Clo)						 	
            new_clone_frequency['T-'+Tu]=Clo2Freq			
        newseq_list_without_redundant, Ignore, new_clone_frequency_without_redundant = Align.CombSimClo(newseq_list, new_clone_frequency, 0.0)		
        return newseq_list_without_redundant, new_clone_frequency_without_redundant	

    def get_other_seq(self, CloTar, CloSeq):
        Dic={}
        for Clo in CloSeq:		 
            if Clo!=CloTar: Dic[Clo]=CloSeq[Clo]
        return Dic
		
    def BranchDecClone(self, seq_list, clone_frequency):		
        Align = MegaAlignment()
        TumorSampleExtract = tsp_information(self.tsp_list)	
        CloFreAna = CloneFrequencyAnalizer()      		
        CloOrder, Clo2Seq = Align.name2seq(seq_list)
        Align.save_mega_alignment_to_file('Test.meg', seq_list)		
        tree_builder = MegaMP()
        tree_builder.mao_file = self.mao_file		
        id = 'branchdec_mega_alignment' 
   
        status = tree_builder.do_mega_mp(seq_list, id)
        if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
        else:
                print 'failed to run megaMP'
        BadPosiLs=[] #multiple mutations
        BadPosi2ChnageCloLs={}		
        for c in Good_posi_info:
           Posi_Inf = Good_posi_info[c]		
           if Posi_Inf!=['Good']:
              if Posi_Inf[0] == 'ToWild':		   
                  BadPosiLs.append(c)
                  BadPosi2ChnageCloLs[c]=Posi_Inf[1][0]				  
        if BadPosiLs!=[]:	 
         NewT2C2F={}
         NewT2Cls={}
         for Tu in clone_frequency:
            NewC2F=	{}
            single_tsp_list = TumorSampleExtract.make_single_tsp_list(Tu)           
            CloFreDic=clone_frequency[Tu]
            Tu=Tu[2:]			
            TuSeq	=self.tumor2seq['#'+Tu]
            NewCloLs=[]	
            NewCloLs1=[]
					
            for Clo in CloFreDic: #original hit clo for the tumor
              ChangeOptions='n'				
              if CloFreDic[Clo]>0:	 
                 CSeq0=Clo2Seq['#'+Clo]
                 ChangePosi=[] #list to fix multiple mutaitons
                 NewBadPosi=[] #remove fixed multiple mutations from BadExtMutPosi				 
                 for Bad in BadPosi2ChnageCloLs:
                     if BadPosi2ChnageCloLs[Bad].count('#'+Clo)!=0:	
                        Change='n'					 
                        for Oth in CloFreDic: #find multiple mutations at the external branch
                            if Oth !=Clo and CloFreDic[Oth]>0:
                               Soth=Clo2Seq['#'+Oth]
                               if Soth[Bad]=='T' and BadPosi2ChnageCloLs[Bad].count('#'+Oth)==0: 
                                         Change='y'
                        if Change=='y':										 
                                         ChangePosi.append(Bad)
                        else: 	NewBadPosi.append(Bad)										 
                 if self.OptionA=='On':										 
                   if ChangePosi!=[]:#fix multiple mutaitons
                              CutCloSeq=Align.ModSeq(CSeq0,ChangePosi,'A',self.Len)					  				  
                              NewCloLs.append(Clo+'Cut'+Tu)
                              NewC2F[Clo+'Cut'+Tu]=CloFreDic[Clo]
                              Clo2Seq['#'+Clo+'Cut'+Tu]=CutCloSeq
                              ChangeOptions='y'							  		               		 
                 if self.OptionC!='Off':
                   OptionCLs=self.OptionC.split(',')		 
                   if len(NewBadPosi)>=int(OptionCLs[1]):#5: #if there are still unfixed multiple mutations, examine if these can be moved to another lineage      		 
                      Code='#'+Clo+'Cut'+Tu in Clo2Seq
                      if Code==True:CSeq0=Clo2Seq['#'+Clo+'Cut'+Tu]			  
                      LarComMutNum=[]       			  
                      TarClo=''			  
                      for OthClo in CloOrder: #clones in other tumors
                         OthClo=OthClo[1:]					  
                         Code= OthClo in CloFreDic
                         if Code!=True or (Code==True and CloFreDic[OthClo]==0):
                              OSeq=Clo2Seq['#'+OthClo]
                              TMP=[]
                              MutCloNum=0
                              for OthClo1 in CloOrder:
                                 Fix=[] #remove from the clone
                                 NewExtra=[] #remove from the other clone						 
                                 for MutP in NewBadPosi: 							  
                                     OthSeq1=Clo2Seq[OthClo1]								  
                                     if OthSeq1[MutP]=='T' :
                                            MutCloNum+=1
                                            Fix.append(MutP)
                                     Len=len(OthSeq1)
                                     cO=0
                                     while cO<Len:
                                         if TuSeq[cO]=='A' and OthSeq1[cO]=='T': NewExtra.append(MutP)
                                         cO+=1										 
                                 UnFixNum=len(NewBadPosi)-len(Fix)                             								 
                                 if len(Fix)>len(NewExtra) and UnFixNum<=len(NewExtra):								 
                                   LarComMutNum=[Fix,NewExtra]
                                   TarClo=OthClo1
                      if TarClo!='': #multiple mutations can be fixed
                            TarClo=TarClo[1:]					  			  					  
                            CutCloSeq=''
                            NewCloSeq=''
                            CutFixPosi=LarComMutNum[0]
                            CutExtraPosi= LarComMutNum[1]                           							
                            OthSeq=Clo2Seq['#'+TarClo]				  
                            c=0
                            UnchangeC=0				
                            while c<Len:
                                  Code=c in CutFixPosi#LarComMutNum                  						  
                                  if Code==True:
                                     CutCloSeq+='A'                        							 
                                  else: CutCloSeq+=CSeq0[c]
                                  Code=c in CutExtraPosi#LarComMutNum                  						  
                                  if Code==True:
                                     NewCloSeq+='A'                        							 
                                  else: NewCloSeq+=OthSeq[c]								  					  
                                  c+=1				                   
                            TestCloLs=[Clo+'Cut1'+Tu,TarClo+'Fix1'+Tu]					
                            Clo2Seq['#'+Clo+'Cut1'+Tu]=CutCloSeq							
                            Clo2Seq['#'+TarClo+'Fix1'+Tu]=NewCloSeq
                            ChangeOptions='y'
                            NewC2F[Clo+'Cut1'+Tu]=1
                            NewC2F[TarClo+'Fix1'+Tu]=1							
              if ChangeOptions=='n':
                   NewC2F[Clo]=1			     
            NewT2C2F[Tu]=NewC2F
         hitseq_align, hitclone_frequency  = CloFreAna.ListHitCloAndSeq(NewT2C2F, Clo2Seq)		 
         outSeqMaj, outSeqAmb, NewT2C2F = Align.CombSimClo(hitseq_align, hitclone_frequency, 0.0)	 
         return  outSeqMaj, NewT2C2F  		 
        else:
         return  seq_list, clone_frequency  			

    def BranchCloneInfer(self, Alignment, clone_frequency, NodeMapLs):
        Align = MegaAlignment()	
        CloOrder, Clo2Seq = Align.name2seq(Alignment)
        TreeAna	= TreeAnalizer() 
		
        print 'Make clusters for predicting ancestral clones'		
        MakeCluster= SNPClusterGenerator(Alignment, self.tsp_list, clone_frequency, self.OptionD[1].strip(), self.CloFreCut)	
        Tumor2CluInfo = MakeCluster.cluster() 
        AllMut=Align.GetSharePosi1(Clo2Seq,'T')
        InitialCloneLs=[]
        for Clo in Clo2Seq:
             if Clo!='#hg19': InitialCloneLs.append(Clo)	   	
        CladeLs=TreeAna.list_clade(NodeMapLs, InitialCloneLs)
	
        NewT2C2F={}
        NewCloLs=[]
        Bra2NewCloneLs={}

        AddSeqID=1		
        for Tu in Tumor2CluInfo:
          CluInfo=Tumor2CluInfo[Tu]		
          NewT2C2F[Tu]=clone_frequency['T-'+Tu]
          CloLs=clone_frequency['T-'+Tu]		  
          CloLsIn=[]
          for Clo0 in CloLs:
                        if CloLs[Clo0]>0: CloLsIn.append(Clo0) 	
          if CluInfo!=[]:
            Clu2Fre=CluInfo[1]['T-'+Tu]
            Clu2EstSNVfre=CluInfo[0]
            Clu_seq=CluInfo[2]		
            HitCluLs0=[]
            for Clu in Clu2Fre: 
                if Clu.find('Clu')!=-1 and Clu2Fre[Clu]>0: HitCluLs0.append(Clu)
			
            if HitCluLs0!=[]:
              CluNameLs, Clu2Seq=Align.name2seq(Clu_seq)			
              BraID=0
              NewCloLs=[]
       			  
              for Clo in CloLsIn: #original hit clo for the tumor		   
                AncLs0=[] #from recent
                for Clade in CladeLs:
                    if Clade.count('#'+Clo)!=0: AncLs0.append(Clade)
                AncLs=self.LineageSort(AncLs0)					

                AncLs=['']+AncLs		
                AncLen=len(AncLs)
                Ac=1
                BraClo=[]	#clusters in anc clone[[Clu1,Clu2],[Clu1, Clu2,Clu3]]

                while Ac<AncLen:
                    #is there a hit cluster in the branch?		
                    Anc0=AncLs[Ac-1] #older anc		
                    Anc1=AncLs[Ac] #younger anc in the branch
                    NewAncLs=AncLs[:Ac][1:]	#ancestors before the branch		                 
                    if Anc0=='': #branch from normal
                       MutPosAnc0=[]
                       MutPosAnc1=AllMut #trunc mutation			   
                    else: #others
                       Seq0ls=Align.UpMeg(Clo2Seq, Anc0)
                       Seq1ls=Align.UpMeg(Clo2Seq, Anc1)					   
                       Clo0Ls0, Dic0=Align.name2seq(Seq0ls)				
                       MutPosAnc0=Align.GetSharePosi1(Dic0,'T')	#Pre
                       Clo1Ls0, Dic=Align.name2seq(Seq1ls)			   
                       MutPosAnc1=Align.GetSharePosi1(Dic,'T')	#Current		
                    AddLs=[] #add clusters by mapping from normal			   
                    for Clu0 in HitCluLs0:
                         MutPosClu=Align.GetMutPos(Clu2Seq['#'+Clu0])		                 			 
                         GoodClu,BadMut=self.CheClu1(MutPosAnc1,MutPosClu)	#GoodClu share mutation with ancestor1			 
                         if GoodClu=='y':#Extra mut in cluster is found in the ancestor0? ancestor0 should not have these extra mutations, but when ancestral seq inferences are difficult, it may happen.
                           Add='y'
                           for Bad in BadMut:
                              Code=Bad in MutPosAnc0
                              if Code!=True: Add='n'               					  
                           if Add=='y':					  
                             AddLs.append(Clu0)                 					 
                    c=0
                    AddLen=len(AddLs)
                    #combine clusters at the branch and previous ancestral clones to make ancestral clone types            			
                    if AddLen>0:			
                       while c<AddLen:
                           AddNewAncLs=[]
                           for Clade in NewAncLs:
                                SeqlsAdd=Align.UpMeg(Clo2Seq, Clade)	
                                ClolsAdd, DicAdd=Align.name2seq(SeqlsAdd)								
                                MutShare = Align.GetSharePosi1(DicAdd,'T')
                                snv_num=len(Clo2Seq[Clade[0]])
                                AddSeq=''
                                cc=0
                                while cc<snv_num:
                                    if MutShare.count(cc)==0: AddSeq+='A'
                                    else: AddSeq+='T'
                                    cc+=1
                                AddNewAncLs.append('AddAnc'+str(AddSeqID))
                                Clu2Seq['#AddAnc'+str(AddSeqID)]=AddSeq								
                                AddSeqID+=1							   
                                                     					   
                           BraClo.append(AddNewAncLs+[AddLs[c]])				   
                           c+=1
                    Ac+=1
                				
                for BcloLs in BraClo:
                    AllMutPosi=[] #get all mutation sites in the ancester		
                    for Bclo in BcloLs:
					 
                         Seq=Clu2Seq['#'+Bclo]						 
                         Dic1={'Anc1':Seq,'Cur':Clo2Seq['#'+Clo]}
                         MutPos=Align.GetSharePosi1(Dic1,'T')				 
                         AllMutPosi+=MutPos
        			######
                    Len=len(Clu2Seq['#'+Bclo])					
                    NormalSeq='A'*Len			
                    BracSeq=Align.ModSeq(NormalSeq,AllMutPosi,'T',Len)
				
                    Redun=self.RedunChe(BracSeq,Clo2Seq,CloLsIn+NewCloLs,int(self.OptionD[1].strip()))		
                    if Redun=='n':	#do not add if ancestral sequence is the same as node-ancestral clone (to define 'same', allow error rate difference)
                        if len(BcloLs)==1: BraID0='Root' 			
                        else:
                            for i in BcloLs:
                                if i.find('Clu')!=-1: break
                                else: BraName=i						
                            BraID0=BraName
                        Code=BraID0 in Bra2NewCloneLs
                        if Code!=True: Bra2NewCloneLs[BraID0]=[]			
                        Clo2Seq['#'+Tu+'BraC'+str(BraID)]=BracSeq
                        NewT2C2F[Tu][Tu+'BraC'+str(BraID)]=1
        
                        NewCloLs.append(Tu+'BraC'+str(BraID))
                        Bra2NewCloneLs[BraID0].append(Tu+'BraC'+str(BraID)) #for each branch, list ancestral clone               				
                        BraID+=1
        				
         #if different tumor sample make ancestral clones at a same branch, these ancestral clones should have ancestor-descendant relationship.
        CombCloLs=[]
        Old2NewName={} #Adjusted ancestral clones
        NewCloLs=CloOrder
        BraComb=0
        for Bra in Bra2NewCloneLs:
            NewCloneLs=Bra2NewCloneLs[Bra]
            if len(NewCloneLs)>1:     
                SpeClo=[]		
                for Clo in NewCloneLs:
                     Seq0=Clo2Seq['#'+Clo]
                
                     for Clo1 in NewCloneLs:
                          if Clo!=Clo1:
                               Seq1=Clo2Seq['#'+Clo1]
                               c=0
                               while c<Len:
                                   if Seq1[c]=='T' and Seq0[c]=='A': SpeClo.append(Clo1)
                                   elif Seq1[c]=='A' and Seq0[c]=='T': SpeClo.append(Clo)
                                   c+=1
                SpeClo=list(set(SpeClo))
                if len(SpeClo)>1:	
                    MaxSeq=self.MakeNewSeq(NewCloneLs,Clo2Seq,'Max',Len,self.tumor2seq)
                    MinSeq=self.MakeNewSeq(NewCloneLs,Clo2Seq,'Min',Len,self.tumor2seq)			
                    Clo2Seq['#BraCombMax'+str(BraComb)]=MaxSeq
                    Clo2Seq['#BraCombMin'+str(BraComb)]=MinSeq
                    CombCloLs+=NewCloneLs
                    NewCloLs+=['BraCombMax'+str(BraComb),'BraCombMin'+str(BraComb)]
                    for Old in NewCloneLs:			
                         Old2NewName[Old]=['BraCombMax'+str(BraComb),'BraCombMin'+str(BraComb)]
                    BraComb+=1
                else: NewCloLs+=NewCloneLs
            else: NewCloLs+=NewCloneLs		
        
        Clo2Seq1={}
        CloLs1=[]
        CloLs11=[]
        Len=len(Clo2Seq[CloOrder[0]])        
        MinNum=self.ErroRate*Len	 
        for Clo in NewCloLs:
          if Clo[0]!='#':Clo='#'+Clo
          Seq=Clo2Seq[Clo]	
          Muts=Align.GetMutPos(Seq)
          MutC=len(Muts)###check if these are not normal
          if MinNum<MutC:
             Clo2Seq1[Clo]=Seq
             CloLs1.append(Clo)	  
             CloLs11.append(Clo[1:])
        AddBraClo1=Align.UpMeg(Clo2Seq1,CloLs1)#Functions.UpMeg(Clo2Seq1,CloLs1,'AA','AddBraClo1.meg') 	
        NewT2C2F1={}
        for Tu in NewT2C2F:
             C2F=NewT2C2F[Tu]
             NewC2F={}
             for C in C2F:
               if C2F[C]>0:	 
                Code=C in Old2NewName
                	
                if Code==True :
                   NewNames=Old2NewName[C]
                   for New in NewNames:
                        Code1='#'+New in CloLs1			   
                        if Code1==True: NewC2F[New]=1
                else:
                        Code1='#'+C in CloLs1		   
                        if Code1==True: 		
                             NewC2F[C]=C2F[C]
             NewT2C2F1[Tu]=NewC2F		   
        #
        AddBraClo1_final, outSeqAmb, NewT2C2F1_final= Align.CombSimClo(AddBraClo1,NewT2C2F1, 0)	
        AddBraClo1_final+=['#hg19','A'*Len]		
        return AddBraClo1_final,NewT2C2F1_final 		
							
    def LineageSort(self, Clades):
        SizeLs=[]
        for Cla in Clades:
            Size = len(Cla)
            if SizeLs.count(Size)==0: SizeLs.append(Size)
        SizeLs.sort(reverse=True)
        Sorted=[]
        for Size in SizeLs:		
            for Cla in Clades:
                if len(Cla) == Size: Sorted.append(Cla)
        return Sorted				
 
    def CheClu1(self, AncCloMut,CluMut):  
        TMP=[]
        C=0   
        for Cmut in CluMut:
           Code=Cmut in AncCloMut
           if Code!=True : TMP.append(Cmut)
           else:
             C+=1 
        if C>0: Good='y'
        else: Good='n'	
        return Good,TMP	
    	
    def RedunChe(self, Seq,N2Seq,NameLs,Min):
        Len=len(Seq)
        Same='n'	
        for N in NameLs:	
           Seq0=N2Seq['#'+N]
           c=0
           Dif='n'
           DifC=0	   
           while c<Len:
               if Seq[c]!=Seq0[c]: 
                     Dif='y'
                     DifC+=1				 
               c+=1 
           if DifC<=Min: Same='y'	   
        return Same		
    	
    def ExtraHit(self, Meg,CloFre):
        CloOrder, Clo2Seq, out2=Functions.ReadMegSeq(Meg)
        Tu2HitCloLs,N2C,T2C2F=Functions.GetCloHitForTu(CloFre,0)
        NewCloLs=[]
        for Tu in Tu2HitCloLs:
           NewCloLs+=Tu2HitCloLs[Tu]
        NewCloLs=list(set(NewCloLs))
        Functions.UpCloFreqTa1('',T2C2F,NewCloLs,CloFre[:-4]+'1.txt')
        Functions.UpMeg(Clo2Seq,NewCloLs,'AA',Meg[:-4]+'1.meg')	
    
    def MakeNewSeq(self, CloneLs,SeqDic,Type,Len,Tu2Seq):
        SelSeq={}
        TuSeqDic={}	
        for Clo in CloneLs: 
             SelSeq['#'+Clo]=SeqDic['#'+Clo]
             Tu=Clo.split('BraC')[0]		 
             TuSeqDic['#'+Tu]=Tu2Seq['#'+Tu]
        TuMut=Functions.GetSharePosi1(TuSeqDic,'T')
        if Type=='Min':	
            MinSeqT=Functions.GetSharePosi1(SelSeq,'T')
            Nseq=''
            c=0
            while c<Len:
               Code=c in MinSeqT
               Code1=c in TuMut		   
               if Code==True and Code1==True: Nseq+='T'
               else: Nseq+='A'
               c+=1		   
        elif Type=='Max':		
            MaxSeqA=Functions.GetSharePosi1(SelSeq,'A')
            Nseq=''
            c=0
            while c<Len:
               Code=c in MaxSeqA
               Code1=c in TuMut		   
               if Code==True: Nseq+='A'
               elif Code1!=True: Nseq+='A'		   
               else: Nseq+='T'
               c+=1
        return Nseq	
    
    def AddTuExtra(self, Alignment, clone_frequency):	
        Align = MegaAlignment()	
        CloOrder, Clo2Seq = Align.name2seq(Alignment)	
        NewT2C2F={}
        AllCloLs=[]
        
        for Tu in clone_frequency:
		
            CloLs=clone_frequency[Tu]
            CloSeqDic={}
            Clo2Chan={}
            print Tu, CloLs			
            for Clo in CloLs:
              if CloLs[Clo]>0:			
                CloSeqDic[Clo]=Clo2Seq['#'+Clo]
                Clo2Chan[Clo]=[]
            print CloSeqDic				
            WildPosi=Align.GetSharePosi1(CloSeqDic,'A')		
            TuSeq=self.tumor2seq['#'+Tu[2:]]
            UniPos=self.GetUni('#'+Tu,self.tumor2seq)	
        
            TuExtPos=[]
            for c in WildPosi:
                Code=c in UniPos
                if Code==True: TuExtPos.append(c)		      
          #  print Tu,TuExtPos        		
            if TuExtPos==[]:
               NewT2C2F[Tu]=clone_frequency[Tu]
               AllCloLs+=CloLs	   
            else:
               NewT2C2F[Tu]={}	
               SNVfreLs=self.Tu2SNV[Tu]
           
               for Pos in TuExtPos:
                     Obs=float(SNVfreLs[Pos])*2
                     Sim=''
                     Small=100
                     for Clo in CloLs:
                         Est=clone_frequency[Tu][Clo]
                         Dif=((Obs-Est)**2)**0.5
                         if Dif<Small: 
                              Small=Dif
                              Sim=Clo			 
                     Clo2Chan[Sim].append(Pos)
        		 
               for C in Clo2Chan:
                  Chan=Clo2Chan[C]
                  if Chan==[]:
                       NewT2C2F[Tu][C]=clone_frequency[Tu][C]
                       AllCloLs.append(C)
                  else:
                     OriSeq=Clo2Seq['#'+C]
                     NewSeq=Align.ModSeq(OriSeq,Chan,'T',Len)						
                     NewName=C+'Add'+Tu
        			 ######
                     Clo2Seq['#'+NewName]=NewSeq			 
                     AllCloLs.append(NewName)
                     NewT2C2F[Tu][NewName]=clone_frequency[Tu][C]
                   
        AllCloLs=list(set(AllCloLs))
        NewSeqListOut = Align.UpMeg(Clo2Seq,AllCloLs)
        return 	NewSeqListOut,NewT2C2F 	

    def GetUni(self, TuN,Dic):
        Align = MegaAlignment()		
        NewDic={}
        for i in Dic:
           if i!=TuN: NewDic[i]=Dic[i]	   
        SharePosi=Align.GetSharePosi1(NewDic,'A')
        TuSeq=Dic['#'+TuN[3:]]
        TMP=[]
        for P in SharePosi:
            if TuSeq[P]=='T': TMP.append(P)
        return TMP		
    


		

	