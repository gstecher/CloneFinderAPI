from alignments.MegaAlignment import MegaAlignment
from parsimony.TreeAnalizer import TreeAnalizer
from parsimony.MegaMP import MegaMP
class DecomposedCloneSimmarize:
    def remove_tumor_and_rename_decomposed(self, tumor2seqs_with_decompose, seqs_with_ancestor, tumor_seqs, REP, clone_frequency_for_seqs_with_ancestor): 
                        Align = MegaAlignment()
                        SeqOrderIni, IniMeg2Seq= Align.name2seq(seqs_with_ancestor)
                        TuLs, TuMeg2Seq= Align.name2seq(tumor_seqs)	
                  #      print 'tuseq',TuMeg2Seq						
                        SNVNum=len(TuMeg2Seq[TuLs[0]])	
                      #  print TuMeg2Seq						
                        IdenLs = Align.identify_similar_seq(tumor_seqs, 0)
                        Tu2IdenTu = Align.make_similar_seq_dic(IdenLs)						
                       # print tumor2seqs_with_decompose						
                       # A, IniMeg2Seq, outAllSeq=Functions.ReadMegSeq('IniMegAT.meg')
                       # A, TuMeg2Seq, B=Functions.ReadMegSeq(filename+'_Seq.meg')
                       # Tu2HitCloLs,N2C,T2C2F=Functions.GetCloHitForTu(IniCloFreq,0)
                        outAllSeq=['MEGA','!Title SNVs;','!Format datatype=dna;',' '] 
                        RmCloLs=[]						
                        for Tu in tumor2seqs_with_decompose:
                               DeComCloLs=	tumor2seqs_with_decompose[Tu]
                               if DeComCloLs!=[]:
                                    IdenTu=Tu2IdenTu['T-'+Tu]
                                    RmCloLs+=IdenTu
                        RmCloLs=list(set(RmCloLs))
                        print RmCloLs
                        Done=[]						
                        for Tu in tumor2seqs_with_decompose:						
                                   # print Tu, tumor2seqs_with_decompose[Tu]						
                                   # Code=Tu in RmCloLs
                                    DeComCloLs=	tumor2seqs_with_decompose[Tu]								   
                                    if RmCloLs.count(Tu)==0:
                                          if Done.count('#'+Tu)==0:									
                                              outAllSeq+=['#'+Tu,TuMeg2Seq['#'+Tu]]
                                              Done.append('#'+Tu)											  
                                    if DeComCloLs!=[]:#os.path.exists(Tu+'Decom.meg')==True:            
                                                DecomCloOrder, Clu2Seq = Align.name2seq(DeComCloLs)#, out2=Functions.ReadMegSeq(Tu+'Decom.meg')
                                                for Clu in Clu2Seq:
                                                            Seq=Clu2Seq[Clu]
                                                            TuClu=Clu[1:].split('Clu')[0]
                                                            Code=TuClu in RmCloLs
                                                            if Code!=True and Clu.find('#Node')==-1: 
                                                                        if Clu.find('#Clu')!=-1: 
                                                                                                Clu='#'+Tu+Clu[1:]+'REP'+str(REP)
                                                                        if Done.count(Clu)==0:																								
                                                                            outAllSeq+=[Clu,Seq]
                                                                            Done.append(Clu)																		
                                    else:
                                 									
                                                HitCloLs=clone_frequency_for_seqs_with_ancestor['T-'+Tu]#Tu2HitCloLs[Tu]
                                                for Clo in HitCloLs:
                                                      if HitCloLs[Clo]>0:												
                                                            TuClo=Clo.split('Clu')[0]
                                                            Code=TuClo in RmCloLs
                                                            if Code!=True and Clo[:4]!='Node':
                                                                if Done.count('#'+Clo)==0:																
                                                                   outAllSeq+=['#'+Clo,IniMeg2Seq['#'+Clo]]
                                                                   Done.append('#'+Clo)																   
                      #  print 'before',outAllSeq															

                       # print 'after',outAllSeq							
                       # Functions.GetOut('Final.meg',outAllSeq)
                        outAllSeq_without_redindant = Align.RmRedunSeq(outAllSeq)
                        outAllSeq_without_redindant+=['#hg19',('A'*SNVNum)] 						
                      #  print outAllSeq_without_redindant						
                        return outAllSeq_without_redindant						
                        #os.system('python '+'RmRedunSeq.py '+'Final.meg')

    def remove_ancestral_decomposed(self, remove_tumor_and_rename_decomposed_seq, Error_rate):
        Align = MegaAlignment()
        SeqOrderIni, Meg2Seq= Align.name2seq(remove_tumor_and_rename_decomposed_seq)	
        good_seq = 	['#MEGA','!Title SNVs;','!Format datatype=dna;',' '] 
        RmCluClo=[]		
        for name1 in SeqOrderIni:
            if name1.find('Clu')!=-1:
                seq1=Meg2Seq[name1]
                for name2 in SeqOrderIni:
                    if name2!='#hg19' and name1!=name2:
                        seq2=Meg2Seq[name2]
                        Additional_mut_num1 = Align.CountAdditionalMut(seq1,seq2)
                        Der=1.0*Additional_mut_num1/len(seq1)
                     #   print name1,name2,Der						
                        if 	name2.find('Clu')!=-1:
                             if Additional_mut_num1==0: RmCluClo.append(name1)
                        else: 
                             if Der<Error_rate: RmCluClo.append(name1)
     #   print RmCluClo							 
        for Name in SeqOrderIni:
             if RmCluClo.count(Name)==0:
                  good_seq+=[Name,Meg2Seq[Name]]
        good_seq	+=['#hg19','A'*len(seq1)]				  
        return good_seq				  
            		
    def fix_back_para_mut_decomposed(self, remove_ancestral_decomposed_seq):                  				  
                        tree_builder=MegaMP()
                        Align = MegaAlignment()						
                     #   DecomClo=filename+'_Decom.meg'
                       # Functions.DoMEGAMPsimple('Final_NoRedun.meg')
                       # print remove_ancestral_decomposed_seq					   
                        status = tree_builder.do_mega_mp(remove_ancestral_decomposed_seq, 'decomposed_seq')
                        if status == True:
                            decom_seqs_with_ancestor, decom_tree, decom_nade_map, decom_mask_seq, decom_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)
                            print 'best alignment',decom_seqs_with_ancestor, decom_tree, decom_nade_map, decom_mask_seq, decom_Good_posi_info
                           # print decom_seqs_with_ancestor
                        else:
                            print 'failed to run megaMP'						
                        adjusted_decomposed_clone_seq = self.AdjDecClo('Final_NoRedun',decom_seqs_with_ancestor, decom_nade_map, decom_Good_posi_info)
                       # adjusted_decomposed_clone_without_redundant_seq =  Align.RmRedunSeq(adjusted_decomposed_clone_seq) 
                       # print 'return',adjusted_decomposed_clone_seq						
                        return adjusted_decomposed_clone_seq 						
                      #  os.system('python '+'RmRedunSeq.py '+DecomClo)
#                        
                      #  os.system('python '+'AddOutGroupHG.py ' + DecomClo[:-4]+'_NoRedun.meg')   
                      #  Functions.Reset(filename,DecomClo[:-4]+'_NoRedunWithOut.meg',IniCloFreq,'_InLoop',filename+'Decom.meg',IniCloFreq,clonetree_dir,dir)	
    def AdjDecClo(self, ID, SeqLs, NodeMap, BackFor):
        Align = MegaAlignment()
      #   = Align.RmUnresolvedPosi(SeqLs)		
       # os.system('python '+'RmUnresolvedPosi2.py '+ID+'.meg')	
        #Sid=open('BestTree.txt','r').readlines()[0]	
        Ao,Anc2Seq = Align.name2seq(SeqLs)
        Len=len(Anc2Seq[Ao[0]])		
        outNew = ['MEGA','!Title SNVs;','!Format datatype=dna;',' ']		
       # Ao,Anc2Seq,outNew=	Functions.ReadMegSeq(ID+'_Anc.meg')
        #T2D2A=NodeMap[0]
        #T2A2D=NodeMap[0]
        #T2I2C=NodeMap[0]
        #T2C2I=NodeMap[0]#Functions.ReadNodeMapMP1(ID)
        Dec2Anc=NodeMap[0]#T2D2A[Sid]
        Code2Name=NodeMap[1]#T2I2C[Sid]
        Anc2Dec=NodeMap[2]#T2A2D[Sid]
        Name2Code=NodeMap[3]#T2C2I[Sid]	
      #  print 'namecode',Name2Code		
        Clu2Change={}
      #  Clu2Index={}
        for Name in Anc2Seq:
             if Name.find('Clu')!=-1: 
               #   Clu2Index[Name]=Name2Code[Name]
                  Clu2Change[Name]={}
      #  BackFor=open(ID+'_resolveMultCount.txt','r').readlines()[1:]
        Posi=0 
        print BackFor		
        while Posi<Len:		
       # for i in BackFor:
            i=BackFor[Posi]	   
            if i!=['Good']:#i.find('Good')==-1:	
               # i=i.strip().split('\t')
                Change=i[0]
                ChangeCloLs=i[1]	#list			
            #    BC=	int(i[1])
            #    MC=int(i[0])			    		
                if Change=='ToMut':#BC>0: #back	
                      for Clo in ChangeCloLs:
                          if Clo.find('Clu')!=-1:
                               Clu2Change[Clo][Posi]='T'						  
                       #   for Clu in Clu2Index:
                        #    if Anc2Seq['#'+Clu][Posi]=='A':					  	
                         #         Anc=Dec2Anc[Clu2Index[Clu]]
                          #        AncS=Anc2Seq['#'+Code2Name[Anc].replace('VerID','0')][Posi]
                           #       if AncS=='T': Clu2Change[Clu][Posi]='T'
                if Change=='ToWild':#MC>0: #multi		
                      for Clo in ChangeCloLs:
                          if Clo.find('Clu')!=-1:
                               Clu2Change[Clo][Posi]='A'					
                     #     for Clu in Clu2Index:				  
                      #      if Anc2Seq['#'+Clu][Posi]=='T':					  
                       #           Anc=Dec2Anc[Clu2Index[Clu]]
                        #          Code=Anc in Dec2Anc					  
                         #         if Code!=True: DecS='A'
                                   
                          #        else:							   
                           #          Decs=Anc2Dec[Anc]
                            #         for Dec in Decs:
                             #           if Dec!=Clu2Index[Clu]:	DecS=Anc2Seq['#'+Code2Name[Dec].replace('VerID','0')][Posi]							  
                                     							 
                              #    AncS=Anc2Seq['#'+Code2Name[Anc].replace('VerID','0')][Posi]						  
                              #    if AncS=='A' or DecS=='A': Clu2Change[Clu][Posi]='A'					    
            Posi+=1
        print Clu2Change			
   #     for Clu in Clu2Change:
    #        Change=Clu2Change[Clu]	
     #       CluSeq=Anc2Seq['#'+Clu]
      #      Len=len(CluSeq)
       #     c=0
        #    NewSeq=''		
         #   while c<Len:
          #       Code=c in Change
           #      if Code==True: NewSeq+=Change[c]
            #     else: NewSeq+=CluSeq[c]
             #    c+=1			 
           # outNew+=['#'+Clu,NewSeq]#+'\n'
        for Clo in Anc2Seq:
            if Clo.find('#Node')==-1 : 
                if Clo.find('Clu')!=-1: 
                    Change=Clu2Change[Clo]
                   # if Clo[0]!='#': Name='#'+Clo
                   # else: Name=Clo					
                    CluSeq=Anc2Seq[Clo]
                    Len=len(CluSeq)
                    c=0
                    NewSeq=''		
                    while c<Len:
                      Code=c in Change
                      if Code==True: NewSeq+=Change[c]
                      else: NewSeq+=CluSeq[c]
                      c+=1			 
                    outNew+=[Clo,NewSeq]#+'\n'			
                else: outNew+=[Clo,Anc2Seq[Clo]]#+'\n'
        outNew_without_redundant_seq = Align.RmRedunSeq(outNew)				
        outNew_without_redundant_seq+=['#hg19','A'*len(Anc2Seq[Clo])]  	
      #  Functions.GetOut(Out,outNew)
      #  os.system('python '+'RmRedunSeq.py ' + Out)
      #  print outNew_without_redundant_seq	  
        return outNew_without_redundant_seq	

    def find_decomposed_clone(self, no_back_para_mut_decomposed_seq, REP):
        Align = MegaAlignment()
        CloLs,Clo2Seq = Align.name2seq(no_back_para_mut_decomposed_seq)
        NewDecom='n'
        for Clo in CloLs:
            if Clo.find('Clu')!=-1:
                 ID='REP'+str(REP)			
                 In=-1*len(ID)
                # print ID, In, Clo[In:]				 
                 if Clo[In:]==ID: NewDecom='y'
        return NewDecom	

	 