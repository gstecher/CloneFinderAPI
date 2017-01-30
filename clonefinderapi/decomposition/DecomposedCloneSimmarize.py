from alignments.MegaAlignment import MegaAlignment
from parsimony.TreeAnalizer import TreeAnalizer
from parsimony.MegaMP import MegaMP

class DecomposedCloneSimmarize:

    def remove_tumor_and_rename_decomposed(self, tumor2seqs_with_decompose, seqs_with_ancestor, tumor_seqs, REP, clone_frequency_for_seqs_with_ancestor): 
                        Align = MegaAlignment()
                        SeqOrderIni, IniMeg2Seq= Align.name2seq(seqs_with_ancestor)
                        TuLs, TuMeg2Seq= Align.name2seq(tumor_seqs)							
                        SNVNum=len(TuMeg2Seq[TuLs[0]])							
                        IdenLs = Align.identify_similar_seq(tumor_seqs, 0)
                        Tu2IdenTu = Align.make_similar_seq_dic(IdenLs)						
                        outAllSeq=['MEGA','!Title SNVs;','!Format datatype=dna;',' '] 
                        RmCloLs=[]						
                        for Tu in tumor2seqs_with_decompose:
                               DeComCloLs=	tumor2seqs_with_decompose[Tu]
                               if DeComCloLs!=[]:
                                    IdenTu=Tu2IdenTu['T-'+Tu]
                                    RmCloLs+=IdenTu
                        RmCloLs=list(set(RmCloLs))
                        Done=[]						
                        for Tu in tumor2seqs_with_decompose:						
                                    DeComCloLs=	tumor2seqs_with_decompose[Tu]								   
                                    if RmCloLs.count(Tu)==0:
                                          if Done.count('#'+Tu)==0:									
                                              outAllSeq+=['#'+Tu,TuMeg2Seq['#'+Tu]]
                                              Done.append('#'+Tu)											  
                                    if DeComCloLs!=[]:           
                                                DecomCloOrder, Clu2Seq = Align.name2seq(DeComCloLs)
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
                                 									
                                                HitCloLs=clone_frequency_for_seqs_with_ancestor['T-'+Tu]
                                                for Clo in HitCloLs:
                                                      if HitCloLs[Clo]>0:												
                                                            TuClo=Clo.split('Clu')[0]
                                                            Code=TuClo in RmCloLs
                                                            if Code!=True and Clo[:4]!='Node':
                                                                if Done.count('#'+Clo)==0:																
                                                                   outAllSeq+=['#'+Clo,IniMeg2Seq['#'+Clo]]
                                                                   Done.append('#'+Clo)																   

                        outAllSeq_without_redindant = Align.RmRedunSeq(outAllSeq)
                        outAllSeq_without_redindant+=['#hg19',('A'*SNVNum)] 										
                        return outAllSeq_without_redindant						

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
					
                        if 	name2.find('Clu')!=-1:
                             if Additional_mut_num1==0: RmCluClo.append(name1)
                        else: 
                             if Der<Error_rate: RmCluClo.append(name1)
					 
        for Name in SeqOrderIni:
             if RmCluClo.count(Name)==0:
                  good_seq+=[Name,Meg2Seq[Name]]
        good_seq	+=['#hg19','A'*len(seq1)]				  
        return good_seq				  
            		
    def fix_back_para_mut_decomposed(self, remove_ancestral_decomposed_seq):                  				  
                        tree_builder=MegaMP()
                        Align = MegaAlignment()									   
                        status = tree_builder.do_mega_mp(remove_ancestral_decomposed_seq, 'decomposed_seq')
                        if status == True:
                            decom_seqs_with_ancestor, decom_tree, decom_nade_map, decom_mask_seq, decom_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)
                        else:
                            print 'failed to run megaMP'						
                        adjusted_decomposed_clone_seq = self.AdjDecClo('Final_NoRedun',decom_seqs_with_ancestor, decom_nade_map, decom_Good_posi_info)					
                        return adjusted_decomposed_clone_seq 						

    def AdjDecClo(self, ID, SeqLs, NodeMap, BackFor):
        Align = MegaAlignment()	
        Ao,Anc2Seq = Align.name2seq(SeqLs)
        Len=len(Anc2Seq[Ao[0]])		
        outNew = ['MEGA','!Title SNVs;','!Format datatype=dna;',' ']		
	
        Clu2Change={}
        for Name in Anc2Seq:
             if Name.find('Clu')!=-1: 
                  Clu2Change[Name]={}
        Posi=0 	
        while Posi<Len:		
            i=BackFor[Posi]	   
            if i!=['Good']:
                Change=i[0]
                ChangeCloLs=i[1][0]	#list
                if 	len(i[1][0])==1 and len(i[1][1])==1 and i[1][0][0].find('Clu')==-1 and i[1][1][0].find('Clu')!=-1:
                    ChangeCloLs=i[1][1]	#list				
			    		
                if Change=='ToMut':#BC>0: #back	
                      for Clo in ChangeCloLs:
                          if Clo.find('Clu')!=-1:
                               Clu2Change[Clo][Posi]='T'						  
                if Change=='ToWild':#MC>0: #multi		
                      for Clo in ChangeCloLs:
                          if Clo.find('Clu')!=-1:
                               Clu2Change[Clo][Posi]='A'					
				    
            Posi+=1

        for Clo in Anc2Seq:
            if Clo.find('#Node')==-1 : 
                if Clo.find('Clu')!=-1: 
                    Change=Clu2Change[Clo]
			
                    CluSeq=Anc2Seq[Clo]
                    Len=len(CluSeq)
                    c=0
                    NewSeq=''		
                    while c<Len:
                      Code=c in Change
                      if Code==True: NewSeq+=Change[c]
                      else: NewSeq+=CluSeq[c]
                      c+=1			 
                    outNew+=[Clo,NewSeq]			
                else: outNew+=[Clo,Anc2Seq[Clo]]
        outNew_without_redundant_seq = Align.RmRedunSeq(outNew)				
        outNew_without_redundant_seq+=['#hg19','A'*len(Anc2Seq[Clo])]  	 
        return outNew_without_redundant_seq	

    def find_decomposed_clone(self, no_back_para_mut_decomposed_seq, REP):
        Align = MegaAlignment()
        CloLs,Clo2Seq = Align.name2seq(no_back_para_mut_decomposed_seq)
        NewDecom='n'
        for Clo in CloLs:
            if Clo.find('Clu')!=-1:
                 ID='REP'+str(REP)			
                 In=-1*len(ID)			 
                 if Clo[In:]==ID: NewDecom='y'
        return NewDecom	

	 