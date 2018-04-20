from alignments.MegaAlignment import MegaAlignment
from parsimony.TreeAnalizer import TreeAnalizer
from parsimony.MegaMP import MegaMP
#from PPcomputer.PredictCellGenotype import PredictCellGenotype
from tsp_profiles.tsp_information import tsp_information
class DecomposedCloneSimmarize:
    def add_back_CNVSNV(self, DecomTu2Seq_builder_sub, CNV_information, original_seqs_builder_all, original_Tumor2Clone_frequency,tsp_list): 
        all_tsp = tsp_information(tsp_list)    				
        v_obs = all_tsp.tumor2alt_frequency()		
        Seq_all_dic={}
        Align=MegaAlignment()	
        Original_clols, Original_clodic_all = Align.name2seq(original_seqs_builder_all)		
        for Tumor in DecomTu2Seq_builder_sub:
            Seq_builder_sub=	DecomTu2Seq_builder_sub[Tumor]
            if Seq_builder_sub!=[]:
              SNVfre_list=v_obs[Tumor]			
              CloLs,Clo2Seq=Align.name2seq(Seq_builder_sub)
              CNVinfo=	CNV_information[Tumor]
              Len=len(CNVinfo)
            #  print Tumor, Clo2Seq.keys()			  
              for Clo in Clo2Seq:
                  Seq_sub=Clo2Seq[Clo]
                  c_seq=0
                  c_all=0
                  Seq_all=''				  
                  while c_all<Len:
                       if CNVinfo[c_all]=='normal': 
                             Seq_all+=Seq_sub[c_seq]
                             c_seq+=1	
                       else:
                         if SNVfre_list[c_all]==0: 	Seq_all+='A'				   
                         else: Seq_all+='?'
                       c_all+=1	
                  if Original_clodic_all.has_key(Clo)==True: Seq_all_dic[Clo]=Original_clodic_all[Clo]					   
                  else: Seq_all_dic[Clo]=Seq_all 
            else:
               CloFre=original_Tumor2Clone_frequency['T-'+Tumor]
          		   
               for Clo in CloFre:
                    if CloFre[Clo]>0:
                         if Seq_all_dic.has_key('#'+Clo)!=True: Seq_all_dic['#'+Clo]=Original_clodic_all['#'+Clo]					
		
    
        decom_all_seq_builder=Align.UpMeg(Seq_all_dic,[])	
   
        return 	decom_all_seq_builder
	
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

    def remove_ancestral_decomposed(self, remove_tumor_and_rename_decomposed_seq, Error_rate, tumor_seqs):
      #  print 'Tu',tumor_seqs	
        Align = MegaAlignment()
        SeqOrderIni, Meg2Seq= Align.name2seq(remove_tumor_and_rename_decomposed_seq)	
        TuLs,Tu2Seq= Align.name2seq(tumor_seqs)		
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
        AddedTuLs=[]					 
        for Name in SeqOrderIni:
             if RmCluClo.count(Name)==0:
                  good_seq+=[Name,Meg2Seq[Name]]
                  AddedTuLs.append(Name.split('Clu')[0])
        for Tu in TuLs:
           if AddedTuLs.count(Tu)==0:good_seq+=[Tu,Tu2Seq[Tu]] 		
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
                        adjusted_decomposed_clone_seq = self.AdjDecClo('Final_NoRedun',decom_seqs_with_ancestor, decom_nade_map, decom_Good_posi_info, decom_tree)					
                        return adjusted_decomposed_clone_seq , decom_tree						

    def AdjDecClo(self, ID, SeqLs, NodeMap, BackFor, Tree):
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
              #  print i,Posi		
                A=i[1][0][0].split('Clu')[0]
                B=	i[1][1][0].split('Clu')[0]
               # print A,B				
                if 	(len(i[1][0])==1 and len(i[1][1])==1 and i[1][0][0].find('Clu')==-1 and i[1][1][0].find('Clu')!=-1):
                    ChangeCloLs=i[1][1]	#list				
                    if A==B: ChangeCloLs=i[1][0]	
                else:
                    ChangeCloLs=i[1][0]	#list	
                    if A==B: ChangeCloLs=i[1][1]					
                #print ChangeCloLs			    		
                if Change=='ToMut':#BC>0: #back	
                      for Clo in ChangeCloLs:
                         # if Clo.find('Clu')!=-1:
                               if Clu2Change.has_key(Clo)!=True: Clu2Change[Clo]={}						  
                               Clu2Change[Clo][Posi]='T'						  
                if Change=='ToWild':#MC>0: #multi		
                      for Clo in ChangeCloLs:
                         # if Clo.find('Clu')!=-1:
                               if Clu2Change.has_key(Clo)!=True: Clu2Change[Clo]={}								  
                               Clu2Change[Clo][Posi]='A'					
				    
            Posi+=1
     #   print Clu2Change			
   #     open('AAA','r').readlines()
        for Clo in Anc2Seq:

          		  
            if Clo.find('#Node')==-1 : 
                TreeAna=TreeAnalizer()
                BraLen=TreeAna.Get_branch_lenghth(Tree,Clo[1:])			
              #  if Clo.find('Clu')!=-1 and BraLen>=1: 
                if BraLen>=1 and Clu2Change.has_key(Clo):				
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
     #   open('AA','r').readlines()		
        return outNew_without_redundant_seq	

    def find_decomposed_clone(self, no_back_para_mut_decomposed_seq, REP, Tree):
        Align = MegaAlignment()
        CloLs,Clo2Seq = Align.name2seq(no_back_para_mut_decomposed_seq)
        DecTipLs=[]        
        DecLs=[]
        DecAncLs=[]	
        RmDecom=[]		
        for Clo in CloLs:
         #   ClosestAnc = Align.find_closest_anc(Clo,Clo2Seq)
          #  if ClosestAnc!='' : 
           #      DecAncLs.appned(Clo)		
            if Clo.find('Clu')!=-1:
                 ID='REP'+str(REP)			
                 In=-1*len(ID)			 
                 if Clo[In:]==ID:
                      DecLs.append(Clo)
                      Posi=Tree.find(Clo[1:]+':')+len(Clo)
                   #   print Tree[Posi]
                      Go='y'
                      BraLen=''					  
                      while Go=='y':
                         BraLen+=Tree[Posi]
                         if Tree[Posi]==',' or Tree[Posi]==')': Go='n'
                         Posi+=1	
                   #   print Clo,BraLen 
                      if float(BraLen[:-1])==0: DecAncLs.append(Clo) ######
                      else: DecTipLs.append(Clo)					  
                      					  
      #  print DecLs,DecAncLs
		
        if DecLs==[]: NewDecom='n'
        elif DecTipLs!=[]: 
             NewDecom='y'
             for Tip in DecTipLs:
                  TipSeq=Clo2Seq[Tip]
                  OriTu=Tip.split('Clu')[0]				  
                #  TipMutC				  
                  Anc='n'				  
                  for Clo in Clo2Seq:
                      if Clo!=Tip: # and OriTu!=Clo:
                          UniNum=Align.CountAdditionalMut(TipSeq,Clo2Seq[Clo])	
                          if UniNum==0: Anc='y'	
                  if Anc=='y': 	RmDecom.append(Tip)					  
        else: NewDecom='anc'		
     #      NewDecom='anc'
      #     for Dclo in DecLs:
       #        if DecAncLs.count(Dclo)==0: 	NewDecom='y'	   
			 
     #   print Clo2Seq.keys() 
        if RmDecom==[]:NewClo2Seq_buil=no_back_para_mut_decomposed_seq
        else:
            NewCloDic={}
            for Clo in Clo2Seq:
                if RmDecom.count(Clo)==0: 	NewCloDic[Clo]=Clo2Seq[Clo]
            NewClo2Seq_buil=Align.UpMeg(NewCloDic,[])				
        return NewDecom,RmDecom	,NewClo2Seq_buil
    def add_back_anc(self, Sub_seq_builder, All_seq_builder):
        Align=MegaAlignment()
        Ls,Sub=Align.name2seq(Sub_seq_builder)	
        Ls,All=Align.name2seq(All_seq_builder)	
        Clo2Seq={}
        for Clo in All:
            if Sub.has_key(Clo)==True: Clo2Seq[Clo]=Sub[Clo]
            else: Clo2Seq[Clo]=All[Clo]	
        Seq_Buil=Align.UpMeg(Clo2Seq,[])
        return Seq_Buil		
    def finalize_results(self, decomposed_seq_builder, decomposed_Tumor2Clone_frequency, origianl_seq_builder, original_Tumor2Clone_frequency, REP):
        Align=MegaAlignment()
        Ls,DecomSeqDic=Align.name2seq(decomposed_seq_builder)	
     #   print Ls		
        Ls,OriSeqDic=Align.name2seq(origianl_seq_builder)
        NewCloSeqDic={}
        NewCloFre={}
    #    print decomposed_Tumor2Clone_frequency,original_Tumor2Clone_frequency		
        for Tu in original_Tumor2Clone_frequency:
             if decomposed_Tumor2Clone_frequency.has_key(Tu)!=True: 
			 
                   CloFre=original_Tumor2Clone_frequency[Tu]
                   for Clo in CloFre:
                       if CloFre[Clo]>0:				   
                          NewCloSeqDic['#'+Clo]=OriSeqDic['#'+Clo]				   
             elif decomposed_Tumor2Clone_frequency[Tu]=={}: 
			 
                   CloFre=original_Tumor2Clone_frequency[Tu]
                   for Clo in CloFre:
                       if CloFre[Clo]>0:				   
                          NewCloSeqDic['#'+Clo]=OriSeqDic['#'+Clo]				   
                                                   
             else:
                 CloFre0=decomposed_Tumor2Clone_frequency[Tu]
                 CloFre={}
                 for Clo in CloFre0:
                   Fre=CloFre0[Clo]
                   if Fre>0:				   
                     if (Clo.find('Clu')!=-1 and Clo.find('REP')==-1) or Clo.find('REP'+str(REP-1))!=-1: 
                         CloFre[Clo+'REP'+str(REP)]=Fre
                         NewCloSeqDic['#'+Clo+'REP'+str(REP)]=DecomSeqDic['#'+Clo]						 
                     else: 
                          CloFre[Clo]=Fre
                          if OriSeqDic.has_key('#'+Clo)==True: NewCloSeqDic['#'+Clo]=OriSeqDic['#'+Clo]	
                          else: NewCloSeqDic['#'+Clo]=DecomSeqDic['#'+Clo]							  
             NewCloFre[Tu]=CloFre
        rename_seq_builder = Align.UpMeg(NewCloSeqDic, [])	
      #  open('AA','r').readlines()		
        return rename_seq_builder, NewCloFre 
    def find_new_clone(self,new_seq_buil,old_seq_buil):
       Align=MegaAlignment()
       Ls,New_dic=Align.name2seq(new_seq_buil)
       Ls,Old_dic=Align.name2seq(old_seq_buil)
    #   print 'old list',Old_dic.keys(),	 '\nnew list',New_dic.keys()  
       Iden='y'	   
       for Clo in New_dic:
           if Clo!='#hg19':	
              NewSeq=New_dic[Clo]
              Redun=Align.find_redundant(NewSeq,Old_dic)
              if Redun==[]: 
                    Iden='n'
                 #   print 'new seq',Clo					


       return Iden	