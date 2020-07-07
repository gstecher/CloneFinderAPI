
from alignments.MegaAlignment import MegaAlignment
import scipy
import os
import sys

class cluster_test(): 
    def remove_insignificant_clones(self, v_obs, CloFre_clone, clone_seq_builder, Tu2CNV, Cut):
         Align=MegaAlignment()	
         OutAncAll='SigTest.txt'		 
         outAncAll='tumor\tDecsendant-Ancestor\tSNV posi\tType\tObsFre\n'		 
         Clone_list, clone_seq_dic = Align.name2seq(clone_seq_builder)	
         new_clone_freq={}
         new_clone_seq_dic={}		 
         for tumor in v_obs:
             CNV=Tu2CNV[tumor]
             Clo2Fre=CloFre_clone['T-'+tumor]	
             ObsFre = v_obs[tumor]

             clone_order=[]
             MutNum2Clo={}
             MutNum_ls=[]			 
             for Clo in Clo2Fre:
                 if Clo2Fre[Clo]>0: 
                      MutPosLs=Align.GetMutPos(clone_seq_dic['#'+Clo])		
                      MutNum=len(MutPosLs)
                      if (MutNum in MutNum2Clo)!=True: MutNum2Clo[MutNum]=[]					  
                      MutNum2Clo[MutNum].append(Clo)
                      MutNum_ls.append(MutNum)
             MutNum_ls=list(set(MutNum_ls))
             MutNum_ls.sort(reverse=True)			 
             for MutNum in MutNum_ls:
                    			 
                      clone_order+=MutNum2Clo[MutNum]
          				  
             CloNum=len(clone_order)
             C1Max=CloNum-1
             InsigLs=[]
         	 
             C1=0	
             while C1<C1Max:
                 Clo1=clone_seq_dic['#'+clone_order[C1]]
                 num_sites=len(Clo1)
                 Min_num=0.01*num_sites				 
                 C2=C1+1
                 while C2<CloNum:
                     Clo2=clone_seq_dic['#'+clone_order[C2]]
           
                     Share=[]
                     Unique=[]
                     c=0
                     while c< num_sites:
                        if CNV[c]=='normal':			   
                          if Clo1[c]=='T' and Clo2[c]=='T':
                                Share.append(ObsFre[c])
                                outAncAll+=tumor+'\t'+clone_order[C1]+'-'+clone_order[C2]+'\t'+str(c)+'\tShare\t'+str(ObsFre[c])+'\n'					   
                          elif Clo1[c]=='T' and Clo2[c]=='A':
                                Unique.append(ObsFre[c])
                                outAncAll+=tumor+'\t'+clone_order[C1]+'-'+clone_order[C2]+'\t'+str(c)+'\tUnique\t'+str(ObsFre[c])+'\n'	
                     							
                        c+=1
                     if (len(Share)<3 or len(Unique)<3) or (len(Share)<Min_num or len(Unique)<Min_num): P=1						
                     else: 
                          P=scipy.stats.ttest_ind(Share,Unique, equal_var = False)
                				  
                          P=P[-1]						  
                     if P>Cut:
                          if clone_order[C1].find('Clu')!=-1 and clone_order[C2].find('Clu')==-1:  InsigLs.append(clone_order[C1])					 
                          else: InsigLs.append(clone_order[C2])					 
  
                     C2+=1	
		
                 C1+=1
             InsigLs=list(set(InsigLs))				 
             if InsigLs!=[]: print('insignificant clones', tumor, InsigLs)
             new_clone_fre_in={}			 
             for Clo in Clo2Fre:
                 if Clo2Fre[Clo]>0 and InsigLs.count(Clo)==0:
                      new_clone_fre_in[Clo]=Clo2Fre[Clo]
                      new_clone_seq_dic['#'+Clo]=clone_seq_dic['#'+Clo]					  
             new_clone_freq['T-'+tumor]= new_clone_fre_in                     				 
         new_seq_builder=Align.UpMeg(new_clone_seq_dic,[])				 
 
         return new_seq_builder, new_clone_freq		 

    def remove_insignificant_clones_add(self, v_obs, CloFre_clone, clone_seq_builder, Tu2CNV, Cut):
         Align=MegaAlignment()	
         OutAncAll='SigTest.txt'		 
         outAncAll='tumor\tDecsendant-Ancestor\tSNV posi\tType\tObsFre\n'		 
         Clone_list, clone_seq_dic = Align.name2seq(clone_seq_builder)	
         new_clone_freq={}
         new_clone_seq_dic={}	
		 
         for tumor in v_obs:
             CNV=Tu2CNV[tumor]
             Clo2Fre=CloFre_clone['T-'+tumor]	
             ObsFre = v_obs[tumor]
             add_clone_freq={}
             clone_order=[]
             MutNum2Clo={}
             MutNum_ls=[]			 
             for Clo in Clo2Fre:
                 if Clo2Fre[Clo]>0: 
                      MutPosLs=Align.GetMutPos(clone_seq_dic['#'+Clo])		
                      MutNum=len(MutPosLs)
                      if (MutNum in MutNum2Clo)!=True: MutNum2Clo[MutNum]=[]					  
                      MutNum2Clo[MutNum].append(Clo)
                      MutNum_ls.append(MutNum)
             MutNum_ls=list(set(MutNum_ls))
             MutNum_ls.sort(reverse=True)			 
             for MutNum in MutNum_ls:
                    			 
                      clone_order+=MutNum2Clo[MutNum]
          				  
             CloNum=len(clone_order)
             C1Max=CloNum-1
             InsigLs=[]
             add_clone_freq[tumor]=[]        	 
             C1=0	
             while C1<C1Max:
                 Clo1=clone_seq_dic['#'+clone_order[C1]]
                 num_sites=len(Clo1)
                 Min_num=0.01*num_sites				 
                 C2=C1+1
                 while C2<CloNum:
                     Clo2=clone_seq_dic['#'+clone_order[C2]]
           
                     Share=[]
                     Unique=[]
                     c=0
                     while c< num_sites:
                        if CNV[c]=='normal':			   
                          if Clo1[c]=='T' and Clo2[c]=='T':
                                Share.append(ObsFre[c])
                                outAncAll+=tumor+'\t'+clone_order[C1]+'-'+clone_order[C2]+'\t'+str(c)+'\tShare\t'+str(ObsFre[c])+'\n'					   
                          elif Clo1[c]=='T' and Clo2[c]=='A':
                                Unique.append(ObsFre[c])
                                outAncAll+=tumor+'\t'+clone_order[C1]+'-'+clone_order[C2]+'\t'+str(c)+'\tUnique\t'+str(ObsFre[c])+'\n'	
                     							
                        c+=1
                     if (len(Share)<3 or len(Unique)<3) or (len(Share)<Min_num or len(Unique)<Min_num): P=1						
                     else: 
                          P=scipy.stats.ttest_ind(Share,Unique, equal_var = False)
                				  
                          P=P[-1]						  
                     if P>Cut:
                          if clone_order[C1].find('Clu')!=-1 and clone_order[C2].find('Clu')==-1:  
                                 InsigLs.append(clone_order[clone_order[C1]])

                                  								 
                          else: 
                                 InsigLs.append(clone_order[C2])						  
			  					 
                             
                     C2+=1	
		
                 C1+=1
             InsigLs=list(set(InsigLs))				 
             if InsigLs!=[]: print('insignificant clones', tumor, InsigLs)
             new_clone_fre_in={}			 
             for Clo in Clo2Fre:
                 New_cloe_hit=[]			 
                 if Clo2Fre[Clo]>0 and InsigLs.count(Clo)==0:
                      if (Clo in new_clone_fre_in)!=True: new_clone_fre_in[Clo]=0				 
                      new_clone_fre_in[Clo]+=Clo2Fre[Clo]
    	  
                      new_clone_seq_dic['#'+Clo]=clone_seq_dic['#'+Clo]	
                      New_cloe_hit.append(Clo)
             print(tumor, InsigLs, new_clone_fre_in)					  
             for InsigClo in InsigLs:
                        Index=clone_order.index(InsigClo)
                        Cont='y'
                        while Cont=='y':
                            NextClo=clone_order[Index-1]
                            if InsigLs.count(NextClo)==0 and Clo2Fre[NextClo]>0: 
                                   Cont='n'
                            else: Index=Index-1	
                        print(NextClo,InsigClo,Clo2Fre[InsigClo])						
                        if (NextClo in new_clone_fre_in)!=True: new_clone_fre_in[NextClo]=0	
                        CurTot=	new_clone_fre_in[NextClo]+Clo2Fre[InsigClo]					
                        if CurTot<=1: new_clone_fre_in[NextClo]+=	Clo2Fre[InsigClo]
                        else: new_clone_fre_in[NextClo]=1.0						
             new_clone_freq[tumor]= new_clone_fre_in  
             print(new_clone_fre_in)
         print(new_clone_freq)			 
         new_seq_builder=Align.UpMeg(new_clone_seq_dic,[])				 
 
         return new_seq_builder, new_clone_freq		
