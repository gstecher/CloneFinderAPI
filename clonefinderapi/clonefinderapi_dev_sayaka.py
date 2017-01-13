from parsers.DefaultTSPParser import DefaultTSPParser
from parsers.DefaultCCFParser import DefaultCCFParser
from alignments.FreqToMegaSeq import FreqToMegaSeq
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
#from binomial.BinomialReplicater import BinomialReplicater
from binomial.SignificanceTester import SignificanceTester
from regression.CloneFrequencyComputer import CloneFrequencyComputer
from decomposition.SNPClusterGenerator import SNPClusterGenerator
from decomposition.SNPGroupCombiner import SNPGroupCombiner
from decomposition.DecomposedCloneSimmarize import DecomposedCloneSimmarize
import os
dir = os.getcwd()
mao_file = dir + '/infer_MP_nucleotide.mao'
print mao_file
BinoRepNum=50
significance_cutoff=0.95
ErroRate=0.01
CloFreCut=0.01
OptionA='On' #On or Off (Cut extra mutations from a clone) 
OptionB='5,3,On' #maximum_number_of_clusters_to_be_produced,Minumum_number_of_SNVs_per_cluster,Decomposed clone has trunc mutaion and is not an ancestral clone)
OptionC='On,5' #On or Off,Number_of_SNVs_per_cluster' (decompose clones based on clone phylogeny) 
OptionD='On,3' #On or Off,Number_of_SNVs_per_cluster (make ancestral clones at middle of a branche)
OptionE='On' #On or Off (Filter backward/paralle mutations)
OptionF='On' #On or Off (Add unassigned SNVs to clones)
OptionG='On' #On or Off (Flag bad data)
OptionH='On' #On or Off (Combine similar clones)
tree_builder = MegaMP()
tree_builder.mao_file = mao_file#'/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/infer_MP_nucleotide.mao'
AnalyzeTree = TreeAnalizer()
if __name__ == "__main__":
   # filename = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\ccf.txt'
   # parser = DefaultCCFParser()
    filename = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\EV005Per1Anc.txt'
    parser = DefaultTSPParser()
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
        tsp_list = parser.get_tumor_sample_profile_list()
        		
        original_align_builder = FreqToMegaSeq()
        original_align_builder.initialize(tsp_list, True) # pass in True to remove duplicates, False to keep duplicates
        tumor_seqs = original_align_builder.get_mega_allalignment()
        print 	'all',tumor_seqs	
        initial_seq = original_align_builder.get_mega_alignment()		
        print initial_seq
      #  clone_frequency = CloneFrequencyComputer(mega_seqs, tsp_list, CloFreCut)	
       # significant_clones = SignificanceTester(BinoRepNum, significance_cutoff, ErroRate, tsp_list, mega_seqs, clone_frequency)
        #significant_clone_in_tumor, significant_clone_list = significant_clones.count_detection()
      #  print 	significant_clone_in_tumor, significant_clone_list	
       # seqs, CloFre = clone_frequency.regress()		
        #clusters = SNPClusterGenerator(mega_seqs, tsp_list, CloFre, OptionB[0].strip(), CloFreCut)
     #   print clusters.cluster()	
        print 'decompose and add ancestor'
		
        id = 'mega_alignment' # id will be used internally for file names
        print original_align_builder.get_mega_alignment()			
        status = tree_builder.do_mega_mp(initial_seq, id)
        if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)
	#outseq, tree_list[best_SetID], best_outset[1], best_outset[2], best_outset[3] #NadeMapInfo, mask_seq, Good_posi_info]			
                print 'best alignment'
                print seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info
        else:
                print 'failed to run megaMP'		
		
        RepeatUnion='y'	
        REP = 0		
        while RepeatUnion=='y':	
       # mega_seqs = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/mega_alignment.meg'
         #   tree_builder = MegaSimpleMP()
         #   tree_builder.mao_file = mao_file#'/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/infer_MP_nucleotide.mao'
		
           # open('AAA','r').readlines()
            REP += 1
            REPs = str(REP)
            clone_frequency = CloneFrequencyComputer(seqs_with_ancestor, tsp_list, CloFreCut)	

            seqs, CloFre_seqs_with_ancestor = clone_frequency.regress()		
            clusters = SNPClusterGenerator(seqs_with_ancestor, tsp_list, CloFre_seqs_with_ancestor, OptionB[0].strip(), CloFreCut)
            tumor2cluster = clusters.cluster()
            decomposition = SNPGroupCombiner(tumor2cluster, seqs_with_ancestor, tumor_seqs, tsp_list, CloFreCut, OptionB[1].strip(), OptionB[2].strip())
            tumor2seqs_with_decompose, RmCloLs = decomposition.get_decomposed_seq()		

            print 'hybrid sample genotypes ',RmCloLs            
            if RmCloLs==[]:  
                        RepeatUnion='n'
                       
            else: 
                  DecomposedCloneSummary = DecomposedCloneSimmarize()
                  remove_tumor_and_rename_decomposed_seq = DecomposedCloneSummary.remove_tumor_and_rename_decomposed(tumor2seqs_with_decompose, seqs_with_ancestor, tumor_seqs, REP, CloFre_seqs_with_ancestor)                  				  
                  print 'original',remove_tumor_and_rename_decomposed_seq
                  decom_without_ancdecom_seq = DecomposedCloneSummary.remove_ancestral_decomposed(remove_tumor_and_rename_decomposed_seq, ErroRate)
                  print 'without ancestraldecom', decom_without_ancdecom_seq
				  ####

                  print 'fix multiple/backward mutations on decomposed clones'
                  no_back_para_mut_decomposed_seq = DecomposedCloneSummary.fix_back_para_mut_decomposed(decom_without_ancdecom_seq)                  				  
                  presence_of_new_decomposed_clone = DecomposedCloneSummary.find_decomposed_clone(no_back_para_mut_decomposed_seq, REP)
                  print presence_of_new_decomposed_clone
                  if presence_of_new_decomposed_clone=='n':
                         RepeatUnion='n'
                  else:
                       print 'decomposition is better?'
                       id = 'mega_alignment' # id will be used internally for file names		
                       status = tree_builder.do_mega_mp(no_back_para_mut_decomposed_seq, id)
                       if status == True:
                             decomseqs_with_ancestor, decom_tree, decom_nade_map, decom_mask_seq, decom_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)
	#outseq, tree_list[best_SetID], best_outset[1], best_outset[2], best_outset[3] #NadeMapInfo, mask_seq, Good_posi_info]			
                             print 'best alignment'
                             print decom_Good_posi_info
                       else:
                             print 'failed to run megaMP'
                       print Initial_Good_posi_info, decom_Good_posi_info							 
                       decomposition_is_better = AnalyzeTree.compare_good_posi_number(Initial_Good_posi_info, decom_Good_posi_info)
                       if decomposition_is_better=='n': 
                             RepeatUnion='n'
                       else:
                            print 'clones were decomposed, so repeat the decomposition'
                            seqs_with_ancestor=decomseqs_with_ancestor
                            tree=decom_tree
                            nade_map=decom_nade_map
                            mask_seq=decom_mask_seq
                            Initial_Good_posi_info=decom_Good_posi_info							
                            							
        print 'clone decomposition is complete!',REP,seqs_with_ancestor,CloFre_seqs_with_ancestor,Initial_Good_posi_info,mask_seq,
        print 'rename clone anmes'
     ##   os.system('python '+'ReNameCloFreMeg.py '+filename+'_CloneSig.meg '+filename+'_CloneSig.txt 0 number') 		
        number_of_good_positions = AnalyzeTree.count_good_posi(Initial_Good_posi_info)
        proportion_of_good_positions = 1.0*number_of_good_positions/len(Initial_Good_posi_info)
        if proportion_of_good_positions<0.5:
    #      #add  Functions.GetOut(Result_dir+'\\'+ID+'_summary.txt',Only str(proportion_of_good_positions)+' of SNV sites are not affected by backward/parallel mutations. Please confirm the tumor profiles.\nOptionA, C, and D were changed to Off, if these were On.\n')
            OptionA='Off'
            OptionC='Off'
            OptionD='Off' 
        if OptionA=='Off' and OptionC=='Off' and OptionD=='Off':
            print 'Do not refine clone collection'	
        else:
            print 'Refine clone collection'



#if OptionA=='On':
#  print 'Cut extra mutations in clones' 
#  IniMeg=filename+'_CloneSeqSig.meg'
#  IniFre=filename+'_CloneFreqSig.txt'
#  os.system('python '+'CutExtraMutClo2.py '+IniFre+' '+IniMeg+' '+filename+' '+ErroR) 
#  MegID='Mod_RandSeq.meg'
#  CloFre='Mod_CloneFreqComb.txt'
#  Functions.Reset(filename,MegID,CloFre,'_AdjExtraMutCloAll',filename+'_Best.meg',filename+'_Best.txt',clonetree_dir,dir)
#  os.system('python '+'ReNameCloFreMeg.py ' + filename+'_Best.meg '+ filename+'_Best.txt 0 number') 
#  Functions.Reset(filename,filename+'_CloneSeqSig.meg', filename+'_CloneFreqSig.txt','_CutExtraInClone',filename+'_CloneSig.meg',filename+'_CloneSig.txt',clonetree_dir,dir)
#
#if OptionA=='On' or OptionC!='Off': 
#  print 'Resolve multiple mutations'
#  IniMeg=filename+'_CloneSig.meg'
#  IniFre=filename+'_CloneSig.txt'
#  os.system('python '+'BranchDecCloneInfer7-Short.py '+IniMeg+' '+IniFre+' '+filename+' '+SLEP_dir+' '+dir+' '+ErroR+' '+OptionA+' '+OptionC ) 
# # open('AAA,'r').readlines()
#  Functions.Reset(filename,'New_RandSeq.meg','New_CloneFreqComb.txt','_AddBraDecCloOldName',filename+'_A.meg',filename+'_A.txt',clonetree_dir,dir)
#  os.system('python '+'ReNameCloFreMeg.py ' + filename+'_A.meg '+filename+'_A.txt 0 number') 
#  Meg=filename+'_CloneSeqSig.meg'
#  Fre=filename+'_CloneFreqSig.txt'
#else:
#  Meg=filename+'_CloneSeqSig.meg'
#  Fre=filename+'_CloneFreqSig.txt'
#  copy(filename+'_CloneSig.meg',Meg)  
#  copy(filename+'_CloneSig.txt',Fre)
# 
#if OptionD=='Off': Functions.Reset(filename,Meg,Fre,'_AddMidBraClo',filename+'_CloneSeqSig.meg',filename+'_CloneFreqSig.txt',clonetree_dir,dir)
#else:
#  Functions.Reset(filename,Meg,Fre,'_AddBraDecClo',filename+'_BraDecClo.meg',filename+'_BraDecClo_CloneFreq.txt',clonetree_dir,dir)
#  print 'add ancestral clone at the middle of branches'
#  IniMeg=filename+'_BraDecClo.meg'
#  IniFre=filename+'_BraDecClo_CloneFreq.txt'
#  PreMegFreID=filename+'_AddBraDecClo'
#  os.system('python '+'BranchCloneInfer5.py '+IniMeg+' '+IniFre+' '+filename+' '+SLEP_dir+' '+dir+' ' +ErroR+' '+OptionD)
#  MegID=filename+'_MidBraAnc'
#  os.system('python '+'CombSimClo.py '+'SigBraClo.meg '+'SigBraClo.txt 0.0 0')  
#  Functions.Reset(filename,'SigBraClo_RandSeq.meg','SigBraCloComb.txt','_AddMidBraClo',filename+'_CloneSeqSig.meg',filename+'_CloneFreqSig.txt',clonetree_dir,dir)
  			
			
			
			
			
			
#print 'Statistics'
#MegID=dir+'\\'+ID+'_AddMidBraClo'
#os.chdir(clonetree_dir)
#c=1
#while c<=BinoNum:
#    File=ID+'Rep'+str(c)+'.txt'
#    copy(dir+'\\'+File,File)	
#    copy(MegID+'.meg','Seq.meg')
#    os.system('python '+'DoSLEPNoBoo.py ' + File+' '+'Seq.meg'+ ' '+SLEP_dir+' For_each_replicate')	
#    Functions.Reset(File[:-4],'Seq.meg',File[:-4]+'_CloneFreq.txt','_Recom','CloneSig.meg','CloneSig.txt',clonetree_dir,dir)
#    c+=1
#os.system('python '+'BinoSampAnaEachTuCloHet.py '+MegID+'.txt '+dir+'\\'+ID+'Rep '+str(BinoNum))
#Tu2Clo,Name2Col,T2C2F=Functions.GetCloHitForTu(MegID+'_Cou.txt',BinoCut)
#OriTu2Clo,OriName2Col,OriT2C2F=Functions.GetCloHitForTu(MegID+'_Cou.txt',0)
#for Tu in Tu2Clo:
#    if Tu2Clo[Tu]==[]:
#      Lar=0
#      LarClo=[]
#      OriCloLs=OriTu2Clo[Tu]
#      for Clo in OriCloLs:
#         Fre=OriT2C2F[Tu][Clo]
#         if Fre>Lar:
#             Lar=Fre
#             LarClo=[Clo]
#         elif Fre==Lar: LarClo.append(Clo)
#      Tu2Clo[Tu]=LarClo
#      for Clo in LarClo:
#         T2C2F[Tu][Clo]=Lar	  
#      	  
#NameOrder, Name2Seq, out2	=Functions.ReadMegSeq(MegID+'.meg')
#NameLs=[]
#Tu2Clo1={}
#for Tu in Tu2Clo:
#    Add=Tu2Clo[Tu]
#    if Add==[]:
#        C2F=T2C2F[Tu]	
#        Order,F2Cls=Functions.Sort(NameOrder,C2F)
#        LarClo=Order[0]
#        LarFre=C2F[LarClo]
#        Add=F2Cls[LarFre]				   
#    NameLs+=Add	
#    Tu2Clo1[Tu]=Add	
#NameLs=list(set(NameLs))
#UpCloFreqTa2('AA',T2C2F,Tu2Clo1,'Last0.txt')
#InFile='Last.meg'
#Functions.UpMeg(Name2Seq,NameLs,'AA',InFile)
#copy(dir+'\\'+ID+'_Seq.meg',ID+'_Seq.meg')
#copy(dir+'\\'+ID+'.txt',ID+'.txt')
#os.system('python '+'ReComCloFreForEachTu.py ' + ID+'.txt '+SLEP_dir+' '+clonetree_dir+' '+clonetree_dir+'\\ '+'Last0.txt '+'Last.meg '+ID+'_Seq.meg'+' 0') 
#copy('Last.meg_CloneFreqSig.txt','Last.txt')
#
#MegID=ID+'_clone'
#if OptionE=='Off' and OptionF=='Off':Functions.Reset(ID,InFile,'Last.txt','_AfterBino',MegID+'.meg',MegID+'.txt',clonetree_dir,dir)
#else: 
# if OptionE=='On':
#     print 'Resolve multiple/backward mutations (these mutations are assigned as wild type)'
#
#     os.system('python '+'AddOutGroupHG.py ' + InFile)   
#     InFile=InFile[:-4]+'WithOut.meg'
#     Functions.DoMEGAMPsimple(InFile)
#     os.system('python '+'RmUnresolvedPosi2.py '+InFile) 
#     copy(InFile[:-4]+'_resolveMultCount.txt', dir+'\\Results\\'+ID+'_resolveMultCount.txt')
#     os.system('python '+'CombSimClo.py '+InFile[:-4]+'_resolve.meg '+'Last.txt 0.0 0')  
#     MegID=ID+'_Final'
#   #  open('AAA','r').readlines()	 
#     Functions.Reset(ID,InFile[:-4]+'_resolve_RandSeq.meg','LastComb.txt','_Resolve',MegID+'.meg',MegID+'.txt',clonetree_dir,dir)
# else:
#     MegID=ID+'_Final' 
#     copy('Last.txt',MegID+'.txt')
#     copy(InFile,MegID+'.meg')
# if OptionF=='On':	 
#     print 'add missing mutations in clone'
#     InFileMeg=MegID+'.meg'
#     InFileCloFre=MegID+'.txt'
#     os.system('python '+'AddTuExtra2.py '+InFileCloFre+' '+InFileMeg+' '+ID)
#     MegID=ID+'_clone'
#     Functions.Reset(ID,ID+'_TruAdd.meg',ID+'_TruAdd_CloneFreq.txt','_AdjExtraMutTu',MegID+'.meg',MegID+'.txt',clonetree_dir,dir)
# else:
#     MegID=ID+'_clone'
#     copy(ID+'_Final.meg',MegID+'.meg')
#     copy(ID+'_Final.txt',MegID+'.txt')	 
#
#if Data=='cnv-post':
#   print 'examine if ancestral clones were identified due to CNV'
#   os.system('python '+'ReComCloFreForEachTuCNV.py ' + ID+'.txt '+SLEP_dir+' '+clonetree_dir+' '+clonetree_dir+'\\ '+MegID+'.txt '+MegID+'.meg '+ID+'_Seq.meg'+' '+str(FreqCut)+' '+dir+'\\'+CNV) 
#   copy(ID+'_CloneFreqSig.txt',MegID+'.txt')
#   copy(ID+'_CloneSeqSig.meg',MegID+'.meg')
##open('AAA','r').readlines()   
#print 'Count missing/extra mutations'
#
#copy(MegID+'.txt', ID+'_CloneFreq.txt')
#copy(MegID+'.meg',ID+'_FinalClones.meg')
#os.system('python '+'CloSeqComForTu.py ' + ID+' 0') 
#copy(ID+'ExtraMissCount.txt', dir+'\\Results\\'+ID+'_ExtraMissCountFinal.txt')



#os.system('python '+'ReNameCloFreMeg.py ' + MegID+'.meg '+MegID+'.txt'+ ' 0 list') 
#MegID=ID+'_CloneSeqSig'
#MegIDFre=ID+'_CloneFreqSig'
#MegIDlast=ID+'_CloneSeqSig1'
#if OptionH=='Off':
#  copy(MegIDFre+'.txt', dir+'\\Results\\'+ID+'_Final.txt') 
#  os.system('python '+'AddOutGroupHG.py ' + MegID+'.meg')   
#  copy(MegID+'WithOut.meg', dir+'\\Results\\'+ID+'_Final.meg')  
#  copy(MegID+'WithOut.meg',MegIDlast+'.meg')  
#else:  
# print 'Combine similar clone sequencies '
# copy(ID+'_CloneFreqSig.txt',MegID+'.txt')
# os.system('python '+'CombSimClo.py '+MegID+'.meg '+MegID+'.txt '+ str(ErroRate) + ' 0')   
# copy(MegID+'_FracNucSeq.txt',dir+'\\Results\\'+ID+'_FracNucSeq.txt')
# copy(MegID+'Comb.txt', dir+'\\Results\\'+ID+'_Final.txt') 
# copy(MegID+'_RandSeq.meg', dir+'\\Results\\'+ID+'_Final.meg')
# copy(MegID+'_RandSeq.meg',MegIDlast+'.meg')
# os.system('python '+'EvalClone.py ' + MegID+'_FracNucSeq.txt')
# copy(ID+'_CloneSNVassign.txt',dir+'\\Results\\'+ID+'_CloneSNVassign.txt')
# #open('AAA','r').readlines()
#print 'Make final clone phylogeny'
#os.system('megacc -a infer_MP_nucleotide.mao -d '+MegIDlast+'.meg'+' -o '+MegIDlast+'.nwk')
##Functions.DoMEGAMPsimple(MegIDlast+'.meg')
#copy(MegIDlast+'.nwk', dir+'\\Results\\'+ID+'_Final.nwk')
#
##########
#Functions.Clean(dir+'\\'+ID+'Rep*')
#Functions.Clean(dir+'\\*'+ID+'_*')
#
#Functions.Clean(Result_dir+'\\'+ID+'Rep*')
#os.remove(Result_dir+'\\'+ID+'_CloneSNVassign.txt')
#os.remove(Result_dir+'\\'+ID+'_ExtraMissCountFinal.txt')
#os.remove(Result_dir+'\\'+ID+'_resolveMultCount.txt')
#Functions.Clean(clonetree_dir+'\\'+ID+'*.txt')
#Functions.Clean(clonetree_dir+'\\'+ID+'*.nwk')
#Functions.Clean(clonetree_dir+'\\'+ID+'*.meg')
#os.chdir(dir+'\\Results\\')		
        		