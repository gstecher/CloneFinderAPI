from parsers.DefaultTSPParser import DefaultTSPParser
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
from decomposition.SNPClusterGenerator_cnv1 import SNPClusterGenerator_cnv1
from decomposition.SNPGroupCombiner_cnv1 import SNPGroupCombiner_cnv1
from parsers.DefaultCNVarser import DefaultCNVarser
from config.ParamsLoader import ParamsLoader
from config.FormatInput import FormatInput
from alignments.FreqToMegaSeq import FreqToMegaSeq
from alignments.MegaAlignment import MegaAlignment
from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from decomposition.DecomposedCloneSimmarize import DecomposedCloneSimmarize
from clone_adjustment.MutationFilter import MutationFilter
from tsp_profiles.tsp_information import tsp_information
from output.OutputWrite import OutputWrite
from significance_test.cluster_test import cluster_test
import os
import sys
import datetime
print 'python CloneFinder.py snv [snv_input]'
startTime = datetime.datetime.now()
print startTime
Significant_cutoff=0.05
dir = os.getcwd()
tree_builder = MegaMP()
tree_builder.mao_file = dir + '/infer_MP_nucleotide.mao'
try: 
    loader = ParamsLoader()
    loader.params_file = 'options.ini' #options_file
    params = loader.load_params()
    total_read_cut = params.total_read_cut
    mutant_read_cut = params.mutant_read_cut	
    summary_file= "input data file: "+params.input_data_file+'\n'+"total read count cutoff: "+str(params.total_read_cut)+'\n'"mutant read count cutoff: "+str(params.mutant_read_cut)+'\n'"clone frequency cutoff: "+str(params.freq_cutoff)+'\n\n'
except:
    print 'Errors in options.ini'	
parser = DefaultTSPParser()
parser.input_data_file = params.input_data_file
parser_cnv = DefaultCNVarser() 
AnalyzeTree = TreeAnalizer()
OutFile = OutputWrite()
Align = MegaAlignment()
CloFreAna = CloneFrequencyAnalizer()
Format = FormatInput()
if parser.parse() == False:
        print parser.messages
else:
        tsp_list = parser.get_tumor_sample_profile_list() 			 
        CNV_information = parser_cnv.get_tumor_cnv_profile(params.cnv_data_file)		
        original_align_builder = FreqToMegaSeq()
        original_align_builder.initialize(tsp_list, True) #True to remove duplicates; False to keep duplicates
        num_sites = tsp_list.num_read_counts()
        all_tsp = tsp_information(tsp_list)    				
        v_obs = all_tsp.tumor2alt_frequency()
        total_read, alt_read, ReadCount_table = all_tsp.tumor2read_count()				
        tumor_seqs = original_align_builder.get_mega_allalignment()
        initial_seq = original_align_builder.get_mega_alignment()
        CNV_information_test = Format.add_low_quality_SNV_info(CNV_information,total_read, alt_read,total_read_cut,mutant_read_cut)
	
        print 'add initial ancestor'		
        id = 'mega_alignment' # id will be used internally for file names
        status = tree_builder.do_mega_mp(initial_seq, id)
        if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
                print 'best alignment'
        else:
                print 'failed to run megaMP'				
        RepeatUnion='y'	
        REP = 0		
        while RepeatUnion=='y':	
            REP += 1
            REPs = str(REP)	
            print 'compute clone frequencies with ancestor'		 
            clone_frequency_cnv = CloneFrequencyComputer_cnv1(seqs_with_ancestor, tsp_list, CNV_information_test, params.freq_cutoff, ReadCount_table)
            clone_frequency_cnv.regress_cnv()		
            print 'test and make decomposed calone sequences'			
            clusters = SNPClusterGenerator_cnv1(clone_frequency_cnv.hitclone_seq_builder, tsp_list, clone_frequency_cnv.Tumor2Clone_frequency, CNV_information, params.freq_cutoff)		
            tumor2cluster = clusters.cluster_cnv()				
            decomposition = SNPGroupCombiner_cnv1(tumor2cluster, seqs_with_ancestor, tumor_seqs, tsp_list, params.freq_cutoff, CNV_information, ReadCount_table)
            tumor2seqsubs_with_decompose_builder, decomposed_tumor_2clonefrequency, RmCloLs = decomposition.get_decomposed_seq()						
            print 'hybrid sample genotypes ',RmCloLs          		
            if RmCloLs==[]:  
                        RepeatUnion='n'                      
                        final_seq = clone_frequency_cnv.hitclone_seq_builder
                        final_clofre = clone_frequency_cnv.Tumor2Clone_frequency
                        print 'no more new clones',clone_frequency_cnv.Tumor2Clone_frequency                   
            else: 
                  summary_file += 'hybrid sample genotypes were suggested: '+str(RmCloLs) +'\n'		
                  print 're-generate whole sequences'
                  original_seq_builder=clone_frequency_cnv.hitclone_seq_builder
                  original_clone_freq =	clone_frequency_cnv.Tumor2Clone_frequency
                  original_set=[original_seq_builder,original_clone_freq]				             				  
                  DecomposedCloneSummary = DecomposedCloneSimmarize()                  		  
                  decom_all_seq_builder = DecomposedCloneSummary.add_back_CNVSNV(tumor2seqsubs_with_decompose_builder,CNV_information,seqs_with_ancestor, original_clone_freq,tsp_list)                  
                  rename_seq_builder, rename_clofre=DecomposedCloneSummary.finalize_results(decom_all_seq_builder, decomposed_tumor_2clonefrequency, original_seq_builder, original_clone_freq, REP)				  
                  print 'fix multiple/backward mutations on decomposed clones'		
                  rename_seq_builder+=['#hg19','A'*num_sites]				  
                  no_back_para_mut_decomposed_seq,Tree = DecomposedCloneSummary.fix_back_para_mut_decomposed(rename_seq_builder) 			  
                  presence_of_new_decomposed_clone, RmDecom, good_seq_builder = DecomposedCloneSummary.find_decomposed_clone(no_back_para_mut_decomposed_seq, REP, Tree)
                  print 'recompute clone frequency'	,RmDecom         		  
                  clone_frequency_cnv2 = CloneFrequencyComputer_cnv1(good_seq_builder, tsp_list, CNV_information_test, params.freq_cutoff, ReadCount_table)
                  clone_frequency_cnv2.regress_cnv()	
                  no_back_para_mut_decomposed_clone_freq=clone_frequency_cnv2.Tumor2Clone_frequency		
                  no_back_para_mut_decomposed_seq=	clone_frequency_cnv2.hitclone_seq_builder	  
                  print 'test clone hit and remove insignificant clones'
                  significant_clone=cluster_test()			
                  no_back_para_mut_decomposed_seq, no_back_para_mut_decomposed_clone_freq= significant_clone.remove_insignificant_clones(v_obs, clone_frequency_cnv2.Tumor2Clone_frequency,  clone_frequency_cnv2.hitclone_seq_builder, CNV_information_test, Significant_cutoff)	
                  no_back_para_mut_decomposed_seq += ['#hg19','A'*num_sites]					  
                  print no_back_para_mut_decomposed_clone_freq	
                  Iden=DecomposedCloneSummary.find_new_clone(no_back_para_mut_decomposed_seq,clone_frequency_cnv.hitclone_seq_builder)			  
                  if Iden=='y': presence_of_new_decomposed_clone='n'				  			  
                  if presence_of_new_decomposed_clone=='n':
                         final_seq = original_set[0]
                         final_clofre = original_set[1]			  
                         print 'there are no hybrid sample genotypes'                       						 
                         RepeatUnion='n'
                  elif 	presence_of_new_decomposed_clone=='anc':                 			  
                         final_seq = rename_seq_builder
                         final_clofre = rename_clofre			  
                         print 'hybrid sample genotypes are ancestral'                      			 
                         RepeatUnion='n'				  
                  else:
                       id = 'dec'+REPs+'_mega_alignment' 				   
                       status = tree_builder.do_mega_mp(no_back_para_mut_decomposed_seq, id)
                       if status == True:
                             decomseqs_with_ancestor, decom_tree, decom_nade_map, decom_mask_seq, decom_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)
                       else:
                             print 'failed to run megaMP'              					 
                       decomposition_is_better = AnalyzeTree.compare_good_posi_number(Initial_Good_posi_info, decom_Good_posi_info, seqs_with_ancestor,no_back_para_mut_decomposed_seq)	      		   
                       if decomposition_is_better=='n': 
                             print 'decomposed clone genotypes are not good (create more backward/parallel mutations)'	
                             summary_file+=	'decomposed clone genotypes are not good (create more backward/parallel mutations), so clones were not decomposed.\n'						 
                             RepeatUnion='n'
                             final_seq = original_set[0]
                             final_clofre = original_set[1]							 				 
                       else:
                            print 'clones were decomposed, so repeat the decomposition'					
                            seqs_with_ancestor=decomseqs_with_ancestor								
                            tree=decom_tree
                            nade_map=decom_nade_map
                            mask_seq=decom_mask_seq
                            Initial_Good_posi_info=decom_Good_posi_info
                            summary_file+='Hybrid sample genotype was decomposed\n'							                  							
        print 'clone decomposition is complete!', REP
        clone_order, clone_seq_dic=Align.name2seq(final_seq)	
        seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName, clone_order  = OutFile.ReNameCloFreMeg(final_seq, final_clofre, 'number')	
        print 'refine clone genotypes'		
        number_of_good_positions = AnalyzeTree.count_good_posi(Initial_Good_posi_info)
        proportion_of_good_positions = 1.0*number_of_good_positions/len(Initial_Good_posi_info)
        print 'proportion of good positions',proportion_of_good_positions			
        summary_file += str(proportion_of_good_positions*100)+'% of SNV sites are not affected by backward/parallel mutations.\n'       
        refine_clone = MutationFilter(tumor_seqs, tsp_list, tree_builder.mao_file)			
        print 'Resolve parallel mutations'
        no_parallel_mutation_seq, no_parallel_mutation_clone_frequency = refine_clone.BranchDecClone(seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName, CNV_information_test)				
        seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName, clone_order  = OutFile.ReNameCloFreMeg(no_parallel_mutation_seq, no_parallel_mutation_clone_frequency, 'number')				
        print 'generate output'					
        clone_frequency_cnv = CloneFrequencyComputer_cnv1(seqs_with_ancestor_newName, tsp_list, CNV_information_test, params.freq_cutoff, ReadCount_table)
        clone_frequency_cnv.regress_cnv()	    		
        final_seq1, final_clone_frequency1, final_clone_order1  = OutFile.ReNameCloFreMeg(clone_frequency_cnv.hitclone_seq_builder, clone_frequency_cnv.Tumor2Clone_frequency, 'number')	
        id = 'mega_alignment'		
        status = tree_builder.do_mega_mp(final_seq1, id)
        if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
                print 'best alignment'                      
        else:
                print 'failed to run megaMP'
        AA, mask_seq_comb, clone_freq = Align.CombSimClo(mask_seq, final_clone_frequency1, 0)		
        print 'test clone hit and remove insignificant clones'
        significant_clone=cluster_test()			
        significant_seq, significant_clone_frequency= significant_clone.remove_insignificant_clones(v_obs, clone_freq,  mask_seq_comb, CNV_information_test, Significant_cutoff)		
        Align.save_mega_alignment_to_file(params.input_id + '_CloneFinder.meg', significant_seq)
        CloFreAna.save_frequency_table_to_file(params.input_id + '_CloneFinder.txt',  significant_clone_frequency,  [])	
	
#######################
os.remove(params.input_id + '.txt')
os.remove(params.input_id + '-CNV.txt')
os.remove('Ini.meg')
os.remove('Ini_freq.txt')
os.remove('Test.meg')

timeFile = params.input_id + '_summary.txt'
endTime = datetime.datetime.now()
print endTime
totalTime = (endTime - startTime)
print totalTime
summary_file += 'Run time: ' + str(totalTime) + '\n'
OutTime=open(timeFile,'w')
OutTime.write(summary_file)
OutTime.close()
#######################			
