
#python CloneFinder.py snv [snv_input]
#python CloneFinder.py ccf [ccf_input] [read coverage]
#python CloneFinder.py cnv [snv_input] ##cnv information is 'snv_input'-CNV.txt

from parsers.DefaultTSPParser import DefaultTSPParser
from parsers.DefaultCCFParser import DefaultCCFParser
from parsers.DefaultCNVarser import DefaultCNVarser
from config.ParamsLoader import ParamsLoader
from alignments.FreqToMegaSeq import FreqToMegaSeq
from alignments.MegaAlignment import MegaAlignment
from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from binomial.SignificanceTester import SignificanceTester
from regression.CloneFrequencyComputer import CloneFrequencyComputer
from regression.CloneFrequencyComputer_cnv import CloneFrequencyComputer_cnv
from decomposition.SNPClusterGenerator import SNPClusterGenerator
from decomposition.SNPClusterGenerator_cnv import SNPClusterGenerator_cnv
from decomposition.SNPGroupCombiner import SNPGroupCombiner
from decomposition.SNPGroupCombiner_cnv import SNPGroupCombiner_cnv
from decomposition.DecomposedCloneSimmarize import DecomposedCloneSimmarize
from clone_adjustment.MutationFilter import MutationFilter
from output.OutputWrite import OutputWrite
import os

import datetime
startTime = datetime.datetime.now()
print startTime

dir = os.getcwd()
tree_builder = MegaMP()
tree_builder.mao_file = dir + '/infer_MP_nucleotide.mao'#mao_file  #'infer_MP_nucleotide.mao'
try: 
    loader = ParamsLoader()
    loader.params_file = 'C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\dev\\options.ini'#options_file
    params = loader.load_params()
    summary_file=params.to_string()+'\n'	
    print params.to_string() # just to see what it looks like
except:
    print 'Errors in options.ini'	
parser = DefaultTSPParser()
parser.input_data_file = params.input_data_file
#if params.data_format=='cnv':
parser_cnv = DefaultCNVarser() 
    
AnalyzeTree = TreeAnalizer()
OutFile = OutputWrite()
Align = MegaAlignment()
CloFreAna = CloneFrequencyAnalizer()

if parser.parse() == False:
        print parser.messages
else:
        tsp_list = parser.get_tumor_sample_profile_list() 
        #if params.data_format=='cnv':
             # print params.cnv_data_file	

			 
        CNV_information = parser_cnv.get_tumor_cnv_profile(params.cnv_data_file)	
             # print 'hhhh',CNV_information			  
        original_align_builder = FreqToMegaSeq()
        original_align_builder.initialize(tsp_list, True) # pass in True to remove duplicates, False to keep duplicates
        num_sites = tsp_list.num_read_counts()
        print params.option_b		
        params.option_b = str(1.0*num_sites*float(params.option_b.split(',')[0].strip())/100).split('.')[0] +','+params.option_b.split(',')[1] +','+params.option_b.split(',')[2] 
        params.option_c = params.option_c.split(',')[0] +','+str(1.0*num_sites*float(params.option_c.split(',')[1].strip())/100).split('.')[0]  
        params.option_d = params.option_d.split(',')[0] +','+str(1.0*num_sites*float(params.option_d.split(',')[1].strip())/100).split('.')[0]  		
        tumor_seqs = original_align_builder.get_mega_allalignment()
        initial_seq = original_align_builder.get_mega_alignment()
        #CNV_information = original_align_builder.get_CNV_information()	#{tumor:['neutral','loss',...],tumor:[]...}	
	
        print 'decompose and add ancestor'
		
        id = 'mega_alignment' # id will be used internally for file names
       # print original_align_builder.get_mega_alignment()			
        status = tree_builder.do_mega_mp(initial_seq, id)
        if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
                print 'best alignment'
             #   print seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info
        else:
                print 'failed to run megaMP'		
		
        RepeatUnion='y'	
        REP = 0		
        while RepeatUnion=='y':	
            REP += 1
            REPs = str(REP)
          #  print seqs_with_ancestor
          #  Align.save_mega_alignment_to_file('Initial.meg', seqs_with_ancestor)    

#############################			
            clone_frequency = CloneFrequencyComputer(seqs_with_ancestor, tsp_list, params.freq_cutoff)
			
            seqs_original, CloFre_seqs_with_ancestor = clone_frequency.regress()
          #  CloFreAna.save_frequency_table_to_file('Original_freq.txt', CloFre_seqs_with_ancestor, [])				
			
           # print 'original',CloFre_seqs_with_ancestor			
            #clusters = SNPClusterGenerator(seqs_with_ancestor, tsp_list, CloFre_seqs_with_ancestor, params.option_b.split(',')[0].strip(), params.freq_cutoff)
        #    Obs_SNVfreq, Est_SNVfreq = clusters.estimatedSNV_save_to_file()	#{tumor:[fre, fre, ...],tumor:[fre..]..}	
         #   new_clone_seq = Align.add_cnv_genotype(tsp_list, CloFre_seqs_with_ancestor, seqs_with_ancestor)	

            clone_frequency_cnv = CloneFrequencyComputer_cnv(seqs_with_ancestor, tsp_list, CNV_information, params.option_b.split(',')[0].strip(), params.freq_cutoff)
			
           # seqs, CloFre_seqs_with_ancestor_cnv = clone_frequency_cnv.regress()#before combining clones with different CNV profiles	
            seqs_with_ancestor_CNVadjusted, CloFre_seqs_with_ancestor_CNVadjusted = clone_frequency_cnv.snvGenotype_after_regress()	
         #   Align.save_mega_alignment_to_file('CNV_seq.meg', seqs_with_ancestor)		
          #  CloFreAna.save_frequency_table_to_file('CNV_freq.txt', CloFre_seqs_with_ancestor, [])					
			
###############################Clone freq computer for cnv			
           # clone_frequency.CNV_file = CNV_information
          #  clone_frequency.v_est = Est_SNVfreq	
          #  seqs, CloFre_seqs_with_ancestor_cnv = clone_frequency.regress()	
           # print 'cnv',CloFre_seqs_with_ancestor_cnv
          #  open('AAA','r').readlinse()

			
            print params.option_b.split(',')
            clusters = SNPClusterGenerator_cnv(seqs_with_ancestor, tsp_list, CloFre_seqs_with_ancestor, CNV_information, params.option_b.split(',')[0].strip(), params.freq_cutoff)
			
           # clusters = SNPClusterGenerator_cnv(seqs_with_ancestor, tsp_list, CloFre_seqs_with_ancestor, params.option_b.split(',')[0].strip(), params.freq_cutoff, CNV_information)
           # clusters.estimatedSNV_save_to_file()			
            tumor2cluster = clusters.cluster()
         #   print 'tumor2clusters ',tumor2cluster 
           # open('AAA','r').readlines()					
            decomposition = SNPGroupCombiner_cnv(tumor2cluster, seqs_with_ancestor, tumor_seqs, tsp_list, params.freq_cutoff, params.option_b.split(',')[1].strip(), params.option_b.split(',')[2].strip(), CNV_information, params.option_b.split(',')[0].strip())
            tumor2seqs_with_decompose, RmCloLs = decomposition.get_decomposed_seq()		

            print 'hybrid sample genotypes ',RmCloLs   
           # open('AAA','r').readlines()			
            if RmCloLs==[]:  
                        RepeatUnion='n'                      
            else: 
                  DecomposedCloneSummary = DecomposedCloneSimmarize()
                  remove_tumor_and_rename_decomposed_seq = DecomposedCloneSummary.remove_tumor_and_rename_decomposed(tumor2seqs_with_decompose, seqs_with_ancestor, tumor_seqs, REP, CloFre_seqs_with_ancestor)  
                #  print remove_tumor_and_rename_decomposed_seq				  
                  decom_without_ancdecom_seq = DecomposedCloneSummary.remove_ancestral_decomposed(remove_tumor_and_rename_decomposed_seq, params.error_rate, tumor_seqs)
				  
                  print 'fix multiple/backward mutations on decomposed clones'
                #  print decom_without_ancdecom_seq				  
                  no_back_para_mut_decomposed_seq = DecomposedCloneSummary.fix_back_para_mut_decomposed(decom_without_ancdecom_seq) 
                #  print 'h'				  
                  presence_of_new_decomposed_clone = DecomposedCloneSummary.find_decomposed_clone(no_back_para_mut_decomposed_seq, REP)
              #    print presence_of_new_decomposed_clone				  
               #   open('AAA','r').readlines()				  
                  if presence_of_new_decomposed_clone=='n':
                         print 'there are no hybrid sample genotypes'
                         						 
                         RepeatUnion='n'
                  else:
                       id = 'dec'+REPs+'_mega_alignment' # id will be used internally for file names
                       #print no_back_para_mut_decomposed_seq					   
                       status = tree_builder.do_mega_mp(no_back_para_mut_decomposed_seq, id)
                       if status == True:
                             decomseqs_with_ancestor, decom_tree, decom_nade_map, decom_mask_seq, decom_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)
                       else:
                             print 'failed to run megaMP'
						 
                       decomposition_is_better = AnalyzeTree.compare_good_posi_number(Initial_Good_posi_info, decom_Good_posi_info)
                     #  print decomposition_is_better				  
                    #   open('AAA','r').readlines()						   
                       if decomposition_is_better=='n': 
                             print 'decomposed clone genotypes are not good (create more backward/parallel mutations)'					   
                             RepeatUnion='n'
                       else:
                            print 'clones were decomposed, so repeat the decomposition'
                            seqs_with_ancestor=decomseqs_with_ancestor
                          #  Align.save_mega_alignment_to_file('AfterDec_seq.meg', seqs_with_ancestor)
                           # open('AAA','r').readlines()								
                            tree=decom_tree
                            nade_map=decom_nade_map
                            mask_seq=decom_mask_seq
                            Initial_Good_posi_info=decom_Good_posi_info							
                            							
        print 'clone decomposition is complete!'#,REP,seqs_with_ancestor,CloFre_seqs_with_ancestor,Initial_Good_posi_info,mask_seq
		
        # 'rename clone anmes'
        seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName, clone_order  = OutFile.ReNameCloFreMeg(seqs_with_ancestor_CNVadjusted, CloFre_seqs_with_ancestor_CNVadjusted, 'number')		
       # Align.save_mega_alignment_to_file('Test_dec_rename.meg', seqs_with_ancestor_newName)
      #  Align.save_mega_alignment_to_file('Test_dec.meg', seqs_with_ancestor_newName)	
      #  CloFreAna.save_frequency_table_to_file('test_dec.txt', CloFre_seqs_with_ancestor_newName, [])	
        print 'refine clone genotypes'		
        number_of_good_positions = AnalyzeTree.count_good_posi(Initial_Good_posi_info)
        proportion_of_good_positions = 1.0*number_of_good_positions/len(Initial_Good_posi_info)
        if proportion_of_good_positions<0.5:
            print 'Only' + str(proportion_of_good_positions)+' of SNV sites are not affected by backward/parallel mutations. Please confirm the tumor profiles.\nOptionA, C, and D were changed to Off, if these were On.'		
            summary_file +=  'Only' + str(proportion_of_good_positions)+' of SNV sites are not affected by backward/parallel mutations. Please confirm the tumor profiles.\nOptionA, C, and D were changed to Off, if these were On.\n'
            params.option_a='Off'
            params.option_c='Off'
            params.option_d='Off' 
        if params.option_a=='Off' and params.option_c=='Off' and params.option_d=='Off': pass
        else:
            refine_clone = MutationFilter(params.option_a, params.option_c, params.option_d, tumor_seqs, tsp_list, params.error_rate, params.freq_cutoff, tree_builder.mao_file)			
            if params.option_a=='On':			
                print 'Cut extra mutations in clones'
                no_extra_mutation_seq, no_extra_mutation_clone_frequency = refine_clone.CutExtraMutClo(seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName)			
                seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName, clone_order  = OutFile.ReNameCloFreMeg(no_extra_mutation_seq, no_extra_mutation_clone_frequency, 'number')
              #  Align.save_mega_alignment_to_file('Test_cut.meg', seqs_with_ancestor_newName)				
              #  CloFreAna.save_frequency_table_to_file('test_cut.txt', CloFre_seqs_with_ancestor_newName,  [])					
            if params.option_a=='On' or params.option_c!='Off': 
                print 'Resolve parallel mutations'#,seqs_with_ancestor_newName
                no_parallel_mutation_seq, no_parallel_mutation_clone_frequency = refine_clone.BranchDecClone(seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName)
                seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName, clone_order  = OutFile.ReNameCloFreMeg(no_parallel_mutation_seq, no_parallel_mutation_clone_frequency, 'number')				
               # Align.save_mega_alignment_to_file('Test_resolve.meg', seqs_with_ancestor_newName)
               # CloFreAna.save_frequency_table_to_file('test_resolve.txt', CloFre_seqs_with_ancestor_newName,  [])				
            if params.option_d=='Off': pass
            else:
                print 'add ancestral clone at the middle of branches'            				
                CloFreCom = CloneFrequencyComputer(seqs_with_ancestor_newName, tsp_list, params.freq_cutoff)
                OutAlign, OutCloFre = CloFreCom.regress_hit_clone(CloFre_seqs_with_ancestor_newName)
                OutAlign+=['#hg19','A'*len(Initial_Good_posi_info)]	
				
                id = 'decanc_alignment' # id will be used internally for file names		
                status = tree_builder.do_mega_mp(OutAlign, id)
                if status == True:
                    decanc_seqs_with_ancestor, decanc_tree, decanc_nade_map, decanc_mask_seq, decanc_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
                else:
                    print 'failed to run megaMP'
			
                add_inter_ancester_seq, add_inter_ancester_clone_frequency = refine_clone.BranchCloneInfer(OutAlign, OutCloFre, decanc_nade_map)
                seqs_with_ancestor_newName0, CloFre_seqs_with_ancestor_newName0, clone_order  = OutFile.ReNameCloFreMeg(add_inter_ancester_seq, add_inter_ancester_clone_frequency, 'number')	  
##################
                clone_frequency_cnv = CloneFrequencyComputer_cnv(seqs_with_ancestor_newName0, tsp_list, CNV_information, params.option_b.split(',')[0].strip(), params.freq_cutoff)
			
           # seqs, CloFre_seqs_with_ancestor_cnv = clone_frequency_cnv.regress()#before combining clones with different CNV profiles	
                seqs_with_ancestor_newName, CloFre_seqs_with_ancestor_newName = clone_frequency_cnv.snvGenotype_after_regress()	

            ##no binomial test			
            significant_clone_frequency = {}
            for Tu in CloFre_seqs_with_ancestor_newName:
               significant_clone_frequency[Tu]= CloFre_seqs_with_ancestor_newName[Tu]			
            significant_seq = seqs_with_ancestor_newName
            #print 'h',significant_clone_frequency         			
            ###
    
            if params.option_e=='Off' and params.option_f=='Off':pass#Functions.Reset(ID,InFile,'Last.txt','_AfterBino',MegID+'.meg',MegID+'.txt',clonetree_dir,dir)
            else: 
                if params.option_e=='On':
                     print 'mask sites with many backward/parallel mutations'					 
                     id = 'mask_alignment' # id will be used internally for file names		
                     status = tree_builder.do_mega_mp(significant_seq, id)
                     if status == True:
                          significant_seq_seqs_with_ancestor, significant_seq_tree, significant_seq_nade_map, significant_mask_seq, significant_seq_decanc_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
                     else:
                          print 'failed to run megaMP'
                   #  Align.save_mega_alignment_to_file('Test_mask0.meg', significant_mask_seq)
                   #  CloFreAna.save_frequency_table_to_file('test_mask0.txt', significant_clone_frequency,  [])							  
                     significant_seq, outSeqAmb, significant_clone_frequency=Align.CombSimClo(significant_mask_seq, significant_clone_frequency, 0)
                   #  Align.save_mega_alignment_to_file('Test_mask.meg', significant_seq)
                  #   CloFreAna.save_frequency_table_to_file('test_mask.txt', significant_clone_frequency,  [])	
					 
                if params.option_f=='On':	 
                     print 'add missing mutations into clone'
                     significant_seq, significant_clone_frequency = refine_clone.AddTuExtra(significant_seq, significant_clone_frequency)
                     significant_seq+=['#hg19','A'*len(Initial_Good_posi_info)]
                  #   Align.save_mega_alignment_to_file('Test_add.meg', significant_seq)
                  #   CloFreAna.save_frequency_table_to_file('test_add.txt', significant_clone_frequency,  [])					 

            print 'final estimation of clone frequencies'						

        #    Align.save_mega_alignment_to_file('Test_before.meg', significant_seq)
        #    CloFreAna.save_frequency_table_to_file('test_before.txt', significant_clone_frequency,  [])	
###################			
           # CloFreCom = CloneFrequencyComputer(significant_seq, tsp_list, params.freq_cutoff)			
          #  final_seq, final_clone_frequency = CloFreCom.regress_hit_clone(significant_clone_frequency)
		######	
            clone_frequency_cnv = CloneFrequencyComputer_cnv(significant_seq, tsp_list, CNV_information, params.option_b.split(',')[0].strip(), params.freq_cutoff)
			
           # seqs, CloFre_seqs_with_ancestor_cnv = clone_frequency_cnv.regress()#before combining clones with different CNV profiles	
            final_seq, final_clone_frequency = clone_frequency_cnv.snvGenotype_after_regress_hitclone(significant_clone_frequency)	
        #    CloFreAna.save_frequency_table_to_file('test_after.txt', final_clone_frequency,  [])
         #   Align.save_mega_alignment_to_file('Test_after.meg', final_seq)
          #  open('AAA','r').readlines()			
################
            if params.option_h=='Off': pass
            else:
                  print 'combine similar clones (mutant/wild is determined by majority of combined clones or is annotated as ?, when mutant/wild does not agree with all combined clones)'			
                  final_seq, final_seq_amb, final_clone_frequency=Align.CombSimClo(final_seq, final_clone_frequency, params.error_rate) 
                 # Align.save_mega_alignment_to_file('Final.meg', final_seq1)

            print 'generate output'				  
            final_seq1, final_clone_frequency1, final_clone_order1  = OutFile.ReNameCloFreMeg(final_seq, final_clone_frequency, 'list')	
            Align.save_mega_alignment_to_file(params.input_id + '_Final.meg', final_seq1)
            id = 'final1_alignment' # id will be used internally for file names
       # print original_align_builder.get_mega_alignment()			
            status = tree_builder.do_mega_mp(final_seq1, id)
            if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
              #  print 'best alignment',tree,final_seq1
                OutF=open(params.input_id + '_Final.nwk','w')
                OutF.write(tree)
                OutF.close()				
             #   print seqs_with_ancestor, tree, nade_map, mask_seq, Initial_Good_posi_info
            else:
                print 'failed to run megaMP'				
            CloFreAna.save_frequency_table_to_file(params.input_id + '_Final.txt', final_clone_frequency1,  final_clone_order1)	

            print 'binomial significance test'		
			
            significant_clones = SignificanceTester(params.bino_num, params.bino_cutoff, params.error_rate, tsp_list, final_seq1, final_clone_frequency1,CNV_information, params.option_b.split(',')[0].strip())
         #   significant_clone_in_tumor  = significant_clones.compute_p()
          #  print significant_clone_in_tumor			
          #  CloFreAna.save_frequency_table_to_file(params.input_id + '_p-valueAll.txt', significant_clone_in_tumor,  final_clone_order1)			
            significant_clone_in_tumor  = significant_clones.compute_p_printAnc()			
            CloFreAna.save_frequency_table_to_file(params.input_id + '_p-value.txt', significant_clone_in_tumor,  final_clone_order1)

            if 	params.option_h!='Off':			  
                final_seq_amb, final_clone_frequency2, final_clone_order1  = OutFile.ReNameCloFreMeg(final_seq_amb, final_clone_frequency, 'list')	
                Align.save_mega_alignment_to_file(params.input_id + '_with_Ambiguous.meg', final_seq_amb)

#######################
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
