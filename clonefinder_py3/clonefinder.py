# 'python CloneFinder.py snv [snv_input]'
from parsers.DefaultTSPParser import DefaultTSPParser
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
from decomposition.SNPGroupCombiner_cnv1 import SNPGroupCombiner_cnv1
from parsers.DefaultCNVarser import DefaultCNVarser
from config.ParamsLoader import ParamsLoader
from config.FormatInput import FormatInput
from alignments.FreqToMegaSeq import FreqToMegaSeq
from alignments.MegaAlignment import MegaAlignment
from output.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from parsimony.MegaAncestor import MegaAncestor
from tsp_profiles.tsp_information import tsp_information
from significance_test.cluster_test import cluster_test
from output.OutputWrite import OutputWrite
import os
import sys
import datetime

Significant_cutoff=0.05

startTime = datetime.datetime.now()
print (startTime)
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
    print ('Errors in options.ini')	
parser = DefaultTSPParser()
parser.input_data_file = params.input_data_file
parser_cnv = DefaultCNVarser() 
AnalyzeTree = TreeAnalizer()
OutFile = OutputWrite()
Align = MegaAlignment()
CloFreAna = CloneFrequencyAnalizer()
Format = FormatInput()


if parser.parse() == False:
        print (parser.messages)
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
       initial_seq0 = original_align_builder.get_mega_alignment()
       initial_seq, A1, A2  = OutFile.ReNameCloFreMeg(initial_seq0, {}, 'number')   
       CNV_information_test = Format.add_low_quality_SNV_info(CNV_information,total_read, alt_read,total_read_cut,mutant_read_cut)          	
       print ('add initial ancestor')		
       id = 'mega_alignment' # id will be used internally for file names
       status = tree_builder.do_mega_mp(initial_seq, id)
       if status != True:
                print ('CloneFinder cannot make inferences, because samples do not have evolutionary structure.')		
       else:
        seqs_with_ancestor, A1 = tree_builder.alignment_least_back_parallel_muts() # True will execute code to remove redundant seqs (True is default)		              				
        Repeat='y'	
        ItC=0		
        while Repeat=='y' and ItC<=5:	
            print ('compute clone frequencies with ancestor')		 
            clone_frequency_cnv = CloneFrequencyComputer_cnv1(seqs_with_ancestor, v_obs, CNV_information_test, params.freq_cutoff)
            clone_frequency_cnv.regress_cnv()	            			
            print ('decompose incorrect candidate clone genotypes')
            decomposition = SNPGroupCombiner_cnv1(clone_frequency_cnv.hitclone_seq_builder, v_obs, clone_frequency_cnv.Tumor2Clone_frequency, CNV_information, params.freq_cutoff, tumor_seqs)			
            hit_seq_dic, Message = decomposition.get_decomposed_seq()			
            print (Message) 
            summary_file += Message+'\n'
            DecomTuNum=len(Message.split(','))
           			
            if Message[:10]!='decomposed' or DecomTuNum>(1.0*len(v_obs)/2):  
                        Repeat='n'                      
                        final_seq = clone_frequency_cnv.hitclone_seq_builder
                        final_clofre = clone_frequency_cnv.Tumor2Clone_frequency
                        print ('no more new clones',clone_frequency_cnv.Tumor2Clone_frequency)
				
            else: 
                  				  
                       hit_seq_dic['#hg19']='A'*num_sites
                       hit_seq_build = Align.UpMeg(hit_seq_dic,list(hit_seq_dic.keys()))					   
                       id = 'dec_mega_alignment' 				   
                       status = tree_builder.do_mega_mp(hit_seq_build, id)
                       if status == True and ItC<5:
                             decomseqs_with_ancestor, A1 = tree_builder.alignment_least_back_parallel_muts() # True will execute code to remove redundant seqs (True is default)
							# self._cleanup_temp_files()	
                             clone_frequency_cnv2 = CloneFrequencyComputer_cnv1(decomseqs_with_ancestor, v_obs, CNV_information_test, params.freq_cutoff)
                             clone_frequency_cnv2.regress_cnv()	
                             decomseqs_with_ancestor_clone_freq=clone_frequency_cnv2.Tumor2Clone_frequency		
                             decomseqs_with_ancestor=	clone_frequency_cnv2.hitclone_seq_builder	
                             seqs_with_ancestor, A1, A2  = OutFile.ReNameCloFreMeg(decomseqs_with_ancestor, decomseqs_with_ancestor_clone_freq, 'number')
                             ItC+=1	
                             
							 
                             print ('clones were decomposed, so repeat the decomposition')									 
                       else:
                             print ('decomposed clones do not have evolutionary structure.. So, we will not decompose them.')  
                             summary_file+=	'decomposed clone genotypes are not good (no evolutionary structure), so clones were not decomposed.\n'								 					 
                             Repeat='n'
                             final_seq = clone_frequency_cnv.hitclone_seq_builder
                             final_clofre = clone_frequency_cnv.Tumor2Clone_frequency					 				 
					                  							
        print ('clone decomposition is complete!')	

        final_seq1, final_clone_frequency1, final_clone_order1  = OutFile.ReNameCloFreMeg(final_seq,final_clofre,   'number') ###	

        print ('test clone hit and remove insignificant clones')
        significant_clone=cluster_test()		
        final_seq, final_clofre= significant_clone.remove_insignificant_clones_add(v_obs, final_clone_frequency1,  final_seq1, CNV_information_test, Significant_cutoff)
        print ('making output files')	

        final_seq1, final_clone_frequency1, final_clone_order1  = OutFile.ReNameCloFreMeg(final_seq,final_clofre,   'number') 		
		
        Align.save_mega_alignment_to_file(params.input_id + '_CloneFinder.meg', final_seq1)
        CloFreAna.save_frequency_table_to_file(params.input_id + '_CloneFinder.txt',  final_clone_frequency1,  [])	

        id = 'final'		
        status = tree_builder.do_mega_mp(final_seq1, id)
        if status == True:
                A1, tree = tree_builder.alignment_least_back_parallel_muts() 
                Rooted=AnalyzeTree.RootTree(tree)
                InferAncestor = MegaAncestor()
                InferAncestor.alignment_file = final_seq1
                InferAncestor.input_tree_file = Rooted
     			
                ancestor_states, offspring2ancestor, cell2code, code2cell = InferAncestor.retrieve_ancestor_states()
                RescaledTree=InferAncestor.get_scaledNWK()					
                OutFile.GetOut(params.input_id + '_CloneFinder.nwk', RescaledTree.replace('hg19:','Normal:'))				                     
	
		

		
	
#######################
os.remove(params.input_id + '.txt')
os.remove(params.input_id + '-CNV.txt')

timeFile = params.input_id + '_summary.txt'
endTime = datetime.datetime.now()
print (endTime)
totalTime = (endTime - startTime)
print (totalTime)
summary_file += 'Run time: ' + str(totalTime) + '\n'
OutTime=open(timeFile,'w')
OutTime.write(summary_file)
OutTime.close()
#######################			
