import sys
import os
from alignments.MegaAlignment import MegaAlignment
from output.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer

CFin=sys.argv[1]
CFmeg=CFin[:-4]+'snv_CloneFinder.meg'
ID=CFin.split('\\')[-1][:-4]

cwd=os.getcwd()
OutFol=CFin.replace(CFin.split('\\')[-1],'')+ID
if os.path.exists(OutFol)!=True: os.mkdir(OutFol)
OutFol+='\\'
print (OutFol)
#open('A','r').readlines()
tree_builder=MegaMP()
tree_analyze=TreeAnalizer()
initial_seq=open(CFmeg,'r').readlines()
status = tree_builder.do_mega_mp(initial_seq, 'initial')
if status != True:
                print ('Tree cannot be inferred.')		
else:
        Nwk=tree_builder.get_trees()
       # print (Nwk)		
        #seqs_with_ancestor, A1 = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)	
        RootedNwk=tree_analyze.RootTree(Nwk[0])
        print (RootedNwk)
        #tree_analyze.MLancetor(initial_seq,RootedNwk)
        tree_analyze.output_PhyloSigFinderIn(initial_seq,RootedNwk,CFin,OutFol,ID+'\\')
        tree_analyze.get_tree_with_branchLen(CFin[:-4])		
     	
        		

		
