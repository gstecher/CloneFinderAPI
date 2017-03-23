
#python CloneFinder.py snv [snv_input]
#python CloneFinder.py ccf [ccf_input] [read coverage]
#python CloneFinder.py cnv [snv_input] ##cnv information is 'snv_input'-CNV.txt

from parsers.DefaultTSPParser import DefaultTSPParser
from parsers.DefaultCCFParser import DefaultCCFParser
from config.ParamsLoader import ParamsLoader
from alignments.FreqToMegaSeq import FreqToMegaSeq
from alignments.MegaAlignment import MegaAlignment
from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from binomial.SignificanceTester import SignificanceTester
from regression.CloneFrequencyComputer import CloneFrequencyComputer
from decomposition.SNPClusterGenerator import SNPClusterGenerator
from decomposition.SNPGroupCombiner import SNPGroupCombiner
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
AnalyzeTree = TreeAnalizer()
OutFile = OutputWrite()
Align = MegaAlignment()
CloFreAna = CloneFrequencyAnalizer()

if parser.parse() == False:
        print parser.messages
else:
        tsp_list = parser.get_tumor_sample_profile_list()      		
        original_align_builder = FreqToMegaSeq()
        original_align_builder.initialize(tsp_list, True) # pass in True to remove duplicates, False to keep duplicates
        tumor_seqs = original_align_builder.get_mega_allalignment()
        initial_seq = original_align_builder.get_mega_alignment()	

def GetHead(Head):
    Head=Head.strip().split('\t')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col 
def GetCloHitForTu(File,Cut):
    File=open(File,'r').readlines()
    AA,Name2Col=GetHead(File[0])
    File=File[1:]
    Tu2Clo={}
    T2C2F={}	
    for i in File:
        i=i.strip().split('\t')	
        Tu=i[0]
        Tu2Clo[Tu]=[]
        T2C2F[Tu]={}		
        for Name in Name2Col:
          if Name!='Tumor':		
            Fre0=i[Name2Col[Name]]		
            Fre=float(Fre0)
            if Fre0.find('-')==-1 and Fre>Cut:
                 T2C2F[Tu][Name]=Fre
                 Tu2Clo[Tu].append(Name)				 
            else: T2C2F[Tu][Name]=0			
    return Tu2Clo,Name2Col,T2C2F

refine_clone = MutationFilter(params.option_a, params.option_c, params.option_d, tumor_seqs, tsp_list, params.error_rate, params.freq_cutoff, tree_builder.mao_file)			
significant_seq0=open('C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\clonefinderapi\\Test_mask.meg','r').readlines()
significant_seq=[]
for i in significant_seq0:
    significant_seq.append(i.strip())
CloFre='C:\\Users\\tuf78332\\Desktop\\CloneFinderAPI\\clonefinderapi\\test_mask.txt'
print significant_seq
#CloFreCom = CloneFrequencyComputer(significant_seq, tsp_list, 0.01)
Tu2Clo,Name2Col,significant_clone_frequency=GetCloHitForTu(CloFre ,0) 
print 'add missing mutations into clone'
significant_seq, significant_clone_frequency = refine_clone.AddTuExtra(significant_seq, significant_clone_frequency)
significant_seq+=['#hg19','A'*len(Initial_Good_posi_info)]
Align.save_mega_alignment_to_file('Test_add.meg', significant_seq)
CloFreAna.save_frequency_table_to_file('test_add.txt', significant_clone_frequency,  [])		
print final_seq, final_clone_frequency			
			