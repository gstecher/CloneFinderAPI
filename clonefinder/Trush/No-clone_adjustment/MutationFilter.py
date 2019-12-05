from alignments.MegaAlignment import MegaAlignment
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from estimated_clone_frequency.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from tsp_profiles.tsp_information import tsp_information
from output.CloneFrequencyWriter import CloneFrequencyWriter
#from decomposition.SNPClusterGenerator import SNPClusterGenerator

class MutationFilter:

    def __init__(self, tumor_seq, tsp_list, mao_file):
        Align = MegaAlignment()
        self.tumor_list, self.tumor2seq = Align.name2seq(tumor_seq)
        self.Len= len(self.tumor2seq[self.tumor_list[0]])
        self.mao_file = mao_file
        self.tsp_list = tsp_list
        TSPinfo=tsp_information(tsp_list)		
        self.Tu2SNV=TSPinfo.tumor2alt_frequency()		
       	

    def BranchDecClone(self, seq_list, clone_frequency, Tu2CNV):		
        Align = MegaAlignment()
        TumorSampleExtract = tsp_information(self.tsp_list)	
        CloFreAna = CloneFrequencyAnalizer()      		
        CloOrder, Clo2Seq = Align.name2seq(seq_list)
        Align.save_mega_alignment_to_file('Test.meg', seq_list)		
        tree_builder = MegaMP()
        tree_builder.mao_file = self.mao_file		
        id = 'branchdec_mega_alignment' 
   
        status = tree_builder.do_mega_mp(seq_list, id)
        if status == True:
                seqs_with_ancestor, tree, nade_map, mask_seq, Good_posi_info = tree_builder.alignment_least_back_parallel_muts(True) # True will execute code to remove redundant seqs (True is default)		
        else:
                print 'failed to run megaMP'
        BadPosiLs=[] #multiple mutations
        BadPosi2ChnageCloLs={}		
        for c in Good_posi_info:
           Posi_Inf = Good_posi_info[c]		
           if Posi_Inf!=['Good']:
              if Posi_Inf[0] == 'ToWild':		   
                  BadPosiLs.append(c)
                  BadPosi2ChnageCloLs[c]=Posi_Inf[1][0]	
        print 'bad positions',BadPosiLs	#,BadPosi2ChnageCloLs			  
        if BadPosiLs!=[]:	 
         NewT2C2F={}
         NewT2Cls={}
         for Tu in clone_frequency:
            NewC2F=	{}
            single_tsp_list = TumorSampleExtract.make_single_tsp_list(Tu)           
            CloFreDic=clone_frequency[Tu]
            CNV=Tu2CNV[Tu[2:]]			
            Tu=Tu[2:]			
            TuSeq	=self.tumor2seq['#'+Tu]
            NewCloLs=[]	
            NewCloLs1=[]
					
            for Clo in CloFreDic: #original hit clo for the tumor
              ChangeOptions='n'
           #   print Tu,CloFreDic			  
              if CloFreDic[Clo]>0:	 
                 CSeq0=Clo2Seq['#'+Clo]
                 ChangePosi=[] #list to fix multiple mutaitons
                 NewBadPosi=[] #remove fixed multiple mutations from BadExtMutPosi				 
                 for Bad in BadPosi2ChnageCloLs:
                     if BadPosi2ChnageCloLs[Bad].count('#'+Clo)!=0 and (CNV[Bad]=='normal' or CNV[Bad]=='Bad-normal'):	
                        Change='n'					 
                        for Oth in CloFreDic: #find multiple mutations at the external branch
                            if Oth !=Clo and CloFreDic[Oth]>0:
                               Soth=Clo2Seq['#'+Oth]
                               if Soth[Bad]=='T' and BadPosi2ChnageCloLs[Bad].count('#'+Oth)==0: 
                                         Change='y'
                        if Change=='y':										 
                                         ChangePosi.append(Bad)
                        else: 	NewBadPosi.append(Bad)										 
                 print 	'change positions',Tu,ChangePosi									 
                 if ChangePosi!=[]:#fix multiple mutaitons
                            #  print 'hhh'				 
                              CutCloSeq=Align.ModSeq(CSeq0,ChangePosi,'A',self.Len)					  				  
                              NewCloLs.append(Clo+'Cut'+Tu)
                              NewC2F[Clo+'Cut'+Tu]=CloFreDic[Clo]
                              Clo2Seq['#'+Clo+'Cut'+Tu]=CutCloSeq
                              ChangeOptions='y'							  		               		 
 						
              if ChangeOptions=='n':
                   NewC2F[Clo]=1			     
            NewT2C2F[Tu]=NewC2F
       #  print Clo2Seq			
         hitseq_align, hitclone_frequency  = CloFreAna.ListHitCloAndSeq(NewT2C2F, Clo2Seq)		 
         outSeqMaj, outSeqAmb, NewT2C2F = Align.CombSimClo(hitseq_align, hitclone_frequency, 0.0)	
      #   print outSeqMaj, NewT2C2F 		 
         return  outSeqMaj, NewT2C2F  		 
        else:
         return  seq_list, clone_frequency  			

   	
							

	