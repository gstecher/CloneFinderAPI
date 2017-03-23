from output.CloneFrequencyWriter import CloneFrequencyWriter
from alignments.MegaAlignment import MegaAlignment
from regression.CloneFrequencyComputer import CloneFrequencyComputer
from decomposition.SNPClusterGenerator import SNPClusterGenerator
import scipy.optimize
import numpy
import scipy

class CloneFrequencyComputer_cnv(object):
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, seqs_with_ancestor, tsp_list, CNV_info, option_b, freq_cutoff):
        Align=MegaAlignment()		
        clone_frequency = CloneFrequencyComputer(seqs_with_ancestor, tsp_list, freq_cutoff)
			
        seqs, CloneFreq = clone_frequency.regress()
     #   CloFreAna.save_frequency_table_to_file('Original_freq.txt', CloneFreq, [])				
			
     #   print 'original',CloneFreq			
        clusters = SNPClusterGenerator(seqs_with_ancestor, tsp_list, CloneFreq, option_b, freq_cutoff)
        Obs_SNVfreq, Est_SNVfreq = clusters.estimatedSNV_save_to_file()	#{tumor:[fre, fre, ...],tumor:[fre..]..}	
        self.M = Align.add_cnv_genotype(tsp_list, CloneFreq, seqs_with_ancestor)	

       # clone_frequency_cnv = CloneFrequencyComputer_cnv(new_clone_seq, tsp_list, params.freq_cutoff,CNV_information,CloFre_seqs_with_ancestor,Obs_SNVfreq, Est_SNVfreq)
			
###########	
       # self.M = MEGAalignment 
        self.tsp_list = tsp_list
        self.CutOff = freq_cutoff
     #   print 'clone frequencies wihtout CNV adjustment',CloneFreq		
        self._Clone_frequency_file=CloneFreq
        self.v_est=Est_SNVfreq		
        self.v_obs = {}
        self.readCount = {}		
        self._CNV_file = CNV_info		
        for profile in self.tsp_list: 
            tumor = profile.name			
            v_list=[]	
            ReadCount_list=[]			
            for read_count in profile:
                 v_list.append(read_count.alt_frequency())
                 ReadCount_list.append([read_count.num_ref,read_count.num_alt])				 
            self.v_obs[tumor]=v_list
            self.readCount[tumor]=ReadCount_list			
#############
        		
        self.out = 'Tumor'			
        self.clone_order=[]
        self.clone_seq = {}
        Name=''		
        for Seq in self.M:
              if Seq!='':		
                if Seq[0]=='#' and Seq!='#MEGA':
                    if Seq!='#hg19' and Seq!='#Normal':
                        self.clone_order.append(Seq)
                        self.out += '\t' + Seq[1:]					
                        self.clone_seq[Seq]=''
                        Name=Seq
                    else: Name=''				
                elif Name!='':
                    self.clone_seq[Name] += Seq	
        self.Len=len(Seq)			
        self.out += '\n'
        self.Min = 	self.make_Min(self.clone_order, self.clone_seq)			
		
    def make_Min(self, clone_order, clone_seq):	
        Min=''		
        snv_num = len(clone_seq[clone_order[0]])
        snv=0
      #  print clone_order		
        while snv < snv_num:
            for Name in clone_order:
               # if snv==0:print Name, clone_seq[Name]			
                if Name!='#hg19' and Name!='#Normal':
                     if Name.find('CNV')==-1:				
                        if clone_seq[Name][snv]=='A' or clone_seq[Name][snv]=='a': value='0'
                        elif clone_seq[Name][snv]=='T' or clone_seq[Name][snv]=='t': value='1'
                        else: value='0.5'					
                        Min+=value+' '
                     else:
                        tumor=Name.split('CNV')[-1]
                        cloneFreq = self._Clone_frequency_file['T-'+tumor][Name.split('CNV')[0][1:]]/2						
                        CNV=self._CNV_file[tumor][snv]
                        Obs=self.v_obs[tumor][snv]
                        Est=self.v_est[tumor][snv]	
                        if clone_seq[Name][snv]=='A' or clone_seq[Name][snv]=='a': value='0'
                        elif clone_seq[Name][snv]=='T' or clone_seq[Name][snv]=='t':
                             Mutant=self.readCount[tumor][snv][1]
                             Total=	self.readCount[tumor][snv][0]+self.readCount[tumor][snv][1]						 
                             BinomialP=scipy.stats.binom_test(Mutant, n=Total, p=cloneFreq)
                           #  print Mutant,Total,cloneFreq,BinomialP							 
                             if CNV=='undecided' or CNV=='normal' or Obs==cloneFreq: value='1'
                             elif Obs>cloneFreq and BinomialP<0.05: ##add more mutant copy
                                 if CNV=='loss': value='2'
                                 elif CNV=='gain':  value='1.333'
                                 elif CNV=='LOH':  value='2'
                                 else: print '1??'
							 
                             elif cloneFreq>Obs and BinomialP<0.05: ##add more wild copy
                                 if CNV=='loss': value='0'
                                 elif CNV=='gain':  value='0.666'
                                 elif CNV=='LOH':  value='0'
                                 else: print '1??'
                             else: 
                                value='1'							 
                              #  print '???'	,CNV,Obs,cloneFreq						 
                        else: value='0.5'					
                        Min+=value+' '
						
			
            Min=Min[:-1]+'; '
            snv+=1	     		
        Min=Min[:-2]
       # print Min
       # print Min	   
        Min=numpy.matrix(Min)
        #print Min		
        return Min		

    def do_ScipyNnls(self,tumor):
        #Res=[]
        Rep=30
        c=0
        while c<Rep:		
            clone_frequency = scipy.optimize.nnls(self.Min,self.v_obs[tumor])[0] #float list
           # print tumor,clone_frequency			
            if c==0 :
              Res=clone_frequency
              CloNum=len(Res)			  
            else:
              Clo=0
              while Clo<CloNum:
                  Res[Clo]+=clone_frequency[Clo]
                  Clo+=1
            c+=1
        Res_up=[]
        for Fre in Res:
            Res_up.append(Fre/Rep)
        return Res_up			
    def regress(self):
        Tumor2Clone2Freq={}	   
        for tumor in self.v_obs:
	
            clone_frequency = scipy.optimize.nnls(self.Min,self.v_obs[tumor])[0] #self.do_ScipyNnls(tumor)#float list	
            Tumor2Clone2Freq['T-'+tumor]={}
            clone_id = 0
            for clone in self.clone_order:
                if 	clone_frequency[clone_id] > self.CutOff	: Fre = clone_frequency[clone_id]
                else: Fre = 0				
                Tumor2Clone2Freq['T-'+tumor][clone[1:]] = Fre*2
                clone_id += 1
        Out = CloneFrequencyWriter(Tumor2Clone2Freq, self.M, 0)	
        OutAlign, OutCloFre = Out.get_hitclone()	
        return OutAlign, OutCloFre

    def snvGenotype_after_regress(self):
        MegSeq, Tu2CloFre = self.regress()
        Align = MegaAlignment()
        NameOrder, Clo2Seq = Align.name2seq(MegSeq)
        NewMeg=[]
        NewMeg_seq=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ','#hg19','A'*self.Len]
        NewCloFre_dic={}		
        for Name in NameOrder:
            Seq=Clo2Seq[Name]		 
            Name=Name.split('CNV')[0]
            if NewMeg.count(Name)==0:
                NewMeg.append(Name)
                NewMeg_seq+=[Name,Seq]
        for Tumor in Tu2CloFre:
            NewCloFre={}
            CloFre=Tu2CloFre[Tumor]
            for Clo in CloFre:
                Fre=CloFre[Clo]
                Clo=Clo.split('CNV')[0]
                if NewCloFre.has_key(Clo)!=True: NewCloFre[Clo]=0
                NewCloFre[Clo]+=Fre
            NewCloFre_dic[Tumor]=NewCloFre
        return 	NewMeg_seq, NewCloFre_dic		

		
    def regress_hit_clone(self, hit_clone_table):
        Tumor2Clone2Freq={}	   
        for tumor in self.v_obs:
            HitCloLs=[]
            if hit_clone_table.has_key(tumor)==True: Tu=tumor
            else: Tu='T-'+tumor
            CloFre=hit_clone_table[Tu]
            for Clo in CloFre:
                 if CloFre[Clo]>0: 
                     HitCloLs.append('#'+Clo)
                     if self.M.count('#'+Clo+'CNV'+tumor+'\n')!=0:
                            print 'tumor clone'
                            HitCloLs.append('#'+Clo+'CNV'+tumor)							
            Min=self.make_Min(HitCloLs, self.clone_seq)		
           # print 	len(Min),'\n',len(self.v_obs[tumor])		
            clone_frequency = scipy.optimize.nnls(Min,self.v_obs[tumor])[0] #self.do_ScipyNnls(tumor)#float list	
            Tumor2Clone2Freq['T-'+tumor]={}

            TotCloFre=0
            while TotCloFre==0:	
              clone_id = 0			
              for clone in HitCloLs:
                if 	clone_frequency[clone_id] > self.CutOff	: Fre = clone_frequency[clone_id]
                else: Fre = 0				
                Tumor2Clone2Freq['T-'+tumor][clone[1:]] = Fre*2
                TotCloFre+=Fre				
                clone_id += 1
        Out = CloneFrequencyWriter(Tumor2Clone2Freq, self.M, 0)	
        OutAlign, OutCloFre = Out.get_hitclone()	
        return OutAlign, OutCloFre    	


    def snvGenotype_after_regress_hitclone(self,hit_clone_table):
        MegSeq, Tu2CloFre = self.regress_hit_clone(hit_clone_table)
        Align = MegaAlignment()
        NameOrder, Clo2Seq = Align.name2seq(MegSeq)
        NewMeg=[]
        NewMeg_seq=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ','#hg19','A'*self.Len]
        NewCloFre_dic={}		
        for Name in NameOrder:
            Seq=Clo2Seq[Name]		 
            Name=Name.split('CNV')[0]
            if NewMeg.count(Name)==0:
                NewMeg.append(Name)
                NewMeg_seq+=[Name,Seq]
        for Tumor in Tu2CloFre:
            NewCloFre={}
            CloFre=Tu2CloFre[Tumor]
            for Clo in CloFre:
                Fre=CloFre[Clo]
                Clo=Clo.split('CNV')[0]
                if NewCloFre.has_key(Clo)!=True: NewCloFre[Clo]=0
                NewCloFre[Clo]+=Fre
            NewCloFre_dic[Tumor]=NewCloFre
        return 	NewMeg_seq, NewCloFre_dic	        	
            			
