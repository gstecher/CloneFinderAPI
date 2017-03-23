from output.CloneFrequencyWriter import CloneFrequencyWriter
import scipy.optimize
import numpy

class CloneFrequencyComputer(object):
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, MEGAalignment, tsp_list, cutoff_clone_frequency):
       # print MEGAalignment	
        self.M = MEGAalignment 
        self.tsp_list = tsp_list
        self.CutOff = cutoff_clone_frequency			
        self.v_obs = {}
        for profile in self.tsp_list: 
            tumor = profile.name			
            v_list=[]				
            for read_count in profile:
                 v_list.append(read_count.alt_frequency())
            self.v_obs[tumor]=v_list
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
			
        self.out += '\n'
        self.Min = 	self.make_Min(self.clone_order, self.clone_seq)	
  #  @property 
   # def CNV_file(self):
    #    return self.CNV_file
    
    #@CNV_file.setter
    #def CNV_file(self, value):
     #   self.CNV_file = value

 #   @property 
  #  def Clone_frequency_file(self):
   #     return self._Clone_frequency_file
    
#    @Clone_frequency_file.setter
 #   def Clone_frequency_file(self, value):
  #      self._Clone_frequency_file = value	

   # @property 
    #def v_est(self):
     #   return self.v_est
    
 #   @v_est.setter
  #  def v_est(self, value):
   #     self.v_est = value	

	
		
	
		
		
    def make_Min(self, clone_order, clone_seq):	
        Min=''		
        snv_num = len(clone_seq[clone_order[0]])
        snv=0
      #  print clone_order		
        while snv < snv_num:
            for Name in clone_order:
               # if snv==0:print Name, clone_seq[Name]			
                if Name!='#hg19' and Name!='#Normal':
                 #    if Name.find('CNV')==-1:				
                        if clone_seq[Name][snv]=='A' or clone_seq[Name][snv]=='a': value='0'
                        elif clone_seq[Name][snv]=='T' or clone_seq[Name][snv]=='t': value='1'
                        else: value='0.5'					
                        Min+=value+' '
                #     else: pass
          #              tumor=Name.split('CNV')[-1]
           #             cloneFreq = self._Clone_frequency_file[tumor][Name.split('CNV')[0]][1:]/2						
            #            CNV=self._CNV_file[tumor][snv]
             #           Obs=self.v_obs[tumor][snv]
              #          Est=self.v_est[tumor][snv]	
               #         if clone_seq[Name][snv]=='A' or clone_seq[Name][snv]=='a': value='0'
                #        elif clone_seq[Name][snv]=='T' or clone_seq[Name][snv]=='t': 
                 #            if CNV=='undecided' or CNV=='normal' or Obs==Est: value='1'
                  #           elif Obs>cloneFreq: ##add more mutant copy
                   #              if CNV=='loss': value='2'
                    #             elif CNV=='gain':  value='1.333'
                     #            elif CNV=='LOH':  value='2'
                      #           else: print '1??'
							 
                       #      elif cloneFreq<Obs: ##add more wild copy
                        #         if CNV=='loss': value='0'
                         #        elif CNV=='gain':  value='0.666'
                          #       elif CNV=='LOH':  value='0'
                           #      else: print '1??'
                            # else: 
                             #   print '???'							 
                    #    else: value='0.5'					
                     #   Min+=value+' '
						
			
            Min=Min[:-1]+'; '
            snv+=1	     		
        Min=Min[:-2]
       # print Min	
        Min=numpy.matrix(Min)
        return Min		

  #  def do_ScipyNnls(self,tumor):
   #     #Res=[]
    #    Rep=30
     #   c=0
      #  while c<Rep:		
       #     clone_frequency = scipy.optimize.nnls(self.Min,self.v_obs[tumor])[0] #float list
        #    print tumor,clone_frequency			
         #   if c==0 :
          #    Res=clone_frequency
           #   CloNum=len(Res)			  
            #else:
             # Clo=0
              #while Clo<CloNum:
               #   Res[Clo]+=clone_frequency[Clo]
                #  Clo+=1
 #           c+=1
  #      Res_up=[]
   #     for Fre in Res:
    #        Res_up.append(Fre/Rep)
     #   return Res_up			
    def regress(self):
        Tumor2Clone2Freq={}	   
        for tumor in self.v_obs:
	
            clone_frequency = scipy.optimize.nnls(self.Min,self.v_obs[tumor])[0] #self.do_ScipyNnls(tumor)#float list	
            Tumor2Clone2Freq['T-'+tumor]={}
            clone_id = 0
            for clone in self.clone_order:
             # if clone.find('CNV')==-1:
               if clone!='#hg19' and clone!='#Normal':			 
                if 	clone_frequency[clone_id] > self.CutOff	: Fre = clone_frequency[clone_id]
                else: Fre = 0				
                Tumor2Clone2Freq['T-'+tumor][clone[1:]] = Fre*2
                clone_id += 1
        Out = CloneFrequencyWriter(Tumor2Clone2Freq, self.M, 0)	
        OutAlign, OutCloFre = Out.get_hitclone()	
        return OutAlign, OutCloFre

    def regress_withCNV(self, CNV_information):
        self._Clone_frequency_file=CNV_information
      #  print 'hhhhhh',self._Clone_frequency_file		
        self.regress()		

		
    def regress_hit_clone(self, hit_clone_table):
        Tumor2Clone2Freq={}	   
        for tumor in self.v_obs:
            HitCloLs=[]
            if hit_clone_table.has_key(tumor)==True: Tu=tumor
            else: Tu='T-'+tumor
            CloFre=hit_clone_table[Tu]
            for Clo in CloFre:
                 if CloFre[Clo]>0: HitCloLs.append('#'+Clo)	
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


        	
            			
