import scipy.optimize
import numpy

class CloneFrequencyComputer(object):
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, MEGAalignment, tsp_list):
        self.M = MEGAalignment 
        self.tsp_list = tsp_list
######		
        self.v_obs = {}
        for profile in self.tsp_list: 
            tumor = profile.name			
            v_list=[]				
            for read_count in profile:
                 v_list.append(read_count.alt_frequency())
            self.v_obs[tumor]=v_list	
#####			
    #    print self.v_obs		
        self.Min=''
        self.out = 'Tumor'		
######		
        self.clone_order=[]
        clone_seq = {}
        Name=''		
        for Seq in self.M:
            if Seq[0]=='#' and Seq!='#MEGA':
                if Seq!='#hg19' and Seq!='#Normal':
                    self.clone_order.append(Seq)
                    self.out += '\t' + Seq[1:]					
                    clone_seq[Seq]=''
                    Name=Seq					
            elif Name!='':
                clone_seq[Name] += Seq	
######				
       # print clone_seq
        self.out += '\n'		
        snv_num = len(clone_seq[Name])
        snv=0
        while snv < snv_num:
            for Name in self.clone_order:
                if Name!='#hg19' and Name!='#Normal':
                    if clone_seq[Name][snv]=='A' or clone_seq[Name][snv]=='a': value='0'
                    elif clone_seq[Name][snv]=='T' or clone_seq[Name][snv]=='t': value='1'
                    else: value='0.5'					
                    self.Min+=value+' '
                  #  print self.Min					
            self.Min=self.Min[:-1]+'; '
            snv+=1	     		
        self.Min=self.Min[:-2]
      #  print self.Min		
        self.Min=numpy.matrix(self.Min)
     #   print self.Min		
    def regress(self):
       # A = numpy.matrix(self.Min)
        Tumor2Clone2Freq={}	   
        for tumor in self.v_obs:	
            clone_frequency = scipy.optimize.nnls(self.Min,self.v_obs[tumor])[0] #float list
           # print tumor, clone_frequency	,self.clone_order		
            Tumor2Clone2Freq['T-'+tumor]={}
            clone_id = 0
            for clone in self.clone_order:
                Tumor2Clone2Freq['T-'+tumor][clone[1:]] = clone_frequency[clone_id]*2
                clone_id += 1
        return Tumor2Clone2Freq				
            			
