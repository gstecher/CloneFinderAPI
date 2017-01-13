from Bio import Phylo

class MegaAlignment():
    """
        Generate TumorSampleProfile replicates
        
        Total read count and SNV frequency is obtained from Original TumorSampleProfile.
        Binomial distribution is used.		
    """
 #   def __init__(self, align_list):
  #      self.M = align_list

    def name2seq(self, align_list):
        self.M = align_list	
        self.clone_order=[]
        self.clone_seq = {}
        Name=''		
        for Seq in self.M:
            if Seq[0]=='#' and Seq!='#MEGA':
                if Seq!='#hg19' and Seq!='#Normal':
                    self.clone_order.append(Seq)
                   				
                    self.clone_seq[Seq]=''
                    Name=Seq
                else: Name=''  					
            elif Name!='':
                self.clone_seq[Name] += Seq	
        return self.clone_order, self.clone_seq	

    def GetSharePosi1(self, sequence, ShareNuc):
        self.name2seq0 = sequence	
        for i in self.name2seq0:
           Name=i
        SharePosi=[]
        Len=len(self.name2seq0[Name])
        c=0
        while c<Len:
            AllMut='y'
            for Ori in self.name2seq0:
              if Ori!='#hg19' and Ori!='#Normal':
               Nuc=self.name2seq0[Ori][c]
               if Nuc!=ShareNuc: AllMut='n'
            if AllMut=='y': 
                 SharePosi.append(c)
            c+=1
        return SharePosi		
    def UpMeg(self, Name2Seq0, NameLs):
            out=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
            for Name in NameLs:
                if Name[0]!='#': Name='#'+Name
                out+=[Name,Name2Seq0[Name]]
            return out #GetOut(Out,out)	
    def RmRedunSeq(self, OriAlign):
        NameOrder, N2Seq = self.name2seq(OriAlign)
     #   print N2Seq		
        out2 = 	['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']	
     #   NameOrder, Name2Seq, out2=Functions.ReadMegSeq(Meg)
        c=0
        Name2IdenLs={}
        SeqNum=len(NameOrder)
        Len=len(N2Seq[NameOrder[0]])
        while c<SeqNum:
            Ref=NameOrder[c]
            RefSeq=N2Seq[Ref]
            Name2IdenLs[Ref]=[Ref]	
            Tc=0	
            while Tc<SeqNum:
                Tar=NameOrder[Tc]
                TarSeq=N2Seq[Tar]	
              #  print Ref, RefSeq				
                DifC=self.CountDifNum(RefSeq,TarSeq)		
                if DifC==0: 
                        Name2IdenLs[Ref].append(Tar)			
                Tc+=1
            c+=1
        Done=[]
        for Name in N2Seq:
          Code=Name in Done
          if Code!=True:
            IdenLs=Name2IdenLs[Name]
            GoodName=''
            for i in IdenLs:
                if GoodName=='': GoodName=i 
                elif i=='#hg19': GoodName=i
                elif GoodName=='#hg19': pass		
                elif i.find('Node')==-1 and i.find('Clu')==-1: GoodName=i	
                elif GoodName.find('Node')==-1 and GoodName.find('Clu')==-1: pass 
                elif GoodName.find('Node')!=-1: GoodName=i 
                elif i.find('Node')!=-1: pass 
                elif i.find('Clu')!=-1 and GoodName.find('Clu')!=-1 and i.find('REP')!=-1 and GoodName.find('REP')==-1: GoodName=i 		
                else: pass
            out2+=[GoodName,N2Seq[Name]]		
            Done+=IdenLs
        return out2#Functions.GetOut(Out2,out2)




		
    def identify_similar_seq(self, Clo2Seq0, CutPro):
            CloOrder, Clo2Seq = self.name2seq(Clo2Seq0)	     
          #  CloOrder=[]
           # for Clo in Clo2Seq:
            #   CloOrder.append(Clo)			
            CombClo=[]
            Done=[]
            for Clo0 in CloOrder:
              Code=Clo0 in Done
              if Code!=True and Clo0!='#hg19':
                Done.append(Clo0)  
                Group=[Clo0]
                Add='y'
                while Add=='y':
                  Add='n'	 
                  for Clo01 in Group:
                   Seq0=Clo2Seq[Clo01]
                   for Clo1 in CloOrder:
                    Code=Clo1 in Done
                    if Code!=True and  Clo1!='#hg19':
                        Seq1=Clo2Seq[Clo1]
                        Len=len(Seq1)
                       # print Len,Seq1,Clo1,Clo2Seq						
                        Dif=self.CountDifNum(Seq0,Seq1)				
                        DifPro=1.0*Dif/Len         			
                        if DifPro<=CutPro:# or Dif<=CutNum: 
                              Group.append(Clo1)
                              Done.append(Clo1)
                              Add='y'  
                CombClo.append(Group)
            return CombClo 
    def make_similar_seq_dic(self, similar_seq_dic):
        Clo2SimCloLs = {}
        for SimCloLs in similar_seq_dic:
            Len1=len(SimCloLs)
           # Len0=Len-1
            c0=0
            while c0 <Len1:
                c1=0
                SimCloLsSub=[]
                while c1<Len1:
                   # if c1!=c0:					
                    SimCloLsSub.append(SimCloLs[c1][1:])
                    c1+=1
                Clo2SimCloLs['T-'+SimCloLs[c0][1:]]=SimCloLsSub						
                c0+=1					
        return 	Clo2SimCloLs	
    def CountDifNum(self, Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif	
    def CountAdditionalMut(self, seq1,seq2):
        Len=len(seq1)
        c=0
        DerMut=0		
        while c<Len:
            if seq1[c]=='T' and seq2[c]=='A': DerMut+=1
            c+=1
        return DerMut			
    def ModSeq(self, CSeq0,ChangePosi,ChanNuc,Len):					  
                      c=0
                      CutCloSeq=''					  
                      while c<Len:
                          Code1=c in ChangePosi						  
                          if Code1==True: CutCloSeq+=ChanNuc
                          else: CutCloSeq+=CSeq0[c]
                          c+=1
                      return CutCloSeq	
    def GetMutPos(self, Seq):
      TMP=[]
      Len=len(Seq)
      c=0
      while c<Len:
        if Seq[c]=='T': TMP.append(c)
        c+=1
      return TMP 


	
    def get_mega_alignment_string(self, SeqLs):
        result = ''
        for item in SeqLs:#self._mega_seqs:
            result += item + "\n"
        return result
    
    def save_mega_alignment_to_file(self, filename, SeqLs):
        destination = open(filename,'w')
        destination.write(self.get_mega_alignment_string(SeqLs))
        destination.close()        	  