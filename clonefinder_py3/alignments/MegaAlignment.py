
class MegaAlignment():

    def name2seq(self, align_list):
        self.M = align_list	
        self.clone_order=[]
        self.clone_seq = {}
        Name=''		
        for Seq in self.M:
          Seq=Seq.strip()		
          if Seq!='':		
            if Seq[0]=='#' and Seq!='#MEGA':
                if Seq!='#hg19' and Seq!='#Normal':
                    self.clone_order.append(Seq)
                   				
                    self.clone_seq[Seq]=''
                    Name=Seq
                else: Name=''  					
            elif Name!='':
                self.clone_seq[Name] += Seq	
        return self.clone_order, self.clone_seq	

		
    def UpMeg(self, Name2Seq0, NameLs):
            if NameLs==[]:
                 for Name in Name2Seq0:
                     NameLs.append(Name)			 
            out=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
            for Name in NameLs:
                if Name in Name2Seq0: Seq=Name2Seq0[Name]
                else: Seq=Name2Seq0['#'+Name]				
                if Name[0]!='#': Name='#'+Name
                				
                out+=[Name,Seq]#Name2Seq0[Name]]
            return out 
						
    def remove_redund_seqs(self, Meg):
        	
        print('removing redundant seqs...')
        	
        NameOrder, Name2Seq=self.name2seq(Meg)
        Name2IdenLs={}
        for Ref in NameOrder:		
      
            RefSeq=Name2Seq[Ref]
            Name2IdenLs[Ref]=self.find_redundant(RefSeq,Name2Seq)			

        Done=[]
        Nonredun_seqdic={}		
        for Name in NameOrder:
         
          if Done.count(Name)==0:
            IdenLs=Name2IdenLs[Name]
            Nonredun_seqdic[Name]=Name2Seq[Name]				
            Done+=IdenLs
        out2=self.UpMeg(Nonredun_seqdic,[])        
        return out2 
		
    def CountDifNum(self, Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif	
			
    def find_redundant(self,Seq,Seq_dic):
           Redun=[]
           for Name in Seq_dic:
                Seq1=Seq_dic[Name]
                DifNum=self.CountDifNum(Seq,Seq1)
                if DifNum==0: Redun.append(Name)
           return Redun				
			

		
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



    def AddNormal(self, seqs): 
        CellOrder, Cell2Seq = self.name2seq(seqs)	
        Len=len(Cell2Seq[CellOrder[0]])	
        seqs+=['#Normal\n',('A'*Len)+'\n']
        return seqs		

