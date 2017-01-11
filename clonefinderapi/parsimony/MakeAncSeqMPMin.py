class MakeAncSeqMPMin(object):
    
    def __init__(self):
        self._filenames = []
        self.ID = ''
        
    def GetRel(self, Ls,SetN):
        A2D={}
        D2A={}
        Index2Lab={}	
        SeqLs={}
        Anc2Seq={}
        TaxaLs=[]	
        outNodeMap='Index\tLabel\tDec_1\tDec_2\n'	
        for Li in Ls:
           if Li!=[]:        	   
            Anc=Li[0]
            NodeN=Li[1]
            Dec1=Li[2]
            Dec2=Li[3]
            A2D[Anc]=[Dec1,Dec2]
            D2A[Dec1]=Anc
            D2A[Dec2]=Anc
            Index2Lab[Anc]=NodeN
            if NodeN=='hg19': NormalIndex=Anc

            if NodeN=='-':
                SeqLs[Anc]=Li[4:]
                Anc2Seq[Anc+'.0']=''
                Len=len(Li[4:])
            elif NodeN!='hg19':	
                SeqIn=''
                Nls=Li[4:]			
                for N in Nls: SeqIn+=N
                Anc2Seq[NodeN]=SeqIn
                TaxaLs.append(NodeN)			
        Root=self._get_root(D2A,Anc)
        RootDec=A2D[Root]	
        Done=[D2A[NormalIndex],NormalIndex]	
        DecLs=[D2A[NormalIndex]]
        Down=[]	
        while DecLs!=[]:
          NewDecLs=[]	
          for Dec in DecLs:	
           InIndex=Dec	
           Code=Dec in D2A
           if Code==True:        	   
            Dec1=D2A[Dec]
            if Dec1!=Root:		  
              Decs=A2D[Dec]
              for Dec0 in Decs:
                  Code=Dec0 in Done		  
                  if Dec0!=Dec and Code!=True: Dec2=Dec0
              if InIndex!=NormalIndex and Dec1!=NormalIndex and Dec2!=NormalIndex: outNodeMap+=InIndex+'\t'+Index2Lab[InIndex]+'\t'+Dec1+'\t'+Dec2+'\n'
              Done.append(InIndex)		  
              NewDecLs.append(Dec1)
              Down.append(Dec2)
            else: #root
              for Dec0 in RootDec:
                  if Dec0!=Dec:  Dec1=Dec0	

              Decs=A2D[Dec]
              for Dec0 in Decs:
                  Code=Dec0 in Done		  
                  if Code!=True: Dec2=Dec0		  
              if InIndex!=NormalIndex and Dec1!=NormalIndex and Dec2!=NormalIndex: outNodeMap+=InIndex+'\t'+Index2Lab[InIndex]+'\t'+Dec1+'\t'+Dec2+'\n'
              Done.append(InIndex)
              Down+=[Dec1,Dec2]
          DecLs=NewDecLs
        if Root==D2A[NormalIndex]: #position of root is correct
            for Dec in RootDec:
                if Dec!=NormalIndex: Start=[Dec]   
            Down+=Start   		
            while Start!=[]:
                  NewStart=[]
                  for Item in Start:
                       Code=Item in A2D
                       if Code==True:NewStart+=A2D[Item]

                  Down+=NewStart
                  Start=NewStart			  
        while Down!=[]:
           NewDown=[]
           for Anc in Down:
              Code=Anc in A2D
              if Code==True:
                 Decs=A2D[Anc]	 
                 if Anc!=NormalIndex and Decs[0]!=NormalIndex and Decs[1]!=NormalIndex: outNodeMap+=Anc+'\t'+Index2Lab[Anc]+'\t'+Decs[0]+'\t'+Decs[1]+'\n'
                 NewDown+=Decs
           Down=NewDown			 			
        self._get_out(self.ID+'_NodeMap'+SetN+'.txt',outNodeMap)
        self._get_out(self.ID+'_NodeMap.txt',outNodeMap)	
        return TaxaLs,'AA',SeqLs,Anc2Seq,Len,Root    
        
    def get_best_alignment(self, files, ID, remove_redundant):        
        SetID=0
        self.ID = ID
        out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
        file_index = 0
        while file_index < len(files):
            AncFile = files[file_index]
            print 'processing file: ' + AncFile
            AncFile=open(AncFile,'r').readlines()
            Lines=[]
            Add='n'	
            for i in AncFile:
                if i.find('Index ')!=-1: Add='y'
                elif Add=='y': 
                    i=i.strip().split(' ')
                    In=[]			
                    for ii in i:
                        if ii!='': In.append(ii)
                    Lines.append(In)				
            TaxaLs,Anc2Cou,AncSeqLs,Anc2Seq,Len,Root=self.GetRel(Lines,str(SetID))
            c=0  	
            while c<Len: #add node sequences
                for Anc in AncSeqLs:
                  if Anc!=Root:		
                    Nuc=AncSeqLs[Anc][c]           			
                    if Nuc=='T':Anc2Seq[Anc+'.0']+='T'			
                    else:Anc2Seq[Anc+'.0']+='A'		               
                c+=1							
            for Anc in Anc2Seq:
              if Anc.split('.')[0]!=Root:
               Code=Anc in TaxaLs	  
               if Code==True:out+='#'+Anc+'\n'+Anc2Seq[Anc]+'\n'	
               else:	   
                out+='#Node'+Anc.replace('.','S')+'T'+str(SetID)+'\n'+Anc2Seq[Anc]+'\n'	
            SetID+=1
            file_index += 1
        if remove_redundant:
            out = self._remove_redund_seqs(out)
        return out  
    
    def _read_mega_seq(self, MegStr): 
        print 'reading mega file...'
        Meg = MegStr.splitlines()
        Read='s'
        out2=''
        NameOrder=[]
        Name2Seq={}
        for i in Meg:
            if i.strip() == '':
                continue
            if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
                Read='n'
                Name=i.strip()
                NameOrder.append(Name)
                Name2Seq[Name]=''
            elif Read=='n': Name2Seq[Name]+=i.strip()
            elif Read=='s': out2+=i + "\n"
        out2 += "\n"
        return NameOrder, Name2Seq, out2    
    
    def _remove_redund_seqs(self, Meg):
        print 'removing redundant seqs...'
        NameOrder, Name2Seq, out2=self._read_mega_seq(Meg)
        c=0
        Name2IdenLs={}
        SeqNum=len(NameOrder)
        Len=len(Name2Seq[NameOrder[0]])
        while c<SeqNum:
            Ref=NameOrder[c]
            RefSeq=Name2Seq[Ref]
            Name2IdenLs[Ref]=[Ref]	
            Tc=0	
            while Tc<SeqNum:
                Tar=NameOrder[Tc]
                TarSeq=Name2Seq[Tar]	
                DifC=self._count_diff_num(RefSeq,TarSeq)		
                if DifC==0: 
                        Name2IdenLs[Ref].append(Tar)			
                Tc+=1
            c+=1
        Done=[]
        for Name in Name2Seq:
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
            out2+=GoodName+'\n'+Name2Seq[Name]+'\n'		
            Done+=IdenLs
        
        return out2 

    def _count_diff_num(self, Seq0,Seq1):
        Len=len(Seq0)		
        Dif=0
        c=0
        while c<Len:
            if Seq0[c]!=Seq1[c]: Dif+=1
            c+=1
        return Dif    
    
    def _get_root(self, Dec2Anc,Node):			  
        Code=Node in Dec2Anc				
        while Code==True: 
            Node=Dec2Anc[Node]
            Code=Node in Dec2Anc				
        return Node	    
    
    def _get_out(self, File,In):
        OutF=open(File,'w')
        OutF.write(In)
        OutF.close()    