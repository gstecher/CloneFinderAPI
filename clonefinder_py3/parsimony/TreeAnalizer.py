from parsimony.MegaMP import MegaMP
from Bio import Phylo
import os
from parsers.nodeMapParser import nodeMapParser
from alignments.MegaAlignment import MegaAlignment
from parsimony.MegaAncestor import MegaAncestor
import os
import tempfile

class TreeAnalizer:
		
    def RootTree(self, OriNwk):
 
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': 'hg19'})
        Phylo.write(trees, 'test1.nwk', "newick")	
        Rooted=open('test1.nwk','r').readlines()[0]	
        os.remove('test1.nwk')	
        os.remove('test.nwk')			
        return Rooted
		
    def MLancetor(self,seqs,Nwk):

        Align = MegaAlignment()
      #  seqs = open(seqsFile,'r').readlines()
        self.CellLs, self.Cell2Seq = Align.name2seq(seqs)
       # print ('h',self.CellLs)
       # print (self.Cell2Seq) 
       # open('A','r').readlines()			
        self.SNVnum = len(self.Cell2Seq[self.CellLs[0]])
      
      #  self.InMeg = Align.AddNormal(seqs)
       
	
        InferAncestor = MegaAncestor()
        InferAncestor.alignment_file = seqs
        InferAncestor.input_tree_file = Nwk	
     			
        self.ancestor_states, self.offspring2ancestor, cell2code, self.code2cell = InferAncestor.retrieve_ancestor_states()
        self.RescaledTree=InferAncestor.get_scaledNWK()		
        self.nodeid2seq={}
       # print (self.ancestor_states)	
        print ('SNV count',self.SNVnum)
     #   open('A','r').readlines()		
        for node in self.ancestor_states:
            Seq=''
            States=	self.ancestor_states[node]			
            c=0
            while c<self.SNVnum:
               # print (node)			
               # print (c)
               # print (States)			
                Nuc=States[c].split('\t')[0]
                Seq+=Nuc
                c+=1 
            self.nodeid2seq[node]=Seq
       # print (RescaledTree)			
        print (self.nodeid2seq)	
        print (self.offspring2ancestor)	
        print (self.code2cell)	
        self.offspring2ancestor_withou_redunSeq=self.RemoveRedun()
    def output_PhyloSigFinderIn(self,initial_seq,RootedNwk,CFin,OutFol,ResDir):
             self.InputFile=CFin[:-4]+'_PSF.input'	
             self.MLancetor(initial_seq,RootedNwk)
             self.sum_posmutInfo(CFin)			 
             self.get_branch_mutations(OutFol,ResDir)
			 

    def sum_posmutInfo(self,CFin):
       	self.Pos2MutInfo={}
        CFin=open(CFin,'r').readlines()[1:]
        c=0
        while c<self.SNVnum:
            In=CFin[c].split('\t')
            if In[2]!='-' and In[3]!='-' and In[4]!='NA':			
           	   self.Pos2MutInfo[c]=In[0]+'\t'+In[1]+'\t1\t'+In[2]+'\t'+In[3]+'\t'+In[4]
            else: self.Pos2MutInfo[c]=''			   
            c+=1			   
    def get_branch_mutations(self,OutFol,ResDir):
       Bra2MutPos={}
       TreeIn='#Tree\n'
       BranchIn='#BranchID\tFile\n'
       AllAncLs=[]	   
       for D in self.offspring2ancestor_withou_redunSeq:
         	   
           Dseq=self.nodeid2seq['Node_'+D]
           A=self.offspring2ancestor_withou_redunSeq[D]
           AllAncLs.append(A)		   
           Aseq=self.nodeid2seq['Node_'+A]
           MutPos=self.get_mutgainPos(Aseq,Dseq)
           DbraID='B'+self.code2cell[D].replace('_','')
           AbraID='B'+self.code2cell[A].replace('_','')
           TreeIn+=	AbraID+'->'+DbraID+'\n'	   
           #print (A+'->'+D)		   
           self.make_MutCountIn(MutPos,OutFol+DbraID+'_MutCount.txt')
           BranchIn+=DbraID+'\t'+ResDir+DbraID+'_MutCount.csv\n'		   
           #open('A','r').readlines()
       AllAncLs=list(set(AllAncLs))
       Root=''	   
       for A in AllAncLs:
          if A not in self.offspring2ancestor_withou_redunSeq:
             if Root!='':
                 print ('more than one root')			 
                 open('A','r').readlines()
             else:	Root=A
       if Root=='':
                 print ('No root')			 
                 open('A','r').readlines()
       else: 
           Rseq=self.nodeid2seq['Node_'+Root]		   
           Aseq='A'*self.SNVnum
           MutPos=self.get_mutgainPos(Aseq,Dseq)
#           DbraID=self.code2cell[D].replace('_','')+'Branch'
           RbraID='B'+self.code2cell[Root].replace('_','')
          # TreeIn+=	AbraID+'->'+DbraID+'\n'	   
           #print (A+'->'+D)		   
           self.make_MutCountIn(MutPos,OutFol+RbraID+'_MutCount.txt')
           BranchIn+=RbraID+'\t'+ResDir+RbraID+'_MutCount.csv\n'
       self.GetOut(self.InputFile,BranchIn+TreeIn)		   
    def make_MutCountIn(self,MutPos,OutMutCountFile):
        out=''
        for Mut in MutPos:
          In=self.Pos2MutInfo[Mut]		
          if In!='':		
            out+=In+'\n'
        self.GetOut(OutMutCountFile,out) 
        self.translate_into_MutCountCSV(OutMutCountFile)
    def translate_into_MutCountCSV(self,Ta):
       # Ta=sys.argv[1] #_MutCount.txt
        Out=Ta[:-4]+'.csv'
        SigOrder='cosmicSig.csv' #please change the path
        SigOrder=open(SigOrder,'r').readlines()[1:]
        Order=[]
        for i in SigOrder:
            Order.append(i.split(',')[0])
        print (Order, len(Order))
        
        Ta=open(Ta,'r').readlines()
        Sig2Count={}
        Sig2PosLs={}		
        for i in Ta:
            i=i.strip().split('\t')
            Pos=i[0]+':'+i[1]			
            Ref=i[3]
            Mut=i[4]
            Tri=i[5]
            Sig=Tri[0]+'['+Ref+'>'+Mut+']'+Tri[2]
            if Sig not in Sig2Count:#.has_key(Sig)!=True: 
                Sig2Count[Sig]=0
                Sig2PosLs[Sig]=[]				
            Sig2Count[Sig]+=1
            Sig2PosLs[Sig].append(Pos)			
        out=''
        for Sig in Order:
            if Sig in Sig2Count:#.has_key(Sig)==True: 
                Count=str(Sig2Count[Sig])
                PosLs=Sig2PosLs[Sig]
                PosLs=';'.join(PosLs)				
            else: 
               Count='0'
               PosLs='NA'			
            out+=Sig+','+Count+','+PosLs+'\n'
        OutF=open(Out,'w')
        OutF.write(out)
        OutF.close()	



	
        			
    def get_mutgainPos(self,Ancseq,Desseq):
       Len=len(Ancseq)
       MutPos=[]
       c=0
       while c<Len:
          if Ancseq[c]=='A' and Desseq[c]=='T': MutPos.append(c)
          c+=1	
       return MutPos		  
       	   
		   
    def RemoveRedun(self):
        Old2New={}
        TipIDls=[]	
        Anc2IdenTip={}		
        for dec in self.offspring2ancestor:
          if dec!='Des1' and dec!='Des2':	
		  
            Dseq=self.nodeid2seq['Node_'+dec]
            if self.code2cell[dec].find('Node_')==-1: 
             # if self.code2cell[dec]=='hg19':Tip='n'
              #else:			  
                Tip='y'
                TipIDls.append(dec)	
			
               # print (dec)				
            else: Tip='n'
            anc=self.offspring2ancestor[dec]			
            Aseq=self.nodeid2seq['Node_'+anc]
            MutC=0
            c=0
            while c<self.SNVnum:
               if Dseq[c]=='T' and Aseq[c]=='A': MutC+=1	
               c+=1	
            if MutC==0: 
               if dec in Old2New:
                   print (dec)
                   print (Old2New[dec])	
                   open('a','r').readlines()				   
               else:			   
                  Old2New[dec]=anc	
                  if Tip=='y': Anc2IdenTip[anc]=dec					  
        NewDec2Anc={}
        print ('Old2New')		
        print (Old2New)   
        Anc2Tip={}		
        for dec in self.offspring2ancestor:
          if dec!='Des1' and dec!='Des2':		
            if self.code2cell[dec].find('Node_')==-1:
         			
                 Anc=self.get_NewID(Old2New,self.offspring2ancestor[dec])
                 NewDec2Anc[dec]=Anc
                 if dec in Old2New:
                   if self.code2cell[dec]!='hg19':				 
                     Anc2Tip[Anc]=dec				 
            else:
               if dec not in Old2New:
                 Anc=self.get_NewID(Old2New,self.offspring2ancestor[dec])
                 NewDec2Anc[dec]=Anc             			   
        print (NewDec2Anc)
        print (Anc2Tip)		
        NewDec2Anc1={}
        for D in NewDec2Anc:
          if self.code2cell[D]!='hg19':		
           A=NewDec2Anc[D]	
           if A in Anc2Tip:
              if Anc2Tip[A]==D: pass
              else: 
                if D in Anc2Tip:			  
                     NewDec2Anc1[Anc2Tip[D]]=Anc2Tip[A]
                else: NewDec2Anc1[D]=Anc2Tip[A]	
           else: 
                if D in Anc2Tip:			  
                     NewDec2Anc1[Anc2Tip[D]]=A
                else: NewDec2Anc1[D]=A
        print (NewDec2Anc1)	
        return NewDec2Anc1		
          			  
        			
    def get_NewID(self,Old2New,Target):
               while Target in Old2New:
                   Target=Old2New[Target]	
               return Target
    def GetOut(self,File,In):
        OutF=open(File,'w')
        OutF.write(In)
        OutF.close()	
    def get_tree_with_branchLen(self,ID):
        Align=MegaAlignment()	
        self.GetOut(ID+'.nwk',self.RescaledTree)
        SeqLs=Align.UpMeg(self.nodeid2seq, [])		
        Align.save_mega_alignment_to_file(ID+'_NodeSeq.meg', SeqLs)