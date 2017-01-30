from alignments.MegaAlignment import MegaAlignment
from Bio import Phylo

class TreeAnalizer:

    def compare_good_posi_number(self, Initial, After):
        IniCou = self.count_good_posi(Initial)
        AftCou = self.count_good_posi(After)
        if IniCou>AftCou: AfterGood='n'
        else: AfterGood='y'

        return AfterGood	
		
    def count_good_posi(self, Posi_info):	
        Good=0
        for Posi in Posi_info:
             if Posi_info[Posi]==['Good']: Good+=1
        return Good	

    def RmUnresolvedPosi(self, SeqLs, NodeMapLs, Nwk):
         Align = MegaAlignment()	
         CloOrder0, Clo2Seq = Align.name2seq(SeqLs)
         CloOrder=[]
         MaxClo=0		 
         for Clo in CloOrder0:
           if Clo!='#hg19' and Clo.find('#Node')==-1: CloOrder.append(Clo)

         Tree=Nwk.strip()
         Trunk=self.ReadNwkExtBra(Tree)

         Pattern = self.list_clade(NodeMapLs, CloOrder)	   
         GoodC=0
         BacC=0
         MulC=0
	
         ChangePos=[]
         Len=len(Clo2Seq[CloOrder0[0]])
         c=0
         out={}
       
         while c<Len:
            Pat=[]
            Num=0	
            for Clo in CloOrder:       
                if Clo2Seq[Clo][c]=='T': 
                   Pat.append(Clo)
                   Num+=1		   
            Change='n'
            SmallestDiff=9999
            DifNum2ChangeClo2Mut={}
            DifNum2ChangeClo2Wild={}
			
            for Tru in Pattern:
                Miss=[] #back
                for Tclo in Tru:
                    if Pat.count(Tclo)==0: Miss.append(Tclo)
                Extra=[] #parallel
                for Pclo in Pat:
                    if Tru.count(Pclo)==0: Extra.append(Pclo)
    
                if Miss==[]: 
                     if len(Pat)==len(Tru):
                         SmallestDiff=0
                         DifNum2ChangeClo2Mut[len(Miss)]=[Miss,Tru]
                         Change='ToMut'							 
                     else:
                      if len(Extra)<=SmallestDiff:					 
                         SmallestDiff=len(Extra)
                         DifNum2ChangeClo2Wild[len(Extra)]=[Extra,Tru]	
                         Change='ToWild'
                elif Extra==[]: 
                     if len(Pat)==len(Tru):
                         SmallestDiff=0
                         DifNum2ChangeClo2Wild[len(Extra)]=[Extra,Tru]
                         Change='ToWild'						 
                     else:
                      if len(Miss)<=SmallestDiff:					 
                         SmallestDiff=len(Miss)
                         DifNum2ChangeClo2Mut[len(Miss)]=[Miss,Tru]	
                         Change='ToMut'							 
                else:						 
                  if len(Miss)>len(Extra) :
                    if len(Miss)<SmallestDiff:
                         SmallestDiff=len(Miss)
                         DifNum2ChangeClo2Mut[len(Miss)]=[Miss,Tru]
                         Change='ToMut'						 
                  else:
                    if len(Extra)<SmallestDiff:
                         SmallestDiff=len(Extra)
                         DifNum2ChangeClo2Wild[len(Extra)]=[Extra,Tru]
                         Change='ToWild'						 
            if SmallestDiff==0: 
               out[c]=['Good']
               GoodC+=1			   
            else:
		
                if SmallestDiff>2 : ChangePos.append(c)			
                if Change=='ToMut':
                        BacC+=1				
                      	out[c]=['ToMut',DifNum2ChangeClo2Mut[SmallestDiff]]
                if Change=='ToWild':
                        MulC+=1				
                      	out[c]=['ToWild',DifNum2ChangeClo2Wild[SmallestDiff]]				

            c+=1
         out2=[]        	
         for Clo in Clo2Seq:
           if Clo[:5]!='#Node':		 
            Seq=Clo2Seq[Clo]
            NewSeq=Align.ModSeq(Seq,ChangePos,'A',Len)   
            out2+=[Clo, NewSeq]#+'\n'

         GoodC={'Good':GoodC,'Mul':MulC,'Bac':BacC,'TrunkMut':Trunk}

         return GoodC, out2, out

    def find_best_set(self, ID2Good):       
        Best=''
        BestC=0
        BestTrunk=0
        BestBack=999999999
	
        for ID in ID2Good:
          GoodC=ID2Good[ID]['Good']
		
          if GoodC>BestC:
              Best=str(ID)
              BestC=GoodC

        return int(Best)			  
  
    def ReadNwkExtBra(self, OriNwk):

        Nwk=self.RootTree(OriNwk)	
        Len=len(Nwk)
        c=0
        Read='s'
        Name=''	
        Dist=''	
        while c<Len:
           if Nwk[c]=='('or Nwk[c]==',':			
                Read='name'
                Name=''			
           elif Nwk[c]==':': 
                 Read='dist'
                 if Name=='hg19':		 
                    TrunkMut=float(Dist)
                    break				
                 Dist=''
           elif Read=='name': Name+=Nwk[c]
           elif Read=='dist':
              if Nwk[c]==')': Read='s'	   
              else: Dist+=Nwk[c]
           c+=1
        return TrunkMut	
		
    def RootTree(self, OriNwk):
 
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': 'hg19'})
        Phylo.write(trees, 'test1.nwk', "newick")	

        return open('test1.nwk','r').readlines()[0]	
    
    def GetDec(self, Anc,A2D,C2N):
      AncLs=[Anc]
	  
      Dec=[]
      while AncLs!=[]:
        TMP=[]
        for Anc in AncLs:
          Code=Anc in A2D 
          if Code==True:
            if A2D[Anc]==['-','-']: Dec.append(C2N[Anc])		  
            else:  TMP+=A2D[Anc] 

        AncLs=TMP	
      return Dec
  
    def list_clade(self, NodeMapLs, CloOrder):
         Pattern=[]
         D2A=NodeMapLs[0]
         A2D=NodeMapLs[1]
         N2C=NodeMapLs[2]
         C2N=NodeMapLs[3]	 
         for MutAnc in N2C:
            Label=N2C[MutAnc]	
            DecLs=self.GetDec(Label,A2D,C2N) 
            Pat=[]	
            for Clo in CloOrder:
                if DecLs.count(Clo[1:])!=0: 
                     Pat.append(Clo)			 
            Pattern.append(Pat)	
         return Pattern			