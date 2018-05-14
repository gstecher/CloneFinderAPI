from alignments.MegaAlignment import MegaAlignment
from Bio import Phylo
#from ete2 import Tree
import os

class TreeAnalizer:
    def Get_branch_lenghth(self, Tree, Clo): #withou #
                      print Tree	
                      Posi=Tree.find(Clo+':')+len(Clo)+1
                      print Tree[Posi]
                      Go='y'
                      BraLen=''					  
                      while Go=='y':
                         BraLen+=Tree[Posi]
                         if Tree[Posi]==',' or Tree[Posi]==')': Go='n'
                         Posi+=1	
                      print Clo,BraLen 
                      return float(BraLen[:-1])		
    def compare_good_posi_number(self, Initial, After, IniSeq_buil, AftSeq_buil):
        IniCou = self.count_good_posi(Initial)
        AftCou = self.count_good_posi(After)
        Align=MegaAlignment()
        CloLs,IniSeq_dic=Align.name2seq(IniSeq_buil)
        CloLs,AftSeq_dic=Align.name2seq(AftSeq_buil)		
        ShareIni=Align.GetSharePosi1(IniSeq_dic,'A')
        ShareAft=Align.GetSharePosi1(AftSeq_dic,'A')		
        print IniCou,AftCou	,len(ShareIni),len(ShareAft)	
        if IniCou>AftCou or len(ShareAft)>len(ShareIni): AfterGood='n'
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
         out2=['#MEGA','!Title SNVs;','!Format datatype=dna;']        	
         for Clo in Clo2Seq:
           if Clo[:5]!='#Node':
              Seq=Clo2Seq[Clo]
              NewSeq=''
              c=0
              while c<Len:
                  if out[c]==['Good']: NewSeq+=Seq[c]
                  else: NewSeq+='A'
                  c+=1
              out2+=[Clo,NewSeq]				  
		 
         #  if Clo[:5]!='#Node':		 
          #  Seq=Clo2Seq[Clo]
           # NewSeq=Align.ModSeq(Seq,ChangePos,'A',Len)   
            #out2+=[Clo, NewSeq]#+'\n'

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
        Rooted=open('test1.nwk','r').readlines()[0]	
        os.remove('test1.nwk')	
        os.remove('test.nwk')			
        return Rooted
    def RootTree_cnv(self, OriNwk, Root):

        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()

        trees = list(Phylo.parse('test.nwk', 'newick'))
		
        for tree in trees:
             tree = tree.root_with_outgroup({'name': Root})
		 
        Phylo.write(trees, 'newtree.nwk', "newick")

        Tree=open('newtree.nwk','r').readlines()[0].strip()
        Len=len(Tree)
        Posi=Tree.find(','+Root+':0.00000')
        PosRev=-1*(Len-Posi)
        
        LastBraLen=''
        Rm=''
        while Tree[PosRev]!=':':
            LastBraLen+=Tree[PosRev]
            PosRev=PosRev-1
        BraLen= LastBraLen[::-1]
       # print BraLen	
        NewTree='('+Root+':'+BraLen+Tree[2:].replace('):'+BraLen+Root+':0.00000','')+'\n'
      #  print NewTree		
        return NewTree	
    def PruneTree(self, OriginalNwk, RmLs):
        t=Tree(OriginalNwk)#("(((A:5,B:5)90:2,C:4)25:3,D:10);")
        t.prune(RmLs)#(["A","C","D"])	
        return t.write(format=1)

    def clone_to_tumor_phylogeny(self, OriginalNwk, Tu2CloFre, CloSeqLs):
        KeepLs=['hg19']
        Keep2TuLs={'hg19':[]}		
        Align=MegaAlignment()
        CloOr, CloSeq = Align.name2seq(CloSeqLs)	
        print Tu2CloFre		
        for Tu in Tu2CloFre:
            CloFre=Tu2CloFre[Tu]
            CloLs=[]			
            for Clo in CloFre:
                 if CloFre[Clo]>0: CloLs.append(Clo)
            LarClo=''
            LarMut=0			
            for Clo0 in CloLs:
                Seq0=CloSeq['#'+Clo0]
                MutC=len(Align.GetMutPos(Seq0))
                if MutC>LarMut:
                      LarMut=MutC
                      LarClo=Clo0	
                Keep='y'					  
                for Clo1 in CloLs:
                    if Clo0!=Clo1:
                        Seq1=CloSeq['#'+Clo1]
                        UniMutNum=0
                        Len=len(Seq1)
                        c=0
                        while c<Len:
                             if Seq0[c]=='T' and Seq1[c]=='A': UniMutNum+=1
                             c+=1
                        Pro=1.0*UniMutNum/Len
                        if Pro<0.05: Keep='n'
                if Keep=='y':
                            if KeepLs.count(Clo0)==0:						
                              KeepLs.append(Clo0)
                              Keep2TuLs[Clo0]=[]
                            Keep2TuLs[Clo0].append(Tu)	
            #KeepLs.append(LarClo)
            if KeepLs.count(LarClo)==0:						
                              KeepLs.append(LarClo)
                              Keep2TuLs[LarClo]=[]
            Keep2TuLs[LarClo].append(Tu)	           			
        RmLs=[]
        for Clo in CloSeq:
            if KeepLs.count(Clo[1:])==0: RmLs.append(Clo[1:])
        print 'remove ancestral clones',RmLs
        print 'tumor ls for each clone', Keep2TuLs		
        Pruned = self.PruneTree(OriginalNwk, KeepLs)
        Pruned_Root =  self.RootTree(Pruned)		
        return Pruned_Root, Keep2TuLs		
                    				 
    def MakeCNVseq(self, Clo2TuLs, PreAbsCNV):
        Len=len(PreAbsCNV['CNVID'])	
        Seq=['#MEGA','!Title SNVs;','!Format datatype=dna;\n']
        for Clo in Clo2TuLs:
             TuLs=Clo2TuLs[Clo]
             CloSeq=''			 
             c=0
             while c<Len:
                Mut='n'			 
                for Tu in TuLs:			 
                     if PreAbsCNV[Tu.split('-')[-1]][c].strip()!='0': Mut='y'
                c+=1	
                if Mut=='y': CloSeq+='T'
                else: CloSeq+='A'
             Seq+=['#'+Clo,CloSeq]				
        print 'CNV seq',Seq		
        return Seq	
    def GetNode2CNVevent(self, Mutorder,Dec2Anc,Index2Seq,Index2Label,Label2Index,C1,C2,C0):
      #  CNVorder=open(CNVorder,'r').readlines()[1:]
       # CNVchr2StartEnd={}	
        Node2Event={}
        for Ind in Index2Seq:
            Node2Event[Ind]=[]		
        Event2Node={}	
        ChrPosi2MutID={}
        Mut2EventChrPosi={}		
        Len=len(Mutorder[C1])
        c=0
        while c<Len:
            Event2Node['M'+str(c)]=[]		
           # i=CNVorder[c].strip().split('\t')
            Chr=Mutorder[C1][c]
            StartEnd=Mutorder[C2][c].split(':')
            Start=int(StartEnd[0])
            End=StartEnd[-1]
            Event=Mutorder[C0][c].replace('1','').replace('2','').replace('3','').replace('4','').replace('5','').replace('6','').replace('7','').replace('8','').replace('9','').replace('0','')
            ChrPosi='Chr'+Chr+':'+str(Start)+'-'+End
            Mut2EventChrPosi['M'+str(c)]=Event+'\t'+ChrPosi			
            if ChrPosi2MutID.has_key(ChrPosi)!=True: ChrPosi2MutID[ChrPosi]=[]
            ChrPosi2MutID[ChrPosi].append('M'+str(c))			
            TimeLs=[]
            for Dec in Dec2Anc:
                DecNuc=Index2Seq[Dec][c]
                AncNuc=Index2Seq[Dec2Anc[Dec]][c]
             #   print Chr,StartEnd,Dec,Dec2Anc[Dec],DecNuc,AncNuc				
                if DecNuc=='T' and AncNuc=='A': TimeLs.append(Dec)
          #  print Chr,StartEnd,TimeLs				
            TimeLs=list(set(TimeLs))
            if TimeLs==[]:
                Normal=Label2Index['hg19']
                SuperRoot=Dec2Anc[Normal]
                Decs=Anc2Dec[SuperRoot]
                Root=''
                for Dec in Decs:
                     if Dec!=Normal: Root=Dec
                TimeLs=[Root]					 
            for Time in TimeLs:          			
                if Node2Event.has_key(Time)!=True: Node2Event[Time]=[]
                Node2Event[Time].append(Event+'\t'+'Chr'+Chr+':'+str(Start)+'-'+End)	
                Event2Node['M'+str(c)].append(Time)			
          #  CNVchr2StartEnd['Chr'+Chr]=[Start,End,TimeLs,Event]
            c+=1	 
      #  print Node2Event,'\n',Event2Node
        return Node2Event,Event2Node,ChrPosi2MutID,Mut2EventChrPosi	

    def RenameNodeID(self, GoodClo2Index,GoodDec2Anc,BadClo2Index,BadDec2Anc,BadSeq):
        NewSeq={}
        for Clo in GoodClo2Index:
             GoodIndex=GoodClo2Index[Clo]
             BadIndex=BadClo2Index[Clo]
             NewSeq[GoodIndex]=BadSeq[BadIndex]
             while GoodDec2Anc.has_key(GoodIndex)==True:
                   GoodIndex=GoodDec2Anc[GoodIndex]
                   BadIndex=BadDec2Anc[BadIndex]
                   NewSeq[GoodIndex]=BadSeq[BadIndex]
        return NewSeq	

    def FindCNVoverlappingSNV(self, Node2CNVEvent,ChrPosi2SNV):
      CNV2SNV={}
      ChrPosi2CNV=[]	
      for node in Node2CNVEvent:	
        ChrPosi2CNV+=Node2CNVEvent[node]
      ChrPosi2CNV=list(set(ChrPosi2CNV))
      for i in ChrPosi2CNV:	
        CNV2SNV[i]=[]		
      for ChrPosi in ChrPosi2SNV:
         SNV=ChrPosi2SNV[ChrPosi][0]
         ChrPosi=ChrPosi.split('\t')[-1].split(':')
         Chr=ChrPosi[0]
         Posi=int(ChrPosi[1].split('-')[0])
         for CNVChrPosi0 in ChrPosi2CNV:
              CNVChroPosi=CNVChrPosi0.split('\t')
              CNV=CNVChroPosi[0]
              ChrRegion=CNVChroPosi[-1].split(':')
              if ChrRegion[0]==Chr: 
                     Start=int(ChrRegion[1].split('-')[0])
                     End=int(ChrRegion[1].split('-')[1])
                     if Start<=Posi and Posi<=End: CNV2SNV[CNVChrPosi0].append(SNV)					 
    		 
	
      return CNV2SNV		
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