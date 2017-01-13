from alignments.MegaAlignment import MegaAlignment
from Bio import Phylo
class TreeAnalizer:
    def compare_good_posi_number(self, Initial, After):
        IniCou = self.count_good_posi(Initial)
        AftCou = self.count_good_posi(After)
        if IniCou>AftCou: AfterGood='n'
        else: AfterGood='y'
       # print IniCou,AftCou,AfterGood
        return AfterGood		
    def count_good_posi(self, Posi_info):	
        Good=0
        for Posi in Posi_info:
             if Posi_info[Posi]==['Good']: Good+=1
        return Good	

    def RmUnresolvedPosi(self, SeqLs, NodeMapLs, Nwk):
         Align = MegaAlignment()	
         CloOrder0, Clo2Seq = Align.name2seq(SeqLs)
		 
      #  Meg=sys.argv[1]
      #  Nwk=Meg[:-4]+'.nwk'
      #  Nwk=open(Nwk,'r').readlines()
      #  TreNum=len(Nwk)
      #  CloOrder0, Clo2Seq, out2=Functions.ReadMegSeq(Meg)
         CloOrder=[]
         MaxClo=0		 
         for Clo in CloOrder0:
           if Clo!='#hg19' and Clo.find('#Node')==-1: CloOrder.append(Clo)
          # if Clo!='#hg19' and Clo.find('#Node')==-1: MaxClo+=1		   
         print 'cloorder',CloOrder ,MaxClo          
     #   c=1	
      #  ID0=0
      #  ID2Good={}
    #    while ID0<TreNum:
         Tree=Nwk.strip()
         Trunk=self.ReadNwkExtBra(Tree)
       #  ID=str(ID0)
        # ID2Pattern={}
         Pattern=[]
       #  NodeMap=Meg[:-4]+'_NodeMap'+ID+'.txt'
         D2A=NodeMapLs[0]
         A2D=NodeMapLs[1]
         N2C=NodeMapLs[2]
         C2N=NodeMapLs[3]#=Functions.ReadNodeMapMPtext(NodeMap)
       #  print N2C		 
         for MutAnc in N2C:
            Label=N2C[MutAnc]
         #   print 'N2C',N2C			
            DecLs=self.GetDec(Label,A2D,C2N) 
         #   print DecLs			
          #  In=''
            Pat=[]	
            for Clo in CloOrder:
               # ID=N2C[Clo[1:]]	
               # Code=ID in DecLs
                if DecLs.count(Clo[1:])!=0: #Code==True:
                   #  In+='\t'+'T'
                     Pat.append(Clo)			 
             #   else: In+='\t'+'A'
           # ID2Pattern[str(c)]=In
            Pattern.append(Pat)	
           # c+=1
         print Pattern		   
         GoodC=0
         BacC=0
         MulC=0
         #Out=Meg[:-4]+'_resolve'+str(ID0)+'.meg'	
         ChangePos=[]
         Len=len(Clo2Seq[CloOrder0[0]])
         c=0
         out={}#'Pos\tMult mut num\tBack mut num\n'
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
            print Pat			
            for Tru in Pattern:
                Miss=[] #back
                for Tclo in Tru:
                    if Pat.count(Tclo)==0: Miss.append(Tclo)
                Extra=[] #parallel
                for Pclo in Pat:
                    if Tru.count(Pclo)==0: Extra.append(Pclo)
                if len(Miss)<len(Extra):
                    if len(Miss)<SmallestDiff:
                         SmallestDiff=len(Miss)
                         DifNum2ChangeClo2Mut[len(Miss)]=Miss
                         Change='ToMut'						 
                else:
                    if len(Extra)<SmallestDiff:
                         SmallestDiff=len(Extra)
                         DifNum2ChangeClo2Wild[len(Extra)]=Extra
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
#                if len(Tru)==Num:
 #                   P=0
  #                  Dif='n'			
   #                 while P<Num:
    #                    if Tru[P]!=Pat[P]: Dif='y'
     #                   P+=1				
      #              if Dif=='n': Change='n'          			
       #     if Change=='y':
        #          SmallMult=9999
         #         SmallBack=9999		  
          #        for Tru in Pattern:
           #           TruNum=len(Tru)		  
            #          if TruNum<Num:               
             #             Miss='n'
              #            P=0
               #           while P<TruNum:
                #              Code=Tru[P] in Pat
                 #             if Code!=True: Miss='y'
                  #            P+=1
                   #       if Miss=='n': 
                    #          Add=Num-TruNum
                     #         if Add<SmallMult: SmallMult=Add                     					  
#                      elif TruNum>Num:                 
 #                         Miss='n'
  #                        P=0
   #                       while P<Num:
    #                          Code=Pat[P] in Tru
     #                         if Code!=True: Miss='y'
      #                        P+=1
       #                   if Miss=='n': 
        #                      Back=TruNum-Num
         #                     if Back<SmallBack: SmallBack=Back	
          #        if SmallMult>2 and SmallBack>2: ChangePos.append(c) 					         
           #       if SmallMult>SmallBack: BacC+=1
            #      else: MulC+=1
             #   #  if MaxClo==SmallMult+1:out[c]=['Good']				  
              #    else:out[c]=[SmallMult, SmallBack]#+=str(c)+'\t'+str(SmallMult)+'\t'+str(SmallBack)+'\n'			   
       #     else: 
        #      out[c]=['Good']#+=str(c)+'\tGood\n' 
         #     GoodC+=1	
            c+=1
         out2=[]        	
         for Clo in Clo2Seq:
            Seq=Clo2Seq[Clo]
            NewSeq=Align.ModSeq(Seq,ChangePos,'A',Len)   
            out2+=[Clo, NewSeq]#+'\n'
       #  Functions.GetOut(Out,out2)	
       #  Functions.GetOut(Out[:-4]+'MultCount.txt',out)
         GoodC={'Good':GoodC,'Mul':MulC,'Bac':BacC,'TrunkMut':Trunk}
        # ID0+=1
         return GoodC, out2, out

    def find_best_set(self, ID2Good):       
        Best=''
        BestC=0
        BestTrunk=0
        BestBack=999999999
        for ID in ID2Good:
          GoodC=ID2Good[ID]['Good']
          if GoodC==BestC:
             if BestTrunk<ID2Good[ID]['TrunkMut']:
               BestBack=ID2Good[ID]['Bac']
               Best=str(ID)
               BestTrunk=ID2Good[ID]['TrunkMut']   
             elif BestTrunk==ID2Good[ID]['TrunkMut']:	 
               BackC=ID2Good[ID]['Bac']
               if BackC<BestBack:
                    BestBack=BackC
                    Best=str(ID)			
          elif GoodC>BestC:
              Best=str(ID)
              BestC=GoodC
              BestTrunk=ID2Good[ID]['TrunkMut']
              BestBack=ID2Good[ID]['Bac']
        return int(Best)			  
       # copy(Meg[:-4]+'_resolve'+Best+'.meg',Meg[:-4]+'_resolve.meg')
     #   copy(Meg[:-4]+'_resolve'+Best+'MultCount.txt',Meg[:-4]+'_resolveMultCount.txt')
      #  Functions.GetOut('BestTree.txt',Best)	  


    def ReadNwkExtBra(self, OriNwk):
        #Functions.GetOut('test.nwk',Nwk)
        #os.system('python '+'root.py test.nwk')
        Nwk=self.RootTree(OriNwk)#open('test_Rooted.nwk','r').readlines()[0]	
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
        print OriNwk	
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': 'hg19'})
        Phylo.write(trees, 'test1.nwk', "newick")	
      #  Rooted_tree = OriNwk.root_with_outgroup({'name': 'hg19'})
		
      #  print Rooted_tree
        print open('test1.nwk','r').readlines()  
        return open('test1.nwk','r').readlines()[0]	
#Input=sys.argv[1]
#Out=Input.replace('.nwk','_Rooted.nwk')
#trees = list(Phylo.parse(Input, 'newick'))
#for tree in trees:
#    tree = tree.root_with_outgroup({'name': 'hg19'})
#Phylo.write(trees, Out, "newick")
    
    def GetDec(self, Anc,A2D,C2N):
      AncLs=[Anc]
      print AncLs,A2D	  
      Dec=[]
      while AncLs!=[]:
        TMP=[]
        for Anc in AncLs:
          Code=Anc in A2D 
          if Code==True:
            if A2D[Anc]==['-','-']: Dec.append(C2N[Anc])		  
            else:  TMP+=A2D[Anc]
          else:	
              print '?'		  
              #Dec.append(C2N[Anc])
        AncLs=TMP
      print 'get pattern',Dec		
      return Dec
##  
	    