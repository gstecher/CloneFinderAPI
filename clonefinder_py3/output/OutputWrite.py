from alignments.MegaAlignment import MegaAlignment
from output.CloneFrequencyAnalizer import CloneFrequencyAnalizer

class OutputWrite():
    def GetOut(self, OutFileName, In):
        OutF=open(OutFileName,'w')
        OutF.write(In)
        OutF.close()		
    def ReNameCloFreMeg(self, seqs, CloFre, Name):
        Align = MegaAlignment()
        CloFreAnalize = CloneFrequencyAnalizer()		

        NameOrder, Clo2Seq = Align.name2seq(seqs)
        if CloFre=={}:
            CloFre['T-A']={}		
            for Clo in Clo2Seq:
                CloFre['T-A'][Clo[1:]]=1			
       # print Clo2Seq,seqs		
        Len = len(Clo2Seq[NameOrder[0]])		
        out = ['#MEGA','!Title SNVs;','!Format datatype=dna;',' ','#hg19','A'*Len]			
        TuLs=[]
        for Tu in CloFre:
             TuLs.append(Tu[3:])
        TuLs.sort()
        Old2NewCloLs={}
        Old2NewCloNum={}
        CloOrder=[]
        Num=1
        for Tu in CloFre:
            Clo2Fre=CloFre[Tu]
            HitClo=[]
            for Clo in Clo2Fre:
                if Clo2Fre[Clo]>0: HitClo.append(Clo)			
            Tu=Tu[2:]			
            C=1
            CloLs,Fre2Clo=CloFreAnalize.Sort(HitClo,Clo2Fre)#from large frequency	
            for Clo in CloLs:
                    Code=Clo in Old2NewCloLs
                    if Code!=True: 
                        Old2NewCloLs[Clo]=''
                        Old2NewCloNum[Clo]='Clone'+str(Num)				
                        CloOrder.append(Clo)
                        Num+=1				
                    Old2NewCloLs[Clo]+=Tu+str(C)
                    C+=1
        if Name=='list': Old2NewClo=Old2NewCloLs
        else: Old2NewClo=Old2NewCloNum			
        NewCloOrder=[]
        NewT2C2F={}         			
        for Clo in CloOrder:
            NewCloOrder.append(Old2NewClo[Clo])
            out+=['#'+Old2NewClo[Clo], Clo2Seq['#'+Clo]]#+'\n'	
        for Tu in CloFre:
            C2F=CloFre[Tu]
            NewC2F={}

            for C in C2F:
              if C2F[C]>0:		
               NewC2F[Old2NewClo[C]]=C2F[C]
            NewT2C2F[Tu]=NewC2F
        return out, NewT2C2F, NewCloOrder			
	