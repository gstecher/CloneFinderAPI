

class FormatInput():

    def ccf2snv(self, ccf_file, read_coverage):	
        CCF=ccf_file
        ReadCov=float(read_coverage)
        Out=CCF[:-4]+'snv.txt'
        
        Tu2CCF=self.ListColStr(CCF)
        TuOrder=[]
        out=''
        for Tu in Tu2CCF:
            TuOrder.append(Tu)
            out+=Tu+':ref\t'+Tu+':alt\t'
        out=out[:-1]+'\n'
        Len=len(Tu2CCF[Tu])

        c=0
        while c<Len:
            for Tu in TuOrder:
                Mut=int(ReadCov*float(Tu2CCF[Tu][c])/2)
                Ref=int(ReadCov-Mut)
                out+=str(Ref)+'\t'+str(Mut)+'\t'
            out=out[:-1]+'\n'
            c+=1
        self.save_file(Out,out)

    def ListColStr(self, File):
          File=open(File,'r').readlines()
          NameOrder,Name2Col=self.GetHead(File[0])
          File=File[1:]
          Tu2Freq={}
          for Tu in NameOrder:
            Tu2Freq[Tu]=[]
          for i in File:
            i=i.strip().split('\t')
            for Tu in Name2Col:
                Tu2Freq[Tu].append(i[Name2Col[Tu]])
          return Tu2Freq

    def GetHead(self, Head):
        Head=Head.strip().split('\t')
        Len=len(Head)
        c=0
        Name2Col={}
        NameOrder=[]	
        while c<Len:
            Name2Col[Head[c]]=c
            NameOrder.append(Head[c])		
            c+=1
        return NameOrder,Name2Col
    def save_file(self, Out, out):
        OutF=open(Out,'w')
        OutF.write(out)
        OutF.close()

    def cnv2snv(self, Ta, CNV):
        Out=Ta[:-4]+'snv.txt'
        SOrder, Samp2CNVfreqIn	,SNVnum,AAA=self.ObsFreqTaHead(CNV)
        Tu2TotRead	,Tu2SNVfre=self.GetTotRead(Ta)

        out=open(CNV,'r').readlines()[0]
        c=0
        while c<SNVnum:
            for S in SOrder:
                OriFre=Tu2SNVfre[S][c]
                MutCopyFra=Samp2CNVfreqIn[S][c]
                AdjFre=OriFre/(MutCopyFra*2)
                NewMut=int(round(AdjFre*Tu2TotRead[S][c],0))
                NewRef=Tu2TotRead[S][c]-NewMut
                out+=str(NewRef)+'\t'+str(NewMut)+'\t'
            out=out[:-1]+'\n'
            c+=1
        self.save_file(Out,out)	
	
    def snv2snv(self, Ta):
        Out=Ta[:-4]+'snv.txt'
        Ta=open(Ta,'r').readlines()		
        SOrder, Samp2Col = self.GetHeadObsFreqTaHead(Ta[0].strip())
        out=''
        for Sample in SOrder:
            out+=Sample+':ref\t'+Sample+':alt\t'	
        out=out[:-1]+'\n'			
        Ta=Ta[1:]
        for i in Ta:		
            i=i.strip().split('\t')
            for S in SOrder:
                out+=i[Samp2Col[S+':ref']]+'\t'+i[Samp2Col[S+':alt']]+'\t'
            out=out[:-1]+'\n'
            
        self.save_file(Out,out)	

    def ObsFreqTaHead(self, Freq): 
      Freq=open(Freq,'r').readlines() 
      SampOrder,Name2Col=self.GetHeadObsFreqTaHead(Freq[0])
      Samp2FreqIn={}
      Samp2TotRead={}
      for Samp in SampOrder:
          Samp2FreqIn[Samp]=[]
          Samp2TotRead[Samp]=[]	  
      Freq=Freq[1:]
      SNVnum=0
      for i in Freq:
        i=i.strip().split('\t')
        TMP={}
        for Samp in SampOrder:
            MutC=int(i[Name2Col[Samp+':alt']])
            RefC=int(i[Name2Col[Samp+':ref']])
            MutFreq=1.0*MutC/(MutC+RefC)
            Tot=MutC+RefC
            Samp2FreqIn[Samp].append(MutFreq)
            Samp2TotRead[Samp].append(Tot)		
        SNVnum+=1
      return SampOrder, Samp2FreqIn	,SNVnum,Samp2TotRead

    def GetHeadObsFreqTaHead(self, Head): 
      Head=Head.strip().split('\t')
      SampOrder=[]
      Name2Col={}
      c=0
      Len=len(Head)
      while c<Len:
        i=Head[c]
        if i.find(':')!=-1:
            Samp=i.split(':')[0]
            Code=Samp in SampOrder
            if Code!=True:
                 SampOrder.append(Samp)
            Name2Col[i]=c
        c+=1
      return SampOrder,Name2Col	  

    def GetHeadObsFreqTaHead1(self, Head): #as a string
      Head=Head.strip().split('\t')
      SampOrder=[]
      Name2Col={}
      c=0
      Len=len(Head)
      while c<Len:
        i=Head[c]
        if i.find(':')!=-1:
            Samp=i.split(':')[0].replace(' ','')
            Code=Samp in SampOrder
            if Code!=True:               		
                 SampOrder.append(Samp)
        Name2Col[i.replace(' ','')]=c
        c+=1
      return SampOrder,Name2Col  
	  
    def GetTotRead(self, Obs):
      Obs=open(Obs,'r').readlines()
      
      TuOrder,Tu2Col =self.GetHeadObsFreqTaHead1(Obs[0])
      Tu2SNVfre={}
      Tu2TotRead={}
      Obs=Obs[1:]
      for i in Obs:
        i=i.strip().split('\t')
        for Tu in TuOrder:
            Mut=int(i[Tu2Col[Tu+':'+'alt']])
            Tot=int(i[Tu2Col[Tu+':'+'ref']])+Mut
            Fre=1.0*Mut/Tot
            Code=Tu in Tu2SNVfre
            if Code!=True: 
                 Tu2SNVfre[Tu]=[]
                 Tu2TotRead[Tu]=[]
            Tu2SNVfre[Tu].append(Fre)
            Tu2TotRead[Tu].append(Tot)
      return Tu2TotRead	,Tu2SNVfre	  