import sys
import Functions

Freq=sys.argv[1]
MegName=sys.argv[2]
Fol=sys.argv[3]
PatID=Freq.split('.txt')[0].split('/')[-1]

SampOrder,Samp2FreqIn,SNVnum,Tu2TotRead=Functions.ObsFreqTaHead(Freq)
CloOrder, Clo2Seq, out2=Functions.ReadMegSeq(MegName)
Clo2PreAbs=Functions.Meg2PreAbs(CloOrder,Clo2Seq,MegName[:-4]+'_PreAbs.txt')
In=''
Len=len(Clo2Seq[CloOrder[0]])
c=0
while c<Len:
  for Name in CloOrder:
       if Name=='#hg19' or Name=='#Normal':pass
       else:  In+=Clo2PreAbs[Name][c]+' '
  In=In[:-1]+'; '
  c+=1	
Meg2Mat={MegName[:-4]:In}
Functions.MakeSLEPinNoBoo(PatID,Meg2Mat,SampOrder,Samp2FreqIn,Fol)

