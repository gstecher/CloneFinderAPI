import Functions
import sys
from scipy import stats

Freq=sys.argv[1]
Rep=sys.argv[2]
BinoNum=int(sys.argv[3])
OriTu2Clo,Name2Col,T2C2F=Functions.GetCloHitForTu(Freq,0)

TuCloHitCou={}
AllClo=[]
for Tu in OriTu2Clo:
    CloLs=OriTu2Clo[Tu]
    AllClo+=CloLs	
    for Clo in CloLs:
       TuClo=Tu+'\t'+Clo
       TuCloHitCou[TuClo]=0
AllClo=list(set(AllClo))
Head='Tumor'
for i in AllClo:
    Head+='\t'+i
Head+='\n'
c=1	
while c<=BinoNum:
    RepTu2Clo,Name2Col,RepT2C2F=Functions.GetCloHitForTu(Rep+str(c)+'_Recom.txt',0)
    for Tu in RepTu2Clo:
        CloLs=RepTu2Clo[Tu]		
        for Clo in CloLs:
            TuClo=Tu+'\t'+Clo
            Code=TuClo in TuCloHitCou
            if Code!=True: TuCloHitCou[TuClo]=0	                 
            TuCloHitCou[TuClo]+=1
    c+=1
	
outCou=Head
for Tu in OriTu2Clo:
    InCou=''
    for Clo in AllClo:
      TuClo=Tu+'\t'+Clo
      Code=TuClo in TuCloHitCou
      if Code==True:        
        InCou+='\t'+str(1.0*(TuCloHitCou[TuClo])/BinoNum)		
      else:
        InCou+='\t0'	  
    outCou+=Tu+InCou+'\n'

Functions.GetOut(Freq[:-4]+'_Cou.txt',outCou)
