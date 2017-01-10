import sys
import Functions

TuLsTa=sys.argv[1]
FreqTa=sys.argv[2]
PreAbsTa=sys.argv[3]
PatID=sys.argv[4]

Head=open(TuLsTa,'r').readlines()[0]
TuLs,AA=Functions.GetHead(Head)

PreAbsTa=PreAbsTa.replace('PatID',PatID)
PreAbsTa=open(PreAbsTa,'r').readlines()
Head=PreAbsTa[0]
CloneOrder,AA=Functions.GetHead(Head)
CloNum=len(CloneOrder)
out='Tumor'
for Clo in CloneOrder: out+='\t'+Clo
out+='\n'
Out=PatID+'_CloneFreq.txt'
AllClone=''
for TuID in TuLs:
    File=FreqTa.replace('PatID',PatID)
    File=File.replace('TuID',TuID)
    File=open(File,'r').readlines()
    c=0
    CloIn=''
    while c<CloNum:	
        CloName=CloneOrder[c]
        Line=File[c].split(',')		
        CloIn+='\t'+str(float(Line[0])*2)
        c+=1
    out+=TuID+CloIn+'\n'
Functions.GetOut(Out,out)