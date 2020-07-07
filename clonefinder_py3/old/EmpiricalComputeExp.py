import sys
import os

Ta=sys.argv[1] #_MutCount.txt
Out=Ta[:-4]+'Format.csv'
SigOrder='C:\\Users\\kumarlab\\Desktop\\CancerSoftware\\TrackSig-master\\TrackSig-master\\annotation\\cosmicSig.csv' #please change the path
SigOrder=open(SigOrder,'r').readlines()[1:]
Order=[]
for i in SigOrder:
    Order.append(i.split(',')[0])
print Order, len(Order)

Ta=open(Ta,'r').readlines()
Sig2Count={}
for i in Ta:
    i=i.strip().split('\t')
    Ref=i[3]
    Mut=i[4]
    Tri=i[5]
    Sig=Tri[0]+'['+Ref+'>'+Mut+']'+Tri[2]
    if Sig2Count.has_key(Sig)!=True: Sig2Count[Sig]=0
    Sig2Count[Sig]+=1
out=''
for Sig in Order:
    if Sig2Count.has_key(Sig)==True: Count=str(Sig2Count[Sig])
    else: Count='0'
    out+=Sig+','+Count+'\n'
OutF=open(Out,'w')
OutF.write(out)
OutF.close()	
#os.system('python EstimateExposure.py '+Out)