import sys
import os
import Functions

Ta=sys.argv[1] 
InMeg=sys.argv[2] 
SLEP_dir=sys.argv[3] 
Com=sys.argv[4]
def SLEPcloFreq(PatID,InMeg,Com):
	SLEPID=InMeg[:-4]
	print 'Do nnls:',PatID,InMeg	
	os.system('python '+'SLEP_NonNeg_NoBoo.py ' + PatID+'.txt ' + InMeg+' '+SLEP_dir)
	Functions.DoSLEP(Com,PatID)
	os.system('python '+'SumCsv2CloneFreqTaNoBoo.py ' + PatID+'_SNVfreq.txt ' + 'TuID.csv ' + InMeg[:-4]+'_PreAbs.txt ' + PatID)
SLEPcloFreq(Ta[:-4],InMeg,Com)