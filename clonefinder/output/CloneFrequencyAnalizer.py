

class CloneFrequencyAnalizer:

    def Sort(self, Ls, C2F):
            
            FreLs=[]
            Fre2Clo={}
            for i in Ls:
               if i[0]=='#': i=i[1:]
               Code=i in C2F
               if Code==True:	   
                F=C2F[i]
                Code=F in FreLs
                if Code!=True:
                   FreLs.append(F)
                   Fre2Clo[F]=[]
                Fre2Clo[F].append(i)
            FreLs.sort(reverse=True) 
            TMP=[]
            for F in FreLs:
                TMP+=Fre2Clo[F]
            return TMP,Fre2Clo
			

    def save_frequency_table_to_file(self, file_name, T2C2F, NameLs):
            if NameLs == []:
                for T in T2C2F:
                    C2F=T2C2F[T]
                    for C in C2F:					
                        if C2F[C]>0: NameLs.append(C)
                NameLs=list(set(NameLs))
            TMP='Tumor'	
            for Name in NameLs: #Clone
              	TMP+='\t'+Name				
            TMP+='\n'
            TuLs=[]
            for T in T2C2F:
                TuLs.append(T)
            TuLs.sort()		
            for T in TuLs:
                TMP+=T
                for C in NameLs:
                    Code=C in T2C2F[T]
                    if Code!=True: In='\t0'			
                    else:In='\t'+str(T2C2F[T][C])		
                    TMP+=In
                TMP+='\n'
     
            destination = open(file_name,'w')
            destination.write(TMP)
            destination.close()    