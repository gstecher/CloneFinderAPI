from alignments.MegaAlignment import MegaAlignment

class CloneFrequencyWriter():
    """
        Compute clone frequency (F) from clone sequences matrix (M) and observed variant frequencies (v_obs).
        
        M x F = v_obs	
    """
    def __init__(self, clone_frequency_table, original_clone_sequences, cutoff_clone_frequency):
        self.CloFreTable = clone_frequency_table
        Align = MegaAlignment()
        SeqOrder, self.OriSeqDic = Align.name2seq(original_clone_sequences)		
        self.CutOff = cutoff_clone_frequency

    def get_hitclone(self):
        self.hitseq = []
        self.hitclone_frequency = {}
	
        for tumor in self.CloFreTable:
            CloFre = self.CloFreTable[tumor]
            for clone in CloFre:
                if CloFre[clone] > self.CutOff:
                   self.hitseq.append(clone)
        self.hitseq = list(set(self.hitseq))
     
        for tumor in self.CloFreTable:
            self.hitclone_frequency[tumor] = {} 		
            CloFre = self.CloFreTable[tumor]
            for clone in CloFre:
                if self.hitseq.count(clone) != 0:
                    self.hitclone_frequency[tumor][clone] = CloFre[clone]
        self.hitseq_align = ['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
        for seq in self.hitseq:
            self.hitseq_align += ['#'+seq, self.OriSeqDic['#'+seq]]
			
        return 	self.hitseq_align, self.hitclone_frequency		
        					
                            		