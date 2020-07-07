
from Bio import Phylo
import os

class TreeAnalizer:
		
    def RootTree(self, OriNwk):
 
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': 'hg19'})
        Phylo.write(trees, 'test1.nwk', "newick")	
        Rooted=open('test1.nwk','r').readlines()[0]	
        os.remove('test1.nwk')	
        os.remove('test.nwk')			
        return Rooted
