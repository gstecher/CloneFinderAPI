from parsers.AbstractTSPParser import AbstractTSPParser
from tsp_profiles.TumorSampleProfile import TumorSampleProfile
from tsp_profiles.ReadCount import ReadCount
import os.path

class DefaultCNVarser():
    def _parse_header(self, header):
        print('parsing header...')
        header = header.strip().split()
        Len=len(header)
        c=0
        Name2Col={}
        NameOrder=[]	
        while c<Len:
            Name2Col[header[c]]=c
            NameOrder.append(header[c])	            
            c+=1
        print('header parsed')
        return NameOrder,Name2Col     

    def get_tumor_cnv_profile(self,cnv_data_file):
            self.cnv_data_file=cnv_data_file	
            print('parsing CNV data file...')
            data = open(self.cnv_data_file, 'r').readlines()
            NameOrder,Name2Col=self._parse_header(data[0])
            data=data[1:]
            self.Tu2CNV={}
            for Tu in NameOrder:
              self.Tu2CNV[Tu]=[]
            for i in data:
                if i.strip() == '':
                    continue
                i=i.strip().split()
                for Tu in Name2Col:
                    if 	i[Name2Col[Tu]]=='Normal': In='normal'
                    elif i[Name2Col[Tu]]=='Loss': In='loss'
                    elif i[Name2Col[Tu]]=='Gain': In='gain'
                    else: In=i[Name2Col[Tu]]					
                    self.Tu2CNV[Tu].append(In)
            return self.Tu2CNV
