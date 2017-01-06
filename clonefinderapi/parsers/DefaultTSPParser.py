from parsers.AbstractTSPParser import AbstractTSPParser
from tsp_profiles.TumorSampleProfile import TumorSampleProfile
from tsp_profiles.ReadCount import ReadCount
import os.path

class DefaultTSPParser(AbstractTSPParser):
           
    def _parse_header(self, header):
        self._init_profile_list(header)
        header = header.strip().split('\t')
        Len=len(header)
        c=0
        Name2Col={}
        NameOrder=[]	
        while c<Len:
            Name2Col[header[c]]=c
            NameOrder.append(header[c])		
            c+=1
        return NameOrder,Name2Col         
   
    def _init_profile_list(self, header):
        header = header.strip().split('\t')
        for segment in header:
            sample_name = segment[:-4]            
            if self.tumor_sample_profiles.profile_exists(sample_name) == False:
                profile = TumorSampleProfile(sample_name)
                self.tumor_sample_profiles.add(profile)
        
    def _build_tumor_sample_profile_list(self, profile_data):
        for profile in self.tumor_sample_profiles:
            ref_col_name = profile.name + ':ref'
            alt_col_name = profile.name + ':alt'
            ref_data = profile_data[ref_col_name]
            alt_data = profile_data[alt_col_name]
            index = 0
            while index < len(ref_data):
                read_count = ReadCount()
                read_count.num_ref = int(ref_data[index])
                read_count.num_alt = int(alt_data[index])                
                profile.add(read_count)
                index += 1

    def parse(self):                
        try:
            if os.path.isfile(self._input_data_file) == False:
                IOError('input data file not found')    
            self.tumor_sample_profiles.name = os.path.basename(self._input_data_file)
            data = open(self._input_data_file, 'r').readlines()
            NameOrder,Name2Col=self._parse_header(data[0])
            data=data[1:]
            Tu2Freq={}
            for Tu in NameOrder:
              Tu2Freq[Tu]=[]
            for i in data:
              i=i.strip().split('\t')
              for Tu in Name2Col:
                  Tu2Freq[Tu].append(i[Name2Col[Tu]])
            self._build_tumor_sample_profile_list(Tu2Freq)            
            return True
        except Exception as e:
            self._messages.append(str(e))
            return False
        
    def get_tumor_sample_profile_list(self):
        return self.tumor_sample_profiles
    