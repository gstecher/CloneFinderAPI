from parsers.AbstractTSPParser import AbstractTSPParser
from tsp_profiles.TumorSampleProfile import TumorSampleProfile
from tsp_profiles.ReadCount import ReadCount
import os.path

class AncestreeTSPParser(AbstractTSPParser):
           
    def _parse_header(self, header):
        header = header.strip().split('\t')
        num_columns = len(header)        
        if header[0].lower() != 'gene_id':
            Exception('unexpected token. Expected gene_id but found ' + header[0])
        index = 1
        while index < num_columns:
            ref_column = header[index]
            index += 1
            alt_column = header[index]
            index += 1
            if alt_column != ref_column:
                Exception('invalid column header. Columns must be given in pairs. Found: ' + alt_column + ' and ' + ref_column)
            profile = TumorSampleProfile(alt_column)
            self.tumor_sample_profiles.add(profile)                    
   
    def _parse_row(self, row):
        row = row.strip().split('\t')  
        num_columns = len(row)
        num_counts_found = 0
        gene_id = row[0]
        index = 1
        while index < num_columns:
            ref_count = int(row[index])            
            alt_count = int(row[index + 1])
            index += 2
            read_count = ReadCount()
            read_count.id = gene_id
            read_count.num_ref = ref_count
            read_count.num_alt = alt_count
            profile = self.tumor_sample_profiles.get_profile(num_counts_found)
            profile.add(read_count)            
            num_counts_found += 1
        if num_counts_found != self.tumor_sample_profiles.count():
            Exception('Expected more sample data. Please check the input file formatting.')

    def parse(self):                
        try:
            if os.path.isfile(self._input_data_file) == False:
                IOError('input data file not found')     
            self.tumor_sample_profiles.name = os.path.basename(self._input_data_file)
            data = open(self._input_data_file, 'r').readlines()
            self._parse_header(data[0])
            index = 1
            while index < len(data):
                self._parse_row(data[index])
                index += 1
            return True
        except Exception as e:
            self._messages.append(str(e))
            return False
        
    def get_tumor_sample_profile_list(self):
        return self.tumor_sample_profiles
    