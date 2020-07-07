class TumorSampleProfileList(object):
    
    def __init__(self):
        self.tumor_profiles = []
        self._name = ''
        self._num_read_counts = 0
        
    def __iter__(self):
        return iter(self.tumor_profiles)
        
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = value
    
    def count(self):
        return len(self.tumor_profiles)
    
    def num_read_counts(self):
        if len(self.tumor_profiles) > 0:
            return self.tumor_profiles[0].count
        return 0
    
    def to_string(self):
        result = self._text_header() + "\n"
        num_loci = self.tumor_profiles[0].count
        index = 0        
        while index < num_loci:            
            line = str(self.tumor_profiles[0].count_id(index)) + "\t"            
            for profile in self.tumor_profiles:
                line = line + str(profile.ref_count(index)) + "\t" + str(profile.alt_count(index)) + "\t"
            index += 1            
            result = result + line + "\n"
        return result
            
    def _text_header(self):
        result = 'ID'
        if self.count > 0:
            for profile in self.tumor_profiles:
                name = profile.name
                temp = "\t" + name + ':ref' + "\t" + name + ':alt'
                result += temp
        return result.strip()
    
    def add(self, profile):
        self.tumor_profiles.append(profile)
        
    
    def get_profile(self, index):
        if (index >= 0) and (index < len(self.tumor_profiles)):
            return self.tumor_profiles[index]
        IndexError('index out of bounds')
        
    def profile_exists(self, profile_name):
        for profile in self.tumor_profiles:
            if profile.name == profile_name:
                return True
        return False
            