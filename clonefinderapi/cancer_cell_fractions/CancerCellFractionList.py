class CancerCellFractionList(object):
    
    def __init__(self):
        self._fraction_list = []
        
    def __iter__(self):
        return iter(self._fraction_list)
    
    def __str__(self):
        string_list = self.to_string_list()
        result = ''
        for item in string_list:
            result += item + "\n"
        return result
    
    def add(self, fraction):
        self._fraction_list.append(fraction)
        
    def count(self):
        return len(self._fraction_list)
    
    def get_item(self, index):
        if (index >= 0) and (index < len(self._fraction_list)):
            return self._fraction_list[index]
        IndexError('index out of bounds')
        
    def to_string_list(self):
        result = []
        result.append('ecm1' + "\t" + 'ecm2' + "\t" + 'ecpt')
        for fraction in self._fraction_list:            
            ecm1 = "{0:.4f}".format(fraction.ecm1)
            ecm2 = "{0:.4f}".format(fraction.ecm2)
            ecpt = "{0:.4f}".format(fraction.ecpt)
            result.append(ecm1 + "\t" + ecm2 + "\t" + ecpt)
        return result
            
    def save_to_file(self, filename):
        destination = open(filename, 'w')
        data = self.__str__()
        destination.write(data)
        destination.close()
        
    