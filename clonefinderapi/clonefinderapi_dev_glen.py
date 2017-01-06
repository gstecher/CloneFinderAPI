from parsers.DefaultTSPParser import DefaultTSPParser

if __name__ == "__main__":
    filename = '/home/gstecher/Documents/NetBeansProjects/CloneFinderAPI/dev/copyNumberVariation.txt'
    parser = DefaultTSPParser()
    parser.input_data_file = filename
    if parser.parse() == False:
        print parser.messages
    else:
        print 'success'
    
