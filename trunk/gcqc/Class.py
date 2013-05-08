''' GCQC.gcqc.Class.py
'''
class Compound(object):
    '''
    class to encapsulate a GCMS compound for QC analysis
    '''
    def __init__(self, name, rt, ion):
        '''
        Constructor
        
        @param name:    Compound name
        @type name:    StringType
        @param rt:    Retention time for compound
        @type rt:    FloatType
        @param ion:    Diagnostic ion for compound
        @type ion:    IntType (GCMS data is unit mass resolution)
        '''
        
        self.__name = name
        self.__rt = rt
        self.__ion = ion

    def get_name(self):
        
        return self.__name

    def get_ion(self):
        
        return self.__ion

    def set_rt(self, rt):
        '''
        The retention time for the ion when found in a sample
        
        @param rt:    Retention time for compound
        @type rt:    FloatType
        '''
        
        self.__rt = rt
        
    def get_rt(self):
        
        return self.__rt

    def set_area(self, area):
        '''
        The area of the ion diagnostic for the compound
        
        @param area:    Peak area for diagnostic ion
        @type area:    FloatType
        '''
        
        self.__area = area

    def get_area(self):
        
        return self.__area

    def set_intensity(self, intensity):
        '''
        The intensity of the ion diagnostic for the compound
        
        @param intensity:    Intensity of diagnostic ion
        @type area:    FloatType
        '''
        
        self.__intensity = intensity

    def get_intensity(self):
        
        return self.__intensity
    
    def set_symmetry(self, sym_value):
        '''
        The ratio of left to right borders,
            as measured from the apex of the peak
        
        @param sym_value:    Symmetry score
        @type sym_value:    FloatType
        '''
        
        self.__symmetry = sym_value
    
    def get_symmetry(self):
        
        return self.__symmetry
    
    def set_delta(self, delta_value):
        '''
        The time difference between the found RT and the expected RT
        
        @param delta_value:    Time (s)
        @type delta_value:    FloatType
        '''
        
        self.__delta = delta_value
    
    def get_delta(self):
        
        return self.__delta

# EOF
