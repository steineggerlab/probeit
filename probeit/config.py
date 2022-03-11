
class Config:
    # thermofilter
    minGC = 0.30
    maxGC = 0.70
    maxHomoDimerTm = 60
    maxHairpinTm = 60
    minProbeTm = 40
    
    # file name
    log = 'log.txt'
    window = 'window.fa'
    probe1 = 'probe1.fa'
    probe2 = 'probe2.fa'
    
    @classmethod
    def getMinGC(cls):
        return cls.minGC
    
    @classmethod
    def getMaxGC(cls):
        return cls.maxGC
    
    @classmethod
    def getMaxhomoDimerTm(cls):
        return cls.maxHomoDimerTm
    
    @classmethod
    def getMaxhairpinTm(cls):
        return cls.maxHairpinTm
    
    @classmethod
    def getMinProbeTm(cls):
        return cls.minProbeTm
        

        
    
