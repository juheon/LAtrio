from .main import *
from main import cMLEMgr

class cEXTMLEMgr(cMLEMgr):
    
    def __getChildLOD(self, a_cIndvFather, a_cIndvMother, a_cIndvChild, a_fTheta):
        cBlueHaplo=a_cIndvChild.getHaplo(1)
        cPinkHaplo=a_cIndvChild.getHaplo(2)
        
        #Find parental origins of each haplotype
        nOriginBlueHaplo=self.__getParentalOriginByMarker( a_cIndvFather, a_cIndvMother, cBlueHaplo)
        nOriginPinkHaplo=self.__getParentalOriginByMarker( a_cIndvFather, a_cIndvMother, cPinkHaplo)
        
        if nOriginBlueHaplo<0:
            nOriginBlueHaplo=self.__updateParentalOriginUsingPairHaplo(nOriginBlueHaplo, nOriginPinkHaplo)
        
        if nOriginPinkHaplo<0:
            nOriginPinkHaplo=self.__updateParentalOriginUsingPairHaplo(nOriginPinkHaplo, nOriginBlueHaplo)
            
            
        #Divide two cases where parental origins of haplotypes is ambiguous or not 
        if nOriginBlueHaplo<0 or nOriginPinkHaplo<0:
            print("EXCEPTION")
            return -1
        else:
            #unambiguous parental origin
            cChildHaploFromFather=cBlueHaplo if nOriginBlueHaplo==1 else cPinkHaplo
            cChildHaploFromMother=cBlueHaplo if nOriginBlueHaplo==0 else cPinkHaplo
            
            newTrioRecomb=self.__getTrioRecomb( a_cIndvFather, a_cIndvMother, cChildHaploFromFather, cChildHaploFromMother )
            fNewTrioLOD=self.__calcLODScore(newTrioRecomb.getMeiosis(), newTrioRecomb.getRecomb(), a_fTheta)
    
            return fNewTrioLOD