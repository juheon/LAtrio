import math
from pedigree import *
from haplopair import cHaploPairMgr
from childhaplo import cChildHaploPairMgr
    
    
class cLODTrio:
    cIndvFather=None
    cIndvMother=None
    cIndvChild=None
    
    lChildHaploPairPhaseUnknown=None
    lChildHaploPairPhaseKnown=None
    lParentHaploStatTable=None
    
    hpmgr=None
    chpmgr=None
    
    def __init__(self, a_cIndvFather, a_cIndvMother, a_cIndvChild):
        self.cIndvFather=a_cIndvFather
        self.cIndvMother=a_cIndvMother
        self.cIndvChild=a_cIndvChild
        
        self.hpmgr=cHaploPairMgr()
        self.lParentHaploStatTable=self.hpmgr.getHaploStatTable(a_cIndvFather, a_cIndvMother)
        
        
        self.chpmgr=cChildHaploPairMgr()
        self.lChildHaploPairPhaseKnown=self.chpmgr.getPossibleChildHaploPair(a_cIndvChild, True)
        self.lChildHaploPairPhaseUnknown=self.chpmgr.getPossibleChildHaploPair(a_cIndvChild, False)

    def IsAffectedChild(self):
        return True if self.cIndvChild.getPheno()==2 else False
        
    def calcLikelihood(self, a_fTheta, a_bPhaseKnown):
        lChildHaploPair=None
        if a_bPhaseKnown:
            lChildHaploPair=self.lChildHaploPairPhaseKnown
        else:
            lChildHaploPair=self.lChildHaploPairPhaseUnknown
            
        #Likelihood value may vary by phase of children
        #If the child has known phase, nChildHaploPairCase==1
        #Otherwise, this will calculate every possible phases and weighted sum
        nChildHaploPairCase=len( lChildHaploPair )
        if nChildHaploPairCase < 1:
            print("exception")
            return None
        fWeightChildHaploCase=1 if nChildHaploPairCase==1 else 1/nChildHaploPairCase
        
        
        fLikelihood=0
        for i in lChildHaploPair:
            newFLikelihood=self.__calcLikelihoodForChildHaploPair(self.lParentHaploStatTable, i, a_fTheta)
            fLikelihood+=fWeightChildHaploCase*newFLikelihood
        
        return fLikelihood

    def __calcLikelihoodForChildHaploPair(self, a_lParentHaploStatTable, a_cChildHaploPair, a_fTheta):
        #For given haplotypes of the child and possible haplotype table of parents, we can find out match cases.
        #Calculate weighted sum of each cases
        lFilteredHaploStatTable=self.chpmgr.getFilteredHaploStatTable(a_lParentHaploStatTable, a_cChildHaploPair)
        nPossibleParentCase=len(lFilteredHaploStatTable)
        if nPossibleParentCase < 1:
            #print("exception")
            return 0
        elif nPossibleParentCase >=2:
            print("exception")
            return 0
        
        cMatchHaploStat=lFilteredHaploStatTable[0]
        
        nRecomb=cMatchHaploStat.getRecomb()
        nNonRecomb=cMatchHaploStat.getNonRecomb()
        fLikelihood=math.pow(a_fTheta, nRecomb)*math.pow(1-a_fTheta, nNonRecomb)
        return fLikelihood
    
"""        for i in lFilteredHaploStatTable:
            nMeiosis=i.getTotalMeiosis()
            nRecomb=i.getRecomb()
            nNonRecomb=i.getNonRecomb()
            newLikelihood=math.pow(a_fTheta, nRecomb)*math.pow(1-a_fTheta, nNonRecomb)            
        return fLikelihood
"""

class cLODPed:
    lLODTrio=None
    cPed=None
    def __init__(self, a_cPed):
        self.cPed=a_cPed
        self.lLODTrio=[]
        for i in a_cPed.getChild():
            newPedTrio=cLODTrio( a_cPed.getParent(1), a_cPed.getParent(2), i )
            self.lLODTrio.append( newPedTrio )
            
    def __getInformativeMeiosis(self):
        cFather=self.cPed.getParent(1)
        cMother=self.cPed.getParent(2)
        
        cBlueHFather=cFather.getHaplo(1)
        cPinkHFather=cFather.getHaplo(2)
        
        cBlueHMother=cMother.getHaplo(1)
        cPinkHMother=cMother.getHaplo(2)
        
        bFatherInfo=False if cBlueHFather.getAllele() == cPinkHFather.getAllele() or cBlueHFather.getMarker() == cBlueHMother.getMarker() else True
        bMotherInfo=False if cBlueHMother.getAllele() == cPinkHMother.getAllele() or cBlueHMother.getMarker() == cPinkHMother.getMarker() else True    

        if bFatherInfo==True and bMotherInfo==True:
            return 2
        elif bFatherInfo==True and bMotherInfo==False:
            return 1
        elif bFatherInfo==False and bMotherInfo==True:
            return 1
        else:
            return 0

    def calcLikelihood(self, a_fTheta, a_bPhaseKnown, a_bAffectedOnly):
        fLikelihood=0
        for i in self.lLODTrio:
            if a_bAffectedOnly==True:
                if i.IsAffectedChild()==False:
                    continue
            
            newfLikelihood=i.calcLikelihood( a_fTheta, a_bPhaseKnown )
            #print(newfLikelihood)
            fLikelihood=newfLikelihood if fLikelihood==0 else fLikelihood*newfLikelihood
        
        return fLikelihood
        
    def calcLOD(self, a_fTheta, a_bPhaseKnown, a_bAffectedOnly):
        nInformativeMeiosis=self.__getInformativeMeiosis()
        nMeiosis=nInformativeMeiosis*self.cPed.getChildCount(a_bAffectedOnly)
        fLikelihood=self.calcLikelihood(a_fTheta, a_bPhaseKnown, a_bAffectedOnly)
        
        fLOD=self.__calcLODscore(nMeiosis, fLikelihood)        
        return fLOD
    
    def getPedID(self):
        return self.cPed.getPedID()


    def __calcLODscore(self, a_nMeiosis, a_fLikelihood):
        fOdds=a_fLikelihood / math.pow( 0.5, a_nMeiosis )
        fLOD=math.log10( fOdds )
        return fLOD



class cLODMgr:
    lLODPed=None
    bPhaseKnown=False
    bAffectedOnly=False
    nPedCount=0
    
    def __init__(self, a_lPed):
        self.nPedCount=len(a_lPed)
        self.lLODPed=[]
        for i in a_lPed:
            newLODPed=cLODPed( i )
            self.lLODPed.append( newLODPed )
    
    def getListTheta(self, a_fMin, a_fMax, a_fBinSize):
        fMax=min(0.5, a_fMax)
        
        lTheta=[]
        nBinCount=int((fMax-a_fMin)/a_fBinSize)
        for i in range(0, nBinCount+1):
            fNewTheta=a_fMin+i*a_fBinSize
            lTheta.append(fNewTheta)
        
        #theta value cannot be 0
        if lTheta[0]==0:
            del lTheta[0]
        return lTheta
    
    def calcListLOD(self, a_lfTheta, a_bPhaseKnown, a_bAffectedOnly):
        self.bPhaseKnown=a_bPhaseKnown
        self.bAffectedOnly=a_bAffectedOnly
        
        lLODStat=[]
        for i in a_lfTheta:
            fLOD=self.calcLOD(i, a_bPhaseKnown, a_bAffectedOnly)
            newLODStat=cLODStat(i, fLOD)
            lLODStat.append(newLODStat)
        return lLODStat
    
    def calcLOD(self, a_fTheta, a_bPhaseKnown, a_bAffectedOnly):
        fLOD=0
        for i in self.lLODPed:
            fPedLOD=i.calcLOD( a_fTheta, a_bPhaseKnown, a_bAffectedOnly)
            fLOD+=fPedLOD
            
            #if fPedLOD<0:
            #print(str(i.getPedID()+"\t"+str(fPedLOD) ))
        return fLOD
    
    def print(self, a_lLODStat):
        cMaxLODStat=None
        for i in a_lLODStat:
            if cMaxLODStat==None:
                cMaxLODStat=i
            else:
                if i.getLOD() >= cMaxLODStat.getLOD():
                    cMaxLODStat=i

        strLinkageAnalysisType="Disease variant" if self.bPhaseKnown==True else "Phenoype"
        strChildType="Only affected children" if self.bAffectedOnly==True else "Affected+Unaffected Children"
                
        print("Note====================")
        print("Pedigree Count: "+str(self.nPedCount))
        print("Linkage Analysis: Marker+" +strLinkageAnalysisType)
        print("Children inclusion: " +strChildType)
        print("1st column: theta")
        print("2nd column: LOD")
       
        print("Estimated Theta=========")
        print(cMaxLODStat.getLine())
        
        print("Grid search result======")
        for i in a_lLODStat:
            print(i.getLine())


class cLODStat:
    fTheta=-1
    fLOD=-1

    def __init__(self, a_fTheta, a_fLOD):
        self.fTheta=a_fTheta
        self.fLOD=a_fLOD

    def getTheta(self):
        return self.fTheta
   
    def getLOD(self):
        return self.fLOD
    
    def getLine(self):
        return str(round(self.fTheta, 3))+"\t"+str(round(self.fLOD, 3))
        