from pedigree import *

class cHaploStatPair:
    cBlueHaploStat=None
    cPinkHaploStat=None
    nTotalMeiosis=-1
    nRecomb=-1
    nNonRecomb=-1
    
    def __init__(self, a_cBlueHaploStat, a_cPinkHaploStat):
        self.cBlueHaploStat=a_cBlueHaploStat
        self.cPinkHaploStat=a_cPinkHaploStat
        self.__calcStat()
    
    def __calcStat(self):
        l_nTotal=0
        l_nRecomb=0
        if self.cBlueHaploStat.getRecombStatus() != -1:
            l_nTotal+=1
            if self.cBlueHaploStat.getRecombStatus() == 1:
                l_nRecomb+=1
        
        if self.cPinkHaploStat.getRecombStatus() != -1:
            l_nTotal+=1
            if self.cPinkHaploStat.getRecombStatus() == 1:
                l_nRecomb+=1
        
        self.nTotalMeiosis=l_nTotal
        self.nRecomb=l_nRecomb
        self.nNonRecomb=l_nTotal-l_nRecomb
    
    def IsMatchHaploPair(self, a_cBlueHaplo, a_cPinkHaplo):
        if a_cBlueHaplo.IsSame( self.cBlueHaploStat.getHaplo( ) ) and a_cPinkHaplo.IsSame( self.cPinkHaploStat.getHaplo() ):
            return True
        
        if a_cBlueHaplo.IsSame( self.cPinkHaploStat.getHaplo( ) ) and a_cBlueHaplo.IsSame( self.cBlueHaploStat.getHaplo() ):
            return True
        
        return False
    
    def getTotalMeiosis(self):
        return self.nTotalMeiosis
    
    def getRecomb(self):
        return self.nRecomb
    
    def getNonRecomb(self):
        return self.nNonRecomb
        
    def getLine(self):
        strBlueHaploStat=self.cBlueHaploStat.getLine()
        strPinkHaploStat=self.cPinkHaploStat.getLine()
        
        return strBlueHaploStat+"\t"+strPinkHaploStat+"\t"+str(self.getTotalMeiosis() )+"\t"+str(self.getRecomb() )+"\t"+str(self.getNonRecomb())
    
class cHaploStat:
    cHaplo=None
    nRecombStatus=-1 #1: recomb, 0: nonrecomb, -1: not informative
    def __init__(self, a_cHaplo, a_nRecombStatus):
        self.cHaplo=a_cHaplo
        self.nRecombStatus=a_nRecombStatus
        
    def getHaplo(self):
        return self.cHaplo
    
    def getRecombStatus(self):
        return self.nRecombStatus
    
    def getLine(self):
        if self.cHaplo==None:
            return "-_-"+"\t"+"0"
        else:
            return self.cHaplo.getLine()+"\t"+str( self.nRecombStatus )
        
class cHaploPairMgr:
    def __init__( self ):
        1
        
    def __IsInformativeMeiosis(self, a_cIndvOneParent):
        #deprecated
        if a_cIndvOneParent==None:
            return False
        
        cBlueHaplo=a_cIndvOneParent.getHaplo(1)
        if cBlueHaplo==None:
            return False;
        
        cPinkHaplo=a_cIndvOneParent.getHaplo(2)
        if cPinkHaplo==None:
            return False
        
        if cBlueHaplo.getMarker()==cPinkHaplo.getMarker():
            return False
        else:
            return True

    def __getPossibleHaploStat(self, a_cIndvOneParent):
        cBlueHaplo=a_cIndvOneParent.getHaplo(1)
        if cBlueHaplo==None:
            return False;
        
        cPinkHaplo=a_cIndvOneParent.getHaplo(2)
        if cPinkHaplo==None:
            return False
        
        nDistinctAllele=1 if cBlueHaplo.getAllele() == cPinkHaplo.getAllele() else 2
        nDistinctMarker=1 if cBlueHaplo.getMarker() == cPinkHaplo.getMarker() else 2
        
        lPossibleHaploStat=[]
        for i in range(1, (nDistinctAllele+1) ):
            for j in range(1, (nDistinctMarker+1) ):
                cHaplo1=a_cIndvOneParent.getHaplo(i)
                cHaplo2=a_cIndvOneParent.getHaplo(j)
                nRecomb= 1 if i!=j else 0
                
                newcHaplo=cHaplo( cHaplo1.getAllele(), cHaplo2.getMarker() )
                cNewHaploStat=cHaploStat( newcHaplo, nRecomb )
                lPossibleHaploStat.append( cNewHaploStat ) 
        return lPossibleHaploStat
            
    def getHaploStatTable(self, a_cIndvFather, a_cIndvMother):
        lHaploStatFather=self.__getPossibleHaploStat( a_cIndvFather )
        lHaploStatMother=self.__getPossibleHaploStat( a_cIndvMother )
        lHaploStatPairTable=[]
        for i in lHaploStatFather:
            for j in lHaploStatMother:
                newHaploStatPair=cHaploStatPair(i, j)
                lHaploStatPairTable.append( newHaploStatPair )
        
        return lHaploStatPairTable
    
    def print(self, a_lHaploStatPairTable):
        for i in a_lHaploStatPairTable:
            print( i.getLine() );
    