from haplopair import *

class cHaploPair:
    cBlueHaplo=None
    cPinkHaplo=None
    
    def __init__(self, a_cBlueHaplo, a_cPinkHaplo):
        self.cBlueHaplo=a_cBlueHaplo
        self.cPinkHaplo=a_cPinkHaplo
    
    def IsMatchHaploPair(self, a_cBlueHaplo, a_cPinkHaplo):
        if a_cBlueHaplo.IsSame( self.cBlueHaplo.getHaplo( ) ) and a_cPinkHaplo.IsSame( self.cPinkHaplo.getHaplo() ):
            return True
        
        if a_cBlueHaplo.IsSame( self.cPinkHaplo.getHaplo( ) ) and a_cBlueHaplo.IsSame( self.cBlueHaplo.getHaplo() ):
            return True
        
        return False

    def getHaplo(self, a_nHaplo):
        if a_nHaplo==1:
            return self.cBlueHaplo
        elif a_nHaplo==2:
            return self.cPinkHaplo
        else:
            print("exception")
            return None
        
    def getLine(self):
        strBlueHaplo=self.cBlueHaplo.getLine()
        strPinkHaplo=self.cPinkHaplo.getLine()
        
        return strBlueHaplo+"\t"+strPinkHaplo 
    
class cChildHaploPairMgr:
    def __init__(self):
        1
    
    def getPossibleChildHaploPair(self, a_cIndvChild, a_bChildPhaseKnown):
        #if phase of child is not known, we need to consider all the possible cases        
        cBlueHaplo=a_cIndvChild.getHaplo(1)
        if cBlueHaplo==None:
            return False;
        
        cPinkHaplo=a_cIndvChild.getHaplo(2)
        if cPinkHaplo==None:
            return False
        
        lPossibleHaploPair=[]
        if a_bChildPhaseKnown:
            cNewBlueHaplo=copy.deepcopy( cBlueHaplo )
            cNewPinkHaplo=copy.deepcopy( cPinkHaplo )
            cNewHaploPair=cHaploPair( cNewBlueHaplo, cNewPinkHaplo)
            lPossibleHaploPair.append( cNewHaploPair ) 
        else:
            #do not use genetic variant
            lAffected=[]
            lUnAffected=[]
            
            #3 possible cases, 2 heterozygous disease allele, 1 homozygous normal allele
            for i in range(1, 3):
                for j in range(1, 3):
                    cNewBlueHaplo=cHaplo(i, cBlueHaplo.getMarker() )
                    cNewPinkHaplo=cHaplo(j, cPinkHaplo.getMarker() )
                    cNewHaploPair=cHaploPair( cNewBlueHaplo, cNewPinkHaplo)
                    if i==2 and j==2:
                        lAffected.append( cNewHaploPair )
                    else:
                        lUnAffected.append( cNewHaploPair )
            
            lPossibleHaploPair=lAffected if a_cIndvChild.getPheno()==2 else lUnAffected
                            
        return lPossibleHaploPair
    
    def getFilteredHaploStatTable(self, a_lHaploStatTable, a_cChildHaploPair):
        cBlueHaploChild=a_cChildHaploPair.getHaplo(1)
        cPinkHaploChild=a_cChildHaploPair.getHaplo(2)
        
        lFilteredHaploStatTable=[]
        for i in a_lHaploStatTable:
            if i.IsMatchHaploPair( cBlueHaploChild, cPinkHaploChild ) == False:
                continue
            lFilteredHaploStatTable.append( i )
        return lFilteredHaploStatTable
            