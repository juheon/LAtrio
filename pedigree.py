import copy
class cHaplo:
    nAllele=-1
    nMarker=-1
    
    def __init__(self, a_nAllele, a_nMarker):
        self.nAllele=a_nAllele;
        self.nMarker=a_nMarker;      
    
    def getAllele(self):
        return self.nAllele;
    
    def setAllele(self, a_nAllele):
        self.nAllele=a_nAllele
    
    def getMarker(self):
        return self.nMarker;
    
    def getLine(self):
        return str(self.nAllele)+"_"+str(self.nMarker)

    def IsSame(self, a_cHaplo):
        if a_cHaplo==None:
            return False
        
        if a_cHaplo.getMarker()==self.getMarker() and a_cHaplo.getAllele()==self.getAllele():
            return True
        else:
            return False
    
class cIndv:
    nPheno=-1
    cBlueHaplo=None
    cPinkHaplo=None
    
    def __init__(self, a_nPheno, a_nBAllele, a_nPAllele, a_nBMarker, a_nPMarker):
        self.nPheno=a_nPheno
        self.cBlueHaplo=cHaplo(a_nBAllele, a_nBMarker);
        self.cPinkHaplo=cHaplo(a_nPAllele, a_nPMarker);
        self.cBlueHaploRecomb=cHaplo(a_nBAllele, a_nPMarker);
        self.cPinkHaploRecomb=cHaplo(a_nPAllele, a_nBMarker);
    
    def getPheno(self):
        return self.nPheno

    def getHaplo(self, a_nHaplo):
        if a_nHaplo==1:
            return self.cBlueHaplo
        elif a_nHaplo==2:
            return self.cPinkHaplo
        else:
            return None
        
    #def getHaploRecomb(self, a_nHaplo):
    #    if a_nHaplo==1:
    #        return self.cBlueHaploRecomb
    #    elif a_nHaplo==2:
    #        return self.cPinkHaploRecomb
    #    else:
    #        return None

    
    def getLine(self):
        strLine=str(self.nPheno)+"\t"+self.cBlueHaplo.getLine()+"\t"+self.cPinkHaplo.getLine()
        return strLine
    
    def recombAllele(self):
        nBlueAllele=self.cBlueHaplo.getAllele()
        nPinkAllele=self.cPinkHaplo.getAllele()
        
        self.cBlueHaplo.setAllele(nPinkAllele)
        self.cPinkHaplo.setAllele(nBlueAllele)
              
class cPed:
    nPedId=-1
    cIFather=None
    cIMother=None
    
    lChild=[]
    
    def __init__(self, a_nPedId, a_nFPheno, a_nFFAllele, a_nFMAllele, a_nFFMarker, a_nFMMarker, a_nMPheno, a_nMFAllele, a_nMMAllele, a_nMFMarker, a_nMMMarker):
        self.nPedId=a_nPedId
        self.cIFather=cIndv(a_nFPheno, a_nFFAllele, a_nFMAllele, a_nFFMarker, a_nFMMarker)
        self.cIMother=cIndv(a_nMPheno, a_nMFAllele, a_nMMAllele, a_nMFMarker, a_nMMMarker)
        self.clearChild()
    
    def clearChild(self):
        if self.lChild != None:
            self.lChild.clear()
        self.lChild=[]
        
    def addChild(self, a_nPheno, a_nFAllele, a_nMAllele, a_nFMarker, a_nMMarker):
        self.lChild.append( cIndv(a_nPheno, a_nFAllele, a_nMAllele, a_nFMarker, a_nMMarker) )
        
    def setChild(self, a_lChild):
        self.lChild=a_lChild
    
    def getParent(self, a_nFather): #father=1, mother=2
        if a_nFather==1:
            return self.cIFather
        elif a_nFather==2:
            return self.cIMother
        else:
            return None

    def getChild(self):
        if len(self.lChild) < 1:
            return None
        else:
            return self.lChild;
   
    def getChildByIndex(self, a_nIndex):
        if a_nIndex >= len(self.lChild):
            return None
        else:
            return self.lChild[a_nIndex];

    def getChildCount(self, a_bAffectedOnly):
        if self.lChild==None:
            return -1
        else:
            if a_bAffectedOnly==False:
                return len(self.lChild);
            
            nCount=0
            for i in self.lChild:
                if i.getPheno()==2:
                    nCount+=1
                    
            return nCount
    
    def getPedID(self):
        return self.nPedId;
    
    def getLine(self):
        strLine=str(self.nPedId)+"\n"+self.cIFather.getLine()+"\n"+self.cIMother.getLine()
        for i in self.lChild:
            strLine+="\n"+i.getLine()
        return strLine 
        
class cPedMgr:
    nChild=-1
    lPed=[]
    fPenentrance=1  ##extend this to use any penentrance number
    
    def __init__(self, a_strFPath, a_nChild):
        self.lPed=[]
        self.nChild=a_nChild
        self.loadFam(a_strFPath)
        
    def loadFam(self, a_strFPath):
        dicIndex={}
        f=open(a_strFPath)
        for line in f:
            arr=line.split(',')
            if arr[0]=="seed":
                j=0
                for i in arr:
                    dicIndex[i]=j
                    j+=1
                continue

            nPedId=arr[ dicIndex["seed"] ]
            #create pedigree
            nFPheno=self.__genotype2pheno( int( arr[ dicIndex["dgf"] ] ) )
            nFFAllele=int( arr[ dicIndex["blue_daf"] ])
            nFMAllele=int( arr[ dicIndex["pink_daf"] ])
            nFFMarker=int( arr[ dicIndex["blue_laf"] ])
            nFMMarker=int( arr[ dicIndex["pink_laf"] ])

            nMPheno=self.__genotype2pheno( int( arr[ dicIndex["dgm"] ] ) )
            nMFAllele=int( arr[ dicIndex["blue_dam"] ])
            nMMAllele=int( arr[ dicIndex["pink_dam"] ])
            nMFMarker=int( arr[ dicIndex["blue_lam"] ])
            nMMMarker=int( arr[ dicIndex["pink_lam"] ])


            #print(str(arr[ (nColIdx_genotype+1) ])+","+str(arr[ (nColIdx_genotype+2) ])+" "+str(nFPheno)+" "+str(nMPheno))
            newCPed=cPed(nPedId, nFPheno, nFFAllele, nFMAllele, nFFMarker, nFMMarker, nMPheno, nMFAllele, nMMAllele, nMFMarker, nMMMarker)
            
            #update child
            for i in range(1, self.nChild+1):
                nCPheno=self.__genotype2pheno( int( arr[ dicIndex[ ("dg"+str(i)) ] ] ) )
                nCFAllele=int( arr[ dicIndex[ ("blue_da"+str(i)) ] ])
                nCMAllele=int( arr[ dicIndex[ ("pink_da"+str(i)) ] ])
                nCFMarker=int( arr[ dicIndex[ ("blue_la"+str(i)) ] ])
                nCMMarker=int( arr[ dicIndex[ ("pink_la"+str(i)) ] ])

                #print( str(nCPheno)+" "+str(nCFAllele)+" "+str(nCMAllele)+" "+str(nCFMarker)+" "+str(nCMMarker) )
                newCPed.addChild( nCPheno, nCFAllele, nCMAllele, nCFMarker, nCMMarker )

            self.lPed.append( newCPed )
        f.close()

    def getPed(self):
        return self.lPed
    
    def getPedDoubleHeteroNonAffectedParent(self, a_lPed):
        lFilteredPed=[]
        for i in a_lPed:
            if self.__IsDoubleHeterozygousParent(i) == False:
                continue
            if self.__IsNonAffectedParent(i) == False:
                continue
           
            lFilteredPed.append( i )
        return lFilteredPed
    
    def getPedHavingAtLeastOneAffectedChild(self, a_lPed):
        lFilteredPed=[]
        for i in a_lPed:
            
            if self.__IsHavingAffectedChild(i) == False:
                continue
           
            lFilteredPed.append( i )
        return lFilteredPed
            
    def __IsHavingAffectedChild(self, a_cPed):
        if a_cPed==None:
            return False

        for i in a_cPed.getChild():
            if i.getPheno()==2:
                return True
        
        return False
                

    def __IsNonAffectedParent(self, a_cPed):
        if a_cPed==None:
            return False

        cFather=a_cPed.getParent(1)
        if cFather==None:
            return False

        cMother=a_cPed.getParent(2)         
        if cMother==None:
            return False;
 
        if cFather.getPheno() == 2 or cMother.getPheno() ==2:
            bUnaffectedParent=False
        else:
            bUnaffectedParent=True
        
        return bUnaffectedParent
        
    def __IsDoubleHeterozygousParent(self, a_cPed):
        #both parents are double heterozygous at marker loci
        #do not care the number of affected children
        #only informative when both parents are not affected, but at least one child is affected
        #double heterozygous parents
        
        if a_cPed==None:
            return False

        cFather=a_cPed.getParent(1)
        if cFather==None:
            return False

        cMother=a_cPed.getParent(2)         
        if cMother==None:
            return False;

        lMarker=[cFather.getHaplo(1).getMarker(), cFather.getHaplo(2).getMarker(), cMother.getHaplo(1).getMarker(), cMother.getHaplo(2).getMarker() ]
        if len( set(lMarker) ) == 4:
            return True
        else:
            return False
        
    def __genotype2pheno(self, a_nValue):
        if a_nValue == 3: #genotype 2/2
            return 2
            #phenotype 2 (disease)
        else:
            return 1
        
    def print(self, a_lPed):
        print(len(a_lPed))
        for i in a_lPed:
            print(i.getLine())
    
