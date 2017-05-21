import math, time
from pedigree import *
from trio import *


def main():
    start_time=time.clock()
    
    cpmgr=cPedMgr("fam3000.theta0.2.csv", 4)
    lPedDBHNonAffP=cpmgr.getPedDoubleHeteroNonAffectedParent( cpmgr.getPed() )
    lPedAffectedChild=cpmgr.getPedHavingAtLeastOneAffectedChild( lPedDBHNonAffP )
    
    file_read_time=time.clock()
    fFileReadTime=file_read_time-start_time
    #cpmgr.print(lDHFilteredPed)
    
    
    lodmgr=cLODMgr( lPedAffectedChild )
    lTheta=lodmgr.getListTheta(0.0, 0.5, 0.05)
    
    bChildPhaseKnown=True
    bAffectedOnly=True
    lLODStat=lodmgr.calcListLOD(lTheta, bChildPhaseKnown, bAffectedOnly)
    lodmgr.print(lLODStat)
    analysis_time=time.clock()
    fAnalysisTime=analysis_time-file_read_time
    fTotalTime=fFileReadTime+fAnalysisTime
    print("Time(total)\t"+str(round(fTotalTime, 4)))
    print("Time(FileRead)\t"+str(round(fFileReadTime, 4)))
    print("Time(Analysis)\t"+str(round(fAnalysisTime, 4)))
    
    bChildPhaseKnown=True
    bAffectedOnly=False
    lLODStat=lodmgr.calcListLOD(lTheta, bChildPhaseKnown, bAffectedOnly)
    lodmgr.print(lLODStat)
    analysis_time=time.clock()
    fAnalysisTime=analysis_time-file_read_time
    fTotalTime=fFileReadTime+fAnalysisTime
    print("Time(total)\t"+str(round(fTotalTime, 4)))
    print("Time(FileRead)\t"+str(round(fFileReadTime, 4)))
    print("Time(Analysis)\t"+str(round(fAnalysisTime, 4)))
              
    bChildPhaseKnown=False
    bAffectedOnly=True
    lLODStat=lodmgr.calcListLOD(lTheta, bChildPhaseKnown, bAffectedOnly)
    lodmgr.print(lLODStat)
    analysis_time=time.clock()
    fAnalysisTime=analysis_time-file_read_time
    fTotalTime=fFileReadTime+fAnalysisTime
    print("Time(total)\t"+str(round(fTotalTime, 4)))
    print("Time(FileRead)\t"+str(round(fFileReadTime, 4)))
    print("Time(Analysis)\t"+str(round(fAnalysisTime, 4)))
    
    bChildPhaseKnown=False
    bAffectedOnly=False
    lLODStat=lodmgr.calcListLOD(lTheta, bChildPhaseKnown, bAffectedOnly)
    lodmgr.print(lLODStat)
    analysis_time=time.clock()
    fAnalysisTime=analysis_time-file_read_time
    fTotalTime=fFileReadTime+fAnalysisTime
    print("Time(total)\t"+str(round(fTotalTime, 4)))
    print("Time(FileRead)\t"+str(round(fFileReadTime, 4)))
    print("Time(Analysis)\t"+str(round(fAnalysisTime, 4)))
    
    
    
#    print("========")
#    bChildPhaseKnown=True
#    lLODStat=lodmgr.calcListLOD(lTheta, bChildPhaseKnown)
#    lodmgr.print(lLODStat)
    
    """
    nTotal=2
    fTheta=0.1
    for i in [0, 1, 2]:
        nRecomb=i
        fOdds=math.pow( fTheta, nRecomb )*math.pow( 1-fTheta, nTotal-nRecomb ) / math.pow(0.5, nTotal )
        fHandLOD=math.log10(fOdds)
        print( "hand\t"+str(fHandLOD) )
    """ 
    
    
    ##todo list
    #expand to qualitative trait
    #larger input data set, with more polymorphic marker
main()
