#! /usr/bin/env python

import numpy as np
import matplotlib
import pylab as pl
import os,string

import sys
results_arr=[]
timetable=np.zeros((5,4))

class results():
    def __init__(self,scenarioname):
        '''
        Assume being handed stats file
        '''
        statsfile=''.join([scenarioname])
        assert os.path.isfile(statsfile), "No such file %s"%statsfile
        blob=np.load(statsfile).tolist()
        self.algo=blob['algo']
        self.bundlerate=blob['bundlerate']
        self.duration=blob['duration']
        self.ngpu=blob['gpu']
        self.powers=blob['linepowers']
        self.rates=blob['linerates'].flatten()
        self.lines=blob['lines']
        self.ratesearch=blob['ratesearch']
        self.origname=scenarioname

def algosum(algostring):
    for result in results_arr:
        if result.algo==algostring:
            for n in range(1,6):
                if len(result.lines)==n:
                    if result.ngpu:
                        timetable[result.ngpu,n-2]=result.duration
                        print("%s,g:%d,n:%d, %f"%(result.algo,result.ngpu,n,result.duration))
                    else:
                        timetable[0,n-2]=result.duration
                        print("%s,g:%d,n:%d, %f"%(result.algo,0,n,result.duration))


    print "%s TABLE"%algostring
    pl.close()
    for nline in range(2,6):
        pl.plot(range(1,5),timetable[1:,nline-2],label="%d"%nline);
        print ("%d&%.3f&%.3f&%.3f&%.3f&%.3f\\\\"%(nline,timetable[0,nline-2],timetable[1,nline-2],timetable[2,nline-2],timetable[3,nline-2],timetable[4,nline-2]))

    pl.xlabel("N-Devices")
    pl.ylabel("Execution time (s)")
    pl.yscale("log")
    pl.legend()
    pl.title("Plot of execution time on multiple GPUs for range of bundle sizes")
    pl.savefig("%s-gpucompare.png"%algostring)
    
    pl.close()
    for nline in range(2,6):
        pl.plot(range(1,5),(timetable[1:,nline-2]/(nline)),label="%d"%nline);

    pl.xlabel("N-Devices")
    pl.ylabel("Ratio of Execution time (s) / N-lines")
    pl.yscale("log")
    pl.legend()
    pl.title("Plot of execution time/line count on multiple GPUs for range of bundle sizes")
    pl.savefig("%s-linetimecompare.png"%algostring)

    pl.close()
    for ngpu in range(1,5):
        pl.plot(range(2,6),timetable[ngpu,:],label="%d"%ngpu);

    pl.xlabel("Bundle Size")
    pl.ylabel("Execution time (s)")
    pl.yscale("log")
    pl.title("Plot of execution time on N-lines for range of GPU counts")
    pl.legend()
    pl.savefig("%s-linecompare.png"%algostring)

    pl.close()
    for nline in range(2,6):
        pl.plot(range(1,5),[timetable[0,nline-2]/ Tp for Tp in timetable[1:,nline-2]],label="%d"%nline);
        print [timetable[0,nline-2]/ Tp for Tp in timetable[1:,nline-2]]
    pl.xlabel("N-Devices")
    pl.ylabel("Parallel Speedup")
    pl.yscale("linear")
    pl.title("Plot of parallel speedup on multiple GPUs for range of bundle sizes")
    pl.legend()
    pl.savefig("%s-gpuspeedup.png"%algostring)
    
    pl.close()
    for nline in range(2,6):
        pl.plot(range(1,5),[(timetable[0,nline-2]/ timetable[n,nline-2])/(240*n) for n in range(1,len(timetable[1:,nline-2])+1)],label="%d"%nline);
        print[(timetable[0,nline-2]/ timetable[n,nline-2])/(240*n) for n in range(1,len(timetable[1:,nline-2])+1)]
    pl.xlabel("N-Devices")
    pl.ylabel("Parallel Efficiency")
    pl.yscale("linear")
    pl.title("Plot of parallel efficiency on multiple GPUs for range of bundle sizes")
    pl.legend()
    pl.savefig("%s-gpuefficiency.png"%algostring)




if __name__=='__main__':
    for arg in sys.argv[1:]:
        try:
            newresult=results(arg)
            results_arr.append(results(arg))
        except:
            print "Failed on %s:"%arg
            raise

    algosum('ISB')
    timetable=np.zeros((5,4))
    algosum('OSB')

