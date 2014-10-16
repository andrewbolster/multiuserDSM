    def _workload_calc(self,workload):
        warpcount=((workload/warpsize)+(0 if ((workload%self.warpsize)==0)else 1))
        warpperblock=max(1,min(8,warpcount))
        threadCount=warpsize * warpperblock
        blockCount=min(self.gridmax/threadCount,max(1,(warpcount/warpperblock)+(0 if ((warpcount%warpperblock)==0)else 1))) 
        return (threadCount,blockCount)
