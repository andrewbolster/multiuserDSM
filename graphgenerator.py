'''
Created on 21 Mar 2011

@author: bolster
'''

import numpy as np
import pylab as pl
import utility

class graphgen():
    scenario=""
    noisefloor = -140
    def __init__(self,scenario):
        """
        Set up some data
        """
        self.scenario=scenario
        self.importdata()
    
    def importdata(self):
        self.p=np.load(self.scenario+"-power.npy")
        self.b=np.load(self.scenario+"-bitrate.npy")
        self.cm=np.load(self.scenario+"-channelmatrix.npy")
    
    def graph_p(self):
        pl.close()
        channels=range(self.p.shape[0])
        for line in range(self.p.shape[1]):
            pl.plot(channels,map(utility.watts_to_dbmhz,self.p[:,line])) #this may be the wrong slicing style
        pl.xlabel("Subchannel Index")
        pl.ylabel("Power (dbmhz)")
        pl.title("Plot of per-tone power assignments for %d lines"%self.p.shape[1])
        pl.ylim([self.noisefloor,0])
        pl.savefig(self.scenario+'-power.png')
        
    def graph_b(self):
        pl.close()
        channels=range(self.b.shape[0])
        for line in range(self.b.shape[1]):
            pl.plot(channels,self.b[:,line]) #this may be the wrong slicing style
        pl.xlabel("Subchannel Index")
        pl.ylabel("Bitloading")
        pl.title("Plot of per-tone bit assignments for %d lines"%self.p.shape[1])
        pl.savefig(self.scenario+'-bitrate.png')
        
    def graph_cm(self):
        pl.close()
        channels=range(self.cm.shape[2])
        for x in range(self.cm.shape[0]):
            for v in range(self.cm.shape[1]):
                pl.plot(channels,map(utility.TodB,self.cm[x,v,:]),label="%d,%d"%(x,v)) #this may be the wrong slicing style
        pl.xlabel("Subchannel Index")
        pl.ylabel("Gain")
        pl.yscale("log")
        pl.legend()
        pl.title("Plot of inter-line crosstalk gains for %d lines"%self.p.shape[1])
        pl.savefig(self.scenario+'-channelmatrix.png')      
        

if __name__ == '__main__':
    g=graphgen("OSB-3k_5k-near_far")
    g.graph_p()
    g.graph_b()
    g.graph_cm()
    