'''
Created on 21 Mar 2011

@author: bolster
'''

import numpy as np
import pylab as pl
import utility as util
import os

class graphgen():
    scenario=""
    noisefloor = -140
    def __init__(self,scenario):
        '''
        Set up some data
        '''
        self.scenario=scenario
        self.dest=util.graphdir
        
        self.importdata()
    
    def importdata(self):
        assert os.path.isdir(util.rawdir), "No Path, Dunno wut happened there...%s"%util.graphdir
        self.p=np.load(util.rawdir+self.scenario+"-power.npy")
        self.b=np.load(util.rawdir+self.scenario+"-bitrate.npy")
        self.cm=np.load(util.rawdir+self.scenario+"-channelmatrix.npy")
    
    def graph_p(self):
        pl.close()
        channels=range(self.p.shape[0])
        for line in range(self.p.shape[1]):
            yvals=np.ma.masked_invalid(map(util.watts_to_dbmhz,self.p[:,line]))
            pl.plot(channels,yvals) #this may be the wrong slicing style
        pl.xlabel("Subchannel Index")
        pl.ylabel("Power (dbmhz)")
        pl.title("Plot of per-tone power assignments for %d lines"%self.p.shape[1])
        pl.savefig(self.dest+self.scenario+'-power.png')
        
    def graph_b(self):
        pl.close()
        channels=range(self.b.shape[0])
        for line in range(self.b.shape[1]):
            pl.plot(channels,self.b[:,line]) #this may be the wrong slicing style
        pl.xlabel("Subchannel Index")
        pl.ylabel("Bits/Frame")
        pl.title("Plot of per-tone bit assignments for %d lines"%self.p.shape[1])
        pl.savefig(self.dest+self.scenario+'-bitrate.png')
        
    def graph_cm(self):
        pl.close()
        channels=range(self.cm.shape[2])
        for x in range(self.cm.shape[0]):
            for v in range(self.cm.shape[1]):
                pl.plot(channels,self.cm[x,v,:],label="%d,%d"%(x,v)) #this may be the wrong slicing style
        pl.xlabel("Subchannel Index")
        pl.ylabel("Gain")
        pl.yscale("log")
        pl.legend()
        pl.title("Plot of inter-line crosstalk gains for %d lines"%self.p.shape[1])
        pl.savefig(self.dest+self.scenario+'-channelmatrix.png')
        
    def graph_all(self,graphdir=util.graphdir):
        if not os.path.isdir(graphdir):
            os.makedirs(graphdir)    
        self.graph_p()
        self.graph_b()
        self.graph_cm()
    

if __name__ == '__main__':
    g=graphgen("OSB-3k_5k-near_far")
    g.graph_all()