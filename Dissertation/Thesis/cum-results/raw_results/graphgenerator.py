#!/usr/bin/env python
'''
Created on 21 Mar 2011

@author: bolster
'''

import numpy as np
import matplotlib
import pylab as pl
import utility as util
import os, string

class graphgen():
    scenario=""
    noisefloor = -140
    def __init__(self,scenario=scenario,textfile=False,type='cm'):
        '''
        Set up some data
        '''

        self.dest=util.graphdir
        
        if textfile:
            for text in textfile
                self.scenario=text.split('/')[-1].split('.')[0]
                self.importtext(text,type)
        else:
            self.scenario=scenario
            self.importdata()

    def importtext(self,textfile,type='cm'):
        assert os.path.isfile(textfile), "No Such File %s"%textfile
        #        Assume 224 channels
        import fileinput
        if type == 'cm':
            cma=[]
            pl.close()
            k=0
            for line in fileinput.input(textfile):
                splitline=string.split(line)
                if (line[0]=='#'):
                    pass
                else:
                    cma.append(splitline[1:])
            cm=np.zeros((len(cma),len(splitline[1:])))
            for k in range(224):
                cm[k]=cma[k]
            channels=range(224)
            nlines=int(np.sqrt(cm.shape[1]))
            for x in range(nlines):
                for v in range(nlines):
                    pl.plot(channels,cm[:,x*nlines+v],label="%d,%d"%(v,x))
            print cm[223]
            
            pl.xlabel("Subchannel Index")
            pl.ylabel("Gain (dB)")
            pl.yscale("log")
            pl.legend()
            pl.title("Plot of inter-line crosstalk gains for %d lines"%nlines)
            pl.savefig(self.dest+self.scenario+'-channelmatrix.png')
        elif type == 'bp':
            ba=[]
            pa=[]
            pl.close()
            k=0
            for line in fileinput.input(textfile):
                splitline=string.split(line)
                if (line[0]=='#'):
                    pass
                else:
                    ba.append([splitline[1],splitline[3]])
                    pa.append([splitline[2],splitline[4]])
            self.b=np.zeros((len(ba),2))
            self.p=np.zeros((len(pa),2))
            for k in range(224):
                self.b[k]=ba[k]
                self.p[k]=[util.dbmhz_to_watts(float(pa[k][i])) for i in range(len(pa[k]))]

            self.graph_p()
            self.graph_b()


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
        channels=range(self.cm.shape[0])
        for x in range(self.cm.shape[1]):
            for v in range(self.cm.shape[2]):
                pl.plot(channels,self.cm[:,x,v],label="%d,%d"%(x,v)) #this may be the wrong slicing style
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
    import sys
    g=graphgen(textfile=sys.argv[2:],type=sys.argv[1])
