'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy as np
import scipy.misc
import multiprocessing
from time import time

from copperhead import *

@cu
def something_like_calc_psd(bitload,k,N,XTG):
    A=np.asmatrix(np.zeros((N,N)))
    B=np.asmatrix(np.zeros((N,1)))
    
    for v in range(N):
        B[v,0]=1<<(bitload[v]-1)/XTG[v,v]
        A[v]=-B[v,0]*XTG[:,v]
        A[v,v]=1
    
    P=np.linalg.solve(A,B)
    return P

class GPU(object):
    def update_delta_p(self,tone,N,P,b):
        _b=np.tile(util.mat2arr(b),N)
        delta_p=np.zeros(N,N)
        #create a square matrix of incremented bitloads
        for line in xrange(self.bundle.N):
            _b[line,line]+=1
            

if __name__ == "__main__":
    
    hasGPU = hasattr(places, 'gpu0')
    
    
    
    def run(fn,*args):
        #test on cpu
        cstart=time()
        cpuresult=fn(*args, targetPlace=places.here)
        cdur=time()-cstart
        if hasGPU:
            gstart=time()
            gpuresult=fn(*args, targetPlace=places.gpu0)
            gdur=time()-gstart
        assert list(cpuresult)==list(gpuresult), "Failed"
        print("C:%f,G:%f"*(cdur,gdur))
        
    N=4
    K=224
    
    #host test values
    b=np.random.randint(0,15,(K,N))
    p=np.random.randn(K,N)
    XTG=np.random.random_sample((N,N))
    dp=np.empty((K,N,N))
    A=np.empty((N,N))
    B=np.empty((N))
    P=np.empty((N,N))

    run(calc_ps,d(b[K/2],K/2,N,XTG))
        
            
            
            