'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy as np
import scipy.misc
import multiprocessing
from time import time
import sys

import utility as util

import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.tools as ctools

import numpy
from pycuda.compiler import SourceModule
from pycuda.gpuarray import GPUArray
from jinja2 import Template

t_kernels=Template("""
#include <pycuda-helpers.hpp>
#define MAT1 {{matrixN}}
#define MAT2 MAT1*MAT1
#define FAILVALUE {{failvalue}}
#define TINY 1.0e-40
#define a(i,j) a[(i)*MAT1+(j)]

#define GO 1
#define NOGO 0

texture<float, cudaTextureType2D,cudaReadModeElementType> XTG;

__device__ void d_pivot_decomp(float *a, int *p, int *q){
    int i,j,k;
    int n=MAT1;
    int pi,pj,tmp;
    float max;
    float ftmp;
    for (k=0;k<n;k++){
        pi=-1,pj=-1,max=FAILVALUE;
        //find pivot in submatrix a(k:n,k:n)
        for (i=k;i<n;i++) {
            for (j=k;j<n;j++) {
                if (fabs(a(i,j))>max){
                    max = fabs(a(i,j));
                    pi=i;
                    pj=j;
                }
            }
        }
        //Swap Row
        tmp=p[k];
        p[k]=p[pi];
        p[pi]=tmp;
        for (j=0;j<n;j++){
            ftmp=a(k,j);
            a(k,j)=a(pi,j);
            a(pi,j)=ftmp;
        }
        //Swap Col
        tmp=q[k];
        q[k]=q[pj];
        q[pj]=tmp;
        for (i=0;i<n;i++){
            ftmp=a(i,k);
            a(i,k)=a(i,pj);
            a(i,pj)=ftmp;
        }
        //END PIVOT

        //check pivot size and decompose
        if ((fabs(a(k,k))>TINY)){//should always be true with pivoting
            for (i=k+1;i<n;i++){
                //Column normalisation
                ftmp=a(i,k)/=a(k,k);
                for (j=k+1;j<n;j++){
                    //a(ik)*a(kj) subtracted from lower right submatrix elements
                    a(i,j)-=(ftmp*a(k,j));
                }
            }
        }
        //END DECOMPOSE
    }
}


__device__ void d_solve(float *a, float *x, int *p, int *q){
    //forward substitution; see  Golub, Van Loan 96
    //And see http://www.cs.rutgers.edu/~richter/cs510/completePivoting.pdf
    int i,j,pi;
    float ftmp;
    float xtmp[MAT1];
    //Swap rows (x=Px)
    for (i=0; i<MAT1; i++){
        pi=p[i];
        xtmp[i]=x[pi]; //value that should be here
    }
    //Lx=x
    //partially taken from Sourcebook on Parallel Computing p577
    for (i=0;i<MAT1;i++){
        ftmp=xtmp[i];
        for (j=0;j<i;j++)
            ftmp-=a(i,j)*xtmp[j];
        xtmp[i]=ftmp; //Unit lower triangular so second division unnecessary
    }
    //backward substitution
    //solves Uy=z
    xtmp[MAT1-1]/=a(MAT1-1,MAT1-1);
    for (i=MAT1-2;i>=0;i--){
        ftmp=xtmp[i];
        for (j=i+1;j<MAT1;j++){
            ftmp-=a(i,j)*xtmp[j];
        }
        xtmp[i]=(ftmp)/a(i,i);//non unit upper triangular so this division is necessary
    }

    //Last bit
    //solves x=Qy
    for (i=0;i<MAT1;i++){
        pi=q[i];
        x[i]=xtmp[pi];
    }
}

__global__ void solve(float *A, float *B){
  //Each thread solves the A[id]x[id]=b[id] problem
  int id= blockDim.x*blockIdx.x + threadIdx.x;
  int p_pivot[MAT1],q_pivot[MAT1];
  if ((GO==1)){
    for (int i=0;i<MAT1;i++) {
        p_pivot[i]=q_pivot[i]=i;
    }

    d_pivot_decomp(&A[id*MAT2],&p_pivot[0],&q_pivot[0]);
    d_solve(&A[id*MAT2],&B[id*MAT1],&p_pivot[0],&q_pivot[0]);
  }
}

//===============================================================================
// CALC_PSD ACCESSORY FUNCTIONS for set channel (woo-hoo!)
//===============================================================================
//Generate the A and B for all possible bitloads (in this offset)
//requires grid(MBPT^N,1,1) block(N,1,1)
//where MBPT^(N-1)>65535, use offset to continue
//thread.y's collaboratively populate A and B for their id
//This probably hammers memory...
__global__ void lk_prepare_permutations(float *A, float *B, int offset){
    //Don't need k as its sorted at the host stage for the creation of xtg
    int j=threadIdx.x;
    int myid=blockIdx.x;
    int bitbangval=myid+offset;

    int bitload[MAT1], i;
    
    //rebase myid to base (MBPT)
    //Unfortunately theres no way around every thread working out its own bitload :( 
    for (i=0; i<MAT1; i++){
        bitload[i]=bitbangval%{{maxbitspertone}};
        bitbangval/={{maxbitspertone}};
    }
    //FIXME Make handle out of range id's
    for (i=0; i<MAT1; i++){
        //Generate a row of A for this permutation and victim y
        //A[myid*MAT2+j*MAT1+i]=-({{channelgap}}*((1<<bitload[j])-1)*d_XTG[i*MAT1+j])/d_XTG[j*MAT1+j];
        A[myid*MAT2+j*MAT1+i]=-({{channelgap}}*((1<<bitload[j])-1)*tex2D(XTG,i,j))/tex2D(XTG,j,j);

    }
    //Generate an item of B
    //B[myid*MAT1+j]=({{noise}}*{{channelgap}}*((1<<bitload[j])-1))/d_XTG[j*MAT1+j];    
    B[myid*MAT1+j]=({{noise}}*{{channelgap}}*((1<<bitload[j])-1))/tex2D(XTG,j,j);

    //Repair an item of A
    //__syncthreads(); //Seems to help with memory coalescing
    A[blockIdx.x*MAT2+j*MAT1+j]=1;
}

//Solve all A and B psds together. 
//requires grid(MBPT^N/threadmax,1,1) block(threadmax,1,1)
__global__ void solve_permutations(float *A, float *B, int offset){
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int bitbangval=id+offset;
    int p_pivot[MAT1],q_pivot[MAT1];
    int i;

    //simulate bitload generation for in-place id check, and pivots at the same time
    for (i=0; i<MAT1; i++){
        bitbangval/={{maxbitspertone}};
        p_pivot[i]=q_pivot[i]=i;
    }
    //Stopper for invalid id's (where bitcombinations is less than maximum blocksize}
    if (bitbangval==0){
        //do the magic
        d_pivot_decomp(&A[id*MAT2],&p_pivot[0],&q_pivot[0]);
        d_solve(&A[id*MAT2],&B[id*MAT1],&p_pivot[0],&q_pivot[0]);
    }
}

//Finally Calculate the LK_Max_permutations
__global__ void lk_max_permutations(float *P, float *LK, float *lambdas, float *w){
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int bitbangval=id;
    float lk=0;
    int bitload[MAT1], i, broken=0;
    
    //At this point, B is populated with the P results.
    for (i=0;i<MAT1;i++){
        bitload[i]=bitbangval%{{maxbitspertone}};
        bitbangval/={{maxbitspertone}};
    }
    if (bitbangval==0){//check for out of range id's
        for (i=0;i<MAT1;i++){
           //Need to check for negative B's
            if (P[id*MAT1+i]<0)
                broken++;
            lk+=(bitload[i]*w[i])-(lambdas[i]*P[id*MAT1+i]);
        }
        
        //If anything is broken return a failing value (around -inf)
        if (broken==0)
            LK[id]=lk;
        else
            LK[id]=FAILVALUE;
    }
}
""")

class GPU(object):
    def __init__(self,bundle):
        self.start=time()
        self.bundle=bundle
        self.N=self.bundle.N
        self.K=self.bundle.K
        self.gamma=bundle.get_GAMMA()
        self.noise=bundle.get_NOISE()
        self.mbpt=bundle.get_MBPT()
        r_kernels=t_kernels.render(matrixN=self.N,
                               channelgap=pow(10,(self.gamma+3)/10), #19.7242
                               noise=self.noise, #4.313e-14
                               maxbitspertone=self.mbpt,
                               failvalue=np.float32(-sys.maxint)
                               )
        self.kernels = SourceModule(r_kernels)
        self.init=time()

    #In progress, ignore this.
    def calc_psd(self,bitload,k,gamma,noise):
        A=np.zeros((self.N*self.N))
        B=np.zeros((self.N))
    
    #Arbitrary solver for destructive Ax=x
    def solve(self,a,b,max):
        d_a=cuda.mem_alloc(a.nbytes)
        d_b=cuda.mem_alloc(b.nbytes)
        cuda.memcpy_htod(d_a,a)
        cuda.memcpy_htod(d_b,b)
        h_b=np.empty_like(b)
        self.go=time()
        
        #Go solve
        go=self.kernels.get_function("solve")
        go(d_a,d_b,block=(1,1,1))
        cuda.memcpy_dtoh(h_b,d_b)
        self.done=time()
        return h_b
    
    def lkmax(self,lambdas,w,xtalk_gain,k):
        #lk_prepare_permutations(float *A, float *B, float *xtg, int *offset)
        #lk_max_permutations(float *A, float *B, int *offset, float *LK, float *lambdas, float *w)
        #Number of expected permutations
        Ncombinations=pow(self.mbpt,self.N)
        
        #loop construct here if Ncombinations > 65535
        global_lk_maxid=-1
        gridmax=65535
        gridsize=min(pow(self.mbpt,(self.N)),gridmax)
        monitor=97

        #Check if this is getting hairy
        (free,total)=cuda.mem_get_info()
        mydev=cuda.Context.get_device()
        threadmax=mydev.get_attribute(cuda.device_attribute.MAX_THREADS_PER_BLOCK)
        #Playing Safe
        threadmax=int(np.floor(threadmax/2))
                
        if False:
            util.log.info("Working on %d combinations for K:%d, Mem %d%% Free"%(Ncombinations,k,(free*100/total)))
        for o in range(0,Ncombinations,gridsize):
            #offset 
            offset = np.int32(o);
            
            #Memories that are presistent
            #A[Nc*N*N],
            d_A=cuda.mem_alloc(np.zeros((gridsize*self.N*self.N)).astype(np.float32).nbytes)
            #B[Ncombinations*N]
            d_B=cuda.mem_alloc(np.zeros((gridsize*self.N)).astype(np.float32).nbytes)
            #LK (where final LK values come to rest...)
            d_lk=cuda.mem_alloc(np.empty((gridsize)).astype(np.float32).nbytes)
            #XTG[N*N] (for this k, clear this after prepare)
            #d_XTG=cuda.mem_alloc(np.zeros((self.N*self.N)).astype(np.float32).nbytes)
            #cuda.memcpy_htod(d_XTG,xtalk_gain.astype(np.float32))
            #Do XTG as a shared texture.
            t_XTG=self.kernels.get_texref("XTG");
            cuda.matrix_to_texref(xtalk_gain.astype(np.float32).copy(),t_XTG, order="F") #F indicates FORTRAN matrix addressing(column major)
        
            #Go prepare A and B
            prepare=self.kernels.get_function("lk_prepare_permutations")
            #if (k>monitor): self.meminfo(prepare,k,o)
            prepare(d_A,d_B,offset,texrefs=[t_XTG],grid=(gridsize,1),block=(self.N,1,1))
            try:
                cuda.Context.synchronize()
            except pycuda._driver.LaunchError:
                util.log.error("Failed on Prepare, Tone %d: XTG:%s"%(k,str(xtalk_gain)))
                raise
            
            if k in []:
                #Bring AB results back to host
                A=cuda.from_device(d_A,(gridsize,self.N,self.N),np.float32)
                B=cuda.from_device(d_B,(gridsize,self.N),np.float32)
                for g in []:
                    P=np.linalg.solve(A[g],B[g].T)
                    util.log.info("====G:%d\nA:%s\nB:%s\nP:%s"%(g,str(A[g]),str(B[g]),str(P)))

            #lambdas
            d_lambdas=cuda.mem_alloc(np.empty((self.N)).astype(np.float32).nbytes)
            cuda.memcpy_htod(d_lambdas,lambdas.astype(np.float32))
            #w
            d_w=cuda.mem_alloc(np.empty((self.N)).astype(np.float32).nbytes)
            cuda.memcpy_htod(d_w,w.astype(np.float32))
            
            #Go Solve
            lksolve=self.kernels.get_function("solve_permutations")
            threadshare_gridsize=int(np.floor(gridsize/threadmax))
            #if (k>monitor): self.meminfo(lksolve,k,o,threadmax)
            lksolve(d_A,d_B,offset, grid=(threadshare_gridsize,1), block=(threadmax,1,1))
            try:
                cuda.Context.synchronize()
            except pycuda._driver.LaunchError:
                util.log.error("Failed on Solve, Tone %d: \nXTG:%s"%(k,str(xtalk_gain)))
                raise            
            #Inter Kernel Housekeeping
            d_A.free()
            
            #Go Find the Max
            lkmax=self.kernels.get_function("lk_max_permutations")
            lkmax(d_B,d_lk,d_lambdas,d_w,grid=(threadshare_gridsize,1), block=(threadmax,1,1))
            try:
                cuda.Context.synchronize()
            except pycuda._driver.LaunchError:
                util.log.error("Failed on LKMax,Tone %d: XTG:%s"%(k,str(xtalk_gain)))
                raise

            #Bring LK results and power back to host
            lk=np.empty((gridsize)).astype(np.float32)
            cuda.memcpy_dtoh(lk,d_lk)
            
            #find the max lk
            lk_maxid=np.argmax(lk)
            
            #Hopefully this stuff goes on in the background
            
            d_lk.free()
            d_lambdas.free()
            d_w.free()
            
            if lk_maxid>global_lk_maxid:
                B=np.empty((gridsize,self.N),np.float32)
                cuda.memcpy_dtoh(B,d_B)
                P=B[lk_maxid]
                global_lk_maxid=lk_maxid
                if False:
                    util.log.info("Tone:%d,lkmaxid:%d,P:%s"%(k,lk_maxid,P))
            
            #'Free at last! Free at last! Thank God Almighty, we are free at last!'
            d_B.free()

        #end for
        bitload=self.bitload_from_id(global_lk_maxid)
        
        #If this works I'm gonna cry
        #print "GPU LKmax %d,%s:%s:%s"%(k,str(lk[lk_maxid]),str(bitload),str(P))
        return (P,bitload)
        
    def calc_psd(self,bitload,channelgap,noise,xtalk_gain):
        #all tones
        d_bitload=gpuarray.to_gpu(bitload)
        d_xtg=gpuarray.togpu(xtalk_gain)
        d_a=gpuarray.empty((self.N,self.N,self.K))
        d_b=gpuarray.empty((self.N,self.K))
        a=channelgap*(pow(2,d_bitload)-1)*xtalk_gain

    def optimise_p(self,lambdas,xtalk_gain):
        pass
        
    def meminfo(self,kernel,k=-1,o=-1,threads=1):
        shared=kernel.shared_size_bytes
        regs=kernel.num_regs
        local=kernel.local_size_bytes
        const=kernel.const_size_bytes
        mbpt=kernel.max_threads_per_block
        devdata=ctools.DeviceData()
        occupancy=ctools.OccupancyRecord(devdata,threads, shared_mem=shared,registers=regs)

        util.log.info("""(%03d,%d)=L:%d,S:%d,R:%d,C:%d,MT:%d,OC:%f"""%(k,o,local,shared,regs,const,mbpt,occupancy.occupancy))
        
    def gpu_test(self):
        self.test=time()
        matrixcount=224;
        a=np.tile(np.asarray([1,3,-2,3,5,6,2,4,3]),matrixcount).astype(np.float32)
        b=np.tile(np.asarray([5,7,8]),matrixcount).astype(np.float32)
        b=self.solve(a, b, matrixcount)
        
        print("Times:\nInit:%f\nGo:%f\nExec:%f"%(self.init-self.start,self.go-self.test,self.done-self.go))
        print b
    
    def bitload_from_id(self,id):
        bitload=np.zeros(self.N)
        for i in range(self.N):
            bitload[i]=id%self.mbpt;
            id/=self.mbpt;
        return bitload

if __name__ == "__main__":
    gpu=GPU(3)
    gpu.gpu_test()
