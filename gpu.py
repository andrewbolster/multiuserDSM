'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy as np
import threading, Queue, collections
from time import time
import sys

import utility as util
try:
    import pycuda
    import pycuda.driver as cuda
    import pycuda.tools as ctools
    from pycuda.compiler import SourceModule
    from pycuda.gpuarray import GPUArray
    anythingbroken=False
except ImportError:
    anythingbroken=True

import numpy
from jinja2 import Template

#Control how smart you want to be
adapt=True

t_kernels=Template("""
#include <pycuda-helpers.hpp>
#define MAT1 {{matrixN}}
#define MAT2 MAT1*MAT1
#define FAILVALUE {{failvalue}}
#define FPT {{floatingpointtype}}
#define CHANNELGAP {{channelgap}}
#define NOISE {{noise}}
#define MBPT {{maxbitspertone}}
#define TINY 1.0e-40
#define MAXID {{maxid}}
#define a(i,j) a[(i)*MAT1+(j)]

#define GO 1
#define NOGO 0

__device__ void d_pivot_decomp(FPT *a, int *p, int *q){
    int i,j,k;
    int n=MAT1;
    int pi,pj,tmp;
    FPT max;
    FPT ftmp;
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


__device__ void d_solve(FPT *a, FPT *x, int *p, int *q){
    //forward substitution; see  Golub, Van Loan 96
    //And see http://www.cs.rutgers.edu/~richter/cs510/completePivoting.pdf
    int i,j,pi;
    FPT ftmp;
    FPT xtmp[MAT1];
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

__global__ void solve(FPT *A, FPT *B){
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
__global__ void lk_osbprepare_permutations(FPT *A, FPT *B, FPT *d_XTG, int offset){
    //Don't need k as its sorted at the host stage for the creation of xtg
    int j=threadIdx.x;
    int myid=blockIdx.x;
    int bitbangval=myid+offset;

    int bitload[MAT1], i;
    
    //rebase myid to base (MBPT)
    //Unfortunately theres no way around every thread working out its own bitload :( 
    for (i=0; i<MAT1; i++){
        bitload[i]=bitbangval%MBPT;
        bitbangval/=MBPT;
    }
    if (myid+offset<MAXID){
      for (i=0; i<MAT1; i++){
          //Generate a row of A for this permutation and victim y
          A[myid*MAT2+j*MAT1+i]=-(CHANNELGAP*((1<<bitload[j])-1)*d_XTG[i*MAT1+j])/d_XTG[j*MAT1+j];
          //A[myid*MAT2+j*MAT1+i]=-(CHANNELGAP*((1<<bitload[j])-1)*tex2D(XTG,i,j))/tex2D(XTG,j,j);
      }
    }
    //Generate an item of B
    B[myid*MAT1+j]=(NOISE*CHANNELGAP*((1<<bitload[j])-1))/d_XTG[j*MAT1+j];    
    //B[myid*MAT1+j]=(NOISE*CHANNELGAP*((1<<bitload[j])-1))/tex2D(XTG,j,j);

    //Repair an item of A
    //__syncthreads(); //Seems to help with memory coalescing
    A[blockIdx.x*MAT2+j*MAT1+j]=1;
}

//Solve all A and B psds together. 
//requires grid(MBPT^N/threadmax,1,1) block(threadmax,1,1)
__global__ void solve_permutations(FPT *A, FPT *B, int offset){
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int bitbangval=id+offset;
    int p_pivot[MAT1],q_pivot[MAT1];
    int i;

    //simulate bitload generation for in-place id check, and pivots at the same time
    for (i=0; i<MAT1; i++){
        bitbangval/=MBPT;
        p_pivot[i]=q_pivot[i]=i;
    }
    //Stopper for invalid id's (where bitcombinations is less than maximum blocksize}
    if (id+offset<MAXID){
        //do the magic
        d_pivot_decomp(&A[id*MAT2],&p_pivot[0],&q_pivot[0]);
        d_solve(&A[id*MAT2],&B[id*MAT1],&p_pivot[0],&q_pivot[0]);
    }
}

//Finally Calculate the LK_Max_permutations
__global__ void lk_max_permutations(FPT *P, FPT *LK, FPT *lambdas, FPT *w, int offset){
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int bitbangval=id;
    FPT lk=0;
    int bitload[MAT1], i, broken=0;
    
    //At this point, B is populated with the P results.
    for (i=0;i<MAT1;i++){
        bitload[i]=bitbangval%MBPT;
        bitbangval/=MBPT;
    }
    if (id+offset<MAXID){//check for out of range id's
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
    else
        LK[id]=FAILVALUE;
}
//===============================================================================
// ISB Functions
//===============================================================================
//Populate the B array with the max permutation for each user
__device__ void isb_peruser_lkmax_B(FPT *LK, int *B, int lineid){
    int i, pmax=-1;
    FPT lkmax=FAILVALUE;
    for (i=0;i<MBPT;i++){
        if (lkmax < LK[lineid*MBPT+i]){
            pmax=i;                     //This is the maxlk bitload for this line
            lkmax=LK[lineid*MBPT+i];
        }    
    }
    B[lineid]=pmax;
}
//Do everything at once
__global__ void isb_optimise_pk(FPT *A, FPT *P, FPT *d_XTG, FPT *LK, FPT *lambdas, FPT *w, int *current_b, int offset){
    int lineid=blockIdx.x;                 //The User that we're playing with
    int permutation=threadIdx.x;            //the permutation we're generating
    int otherline=threadIdx.y;              //The other Line we're calculating for row of A/B
    int index=lineid*MBPT+permutation;     //The particular array index we're building
    int g_index=index*MAT1*MBPT+otherline; //used as a threadshare go-between
    
    int bitload[MAT1], i, broken=0;
    int p_pivot[MAT1],q_pivot[MAT1];
    FPT lk=0;


    //copy current bitload into thread memory
    for (i=0; i<MAT1; i++){
        bitload[i]=current_b[i];
        p_pivot[i]=q_pivot[i]=i;
    }
    
    if (index<MAT1*MBPT){
      //make this threads permutation
      bitload[lineid]=permutation;
      for (i=0; i<MAT1; i++){
          //Generate a row of A for this permutation and victim y
          A[index*MAT2+otherline*MAT1+i]=-(CHANNELGAP*((1<<bitload[otherline])-1)*d_XTG[i*MAT1+otherline])/d_XTG[otherline*MAT1+otherline];
      }
      //Generate an item of P
      P[index*MAT1+otherline]=(NOISE*CHANNELGAP*((1<<bitload[otherline])-1))/d_XTG[otherline*MAT1+otherline];    
    
      //Repair an item of A
      //__syncthreads(); //Seems to help with memory coalescing
      A[index*MAT2+otherline*MAT1+otherline]=1;
    }

    //Stopper for invalid id's (gets around multiple kernels)
    if (g_index<MAT1*MBPT){
        //do the magic
        d_pivot_decomp(&A[g_index*MAT2],&p_pivot[0],&q_pivot[0]);
        d_solve(&A[g_index*MAT2],&P[g_index*MAT1],&p_pivot[0],&q_pivot[0]);
    }
    
    if (index<MAT1*MBPT){//check for out of range id's
        for (i=0;i<MAT1;i++){
           //Need to check for negative P's
            if (P[index*MAT1+i]<0)
                broken++;
            lk+=(bitload[i]*w[i])-(lambdas[i]*P[index*MAT1+i]);
        }
        
        //If anything is broken return a failing value (around -inf)
        if (broken==0)
            LK[index]=lk;
        else
            LK[index]=FAILVALUE;
    }
    else
        LK[index]=FAILVALUE;
    
    if (g_index<MAT1){
        isb_peruser_lkmax_B(LK,current_b,g_index);
    }
}


""")

class GPU(object):
    def __init__(self,bundle):
        if anythingbroken:
            util.log.error("GPU imports failed miserably")
        self.bundle=bundle
        self.N=self.bundle.N
        self.K=self.bundle.K
        self.gamma=bundle.get_GAMMA()
        self.noise=bundle.get_NOISE()
        self.mbpt=bundle.get_MBPT()
        self.print_config=True

        #Set up context for initial setup
        cuda.init()
        mydev=cuda.Device(0)
        self.devcount=mydev.count()
        ctx=mydev.make_context()
        
        #Work out some context sensitive runtime parameters (currently assumes homogenous gpus)
        compute=mydev.compute_capability()
        self.threadmax=mydev.get_attribute(cuda.device_attribute.MAX_THREADS_PER_BLOCK)
        self.warpsize=mydev.get_attribute(cuda.device_attribute.WARP_SIZE)
        self.mps=mydev.get_attribute(cuda.device_attribute.MULTIPROCESSOR_COUNT)
        self.blockpermp=32
        self.gridmax=min(self.blockpermp*self.mps,65535)
        
        ctx.pop()
        ctx.detach()
        del ctx
        del mydev

        
        #Some more hardware based intelligence
        if (compute>=(1,3) and adapt):
            self.type=np.double
            typestr="double"
        else:
            self.type=np.float32
            typestr="float"
            
        #Pre-compile kernels
        self.r_kernels=t_kernels.render(matrixN=self.N,
                               channelgap=pow(10,(self.gamma+3)/10), #19.7242
                               noise=self.noise, #4.313e-14
                               maxbitspertone=self.mbpt,
                               failvalue=self.type(-sys.maxint),
                               floatingpointtype=typestr,
                               maxid=pow(self.mbpt,self.N)-1
                               )
               
        #Set up threading queues
        self.argqueue=Queue.Queue()
        self.resqueue=Queue.Queue()
        #spawn a threadpool
        self.threadpool = [None]*self.devcount
        util.log.info("Spawning %d GPU Threads"%self.devcount)
        for dev in range(len(self.threadpool)):
            self.threadpool[dev] = gpu_thread(self.argqueue,self.resqueue,self,device=dev)
            self.threadpool[dev].start()
        
    #destructor (for CUDA tidyness)
    def __del__(self):
        #kill threads
        for dev in range(len(self.threadpool)):
            del self.threadpool[dev]

    #Arbitrary solver for destructive Ax=x
    def solve(self,a,b,max):
        #context and kernel initialisation
        util.log.info("Initialising CUDA device")
        self.ctx = ctools.make_default_context()
        self.ctx.push()
        self.kernels=SourceModule(self.r_kernels)
        
        #Memory
        d_a=cuda.mem_alloc(a.astype(self.type).nbytes)
        d_b=cuda.mem_alloc(b.astype(self.type).nbytes)
        cuda.memcpy_htod(d_a,a.astype(self.type))
        cuda.memcpy_htod(d_b,b.astype(self.type))
        h_b=np.empty_like(b.astype(self.type))
        self.go=time()
        
        #Go solve
        go=self.kernels.get_function("solve")
        go(d_a,d_b,block=(1,1,1),grid=(1,1))
        cuda.memcpy_dtoh(h_b,d_b)
        self.done=time()
        self.ctx.pop()
        self.ctx.detach()
        return h_b

    def meminfo(self,kernel,k=-1,o=-1,threads=[],name=""):
        (free,total)=cuda.mem_get_info()
        shared=kernel.shared_size_bytes
        regs=kernel.num_regs
        local=kernel.local_size_bytes
        const=kernel.const_size_bytes
        mbpt=kernel.max_threads_per_block
        devdata=ctools.DeviceData()
        occupancy=ctools.OccupancyRecord(devdata,threads[0], shared_mem=shared,registers=regs)

        util.log.info("%s(%03d,%d)=L:%d,S:%d,R:%d,C:%d,MT:%d,T:%d,OC:%f,Free:%d"%(name,k,o,local,shared,regs,const,mbpt,threads[0],occupancy.occupancy,(free*100)/total))
        
    def gpu_test(self):
        self.test=time()
        matrixcount=224;
        a=np.tile(np.asarray([1,3,-2,3,5,6,2,4,3]),matrixcount).astype(np.float32)
        b=np.tile(np.asarray([5,7,8]),matrixcount).astype(np.float32)
        b=self.solve(a, b, matrixcount)
        
        print("Times:\nInit:%f\nGo:%f\nExec:%f"%(self.init-self.start,self.go-self.test,self.done-self.go))
        print b
      
    def osb_optimise_p(self,lambdas,w,xtalk_gain):
        funcname='osb_optimise_p'
        try:
            #construct the queue
            for k in range(self.K):
                self.argqueue.put((funcname,lambdas,w,xtalk_gain[k],k))
            
            #Wait for everything to end
            self.argqueue.join()
                  
            p=np.zeros((self.K,self.N)) #per tone per user power in watts
            b=np.asmatrix(np.zeros((self.K,self.N)))        
            while True:
                try:
                    queueitem=self.resqueue.get_nowait()
                    (func,(k,power,bitload))=queueitem
                except Queue.Empty:
                    break
                except ValueError:
                    util.log.error("Invalid Queueitem %s"%(str(queueitem)))
                    continue
                
                if func==funcname:
                    p[k]=power
                    b[k]=bitload
                    #util.log.info("%d:%s:%s"%(k,str(bitload),str(power)))
                else:
                    util.log.error("Invalid Functions %s"%(str(queueitem)))
            return (p,b)
        except(KeyboardInterrupt,SystemExit):
            util.log.error("Suicide Error: Results tainted, quitting incase")            
            raise Exception
        except pycuda._driver.MemoryError:
            util.log.error("Memory Error: Results tainted, quitting incase")
            raise MemoryError
            sys.exit(1)
    
    def isb_optimise_p(self,lambdas,w,xtalk_gain, bitload):
        funcname='isb_optimise_p'
        try:
            #construct the queue
            for k in range(self.K):
                self.argqueue.put((funcname,lambdas,w,xtalk_gain[k],k,bitload[k]))
            
            #Wait for everything to end
            self.argqueue.join()
                  
            p=np.zeros((self.K,self.N)) #per tone per user power in watts
            b=np.asmatrix(np.zeros((self.K,self.N)))        
            while True:
                try:
                    queueitem=self.resqueue.get_nowait()
                    (func,(k,power,bitload))=queueitem
                except Queue.Empty:
                    break
                except ValueError:
                    util.log.error("Invalid Queueitem %s"%(str(queueitem)))
                    continue
                
                if func==funcname:
                    p[k]=power
                    b[k]=bitload
                    #util.log.info("%d:%s:%s"%(k,str(bitload),str(power)))
                else:
                    util.log.error("Invalid Functions %s"%(str(queueitem)))
            return (p,b)
        except(KeyboardInterrupt,SystemExit):
            util.log.error("Suicide Error: Results tainted, quitting incase")            
            raise Exception
        except pycuda._driver.MemoryError:
            util.log.error("Memory Error: Results tainted, quitting incase")
            raise MemoryError
            sys.exit(1)
        
        
class gpu_thread(threading.Thread):
    def __init__(self,argqueue,resqueue,parent,device=0):
        threading.Thread.__init__(self)
        self.device=device
        self.argqueue=argqueue
        self.resqueue=resqueue
        self.warpsize=parent.warpsize
        self.gridmax=parent.gridmax
        self.N=parent.N
        self.mbpt=parent.mbpt
        self.type=parent.type
        self.print_config=parent.print_config
        self.r_kernels = parent.r_kernels
        self.local=threading.local()
        
    def run(self):
        try:
            #Initialise this device
            util.log.info("Initialising CUDA device")
            self.local.dev = cuda.Device(self.device)
            self.local.ctx = self.local.dev.make_context()
            self.local.ctx.push()
        except pycuda._driver.MemoryError:
            util.log.info("Balls")
            raise
            return
        
        #Initialise the kernel
        self.local.kernels=SourceModule(self.r_kernels)
                
        gridmax=65535
        
        
        #Kernels
        self.k_osbprepare=self.local.kernels.get_function("lk_osbprepare_permutations")
        self.k_osbsolve=self.local.kernels.get_function("solve_permutations")
        self.k_osblk=self.local.kernels.get_function("lk_max_permutations")
        self.k_solve=self.local.kernels.get_function("solve")
        self.k_isboptimise=self.local.kernels.get_function("isb_optimise_pk")


        #loop to empty queue
        while True:
            #grab args from queue (block until recieved)
            queueitem=self.argqueue.get()
            func=queueitem[0]
            args=queueitem[1:]
            
            if func=='osb_optimise_p':
                result=self.osb_optimise_p(*args)
                self.resqueue.put((func,result))
            elif func=='isb_optimise_p':
                result=self.isb_optimise_p(*args)
                self.resqueue.put((func,result))
            else:
                self.resqueue.put(None)
            
            self.argqueue.task_done()#nothing seems to get past this

        #end queue loop
        
    #Arbitrary solver for destructive Ax=x (Not finished)
    def solve(self,a,b,max):
                
        #Memory
        d_a=cuda.mem_alloc(a.astype(self.type).nbytes)
        d_b=cuda.mem_alloc(b.astype(self.type).nbytes)
        cuda.memcpy_htod(d_a,a.astype(self.type))
        cuda.memcpy_htod(d_b,b.astype(self.type))
        h_b=np.empty_like(b.astype(self.type))        
        #Go solve
        
        go(d_a,d_b,block=(1,1,1),grid=(1,1))
        cuda.memcpy_dtoh(h_b,d_b)
        d_a.free()
        d_b.free()
        return h_b
    
    def osb_optimise_p(self,lambdas,w,xtalk_gain,k):
        #Number of expected permutations
        Ncombinations=pow(self.mbpt,self.N)-1
        
        gridsize=min(Ncombinations,self.gridmax)

        #Check if this is getting hairy and assign grid/block dimensions
        warpcount=(Ncombinations/self.warpsize)+(0 if ((Ncombinations%self.warpsize)==0)else 1)
        warpperblock=max(1,min(8,warpcount))
        threadCount=self.warpsize * warpperblock
        blockCount=min(self.gridmax,max(1,warpcount/warpperblock))

        #How many individual lk's
        memdim=blockCount*threadCount

        N_grid=((memdim),1)
        N_block=(self.N,1,1)
        threadshare_grid=(blockCount,1)
        threadshare_block=(threadCount,1,1)
        
        monitor=255
        gpudiag=False
        prepdiag=False
        combodiag=True
        
        #Mallocs
        d_A=cuda.mem_alloc(np.zeros((memdim*self.N*self.N)).astype(self.type).nbytes)
        d_B=cuda.mem_alloc(np.zeros((memdim*self.N)).astype(self.type).nbytes)
        d_lk=cuda.mem_alloc(np.empty((memdim)).astype(self.type).nbytes)
        d_XTG=cuda.mem_alloc(np.zeros((self.N*self.N)).astype(self.type).nbytes)
        d_lambdas=cuda.mem_alloc(np.empty((self.N)).astype(self.type).nbytes)
        d_w=cuda.mem_alloc(np.empty((self.N)).astype(self.type).nbytes)
        
        #reset counter
        global_lk_maxid=-1

        #copy arguments to device
        cuda.memcpy_htod(d_XTG,xtalk_gain.astype(self.type))   
        cuda.memcpy_htod(d_lambdas,lambdas.astype(self.type))
        cuda.memcpy_htod(d_w,w.astype(self.type))
           
        #Print some information regarding the thread execution structure
        if (self.print_config):
            #Some info about combinations being run
            for o in range(0,self.Ncombinations,memdim):
                if combodiag:
                    (free,total)=cuda.mem_get_info()
                    util.log.info("Working on %d-%d combinations of %d for K:%d, L:%s, Mem %d%% Free"%(o,o+memdim,self.Ncombinations,k,str(lambdas),(free*100/total)))
            self.print_config=False
        
        #Perform LK Calculation and Maximisation for this channel, however many sectors it takes.
        for o in range(0,self.Ncombinations,memdim):
            #offset 
            offset = np.int32(o);

            #Go prepare A and B
            try:
                #prepare(d_A,d_B,offset,texrefs=[t_XTG],grid=default_grid,block=self.N_block)                                                                                                                                              
                self.k_osbprepare(d_A,d_B,d_XTG,offset,grid=self.N_grid,block=self.N_block)
                cuda.Context.synchronize()
            except (pycuda._driver.LaunchError,pycuda._driver.LogicError):
                util.log.error("Failed on Prepare,Tone %d: XTG:%s\nGridDim:%s,BlockDim:%s"%(k,str(xtalk_gain.flatten()),str(self.N_grid),str(self.N_block)))
                raise

            if prepdiag:
                #Bring AB results back to host
                A=cuda.from_device(d_A,(memdim,self.N,self.N),self.type)
                B=cuda.from_device(d_B,(memdim,self.N),self.type)
                np.save("A",A)
                np.save("B",B)
                for g in [223]:
                    P=np.linalg.solve(A[g],B[g].T)
                    if not (numpy.isfinite(P)).all():
                        util.log.info("====G:%d\nA:%s\nB:%s\nP:%s"%(g,str(A[g]),str(B[g]),str(P)))

            #Go Solve
            try:
                self.k_osbsolve(d_A,d_B,offset, grid=self.threadshare_grid, block=self.threadshare_block)
                cuda.Context.synchronize()
            except:
                util.log.error("Failed on Solve,Tone %d: XTG:%s\nGridDim:%s,BlockDim:%s"%(k,str(xtalk_gain.flatten()),str(self.threadshare_grid),str(self.threadshare_block)))
                raise            
            
            #Go Find the LK Values
            if (k>monitor) and gpudiag: self.meminfo(lkmax,k,o,self.threadshare_block,"Max")
            try:
                self.k_osblk(d_B,d_lk,d_lambdas,d_w,offset,grid=self.threadshare_grid,block=self.threadshare_block)
                cuda.Context.synchronize()
            except:
                util.log.error("Failed on LKMax,Tone %d: XTG:%s\nGridDim:%s,BlockDim:%s"%(k,str(xtalk_gain),str(self.threadshare_grid),str(self.threadshare_block)))
                raise

            #Bring LK results and power back to host
            lk=np.empty((memdim)).astype(self.type)
            cuda.memcpy_dtoh(lk,d_lk)
            
            #find the max lk
            lk_maxid=np.argmax(lk)
            
            if lk_maxid>global_lk_maxid:
                B=np.empty((memdim,self.N),self.type)
                cuda.memcpy_dtoh(B,d_B)
                cuda.Context.synchronize()
                P=B[lk_maxid]
                global_lk_maxid=lk_maxid
                bitload=util.bitload_from_id(global_lk_maxid,self.N,self.mbpt)
                if (k>monitor):
                    util.log.info("GPU LKmax %d,%s:%s:%s"%(k,str(lk[lk_maxid]),str(bitload),str(P)))

        #end for
        d_A.free()
        d_B.free()
        d_lk.free()
        d_lambdas.free()
        d_w.free()
        d_XTG.free()
        return (k,P,bitload)
        
    def isb_optimise_p(self,lambdas,w,xtalk_gain,k,current_bitload):
        '''
        __global__ void isb_bitload_permutations(FPT *A, FPT *B, FPT *d_XTG, FPT *current_b, int offset){
        __global__ void isb_solve_permutations(FPT *A, FPT *B, int offset){
        __global__ void isb_generate_lk(FPT *P, FPT *LK, FPT *lambdas, FPT *w, FPT *current_b, int offset){
        __global__ void isb_peruser_lkmax_B(FPT *LK, int *B){
        '''
        #Number of expected permutations
        Ncombinations=self.mbpt*self.N
        
        gridsize=min(Ncombinations,self.gridmax)

        #Check if this is getting hairy and assign grid/block dimensions
        warpcount=(Ncombinations/self.warpsize)+(0 if ((Ncombinations%self.warpsize)==0)else 1)
        warpperblock=max(1,min(8,warpcount))
        threadCount=self.warpsize * warpperblock
        blockCount=min(self.gridmax,max(1,warpcount/warpperblock))

        #How many individual lk's
        memdim=blockCount*threadCount
        assert memdim>Ncombinations, "Too many combinations for no loop construct"

        memdim=Ncombinations

        N_grid=((memdim),1)
        N_block=(self.N,1,1)
        
        #threadshare_grid=(blockCount,1)
        #threadshare_block=(threadCount,1,1)
        threadshare_grid=(self.N,1)
        threadshare_block=(self.mbpt,self.N,1)
        monitor=223
        gpudiag=False
        prepdiag=False
        combodiag=True
        
        #Mallocs
        d_A=cuda.mem_alloc(np.zeros((memdim*self.N*self.N)).astype(self.type).nbytes)
        d_B=cuda.mem_alloc(np.zeros((memdim*self.N)).astype(self.type).nbytes)
        d_lk=cuda.mem_alloc(np.empty((memdim)).astype(self.type).nbytes)
        d_XTG=cuda.mem_alloc(np.zeros((self.N*self.N)).astype(self.type).nbytes)
        d_lambdas=cuda.mem_alloc(np.empty((self.N)).astype(self.type).nbytes)
        d_bitload=cuda.mem_alloc(np.empty((self.N)).astype(np.int32).nbytes)
        d_w=cuda.mem_alloc(np.empty((self.N)).astype(self.type).nbytes)
        
        #reset counter
        global_lk_maxid=-1

        #copy arguments to device
        cuda.memcpy_htod(d_XTG,xtalk_gain.astype(self.type))   
        cuda.memcpy_htod(d_lambdas,lambdas.astype(self.type))
        cuda.memcpy_htod(d_w,w.astype(self.type))
        
        #Bitload locations
        bitload=np.asarray(current_bitload)[0].astype(np.int32)
        lastload=np.tile(-1,(self.N)).astype(np.int32)
           
        #Print some information regarding the thread execution structure
        o=0
        if (self.print_config):
            if combodiag:
                (free,total)=cuda.mem_get_info()
                util.log.info("Executing %d tests, tone:%d, L:%s, Mem %d%% Free"%(Ncombinations,k,str(lambdas),(free*100/total)))
                util.log.info("Grid:%s,Block:%s"%(str(threadshare_grid),str(threadshare_block)))

            self.print_config=False
        
        #Have to deal with non-convergence; simplest way is to keep a record of the past attempts and when it gets repeated, that'll do
        past_bitloads=[]
        bitloadcounter=collections.Counter()
        stuck=False
        
        #Perform LK Calculation and Maximisation for this channel, however many sectors it takes.
        #for o in range(0,Ncombinations,memdim):
        #offset 
        offset = np.int32(o);
        its=0
        while not (lastload==bitload).all() and not stuck:
            lastload=bitload.copy()
            cuda.memcpy_htod(d_bitload,bitload.astype(np.int32)) #int
            
            #Go prepare A and B
            try:
                #void isb_optimise_pk(FPT *A, FPT *B, FPT *d_XTG, FPT *LK, FPT *lambdas, FPT *w, FPT *current_b, int offset){
                self.k_isboptimise(d_A,d_B,d_XTG,d_lk,d_lambdas,d_w,d_bitload,offset,grid=threadshare_grid, block=threadshare_block)                 
                cuda.Context.synchronize()
            except (pycuda._driver.LaunchError,pycuda._driver.LogicError):
                util.log.error("Failed on Optimise,Tone %d: XTG:%s\nGridDim:%s,BlockDim:%s"%(k,str(xtalk_gain.flatten()),str(threadhshare_grid),str(threadshare_block)))
                raise
    
            #Bring peruser bitload results back to host
            cuda.memcpy_dtoh(bitload,d_bitload)
            cuda.Context.synchronize()
            its+=1
            #util.log.info("Bitload:%s, last:%s, Iteration:%d"%(str(bitload),str(lastload),its))
            past_bitloads.append(bitload.dumps())
            bitloadcounter.update(past_bitloads)
            (mostcommon,count)=bitloadcounter.most_common(1)[0]
            if count > 3:
                util.log.info("Got stuck on %s"%str(np.loads(mostcommon)))
                stuck=True
        
        #bring the power back for final result
        P=np.empty((memdim,self.N),self.type)
        cuda.memcpy_dtoh(P,d_B)
        cuda.Context.synchronize()
        #In theory, each P for the relevant bit permutation for each user should be the same
        for jump in range(self.N):
            jumped=(jump+1)%self.N
            assert (P[jump*self.mbpt+bitload[jump]]==P[jumped*self.mbpt+bitload[jumped]]).all(), "Assumptions wrong:P:%s"%str(P)
        
        P=P[bitload[0]]
        if k > monitor:
            (free,total)=cuda.mem_get_info()
            util.log.info("GPU LKmax %d,%s:%s: %d%% Free"%(k,str(bitload),str(P),(free*100/total)))
        
        if (bitload==0).any():
            util.log.info("Zero is allowed! %d:%s"%(k,str(bitload)))

        #end for
        d_A.free()
        d_B.free()
        d_lk.free()
        d_lambdas.free()
        d_w.free()
        d_XTG.free()
        return (k,P,bitload)
  
    def __del__(self):
        try:
            #kill the device
            self.local.ctx.pop()
            self.local.ctx.detach()
        except AttributeError:
            pass


if __name__ == "__main__":
    gpu=GPU(3)
    gpu.gpu_test()
