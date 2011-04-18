'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy as np
import scipy.misc
import multiprocessing
from time import time

import pycuda.autoinit
import pycuda.driver as cuda
import numpy
from pycuda.compiler import SourceModule
from pycuda.gpuarray import GPUArray
from jinja2 import Template

t_solve=Template("""
#define MAT1 {{matrixN}}
#define MAT2 MAT1*MAT1
#define TINY 1.0e-40
#define a(i,j) a[(i)*MAT1+(j)]

#define GO 1
#define NOGO 0

__device__ void d_pivot_decomp(float *a, int *p, int *q){
    int i,j,k;
    int n=MAT1;
    int pi,pj,tmp;
    float max;
    float ftmp;
    for (k=0;k<n;k++){
        pi=-1,pj=-1,max=0.0;
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
        if ((fabs(a(k,k))>TINY)){
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
    int i,ii=0,j;
    float ftmp;
    float xtmp[MAT1];
    //Swap rows (x=Px)
    for (i=0; i<MAT1; i++){
        xtmp[i]=x[p[i]]; //value that should be here
    }
    //Lx=x
    for (i=0;i<MAT1;i++){
        ftmp=xtmp[i];
        if (ii != 0)
            for (j=ii-1;j<i;j++)
                ftmp-=a(i,j)*xtmp[j];
        else
            if (ftmp!=0.0)
                ii=i+1;
        xtmp[i]=ftmp;
    }
    //backward substitution
    //partially taken from Sourcebook on Parallel Computing p577
    //solves Uy=z
    xtmp[MAT1-1]/=a(MAT1-1,MAT1-1);
    for (i=MAT1-2;i>=0;i--){
        ftmp=xtmp[i];
        for (j=i+1;j<MAT1;j++){
            ftmp-=a(i,j)*xtmp[j];
        }
        xtmp[i]=(ftmp)/a(i,i);
    }
    for (i=0;i<MAT1;i++)

    //Last bit
    //solves x=Qy
    for (i=0;i<MAT1;i++){
        x[i]=xtmp[q[i]];
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

__global__ void optimise_p(float *lambdas ){
    int k=blockIdx.x;
    int x=threadIdx.x;
    int y=threadIdx.y;
}  
""")

class GPU(object):
    def __init__(self,bundle):
        self.start=time()
        self.bundle=bundle
        self.N=self.bundle.N
        self.K=self.bundle.K
        r_solve=t_solve.render(matrixN=self.N, mbpt=15)
        self.g_solve = SourceModule(r_solve)
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
        go=self.g_solve.get_function("solve")
        go(d_a,d_b,block=(1,1,1))
        cuda.memcpy_dtoh(h_b,d_b)
        self.done=time()
        return h_b
    
    def calc_psd(self,bitload,channelgap,noise,xtalk_gain):
        #all tones
        d_bitload=gpuarray.to_gpu(bitload)
        d_xtg=gpuarray.togpu(xtalk_gain)
        d_a=gpuarray.empty((self.N,self.N,self.K))
        d_b=gpuarray.empty((self.N,self.K))
        a=channelgap*(pow(2,d_bitload)-1)*xtalk_gain
        
        

    def optimise_p(self,lambdas,xtalk_gain):
        pass
        
        
    def gpu_test(self):
        self.test=time()
        matrixcount=224;
        a=np.tile(np.asarray([1,3,-2,3,5,6,2,4,3]),matrixcount).astype(np.float32)
        b=np.tile(np.asarray([5,7,8]),matrixcount).astype(np.float32)
        b=self.solve(a, b, matrixcount)
        
        print("Times:\nInit:%f\nGo:%f\nExec:%f"%(self.init-self.start,self.go-self.test,self.done-self.go))
        print b
        

if __name__ == "__main__":
    gpu=GPU(3)
    gpu.gpu_test()
