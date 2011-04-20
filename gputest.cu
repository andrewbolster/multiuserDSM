#include <utility>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sm_11_atomic_functions.h>

#define CUDA_CHK(NAME, ARGS) { \
  cudaError_t cuda_err_code = NAME ARGS; \
  if (cuda_err_code != cudaSuccess) { \
    printf("%s failed with code %d\n", #NAME, cuda_err_code); \
    abort(); \
  } \
}
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
#define MAT1 3
#define MAT2 MAT1*MAT1
#define TINY 1.0e-40
#define a(i,j) a[(i)*MAT1+(j)]

#define channelgap 19.7242
#define noise 4.313e-14
#define maxbitspertone 15

#define GO 1
#define NOGO 0

void Check_Kernel(const char *message){
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess){
    fprintf(stderr,"Error: %s:%s\n",message, cudaGetErrorString(error));
  }
}

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
//===============================================================================
// CALC_PSD ACCESSORY FUNCTIONS for set channel (woo-hoo!)
//===============================================================================
//Generate the A dna B for all possible bitloads (in this offset)
//requires grid(MBPT^N,1,1) block(N,1,1)
//where MBPT^(N-1)>65535, use offset to continue
//thread.y's collaboratively populate A and B for their id
//This probably hammers memory...
__global__ void lk_prepare_permutations(float *A, float *B, float *xtg, int *offset){
    //Don't need k as its sorted at the host stage for the creation of xtg
    int y=threadIdx.x;
    int myid=blockIdx.x+(*offset);
    int bitbangval=myid;

    int bitload[MAT1], i;
    
    //rebase myid to base (MBPT)
    //Unfortunately theres no way around every thread working out its own bitload :( 
    for (i=MAT1-1; i>=0; i--){
        bitload[i]=bitbangval%maxbitspertone;
        bitbangval/=maxbitspertone;
    }
    //Generate a row of A for this permutation and victim y
    for (i=0;i<MAT1;i++){
        A[myid*MAT2+y*MAT1+i]=-(channelgap*((1<<bitload[i])-1)*xtg[i*MAT1+y])/xtg[y*MAT1+y];
    }
    //Generate an item of B
    B[myid*MAT1+y]=(noise*channelgap*((1<<bitload[y])-1))/xtg[y*MAT1+y];
    //Repair an item of A
    __syncthreads();
    A[myid*MAT2+y*MAT1+y]=1;
}

//Solve all A and B psds together. 
//requires grid(MBPT^N,1,1) block(1,1,1)
__global__ void lk_max_permutations(float *A, float *B, int *offset, float *LK, float *lambdas, float *w){
    //MyID is now effectively 3D
    int id=blockIdx.x+(*offset);
    int bitbangval=id;
    int p_pivot[MAT1],q_pivot[MAT1];
    float lk=0;
    int bitload[MAT1], i;

    
    //generate bitload and pivots at the same time
    for (i=MAT1-1; i>=0; i--){
        bitload[i]=bitbangval%maxbitspertone;
        bitbangval/=maxbitspertone;
        p_pivot[i]=q_pivot[i]=i;
    }

    //do the magic
    d_pivot_decomp(&A[id*MAT2],&p_pivot[0],&q_pivot[0]);
    d_solve(&A[id*MAT2],&B[id*MAT1],&p_pivot[0],&q_pivot[0]);
    
    //At this point, B is populated with the P results.
    for (i=0;i<MAT1;i++){
        lk+=(bitload[i]*w[i])-(lambdas[i]*B[id*MAT1+i]);
    }
    LK[id]=lk;
}

int main(){
  //What are you actually trying to do:
  //  generate 2 input matrixes, (NxN,Nx1) and 1 output (1xN)
  //  do this over outerloop length for threadiding
  cudaSetDevice(0);
  const unsigned int matrixcount=4;
  const unsigned int ncombos=pow(maxbitspertone,MAT1);
  const unsigned int outerloop=ncombos;
  const unsigned int matsize=MAT2*outerloop;
  const unsigned int vecsize=MAT1*outerloop;
  //float a[]={1,3,-2,3,5,6,2,4,3};
  //const float exampleA[]={7,3,-11,-6,7,10,-11,2,-2};
  //const float exampleA[]={4,3,6,3};
  //const float b[]={5,7,8}; 
  //const float exampleB[]={4,5};
  const float x[]={5e-10,7e-20,8e-8,7e-20,8e-8,5e-10,8e-20,5e-8,7e-10}; //xtg analogue

  //memory allocations
  int* h_offset=(int*)malloc(sizeof(int));
  (*h_offset)=0;
  int* d_offset;
  CUDA_CHK(cudaMalloc, (&d_offset, sizeof(int)));
  CUDA_CHK(cudaMemcpy, (d_offset, h_offset, sizeof(int), cudaMemcpyHostToDevice));

  float* h_A = (float*)malloc(sizeof(float)*matsize);//empty till after
  float* h_b = (float*)malloc(sizeof(float)*vecsize);//empty till after

  float* h_l = (float*)malloc(sizeof(float)*vecsize);
  float* h_x = (float*)malloc(sizeof(float)*matsize);
  float* d_A;
  float* d_b;
  float* d_x;
  float* d_l;
  CUDA_CHK(cudaMalloc, (&d_A, sizeof(float)*matsize));
  CUDA_CHK(cudaMalloc, (&d_b, sizeof(float)*vecsize));
  CUDA_CHK(cudaMalloc, (&d_l, sizeof(float)*vecsize));
  CUDA_CHK(cudaMalloc, (&d_x, sizeof(float)*matsize));

  //don't need these on device till after permutation generation
  float* h_lk = (float*)malloc(sizeof(float)*outerloop);//empty till after
  float* h_w = (float*)malloc(sizeof(float)*vecsize);
  float* d_lk;
  float* d_w;

  printf("Mallocd\n");
  //fill matrix and vector with stuff
  for (unsigned int i = 0;i<outerloop;i++){
    //printf("\n%d\n",i);
    for (unsigned int j = 0; j < MAT1; j++){
      h_l[(i*MAT1)+j]=1.0;
      h_w[(i*MAT1)+j]=1.0;
      //printf("\n%d:",j);
      for (unsigned int k=0; k < MAT1; k++){
        //printf("%d,",k);
        //h_A[(i*MAT2)+(j*MAT1)+k]=a(j,k);
        h_x[(i*MAT2)+(j*MAT1)+k]=x[j*MAT1+k];
        
      }
    }
  }

  printf("Generated\n");
    //copy values to device
 // CUDA_CHK(cudaMemcpy, (d_A, h_A, sizeof(float)*matsize, cudaMemcpyHostToDevice));
  CUDA_CHK(cudaMemcpy, (d_l, h_l, sizeof(float)*vecsize, cudaMemcpyHostToDevice));
  CUDA_CHK(cudaMemcpy, (d_x, h_x, sizeof(float)*matsize, cudaMemcpyHostToDevice));

  printf("Copied\n");/*
  for (unsigned int i=0; i<outerloop; i++){
    printf("\n%d:x:A|l",i);
    //printf("%.3lf|",h_x[i*MAT1]);
    for (unsigned int j=0; j<MAT1; j++){
      printf("\n%g:",h_x[i*MAT1+j]);
      for (unsigned int k=0;k<MAT1; k++){
        printf("%g,",h_A[(i*MAT2)+(j*MAT1)+k]);
      }
      printf("|%g",h_l[i*MAT1+j]);
    }
  }
  puts("\n");
  */
  //parameters
  //dim3 blocksPerGrid((outerloop + threadsPerBlock.x -1)/threadsPerBlock.x,1,1);
  dim3 blocksPerGrid(pow(maxbitspertone,MAT1),1,1);
  dim3 threadsPerBlock(MAT1,1,1);
  
  printf("TPB:%d,BPG:%d\n",threadsPerBlock.x,blocksPerGrid.x);
  //Execute
  cudaEvent_t evt_start, evt_stop;
  CUDA_CHK(cudaEventCreate, (&evt_start));
  CUDA_CHK(cudaEventCreate, (&evt_stop));
  CUDA_CHK(cudaEventRecord, (evt_start,0));

  lk_prepare_permutations<<<blocksPerGrid,threadsPerBlock>>>(d_A,d_b,d_x,d_offset);
  cudaDeviceSynchronize();
  Check_Kernel("Generate");
  CUDA_CHK(cudaMemcpy, (h_A,d_A, sizeof(float)*matsize, cudaMemcpyDeviceToHost));
  CUDA_CHK(cudaMemcpy, (h_b,d_b, sizeof(float)*vecsize, cudaMemcpyDeviceToHost));
  CUDA_CHK(cudaMemcpy, (h_x,d_x, sizeof(float)*matsize, cudaMemcpyDeviceToHost));

  printf("Ran solve\n");

  

  CUDA_CHK(cudaFree, (d_x));
  CUDA_CHK(cudaMalloc, (&d_lk, sizeof(float)*vecsize));
  CUDA_CHK(cudaMalloc, (&d_w, sizeof(float)*vecsize));
  CUDA_CHK(cudaMemcpy, (d_w, h_w, sizeof(float)*vecsize, cudaMemcpyHostToDevice));
  dim3 threadsPerBlock_lkmax(1,1,1);
  lk_max_permutations<<<blocksPerGrid,threadsPerBlock>>>(d_A,d_b,d_offset,d_lk,d_l,d_w);
  cudaDeviceSynchronize();
  Check_Kernel("LKmax");
  CUDA_CHK(cudaMemcpy, (h_lk, d_lk, sizeof(float)*outerloop, cudaMemcpyDeviceToHost));
  CUDA_CHK(cudaEventRecord, (evt_stop, 0));
  CUDA_CHK(cudaEventSynchronize, (evt_stop));
  float total_time;
  CUDA_CHK(cudaEventElapsedTime, (&total_time, evt_start, evt_stop));
  float one_time = total_time * 1e-3;

  for (unsigned int i=0; i<outerloop; i++){
    printf("\n%d:x:A:LK:%g",i,h_lk[i]);
    //printf("%.3lf|",h_x[i*MAT1]);
    for (unsigned int j=0; j<MAT1; j++){
      printf("\n%g:",h_x[i*MAT1+j]);
      for (unsigned int k=0;k<MAT1; k++){
        printf("%g,",h_A[(i*MAT2)+(j*MAT1)+k]);
      }
    }
  }
  puts("\n");
  printf("time: %g s\n", one_time);

  cudaEventDestroy(evt_start);
  cudaEventDestroy(evt_stop);
  free(h_A);
  free(h_b);
  free(h_x);
  free(h_w);
  free(h_l);
  free(h_lk);
  CUDA_CHK(cudaFree, (d_A)); 
  CUDA_CHK(cudaFree, (d_b)); 
  CUDA_CHK(cudaFree, (d_l)); 
  CUDA_CHK(cudaFree, (d_lk)); 
  CUDA_CHK(cudaFree, (d_w)); 
}

