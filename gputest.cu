#include <utility>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define CUDA_CHK(NAME, ARGS) { \
  cudaError_t cuda_err_code = NAME ARGS; \
  if (cuda_err_code != cudaSuccess) { \
    printf("%s failed with %s\n", #NAME, cudaGetErrorString(cuda_err_code)); \
    abort(); \
  } \
}
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
#define MAT1 4
#define MAT2 MAT1*MAT1
#define TINY 1.0e-40
#define FAILVALUE -1.0e40
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
    int i,ii=0,j,pi;
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
        if (ii!=0){
          for (j=ii-1;j<i;j++)
              ftmp-=a(i,j)*xtmp[j];
        } else{
          if (ftmp !=0.0)
            ii=i+1;
        }
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
__global__ void lk_prepare_permutations(float *A, float *B, float *d_XTG, int offset){
    //Don't need k as its sorted at the host stage for the creation of xtg
    int j=threadIdx.x;
    int myid=blockIdx.x; //offset doesn't matter for this
    int bitbangval=myid;

    int bitload[MAT1], i;
    
    //rebase myid to base (MBPT)
    //Unfortunately theres no way around every thread working out its own bitload :( 
    for (i=0; i<MAT1; i++){
        bitload[i]=bitbangval%maxbitspertone;
        bitbangval/=maxbitspertone;
    }
    for (i=0; i<MAT1; i++){
        //Generate a row of A for this permutation and victim j
        //A[myid*MAT2+j*MAT1+i]=-(channelgap*((1<<bitload[j])-1)*tex2D(XTG,i,j))/tex2D(XTG,j,j);
        A[myid*MAT2+j*MAT1+i]=-(channelgap*((1<<bitload[j])-1)*d_XTG[i*MAT1+j])/d_XTG[j*MAT1+j];
    }
    //Generate an item of B
    B[myid*MAT1+j]=(noise*channelgap*((1<<bitload[j])-1))/d_XTG[j*MAT1+j];
    
    //Repair an item of A
    __syncthreads(); //Seems to help with memory coalescing
    A[myid*MAT2+j*MAT1+j]=1;
}

//Solve all A and B psds together. 
//requires grid(MBPT^N/threadmax,1,1) block(threadmax,1,1)
__global__ void solve_permutations(float *A, float *B, int offset){
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int bitbangval=id;
    int p_pivot[MAT1],q_pivot[MAT1];
    int i;

    //simulate bitload generation for in-place id check, and pivots at the same time
    for (i=0; i<MAT1; i++){
        bitbangval/=maxbitspertone;
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
        bitload[i]=bitbangval%maxbitspertone;
        bitbangval/=maxbitspertone;
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
int main(){
  //What are you actually trying to do:
  //  generate 2 input matrixes, (NxN,Nx1) and 1 output (1xN)
  //  do this over outerloop length for threadiding
  printf("Hello there\n");
  cudaSetDevice(0);
  printf("Set Device\n");
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
  //const float x[]={5e-10,7e-20,8e-8,7e-20,8e-8,5e-10,8e-20,5e-8,7e-10}; //xtg analogue
  const float x[]={//this one breaks
           7.15337843e-09,   9.98540799e-18,   8.27619149e-13,
           9.98540799e-18,   1.79722338e-07,   7.45386129e-06,
           1.79722338e-07,   5.17430336e-10,   8.27619149e-13,
           9.98540799e-18,   7.15337843e-09,   9.98540799e-18,
           1.79722338e-07,   5.17430336e-10,   1.79722338e-07,
           7.45386129e-06};
  const float x2[]={//this one should be fine
           7.66152695e-09,   1.08253155e-17,   8.72877254e-13,
           1.08253155e-17,   1.79434933e-07,   7.76722237e-06,
           1.79434933e-07,   5.30951476e-10,   8.72877254e-13,
           1.08253155e-17,   7.66152695e-09,   1.08253155e-17,
           1.79434933e-07,   5.30951476e-10,   1.79434933e-07,
           7.76722237e-06};
  const float x178[]={
         5.15676578e-11,   2.60163643e-20,   1.55231521e-14,
         2.60163643e-20,   1.74280845e-07,   3.86365544e-07,
         1.74280845e-07,   6.97834034e-11,   1.55231521e-14,
         2.60163643e-20,   5.15676578e-11,   2.60163643e-20,
         1.74280845e-07,   6.97834034e-11,   1.74280845e-07,
         3.86365544e-07};
  const float x179[]={//The 'new' broken one
         4.87663895e-11,   2.42887090e-20,   1.48197082e-14,
         2.42887090e-20,   1.73987318e-07,   3.73631407e-07,
         1.73987318e-07,   6.81261232e-11,   1.48197082e-14,
         2.42887090e-20,   4.87663895e-11,   2.42887090e-20,
         1.73987318e-07,   6.81261232e-11,   1.73987318e-07,
         3.73631407e-07};

  //memory allocations
  int h_offset=0;

  float* h_A = (float*)malloc(sizeof(float)*matsize);//empty till after
  float* h_b = (float*)malloc(sizeof(float)*vecsize);//empty till after

  float* h_l = (float*)malloc(sizeof(float)*vecsize);
  float* h_x = (float*)malloc(sizeof(float)*MAT2);
  float* d_A;
  float* d_b;
  float* d_x;
  float* d_l;
  CUDA_CHK(cudaMalloc, (&d_A, sizeof(float)*matsize));
  CUDA_CHK(cudaMalloc, (&d_b, sizeof(float)*vecsize));
  CUDA_CHK(cudaMalloc, (&d_l, sizeof(float)*vecsize));
  CUDA_CHK(cudaMalloc, (&d_x, sizeof(float)*MAT2));

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
      if (i<MAT1) h_x[(i*MAT1)+j]=x179[i*MAT1+j];
      //printf("\n%d:",j);
      //for (unsigned int k=0; k < MAT1; k++){
        //printf("%d,",k);
        //h_A[(i*MAT2)+(j*MAT1)+k]=a(j,k);
        
      //}
    }
  }
  printf("Generated\n");
    //copy values to device
 // CUDA_CHK(cudaMemcpy, (d_A, h_A, sizeof(float)*matsize, cudaMemcpyHostToDevice));
  CUDA_CHK(cudaMemcpy, (d_l, h_l, sizeof(float)*vecsize, cudaMemcpyHostToDevice));
  CUDA_CHK(cudaMemcpy, (d_x, h_x, sizeof(float)*MAT2, cudaMemcpyHostToDevice));
  CUDA_CHK(cudaBindTexture,(NULL,XTG,d_x,sizeof(float)*MAT2));

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
  dim3 blocksPerGrid(ncombos,1,1);
  dim3 threadsPerBlock(MAT1,1,1);
  
  printf("TPB:%d,BPG:%d\n",threadsPerBlock.x,blocksPerGrid.x);
  //Execute
  cudaEvent_t evt_start, evt_stop;
  CUDA_CHK(cudaEventCreate, (&evt_start));
  CUDA_CHK(cudaEventCreate, (&evt_stop));
  CUDA_CHK(cudaEventRecord, (evt_start,0));

  lk_prepare_permutations<<<blocksPerGrid,threadsPerBlock>>>(d_A,d_b,d_x,h_offset);
  cudaDeviceSynchronize();
  Check_Kernel("Generate");
  CUDA_CHK(cudaMemcpy, (h_A,d_A, sizeof(float)*matsize, cudaMemcpyDeviceToHost));
  CUDA_CHK(cudaMemcpy, (h_b,d_b, sizeof(float)*vecsize, cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();
  for (unsigned int i=10000; i<10100; i++){
    printf("\n%d:A|b\n",i);
    //printf("%.3lf|",h_x[i*MAT1]);
    for (unsigned int j=0; j<MAT1; j++){
      for (unsigned int k=0;k<MAT1; k++){
        printf("%g,",h_A[(i*MAT2)+(j*MAT1)+k]);
      }
      printf("|%g\n",h_b[i*MAT1+j]);
    }
  }
  puts("\n");

  printf("Ran Generate\n");

  //CUDA_CHK(cudaFree, (d_x));
  CUDA_CHK(cudaMalloc, (&d_lk, sizeof(float)*outerloop));
  CUDA_CHK(cudaMalloc, (&d_w, sizeof(float)*vecsize));
  CUDA_CHK(cudaMemcpy, (d_w, h_w, sizeof(float)*vecsize, cudaMemcpyHostToDevice));
  dim3 threadsPerBlock_lksolve(256,1,1);
  dim3 blocksPerGrid_lksolve(ncombos/256,1,1);
  solve_permutations<<<blocksPerGrid_lksolve,threadsPerBlock_lksolve>>>(d_A,d_b,h_offset);
  Check_Kernel("Solve");
  printf("Ran Solve\n");
lk_max_permutations<<<blocksPerGrid_lksolve,threadsPerBlock_lksolve>>>(d_b,d_lk, d_l, d_w);
  Check_Kernel("Max");
  CUDA_CHK(cudaMemcpy, (h_lk, d_lk, sizeof(float)*outerloop, cudaMemcpyDeviceToHost));
  CUDA_CHK(cudaEventRecord, (evt_stop, 0));
  CUDA_CHK(cudaEventSynchronize, (evt_stop));
  float total_time;
  CUDA_CHK(cudaEventElapsedTime, (&total_time, evt_start, evt_stop));
  float one_time = total_time * 1e-3;

  float lk_max=FAILVALUE;
  int lk;
  for (unsigned int i=0; i<outerloop; i++){
    //printf("%d:A:LK:%g\n",i,h_lk[i]);
    if (h_lk[i]>lk_max){
      lk=i;
      lk_max=h_lk[i];
      printf("New LKMax:%d,(%g)\n",i,h_lk[i]);
    
      //printf("%.3lf|",h_x[i*MAT1]);
      for (unsigned int j=0; j<MAT1; j++){
        for (unsigned int k=0;k<MAT1; k++){
          printf("%g,",h_A[(i*MAT2)+(j*MAT1)+k]);
        }
      puts("\n");
      }
    } 
  }
  printf("time: %g s\nlkmax:%g@%d\n", one_time,h_lk[lk],lk);

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

