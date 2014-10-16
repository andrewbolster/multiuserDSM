#include <pycuda-helpers.hpp>
#define MAT1 {{matrixN}}
#define MAT2 MAT1*MAT1
#define FAILVALUE {{failvalue}}
#define FPT {{floatingpointtype}}
#define CHANNELGAP {{channelgap}}
#define NOISE {{noise}}
#define MBPT {{maxbitspertone}}
#define TINY 1.0e-40
#define MAXBITPERM {{maxbitperm}}
#define K {{K}}
#define a(i,j) a[(i)*MAT1+(j)]

//Each thread solves the A[id]x[id]=b[id] problem
__global__ void solve(FPT *A, FPT *B){
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

//Generate A/B for a particular bitload[index]
//Leverages mixed-parallelism using otherline variable
//  i.e can be executed as part of a serial loop
__device__ void generate_AB(FPT *A, FPT *P, FPT *d_XTG, int *bitload, int index, int otherline){
    for (int i=0; i<MAT1; i++){
          //Generate a row of A for this permutation and victim y
          A[index*MAT2+otherline*MAT1+i]=-(CHANNELGAP*((1<<bitload[otherline])-1)*d_XTG[i*MAT1+otherline])/d_XTG[otherline*MAT1+otherline];
      }
      //Generate an item of P
      P[index*MAT1+otherline]=(NOISE*CHANNELGAP*((1<<bitload[otherline])-1))/d_XTG[otherline*MAT1+otherline];    
    
      //Repair an item of A
      A[index*MAT2+otherline*MAT1+otherline]=1;
}

__device__ void lkcalc(int *bitload, FPT *lambdas, FPT *w, FPT *P, FPT *LK, int index){
    FPT lk = 0;
    int broken = 0;
    for (int i=0;i<MAT1;i++){
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

//Given space for A and P, and current_b[N*MAT1] populate P with the psds
//Assume d_XTG is relevant to index
__device__ void d_calc_psd(FPT *A, FPT *P, FPT *d_XTG, int *bitload, int index){
    int i;
    int p_pivot[MAT1], q_pivot[MAT1];
    
    //generate A and B
    for (i=0;i<MAT1;i++){
        generate_AB(A,P,d_XTG,bitload,index,i);
        p_pivot[i]=q_pivot[i]=i;
    }
    __syncthreads();
    d_pivot_decomp(&A[index*MAT2],&p_pivot[0],&q_pivot[0]);
    d_solve(&A[index*MAT2],&P[index*MAT1],&p_pivot[0],&q_pivot[0]);
}

__global__ void calc_psd(FPT *A, FPT *P, FPT *d_XTG, int *current_b, int N){
    //Assume we're doing a full channel range recalculation
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id<N){
        d_calc_psd(A,P,&d_XTG[id*MAT2],&current_b[id*MAT1],id);
    }
}
//===============================================================================
// OSB FUNCTIONS for set channel execution
//===============================================================================
//Generate the A and B for all possible bitloads (in this offset)
//requires grid(MBPT^N,1,1) block(N,1,1)
//thread.y's collaboratively populate A and B for their id
__global__ void lk_osbprepare_permutations(FPT *A, FPT *B, FPT *d_XTG, int offset){
    //Don't need k as its sorted at the host stage for the creation of xtg
    int j=threadIdx.x;
    int myid=blockIdx.x;
    int bitbangval=myid+offset;

    int bitload[MAT1], i;

    //rebase myid to base (MBPT)
    //Unfortunately theres no way around every thread working out its own bitload
    for (i=0; i<MAT1; i++){
        bitload[i]=bitbangval%MBPT;
        bitbangval/=MBPT;
    }
    if (myid+offset<MAXBITPERM){
      generate_AB(A,B,d_XTG,bitload,myid,j);
    }
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
    if (id+offset<MAXBITPERM){
        //do the magic
        d_pivot_decomp(&A[id*MAT2],&p_pivot[0],&q_pivot[0]);
        d_solve(&A[id*MAT2],&B[id*MAT1],&p_pivot[0],&q_pivot[0]);
    }
}

//Finally Calculate the LK_Max_permutations
__global__ void lk_max_permutations(FPT *P, FPT *LK, FPT *lambdas, FPT *w, int offset){
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int bitbangval=id;
    int bitload[MAT1], i;
    
    //At this point, B is populated with the P results.
    for (i=0;i<MAT1;i++){
        bitload[i]=bitbangval%MBPT;
        bitbangval/=MBPT;
    }
    if (id+offset<MAXBITPERM){//check for out of range id's
        lkcalc(bitload,lambdas,w,P,LK,id);
    }
    else
        LK[id]=FAILVALUE;
}

//===============================================================================
// ISB FUNCTIONS for channel range
//===============================================================================

//Do all channels at once for each permutation range
__global__ void isb_optimise_inc(FPT *A, FPT *P, FPT *d_XTG, FPT *LK, FPT *lambdas, FPT *w, int *current_b, int offset){
    const int k=blockIdx.x;                 //The channel that we're playing with
    const int permutation=threadIdx.x;            //the permutation we're generating
    const int index=k*MBPT+permutation;     //The channel-permutation array index we're building
    int line,i;                        //The internal line loop and generic incrementer
    
    __shared__ int bitload[MAT1];
    int threadbit[MAT1];

    //copy initial bitload into block memory
    if (permutation < MAT1){
        //current_b.shape(K,N)
        bitload[permutation]=0;
    }
    __syncthreads();

    // This algorithm swaps the k-range and last=this loops
    for (line=0;line<MAT1;line++){
        //copy shared bitload into thread memory
        for (i=0; i<MAT1; i++){
            threadbit[i]=bitload[i]; //Copy into thread memory
        }
        //For this new user, make him special (line-1 should be optimised)
        threadbit[line]=permutation;    
        
        //Solve!
        d_calc_psd(A,P,&d_XTG[k*MAT2],threadbit,index);
        
        //Calculate lk for each permutation on each channel in parallel
        lkcalc(threadbit,lambdas,w,P,LK,index);
        
        //Return maxlk bitload for this user on this channel (threadbit partially overwritten, no problem
        bitload[line]=isb_perblock_lkmax_B(LK,NULL,k,(int)MBPT);
        
        __syncthreads();
    }

    __syncthreads();
    //For each channel, copy bitload back to current_b
    if (permutation<MAT1){
        current_b[k*MAT1+permutation]=bitload[permutation];
    }
    __syncthreads();
    //At the end of this, current_b will contain the optimal bitloads for all channels, addressible as [k*MAT1+line]
    //    P, addressable as [k*MBPT+bitload] (for the last user)
    //    lk is more or less useless in this case.
}

//Populate the B array with the max bitload permutation for each block of lk's
__device__ int isb_perblock_lkmax_B(FPT *LK, int *B, int id, int blocksize){
    int i, pmax=-1;
    FPT lkmax=FAILVALUE;
    for (i=0;i<blocksize;i++){
        if (lkmax <= LK[id*blocksize+i]){
            pmax=i;                     //This is the maxlk bitload for this line
            lkmax=LK[id*blocksize+i];
        }
    }
    if (B!=NULL){
        B[id]=pmax;
    }
    return pmax;
}

