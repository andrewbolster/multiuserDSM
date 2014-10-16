#define MAT1 4
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
  //solves Uy=z
  xtmp[MAT1-1]/=a(MAT1-1,MAT1-1);
  for (i=MAT1-2;i>=0;i--){
    ftmp=xtmp[i];
    for (j=i+1;j<MAT1;j++){
      ftmp-=a(i,j)*xtmp[j];
    }
    xtmp[i]=(ftmp)/a(i,i);
  }
  //solves x=Qy
  for (i=0;i<MAT1;i++){
    x[i]=xtmp[q[i]];
  }
}

__global__ void solve(float *A, float *B, int max){
  //Each thread solves the A[id]x[id]=b[id] problem
  int id= blockDim.x*blockIdx.x + threadIdx.x;
  int p_pivot[MAT1],q_pivot[MAT1];
  if ((GO==1) && (id < max)){
    for (int i=0;i<MAT1;i++) {
      p_pivot[i]=q_pivot[i]=i;
    }

    d_pivot_decomp(&A[id*MAT2],&p_pivot[0],&q_pivot[0]);
    d_solve(&A[id*MAT2],&B[id*MAT1],&p_pivot[0],&q_pivot[0]);
  }
}
