'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy as np
import scipy.misc
import multiprocessing

def have_pycuda():
    try:
        import pycuda
        return True
    except:
        return False

if have_pycuda():
    from pycuda.tools import make_default_context
    import pycuda.autoinit
    import pycuda.cumath
    import pycuda.driver as drv
    import pycuda.gpuarray as garray
    import pycuda.tools as pytools
    from pycuda.compiler import SourceModule
    pointless=False
else:
    pointless=True

kernels=SourceModule("""
__device__ int Crout_LU_Decomposition (double *A, int n)
{
   int i, j, k, p;
   double *p_k, *p_row, *p_col;

//         For each row and column, k = 0, ..., n-1,
//            find the lower triangular matrix elements for column k
//            and if the matrix is non-singular (nonzero diagonal element).
//            find the upper triangular matrix elements for row k. 
 
   for (k = 0, p_k = A; k < n; p_k += n, k++) {
      for (i = k, p_row = p_k; i < n; p_row += n, i++) {
         for (p = 0, p_col = A; p < k; p_col += n, p++)
            *(p_row + k) -= *(p_row + p) * *(p_col + k);
      }  
      if ( *(p_k + k) == 0.0 ) return -1;
      for (j = k+1; j < n; j++) {
         for (p = 0, p_col = A; p < k; p_col += n,  p++)
            *(p_k + j) -= *(p_k + p) * *(p_col + j);
         *(p_k + j) /= *(p_k + k);
      }
   }
   return 0;
}

__device__ void Lower_Triangular_Solve(double *L, double B[], double x[], int n){
   int i, k;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   for (k = 0; k < n; L += n, k++) {
      if (*(L + k) == 0.0) return;           // The matrix L is singular
      x[k] = B[k];
      for (i = 0; i < k; i++) x[k] -= x[i] * *(L + i);
      x[k] /= *(L + k);
   }
}

__device__ void Unit_Upper_Triangular_Solve (double *U, double B[], double x[], int n){
   int i, k;

//         Solve the linear equation Ux = B for x, where U is an upper
//         triangular matrix.                                      
   x[n-1] = B[n-1]; 
   for (k = n-2, U += n * (n - 2); k >= 0; U -= n, k--) {
      x[k] = B[k];
      for (i = k + 1; i < n; i++) x[k] -= x[i] * *(U + i);
   }
}
                                                                           //
__device__ void Crout_LU_Solve (double *LU, double B[], double x[], int n){
//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   Lower_Triangular_Solve(LU, B, x, n);

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.
//         The diagonal part of the upper triangular part of the matrix is
//         assumed to be 1.0.

   Unit_Upper_Triangular_Solve(LU, x, x, n);
}

//Calc PSD for all channels in parallel (probably fucks memory accesses)
//Do calc_psd in blocks so instead of returning a linear array for one bitload and one channel
//populate the cache with all permutations of all incrememntal bitloads
__device__ void calc_psd (double **A, double *B, double *P, int *bitload_k, int *tone, 
                          float channelgap, float noise, 
                          int N, float **XTG)
                          {
   int v,x;
   float v_xtalk;
   for ( v = 0; v < N; v++ ) {
        v_xtalk=XTG[v][v];
        B[v]=((channelgap*(__powf(2,bitload_k[v])-1)*noise)/v_xtalk);
        for ( x = 0; x < N; x++ ) {
            A[v][x]=(-1.0*channelgap*(__powf(2,bitload_k[v]-1)*noise*XTG[x][v]))/v_xtalk;
        }
        A[v][v]=1;
    
    Crout_LU_Decomposition(&A[0][0],N);
    Crout_LU_Solve(&A[0][0], B, P, N);    
    }
} 

""")

generateA=SourceModule(
"""

    
"""
)
class GPU(object):
    def update_delta_p(self,tone,N,P,b):
        _b=np.tile(util.mat2arr(b),N)
        delta_p=np.zeros(N,N)
        #create a square matrix of incremented bitloads
        for line in xrange(self.bundle.N):
            _b[line,line]+=1
            

if __name__ == "__main__":
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
        
    calc_psd=kernels.get_function("calc_psd")
    calc_psd(A,B,P,b[K/2],K/2,cg,0,N,XTG)
        
            
            
            