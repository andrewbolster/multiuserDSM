'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy
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
    import pycuda.driver as drv
    import pycuda.gpuarray as garray
    import pycuda.tools as pytools
    from pycuda.compiler import SourceModule
    pointless=False
else:
    pointless=True

kernels=SourceModule("""
//Calc PSD for all channels in parallel (probably fucks memory accesses)
//Do calc_psd in blocks so instead of returning a linear array for one bitload and one channel
//populate the cache with all permutations of all incrememntal bitloads
__device__ void _calc_psd(float **A, float *B, int[] *bitload_k, int *tone, float *channelgap, float *noise, int N, float[][] *XTG){
   int v,x;
   float v_xtalk;
   float P[N];
   for ( v = 0; v < N; v++){
        v_xtalk=XTG[v][v];
        B[v]=(channelgap*(pow(2,bitload[v])-1)*noise)/v_xtalk;
        for ( x = 0; x < N; x++){
            A[v][x]=(-1.0*channelgap*(pow(2,bitload_k[v]-1)*noise*XTG[x][v])/v_xtalk;
        }
        A[v][v]=1;
    
    Crout_LU_Decomposition(&A[0][0],N);
    Crout_LU_Solve(&A[0][0], B, P, N);
    
      
    
    
////////////////////////////////////////////////////////////////////////////////
// File: crout.c                                                              //
// Routines:                                                                  //
//    Crout_LU_Decomposition                                                  //
//    Crout_LU_Solve                                                          //
//                                                                            //
// Required Externally Defined Routines:                                      //
//    Lower_Triangular_Solve                                                  //
//    Unit_Upper_Triangular_Solve                                             //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines 
//int  Lower_Triangular_Solve(double *L, double B[], double x[], int n);
//void Unit_Upper_Triangular_Solve(double *U, double B[], double x[], int n);

////////////////////////////////////////////////////////////////////////////////
//  int Crout_LU_Decomposition(double *A, int n)                              //
//                                                                            //
//  Description:                                                              //
//     This routine uses Crout's method to decompose the n x n matrix A       //
//     into a lower triangular matrix L and a unit upper triangular matrix U  //
//     such that A = LU.                                                      //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Crout's method the diagonal elements of U are 1 and are      //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of L.  (det A = det L * det U = det L).                         //
//     This routine is suitable for those classes of matrices which when      //
//     performing Gaussian elimination do not need to undergo partial         //
//     pivoting, e.g. positive definite symmetric matrices, diagonally        //
//     dominant band matrices, etc.                                           //
//     For the more general case in which partial pivoting is needed use      //
//                    Crout_LU_Decomposition_with_Pivoting.                   //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is Crout_LU_Solve          //
//     (see below).                                                           //
//                                                                            //
//     The Crout method is given by evaluating, in order, the following       //
//     pair of expressions for k = 0, ... , n - 1:                            //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                                 for i = k, ... , n-1,                      //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                                                  / L[k][k] //
//                                      for j = k+1, ... , n-1.               //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A   Pointer to the first element of the matrix A[n][n].        //
//     int     n   The number of rows or columns of the matrix A.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//                                                                            //
//     (your code to intialize the matrix A)                                  //
//                                                                            //
//     err = Crout_LU_Decomposition(&A[0][0], N);                             //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
__global__ int Crout_LU_Decomposition(double *A, int n)
{
   int row, i, j, k, p;
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

__global__ int Lower_Triangular_Solve(double *L, double B[], double x[], int n)
{
   int i, k;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   for (k = 0; k < n; L += n, k++) {
      if (*(L + k) == 0.0) return -1;           // The matrix L is singular
      x[k] = B[k];
      for (i = 0; i < k; i++) x[k] -= x[i] * *(L + i);
      x[k] /= *(L + k);
   }

   return 0;
}

__global__ void Unit_Upper_Triangular_Solve(double *U, double B[], double x[], int n)
{
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
__global__ int Crout_LU_Solve(double *LU, double B[], double x[], int n)
{
//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   if ( Lower_Triangular_Solve(LU, B, x, n) < 0 ) return -1;

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.
//         The diagonal part of the upper triangular part of the matrix is
//         assumed to be 1.0.

   Unit_Upper_Triangular_Solve(LU, x, x, n);
  
   return 0;
}

""")

class GPU(object):
    def update_delta_p(self,tone,N,P):
        for line in range(self.bundle.N):
            
            
            