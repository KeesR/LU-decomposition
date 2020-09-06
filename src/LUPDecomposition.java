
public class LUPDecomposition {

    // INPUT:  A - array of pointers to rows of a square matrix having
    //         dimension N
    //         Tol - small tolerance number to detect failure when the
    //         matrix is near degenerate
    // OUTPUT: Matrix A is changed, it contains a copy of both matrices
    //         L-E and U as A=(L-E)+U such that P*A=L*U. The permutation 
    //         matrix is not stored as a matrix, but in an integer vector 
    //         P of size N+1 containing column indexes where the
    //         permutation matrix has "1". The last element P[N]=S+N,
    //         where S is the number of row exchanges needed for
    //         determinant computation, det(P)=(-1)^S
  
    static boolean LUPDecompose(double[][] A, int N, double Tol, int[] P) {
      
        //  Check preconditions
        if (A.length != N || A[0].length != N) {
          throw new IllegalArgumentException("Matrix A must be square");
        }
        if (P.length != N + 1) {
          throw new IllegalArgumentException("P must be a N+1 sixed int array");
        }
        
        for (int i = 0; i <= N; i++)
            P[i] = i; //  Unit permutation matrix, P[N] initialized with N
    
        for (int i = 0; i < N; i++) {
            double maxA = 0.0;
            int imax = i;
    
          for (int k = i; k < N; k++) {
              double absA = Math.abs(A[k][i]);
              if (absA > maxA) { 
                  maxA = absA;
                  imax = k;
              }
          }
      
          if (maxA < Tol) return false; //failure, matrix is degenerate
          
          if (imax != i) {
              //  pivoting P
              int j = P[i];
              P[i] = P[imax];
              P[imax] = j;
  
              //  pivoting rows of A
              double[] ptr = A[i];
              A[i] = A[imax];
              A[imax] = ptr;
  
              //  counting pivots starting from N (for determinant)
              P[N]++;
          }
      
          for (int j = i + 1; j < N; j++) {
              A[j][i] /= A[i][i];
              for (int k = i + 1; k < N; k++)
                  A[j][k] -= A[j][i] * A[i][k];
          }
      }
    
        return true;  //  decomposition done 
    }
  
    //  INPUT:  A,P filled in LUPDecompose; b - rhs vector; N - dimension
    //  OUTPUT: x - solution vector of A*x=b

    static void LUPSolve(double[][] A, int[] P, double[] b, int N, double[] x) {
    
        //  Check preconditions
        if (A.length != N || A[0].length != N) {
          throw new IllegalArgumentException("Matrix A must be square");
        }
        if (P.length != N + 1) {
          throw new IllegalArgumentException("P must be a N+1 sixed integer array");
        }
        if (x.length != N ) {
          throw new IllegalArgumentException("x must be a N sixed double array");
        }

        for (int i = 0; i < N; i++) {
            x[i] = b[P[i]];
    
            for (int k = 0; k < i; k++)
                x[i] -= A[i][k] * x[k];
        }
    
        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                x[i] -= A[i][k] * x[k];
    
            x[i] = x[i] / A[i][i];
        }
    }
    
    //  INPUT:  A,P filled in LUPDecompose; N - dimension
    //  OUTPUT: IA is the inverse of the initial matrix

    static void LUPInvert(double[][] A, int[] P, int N, double[][] IA) {
      
        //  Check preconditions
        if (A.length != N || A[0].length != N) {
          throw new IllegalArgumentException("Matrix A must be square");
        }
        if (P.length != N + 1) {
          throw new IllegalArgumentException("P must be a N+1 sixed integer array");
        }
        if (IA.length != A.length || IA[0].length != A[0].length) {
          throw new IllegalArgumentException("Matrix IA must be same size as matrix A");
        }
  
        //  create unity matrix, but with the rows pivotted acc. to P
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                if (P[i] == j) 
                    IA[i][j] = 1.0;
                else
                    IA[i][j] = 0.0;
    
                for (int k = 0; k < i; k++)
                    IA[i][j] -= A[i][k] * IA[k][j];
            }
            
            for (int i = N - 1; i >= 0; i--) {
                for (int k = i + 1; k < N; k++) {
                    IA[i][j] -= A[i][k] * IA[k][j];
                }
                IA[i][j] = IA[i][j] / A[i][i];
            }
        }
    }
    
    //  INPUT:  A,P filled in LUPDecompose; N - dimension. 
    //  OUTPUT: Function returns the determinant of the initial matrix
  
    static double LUPDeterminant(double[][]A, int[] P, int N) {
    
        //  Check preconditions
        if (A.length != N || A[0].length != N) {
          throw new IllegalArgumentException("Matrix A must be square");
        }
        if (P.length != N + 1) {
          throw new IllegalArgumentException("P must be a N+1 sixed integer array");
        }

        double det = A[0][0];
    
        for (int i = 1; i < N; i++)
            det *= A[i][i];
    
        if ((P[N] - N) % 2 == 0)
            return det; 
        else
            return -det;
    }
}

