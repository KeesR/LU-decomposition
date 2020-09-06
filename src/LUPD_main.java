
public class LUPD_main {

    public static void main(String[] args) {
        test_inversion();
    }
    
    static void test_inversion() {
        double[][] A = { { 2, 1 }, { -1, 2} };
        int[] P = new int[3];

        printMatrix("Matrix", A);

        LUPDecomposition.LUPDecompose( A, 2, 1e-8,  P);

        double[][] IA = new double[2][2];
        LUPDecomposition.LUPInvert(A, P, 2, IA);

        printMatrix("Inverse", IA);

        double det = LUPDecomposition.LUPDeterminant(A, P, 2);
        System.out.println("Determinant" + det);

        double[] b = { 4, 0 };
        printVector("b vector", b);

        double[] x = { 0, 0 };
        LUPDecomposition.LUPSolve(A, P, b, 2, x);
        printVector("Solution", x);
      
    }

    static void printMatrix(String message, double[][] m) {
        System.out.println(message);
        for (int r = 0; r < m.length; r++) {
            for (int c = 0; c < m[0].length; c++) {
                System.out.print("  " + m[r][c]);
            }
            System.out.println();
        }
    }

    static void printVector(String message, double[] v) {
        System.out.println(message);
        for (int c = 0; c < v.length; c++) {
            System.out.println("  " + v[c]);
        }
    }
}
