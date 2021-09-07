package KPC;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

/**
 * Helper static class that performs simple matrix operations useful in MAP fitting.
 */
public class MatrixOps {

    /**
     * I understand this is horrible, but it is far too late to redo the backbone of the project,
     * and EJML does not support fast enough matrix powers. This method isn't altogether
     * too crucial, as it is only used for MAP joints. If you come across this and
     * find it inadequate, please fix. JNF maybe? :)
     * @param input Matrix to be raised to a power.
     * @param exp number of times we multiply the matrix by itself
     * @return new resultant Matrix.
     */
    protected static DMatrixRMaj matrixPower(DMatrixRMaj input, int exp) {

        int size = input.getNumRows();
        double[][] newAr = new double[size][size];
        for (int i = 0; i < size; i++) {
            DMatrixRMaj temp = new DMatrixRMaj(size);
            CommonOps_DDRM.extractRow(input, i, temp);
            newAr[i] = temp.data;
        }
        Array2DRowRealMatrix mat = new Array2DRowRealMatrix(newAr);
        RealMatrix result = mat.power(exp);
        return new DMatrixRMaj(result.getData());
    }

    /**
     * Negates the input matrix and then computes the inverse.
     * @param input Matrix to be manipulated.
     * @return New matrix with desired changes.
     */
    protected static DMatrixRMaj negateInvert(DMatrixRMaj input) {
        DMatrixRMaj A = new DMatrixRMaj();
        CommonOps_DDRM.changeSign(input, A);
        CommonOps_DDRM.invert(A);
        return A;
    }

    /**
     * Generates vector of correct size with all 1s.
     * @param n Number of entries in output.
     * @return vector of desired size of all 1s.
     */
    protected static DMatrixRMaj getOnes(int n) {
        DMatrixRMaj e = new DMatrixRMaj(n, 1);
        CommonOps_DDRM.fill(e, 1);
        return e;
    }
}
