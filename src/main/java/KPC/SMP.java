package KPC;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

/**
 * Represents a Semi Markov Process, initialized by D0, and P.
 * Note: may still be a feasible MAP.
 */
public class SMP extends MAP {
    /**
     *
     * @param d0 D0 Matrix (same as MAP).
     * @param p Transition probability matrix.
     */
    public SMP(DMatrixRMaj d0, DMatrixRMaj p) {
        super();
        this.D0 = d0;
        this.P = p;
        this.numStates = d0.getNumRows();
        D1 = new DMatrixRMaj(numStates, numStates);
        Q = new DMatrixRMaj(numStates, numStates);
        DMatrixRMaj negd0 = new DMatrixRMaj(numStates, numStates);
        CommonOps_DDRM.changeSign(d0, negd0);
        CommonOps_DDRM.mult(negd0, p, D1);
        CommonOps_DDRM.add(d0, D1, Q);
        pi = dtmc();
    }

    /**
     * Checks if the SMP is valid (feasible without D1 constraint).
     * @return whether SMP is valid or not.
     */
    public boolean isValid() {
        return isValid(Constants.FEASIBLE_TOL);
    }

    /**
     * Checks if the SMP is valid (feasible without D1 constraint).
     * @return whether SMP is valid or not.
     */
    public boolean isValid(double tol) {
        return feasiblePi() && stochasticP(tol) && irreducibleP(tol);
    }
}
