package org.qore.KPC;

import org.ejml.data.*;
import org.ejml.dense.row.EigenOps_DDRM;
import org.ejml.dense.row.SingularOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;

import java.util.Arrays;
import java.util.logging.Level;

import static org.apache.commons.math3.util.CombinatoricsUtils.factorialDouble;

/**
 * Class representing a Markov Arrival Process in the {D0, D1} representation.
 */
public class MAP {

    public DMatrixRMaj D0;
    public DMatrixRMaj D1;
    public DMatrixRMaj P;
    public DMatrixRMaj Q;
    public DMatrixRMaj pi;
    public int numStates;

    public MAP() {}

    /**
     * Creates a MAP from the input matrices.
     *
     * @param D0 D0 of the MAP.
     * @param D1 D1 of the MAP.
     */
    public MAP(DMatrixRMaj D0, DMatrixRMaj D1) {
        checkSquareAndEqual(D0, D1);
        numStates = D0.getNumRows();
        this.D0 = D0;
        this.D1 = D1;
        P = new DMatrixRMaj(numStates, numStates);
        Q = new DMatrixRMaj(numStates, numStates);
        construct();
    }

    /**
     * Generates P, Q, and pi vectors for a given D0, D1
     */
    protected void construct() {
        DMatrixRMaj d0_inv = MatrixOps.negateInvert(D0);
        CommonOps_DDRM.mult(d0_inv, D1, P);
        CommonOps_DDRM.add(D0, D1, Q);
        pi = dtmc();
    }



    /**
     * Size checking helper method for MAP and SMP constructors.
     *
     * @param a First matrix.
     * @param b Second matrix.
     */
    protected void checkSquareAndEqual(DMatrixRMaj a, DMatrixRMaj b) {
        if ((a.getNumRows() != b.getNumRows()) || (a.getNumCols() != b.getNumCols()) ||
                (a.getNumRows() != a.getNumCols())) {
            throw new IllegalArgumentException("D0 and D1 Matrices are invalid shapes");
        }
    }

    /**
     * Checks if a MAP is feasible according to relations determined in
     * Casale et. al. 2007.
     *
     * @return if MAP is feasible or not.
     */
    public boolean isFeasible() {
        return isFeasible(Constants.FEASIBLE_TOL);
    }

    /**
     * Checks if a MAP is feasible within specified tolerance
     * according to relations determined in Casale et. al. 2007.
     *
     * @return if MAP is feasible or not.
     */
    public boolean isFeasible(double tol) {
        return feasiblePi() && feasibleD0() && feasibleD1() && stochasticP(tol) && stochasticQ(tol) &&
                irreducibleP(tol) && irreducibleQ(tol);
    }

    /**
     * Checks if Pi vector is feasible.
     */
    protected boolean feasiblePi() {
        return CommonOps_DDRM.elementSum(pi) != 0 && !Double.isNaN(CommonOps_DDRM.elementSum(pi));
    }

    /**
     * Checks if diagonal elements of D0 are negative and off-diagonals are positive.
     */
    protected boolean feasibleD0() {
        boolean feasible = true;
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (i == j && D0.get(i, j) > 0) {
                    Constants.LOGGER.info(String.format("Element (%d, %d) of D0 is not negative", i, j));
                    feasible = false;
                }
                if (i != j && D0.get(i, j) < 0) {
                    Constants.LOGGER.info(String.format("Element (%d, %d) of D0 is negative", i, j));
                    feasible = false;
                }
            }
        }
        return feasible;
    }

    /**
     * Checks if all elements of D1 are postive.
     */
    protected boolean feasibleD1() {
        boolean feasible = true;
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (D1.get(i, j) < 0) {
                    Constants.LOGGER.info(String.format("Element (%d, %d) of D1 is negative", i, j));
                    feasible = false;
                }
            }
        }
        return feasible;
    }

    /**
     * Checks if Q is stochastic.
     */
    protected boolean stochasticQ(double tol) {
        boolean feasible = true;
        DMatrixRMaj sums = new DMatrixRMaj(Q.getNumRows());
        CommonOps_DDRM.sumRows(Q, sums);
        for (int i = 0; i < sums.getNumElements(); i++) {
            if (Math.abs(sums.get(i)) > tol) {
                Constants.LOGGER.info(String.format("Row %d of Q is not stochastic", i));
                feasible = false;
            }
        }
        return feasible;
    }

    /**
     * Checks if P is stochastic.
     */
    protected boolean stochasticP(double tol) {
        boolean feasible = true;
        DMatrixRMaj sums = new DMatrixRMaj(P.getNumRows());
        CommonOps_DDRM.sumRows(P, sums);
        for (int i = 0; i < sums.getNumElements(); i++) {
            if (Double.isNaN(sums.get(i)) || Math.abs(sums.get(i) - 1) > tol) {
                Constants.LOGGER.info(String.format("Row %d of P is not stochastic", i));
                feasible = false;
            }
        }
        return feasible;
    }

    /**
     * Checks if Q is irreducible.
     */
    protected boolean irreducibleP(double tol) {
        boolean feasible = true;
        if (numStates <= Constants.MAX_MAP_SIZE) {
            EigenDecomposition_F64<DMatrixRMaj> eigenP = DecompositionFactory_DDRM.eig(numStates, false);
            DMatrixRMaj evals = new DMatrixRMaj(numStates, 1);
            BMatrixRMaj isOne = new BMatrixRMaj(numStates, 1);
            eigenP.decompose(P);
            CommonOps_DDRM.extractDiag(EigenOps_DDRM.createMatrixD(eigenP), evals);
            CommonOps_DDRM.elementMoreThan(evals, 1 - tol, isOne);
            if (CommonOps_DDRM.countTrue(isOne) > 1) {
                Constants.LOGGER.info("P is not irreducible");
                feasible = false;
            }
        } else {
            Constants.LOGGER.info("Number of states is too high");
            return false;
        }
        return feasible;
    }

    /**
     * Checks if P is irreducible.
     */
    protected boolean irreducibleQ(double tol) {
        boolean feasible = true;
        EigenDecomposition_F64<DMatrixRMaj> eigenQ = DecompositionFactory_DDRM.eig(numStates, false);
        DMatrixRMaj evals = new DMatrixRMaj(numStates, 1);
        BMatrixRMaj isZero = new BMatrixRMaj(numStates, 1);
        eigenQ.decompose(Q);
        CommonOps_DDRM.extractDiag(EigenOps_DDRM.createMatrixD(eigenQ), evals);
        CommonOps_DDRM.elementMoreThan(evals, 0 + tol, isZero);
        if (CommonOps_DDRM.countTrue(isZero) > 1) {
            EigenOps_DDRM.createMatrixD(eigenQ);
            Constants.LOGGER.info("Q is not irreducible");
            feasible = false;
        }
        return feasible;
    }

    /**
     * Calculates the first 3 moments of the MAP according to Casale et. al. 2007.
     *
     * @return First three moments of the MAP
     */
    public double[] getMoments() {
        return getMoments(new int[] {1,2,3});
    }

    /**
     * Calculates the moments of the MAP according to Casale et. al. 2007.
     *
     * @param orders Order of moments desired (e.g. E[X] and E[X^3] => [1, 3])
     * @return specified moments of the MAP.
     */
    public double[] getMoments(int[] orders) {
        DMatrixRMaj e = MatrixOps.getOnes(numStates);
        double[] moments = new double[orders.length];

        int i = 0;
        // TODO: error check MAP is not null
        for (int order : orders) {
            DMatrixRMaj d0_inv = MatrixOps.negateInvert(D0);
            DMatrixRMaj d0_k = MatrixOps.matrixPower(d0_inv, order);
            double fact = factorialDouble(order);
            DMatrixRMaj pi_d0k = new DMatrixRMaj();
            DMatrixRMaj result = new DMatrixRMaj();
            CommonOps_DDRM.mult(pi, d0_k, pi_d0k);
            CommonOps_DDRM.mult(fact, pi_d0k, e, result);
            moments[i++] = result.get(0);
        }
        return moments;
    }

    /**
     * Calculates the autocorrelations for the list of lags specified
     * according to Casale et. al. 2007.
     *
     * @param lags List of lags of which we want the AC values.
     * @return autocorrelation coefficients of the input lags.
     */
    public double[] getAcf(int[] lags) {
        double[] acs = new double[lags.length];
        DMatrixRMaj d0_inv = MatrixOps.negateInvert(D0);
        DMatrixRMaj e = MatrixOps.getOnes(numStates);
        int i = 0;
        for (int lag : lags) {
            DMatrixRMaj t1 = new DMatrixRMaj();
            DMatrixRMaj t2 = new DMatrixRMaj();
            DMatrixRMaj t3 = new DMatrixRMaj();
            DMatrixRMaj corr_k = new DMatrixRMaj();
            DMatrixRMaj p_k = MatrixOps.matrixPower(P, lag);
            CommonOps_DDRM.mult(pi, d0_inv, t1);
            CommonOps_DDRM.mult(t1, p_k, t2);
            CommonOps_DDRM.mult(t2, d0_inv, t3);
            CommonOps_DDRM.mult(t3, e, corr_k);
            double mean_squared = getMean()*getMean();
            acs[i] = (corr_k.get(0) - mean_squared) / (getMoments(new int[]{2})[0]- mean_squared);
            i++;
        }
        return acs;
    }

    // input: a_k =
    //        i_k =

    public double getGamma() {
        double[] acs = getAcf(new int[] {1,2});
        return acs[1] / acs[0];
    }

    public double[] getBCs(int[][] lags) {
        double[] bc = new double[lags.length];
        int i =0;
        int[] powers = {1,1,1};
        for (int[] joint : lags) {
            bc[i] = getJoint(joint, powers);
            i++;
        }
        return bc;
    }

    /**
     * Calculates the joint characteristics for the list of
     * lags and powers according to Casale et. al. 2007.
     *
     * @param a_k Ordered, increasing indices of ACs.
     * @param i_k Powers corresponding to the indices of AC.
     * @return Joints moments of the MAP corresponding to input.
     */
    public double getJoint(int[] a_k, int[] i_k) {
        a_k = Equations.cumSum(a_k);
        DMatrixRMaj d0_inv = MatrixOps.negateInvert(D0);
        if (!jointVectorsLegal(a_k, i_k)) {
            throw new IllegalArgumentException("Joint Moment vectors do not match");
        }

        DMatrixRMaj previous = pi;
        DMatrixRMaj intermediate = new DMatrixRMaj(1, numStates);
        DMatrixRMaj d_p = new DMatrixRMaj();
        DMatrixRMaj d_kl;
        DMatrixRMaj p_kl;
        for (int l = 0; l < i_k.length; l++) {
            d_kl = MatrixOps.matrixPower(d0_inv, i_k[l]); // -D0^-k
            if (l < i_k.length - 1) {
                p_kl = MatrixOps.matrixPower(P, a_k[l+1]- a_k[l]); // first ones are as in eqn
            } else {
                p_kl = new DMatrixRMaj(numStates, 1); // col vector for last P^k
                CommonOps_DDRM.fill(p_kl, 1);
            }
            CommonOps_DDRM.mult(factorialDouble(i_k[l]), d_kl, p_kl, d_p);
            CommonOps_DDRM.mult(previous, d_p, intermediate);
            previous = intermediate.copy();
        }
        return intermediate.get(0);
    }

    /**
     * Helper method determines if the inputs to the joint method are legal.
     */
    private boolean jointVectorsLegal(int[] a_k, int[] i_k) {
        if (a_k.length != i_k.length) {
            Constants.LOGGER.log(Level.SEVERE, "Vector lengths do not match");
            return false;
        }

        for (int j = 1; j < a_k.length; j++) {
            if (a_k[j] <= a_k[j - 1]) {
                Constants.LOGGER.severe("vector of 'a's not sorted");
                return false;
            }
            if (i_k[j] <= 0) {
                Constants.LOGGER.severe("vector of 'i's must be positive");
                return false;
            }
        }
        return true;
    }

    /**
     * Solves the equilibrium distribution of the Markov chain.
     *
     * @param continuous Solve continuous equilibrium distribtution.
     * @return pi vector s.t. pi*P = pi (if false) or Q*pi = 0 (if true).
     */
    public DMatrixRMaj ctmc(boolean continuous) {
        DMatrixRMaj pi = new DMatrixRMaj(numStates);
        DMatrixRMaj aT = new DMatrixRMaj();
        DMatrixRMaj A;

        if (continuous) {
            A = new DMatrixRMaj();
            CommonOps_DDRM.subtract(P, CommonOps_DDRM.identity(P.getNumRows()), A);
        } else {
            A = Q;
        }

        // TODO: handle case where A is reducible?
        // assuming A irreducible, and therefore rank == numStates
        CommonOps_DDRM.transpose(A, aT);
        SingularValueDecomposition_F64<DMatrixRMaj> svd =
                DecompositionFactory_DDRM.svd(A.getNumRows(), A.getNumCols(),true, true, false);
        svd.decompose(aT);
        try {
            SingularOps_DDRM.nullVector(svd, true, pi);
        } catch (Exception e) {
  //          Constants.LOGGER.warning("Cannot calculate pi vector");
        }
        CommonOps_DDRM.transpose(pi);
        CommonOps_DDRM.divide(pi, CommonOps_DDRM.elementSum(pi));
        return pi;
    }

    /**
     * Solves the equilibrium distribution of the Continuous Markov chain.
     *
     * @return pi vector s.t. Q*pi = 0.
     */
    public DMatrixRMaj dtmc() {
        return ctmc(true);
    }

    /**
     * @return SCV of the MAP.
     */
    public double getSCV() {
        return getVariance() / (getMean() * getMean());
    }

    /**
     * @return Variance of the MAP.
     */
    public double getVariance() {
        int[] moments = {2};
        return getMoments(moments)[0] - (getMean() * getMean());
    }

    /**
     * @return Mean of the MAP.
     */
    public double getMean() {
        return 1 / getLambda();
    }

    /**
     * @return Mean arrival rate of the MAP.
     */
    public double getLambda() {
        DMatrixRMaj e = MatrixOps.getOnes(numStates);
        DMatrixRMaj x = ctmc(false); // xQ = 0
        DMatrixRMaj result = new DMatrixRMaj();
        DMatrixRMaj x_d1 = new DMatrixRMaj();
        CommonOps_DDRM.mult(x, D1, x_d1);
        CommonOps_DDRM.mult(x_d1, e, result);
        return result.get(0); // should be only value
    }

    /**
     * Scales MAP to specified mean.
     *
     * @param newMean New mean of MAP.
     * @param allowSM Allow SMP (infeasible D1).
     */
    public void scale(double newMean, boolean allowSM) {
        double ratio = getMean() / newMean;
        CommonOps_DDRM.scale(ratio, D0);
        CommonOps_DDRM.scale(ratio, D1);
        if (!allowSM) {
            normalize();
        }

    }

    /**
     * Normalizes the MAP so that it is feasible.
     */
    public void normalize() {
        // set all neg entries to zero
        for (int i = 0; i < D0.getNumElements(); i++) {
            D0.set(i, Math.max(D0.get(i), 0));
            D1.set(i, Math.max(D1.get(i), 0));
        }
        for (int i = 0; i < D0.getNumRows(); i++) {
            D0.set(i,i,0);
        }
        DMatrixRMaj d0RowSums = new DMatrixRMaj(D0.getNumRows());
        DMatrixRMaj d1RowSums = new DMatrixRMaj(D0.getNumRows());
        CommonOps_DDRM.sumRows(D0, d0RowSums);
        CommonOps_DDRM.sumRows(D1, d1RowSums);
        CommonOps_DDRM.addEquals(d0RowSums, d1RowSums);
        for (int i = 0; i < D0.getNumRows(); i++) {
            D0.set(i,i, -1*d0RowSums.get(i));
        }
        construct();
    }

    /**
     * Prints the MAP to the screen. Useful for debugging.
     */
    public void print() {
        System.out.println("D0: " + D0.toString());
        System.out.println("D1: " + D1.toString());
    }


}