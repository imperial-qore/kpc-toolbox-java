package org.qore.KPC;

import de.xypron.jcobyla.Calcfc;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;

import java.util.Arrays;

/**
 * BC Objective function class used by JCOBYLA nonlinear
 * constrained optimizer.
 */
public class BCObjFunc implements Calcfc {
    Trace t;
    int J;
    boolean verbose;
    double[] scv;
    double[] gamma;
    boolean allowSM;

    /**
     * Generates BC Objective Function instance with Trace member and existing
     * SCV and gamma values of MAP assumed to give good AC fit.
     *
     * @param t Trace object MAP aims to fit.
     * @param J Number of composing MAPs.
     * @param resSCV SCV values from a successful iteration of AC fit.
     * @param resGamma Gamma values from a successful iteration of AC fit.
     * @param sm Does MAP generator allow semi-markov fits as well.
     */
    public BCObjFunc(Trace t, int J, double[] resSCV, double[] resGamma, boolean sm, boolean verbose) {
        this.t = t;
        this.J = J;
        this.verbose = verbose;
        scv = resSCV;
        gamma = resGamma;
        allowSM = sm;
    }

    /**
     * The objective and constraints function evaluation method used in COBYLA2 minimization.
     *
     * @param n Number of variables. 2 for each MAP.
     * @param m Number of constraints. 3 for each MAP.
     * @param x Variable values to be employed in function and constraints calculation.
     * @param con Calculated function values of the constraints. Enforced as non-negative.
     * @return Calculated objective function value (closeness of BC fit).
     */
    @Override
    public double compute(int n, int m, double[] x, double[] con) {
        double[] e1j = Arrays.copyOfRange(x, 0, J); // e1 of the j composing maps
        double[] e3j = Arrays.copyOfRange(x, J, x.length); // e3 of the j composing maps
        double[] e2j = new double[J];
        for (int i = 0; i < J; i++) {
            e2j[i] = (1 + scv[i])*e1j[i]*e1j[i];
        }
        int index = 0;

        // if SCV for the first arbitrary MAP is > 1, then we constrain on that MAP too
        int start = scv[0] > 1 ? 0 : 1;
        for (int j = start; j < J; j++) { // start at 2, first MAP is arbitrary
            con[index++] = e2j[j] - (2 + Constants.CONSTRAINT_TOL) * e1j[j] * e1j[j];            // E2j(j) > 2*E1j(j)^2
            con[index++] = e3j[j] - (1.5 + Constants.CONSTRAINT_TOL) * e2j[j] * e2j[j] / e1j[j]; // E3j(j) > 3*E2j(j)^2/(2*E1j(j))
        }
        double temp = 1;
        for (double v : e3j) {
            temp *= v;
        }
        double denom = Math.pow(6, J-1);
        temp /= denom;
        con[index++] = 2 - temp / t.moments[2];
        con[index++] = temp / t.moments[2] - 0.5;

        for (int j = 1; j < J; j++) {
            double t = e2j[j] - 2*e1j[j]*e1j[j];
            double a = (1.0/3)*e1j[j] / t * e3j[j];
            double b = 0.5*e2j[j]*e2j[j] / t;
            con[index++] = 1.0e-16 - (a - b);
            con[index++] = a + b - 1.0e-16 ;
            // 1/3*E1j(j)/(E2j(j)-2*E1j(j)^2)*E3j(j)-1/2*E2j(j)^2/(E2j(j)-2*E1j(j)^2) - 1/(1e+16);
            // 1e-16 - 1/3*E1j(j)/(E2j(j)-2*E1j(j)^2)*E3j(j)-1/2*E2j(j)^2/(E2j(j)-2*E1j(j)^2);
        }
        return value(e1j, e3j);
    }

    /** Helper function for compute that handles the evaluation of the obj function
     * as its closeness of fit to the Trace member.
     *
     * @param e1 First moments of MAP(2)s.
     * @param e3 Third moments of MAP(2)s.
     */
    double value(double[] e1, double[] e3) {
        MAP composedMAP;
        try {
            if (allowSM) {
                composedMAP = KPC.composeSMMAP(e1, e3, scv, gamma, J, verbose);
            } else {
                composedMAP = KPC.composeMAP(e1, e3, scv, gamma, J, verbose);
            }
        } catch (InfeasibleMAPException e) {
            Constants.LOGGER.warning(e.getMessage());
            Constants.LOGGER.warning("E1: " + Arrays.toString(e1) + "; E3: " + Arrays.toString(e3));
            return Double.MAX_VALUE;
        }
        composedMAP.scale(t.moments[0], allowSM);
        DMatrixRMaj trace_bc = new DMatrixRMaj(t.bc);
        double nbc = NormOps_DDRM.normP2(trace_bc);
        DMatrixRMaj map_bc = new DMatrixRMaj(t.bc.length, 1);
        for (int i = 0; i < t.bc.length; i++) {
            double joint = composedMAP.getJoint(t.bcLags[i], new int[]{1,1,1});
            map_bc.set(i,0, joint);
        }
        DMatrixRMaj diff = new DMatrixRMaj(t.bc.length);
        CommonOps_DDRM.subtract(trace_bc, map_bc, diff);
        return NormOps_DDRM.normP2(diff) / nbc;
    }
}
