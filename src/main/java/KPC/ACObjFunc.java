package org.qore.KPC;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;

import java.util.Arrays;

/**
 * Autocorrelation fit objective function usable by Apache Commons
 * multivariate optimizer, representing the closeness of fit between
 * a Trace and a MAP generated to fit it.
 */
public class ACObjFunc implements MultivariateFunction {
    Trace t;
    int J;

    /** Generates new instance of value the value function to be used by optimizer.
     *
     * @param t Trace object MAP aims to fit.
     * @param J Number of composing MAPs.
     */
    public ACObjFunc(Trace t, int J) {
        this.t = t;
        this.J = J;
    }

    /**
     * Computes the objective value of the data.
     *
     * @param x input data to be optimized with length 2*J,
     *          where first J instances are scv, last J are gamma.
     *
     * @return the closeness of the fit for a MAP generated
     * with the given parameters.
     */
    @Override
    public double value(double[] x) {
        double[] scv = Arrays.copyOfRange(x, 0, J);
        double[] gamma = Arrays.copyOfRange(x, J, x.length);

//        DMatrixRMaj lags = new DMatrixRMaj(Arrays.stream(t.acLags).asDoubleStream().toArray());
        DMatrixRMaj lags = new DMatrixRMaj(t.acLags.length,1);
        for (int i = 0; i < t.acLags.length; i++) {
            lags.set(i, t.acLags[i]);
        }
        int nLags = lags.getNumElements();
        DMatrixRMaj acfCoeff = new DMatrixRMaj(nLags);
        double scv_j = scv[0];
        double scv_prev;
        double scalar = 0.5 * (1 - 1 / scv_j);
        CommonOps_DDRM.elementPower(gamma[0], lags, acfCoeff);
        CommonOps_DDRM.scale(scalar, acfCoeff);
        DMatrixRMaj X = new DMatrixRMaj(nLags);
        DMatrixRMaj Xp1 = new DMatrixRMaj(nLags);
        DMatrixRMaj temp = new DMatrixRMaj(nLags);
        DMatrixRMaj result = new DMatrixRMaj(nLags);
        for (int j = 1; j < J; j++) {
            scv_prev = scv_j;
            scv_j = (1 + scv_j) * ((1 + scv[j])) / 2 - 1;
            double r0j = 0.5 * (1 - 1 / scv[j]);
            CommonOps_DDRM.elementPower(gamma[j], lags, X);
            CommonOps_DDRM.scale(scv[j] * r0j, X); // create X
            CommonOps_DDRM.add(X, 1, Xp1); // X + 1
            CommonOps_DDRM.elementMult(acfCoeff, Xp1, temp);
            CommonOps_DDRM.scale(scv_prev, temp);
            CommonOps_DDRM.add(temp, X, result);
            CommonOps_DDRM.divide(result, scv_j, acfCoeff);
        }
        double e0_sq = t.moments[0] * t.moments[0];
        double trace_scv = (t.moments[1] - e0_sq) / e0_sq;
        DMatrixRMaj trace_ac = new DMatrixRMaj(t.ac);
        double acNorm = NormOps_DDRM.normP2(trace_ac);
        DMatrixRMaj ac_difference = new DMatrixRMaj();
        CommonOps_DDRM.subtract(trace_ac, acfCoeff, ac_difference);
        double ac_output = NormOps_DDRM.normP1(ac_difference) / acNorm;
        double scv_output = ((scv_j - trace_scv) * (scv_j - trace_scv)) / (trace_scv * trace_scv);
        return ac_output + scv_output;
    }
}
