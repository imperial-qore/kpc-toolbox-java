package org.qore.KPC;

import de.xypron.jcobyla.Calcfc;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;

import java.util.Arrays;

@Deprecated
public class OldACObjFunc implements Calcfc {
    Trace t;
    int J;

    public OldACObjFunc(Trace t, int J) {
        this.t = t;
        this.J = J;
    }


    /**
     * The objective and constraints function evaluation method used in COBYLA2 minimization.
     * @param n Number of variables. 2 for each MAP
     * @param m Number of constraints. 3 for each MAP
     * @param x Variable values to be employed in function and constraints calculation.
     * @param con Calculated function values of the constraints. Enforced as non-negative
     * @return Calculated objective function value.
     */
    @Override
    public double compute(int n, int m, double[] x, double[] con) {
        //System.out.println("data: " + Arrays.toString(x));
        double[] scv = Arrays.copyOfRange(x, 0, J); // scv of the j composing maps
        double[] gamma = Arrays.copyOfRange(x, J, x.length); // gamma of the j composing maps
        con[0] = scv[0] - (0.5 - Constants.CONSTRAINT_TOL); // scv(0) > 0.5
        int index = 1;
        for (int j = 1; j < J; j++) {
            con[index++] = scv[j] - (1 + Constants.CONSTRAINT_TOL); // scv > 1
        }
        for (int i = 0; i < J; i++) {
            con[index++] = (1 - Constants.CONSTRAINT_TOL) - gamma[i]; // gamma < 1
            con[index++] = gamma[i] - Constants.CONSTRAINT_TOL; // gamma > 0
        }
        //System.out.println("Constraints: " + Arrays.toString(con));
        return value(scv, gamma);
    }

    // J = nummaps. 2*j = number of variables
    public double value(double[] scv, double[] gamma) { // x = [scv(1:J), gamma(1:J)]
        DMatrixRMaj acfCoeff = new DMatrixRMaj();
        DMatrixRMaj lags = new DMatrixRMaj(t.acLags.length);
        for (int i = 0; i < t.acLags.length; i++) {
            lags.set(i, t.acLags[i]);
        }
        double scv_j = scv[0];
        double scv_prev;
        double scalar = 0.5 * (1 - 1 / scv_j);
        //System.out.println("base: " + base); // correct
        CommonOps_DDRM.elementPower(gamma[0], lags, acfCoeff);
        CommonOps_DDRM.scale(scalar, acfCoeff);
        //System.out.println(acfCoeff.toString()); // correct
        DMatrixRMaj X = new DMatrixRMaj();
        DMatrixRMaj Xp1 = new DMatrixRMaj();
        DMatrixRMaj temp = new DMatrixRMaj();
        DMatrixRMaj result = new DMatrixRMaj();
        for (int j = 1; j < J; j++) {
            scv_prev = scv_j;
            scv_j = (1 + scv_j) * ((1 + scv[j])) / 2 - 1;
            double r0j = 0.5 * (1 - 1 / scv[j]);
            CommonOps_DDRM.elementPower(gamma[j], lags, X);
            CommonOps_DDRM.scale(scv[j] * r0j, X); // create X
            CommonOps_DDRM.add(X, 1, Xp1); // X + 1
//            System.out.println("Existing" + acfCoeff.toString());
//            System.out.println("Multiplier" + Xp1.toString());
            CommonOps_DDRM.elementMult(acfCoeff, Xp1, temp);
            CommonOps_DDRM.scale(scv_prev, temp);
            CommonOps_DDRM.add(temp, X, result);
            CommonOps_DDRM.divide(result, scv_j, acfCoeff);
        }
//        System.out.println("scvj: " + scv_j);  // correct
        double e0_sq = t.moments[0] * t.moments[0];
        double trace_scv = (t.moments[1] - e0_sq) / e0_sq;
//        System.out.println("trace scv: " + trace_scv); // correct
        DMatrixRMaj trace_ac = new DMatrixRMaj(t.ac);
        double acNorm = NormOps_DDRM.normP2(trace_ac);
        DMatrixRMaj ac_difference = new DMatrixRMaj();
        CommonOps_DDRM.subtract(trace_ac, acfCoeff, ac_difference);
        double ac_output = NormOps_DDRM.normP1(ac_difference) / acNorm;
        double scv_output = ((scv_j - trace_scv) * (scv_j - trace_scv)) / (trace_scv * trace_scv);
        double value = ac_output + scv_output;
//        System.out.println("value: " + value);
        return value;
    }
}
