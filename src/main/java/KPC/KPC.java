package KPC;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;

import java.util.Arrays;

/**
 * Class handles the composition logic of MAP(2)s to form MAPs
 * with larger numbers of states.
 */
public class KPC {

    /**
     * Kronecker composes two MAP objects according to the relation defined
     * in Casale et. al. (2007).
     *
     * @param a First MAP.
     * @param b Second MAP.
     * @return Kronecker composed MAP of inputs.
     */
    public static MAP kroneckerCompose(MAP a, MAP b) {
        if (b == null) {
            return a;
        }
        DMatrixRMaj d0_c = new DMatrixRMaj();
        DMatrixRMaj d1_c = new DMatrixRMaj();
        DMatrixRMaj neg_a_d0 = new DMatrixRMaj();
        CommonOps_DDRM.changeSign(a.D0, neg_a_d0);
        CommonOps_DDRM.kron(neg_a_d0, b.D0, d0_c);
        CommonOps_DDRM.kron(a.D1, b.D1, d1_c);

        return new MAP(d0_c, d1_c);
    }

    /**
     * Kronecker composes two MAP objects according to the relation defined
     * in Casale et. al. (2007).
     *
     * @param a First MAP.
     * @param b Second MAP.
     * @return Kronecker composed MAP of inputs.
     */
    public static SMP kroneckerCompose(SMP a, MAP b) {
        DMatrixRMaj d0_c = new DMatrixRMaj();
        DMatrixRMaj p_c = new DMatrixRMaj();
        DMatrixRMaj neg_a_d0 = new DMatrixRMaj();
        CommonOps_DDRM.changeSign(a.D0, neg_a_d0);
        CommonOps_DDRM.kron(neg_a_d0, b.D0, d0_c);
        CommonOps_DDRM.kron(a.P, b.P, p_c);
        return new SMP(d0_c, p_c);
    }

    /**
     * Function creates a MAP(J) of chosen order based
     * by constructing the MAP from the characteristics given
     *
     * @param e1j vector of values of first moment to be fit
     * @param e3j vector of values of 3rd moment to be fit
     * @param scv vector of values of SCV
     * @param gamma vector of values of Gamma
     * @param J numStates of the composed MAP
     * @return composed MAP of the proposed values, or null if
     * we can not find a suitable MAP for those characteristics
     */
    public static MAP composeMAP(double[] e1j, double[] e3j, double[] scv, double[] gamma, int J, boolean verbose) {
        double e2_0 = (1 + scv[0])*e1j[0]*e1j[0];
        MAP arbitraryMAP = Fittings2.mmpp(e1j[0], e2_0, e3j[0], scv[0], gamma[0]);
        if (!arbitraryMAP.isFeasible()) {
            if (scv[0] < 0.5) {
                arbitraryMAP = Fittings2.erlang(e1j[0], 2);
                if (verbose) {
                    Constants.LOGGER.info("MAP 1 is erlang");
                }
            } else {
                if (verbose) {
                    Constants.LOGGER.info("MAP 1 presumably infeasible E3");
                }
                arbitraryMAP = Fittings2.map(e1j[0], e2_0, -1, scv[0], gamma[0]);
                if (arbitraryMAP == null || !arbitraryMAP.isFeasible()) {
                    if (verbose) {
                        Constants.LOGGER.info("MAP 1 presumably infeasible Gamma");
                    }
                    arbitraryMAP = Fittings2.map(e1j[0], e2_0, -1, scv[0], 0);
                    if (arbitraryMAP == null || !arbitraryMAP.isFeasible()) {
                        throw new InfeasibleMAPException("Can not generate arbitrary MAP with selected values");
                    }
                }
            }
        } else {
            Constants.LOGGER.info("MAP 1 generated with map2");
        }
        if (arbitraryMAP == null || !arbitraryMAP.isFeasible()) {
            throw new InfeasibleMAPException("Arbitrary MAP 1 can not be formed");
        }
        MAP kpc = composeDiagonalMAPs(e1j, e3j, scv, gamma, J, verbose);
        kpc = kroneckerCompose(arbitraryMAP, kpc);
        if (!kpc.isFeasible()) {
            throw new InfeasibleMAPException("MAP composition results in Infeasible MAP");
        }
        kpc.normalize();
        return kpc;
    }


    /**
     * Function creates a SMP(J) of chosen order based
     * by constructing the MAP from the characteristics given.
     *
     * @param e1j vector of values of first moment to be fit.
     * @param e3j vector of values of 3rd moment to be fit.
     * @param scv vector of values of SCV.
     * @param gamma vector of values of Gamma.
     * @param J numStates of the composed MAP.
     * @return composed MAP of the proposed values, or null if
     * we can not find a suitable MAP for those characteristics.
     */
    public static SMP composeSMMAP(double[] e1j, double[] e3j, double[] scv, double[] gamma, int J, boolean verbose) {
        double e2_0 = (1 + scv[0])*e1j[0]*e1j[0];
        SMP sm = Fittings2.sm(e1j[0], e2_0, e3j[0], gamma[0]); // arbitrary MAP
        MAP kpc = composeDiagonalMAPs(e1j, e3j, scv, gamma, J, verbose); // Rest of MAPs

        return kroneckerCompose(sm, kpc);
    }

    /**
     * Given the first arbitrary MAP, composes the remaining J-1 diagonal D0 maps to fit the characteristics
     * given in the input parameters.
     *
     * @param e1j vector of values of first moment to be fit.
     * @param e3j vector of values of 3rd moment to be fit.
     * @param scv vector of values of SCV.
     * @param gamma vector of values of Gamma.
     * @param J numStates of the composed MAP.
     * @return composed MAP of the proposed values.
     */
    public static MAP composeDiagonalMAPs(double[] e1j, double[] e3j, double[] scv, double[] gamma, int J, boolean verbose) {
        MAP kpc = null;
        for (int j = 1; j < J; j++) {
            double e2_j = (1 + scv[j]) * Math.pow(e1j[j], 2);
            // Solves inconsistent rounding error
            String sValue = String.format("%.15f", e2_j);
            e2_j = Double.parseDouble(sValue);
            MAP mapj = Fittings2.feasblock(e1j[j], e2_j, e3j[j], gamma[j]);
            if (mapj == null || !mapj.isFeasible()) {
                if (scv[j] < 1) {
                    if (verbose) {
                        Constants.LOGGER.info(String.format("MAP %d has low variability", j + 1));
                    }
                    mapj = Fittings2.mapExp(e1j[j]);
                } else {
                    if (verbose) {
                        Constants.LOGGER.info(String.format("MAP %d likely has infeasible E3", j + 1));
                    }
                    mapj = Fittings2.map(e1j[j], e2_j, -1, scv[0], gamma[j]);
                    if (mapj == null) {
                        if (verbose) {
                            Constants.LOGGER.info(String.format("MAP %d has infeasible Gamma", j + 1));
                        }
                        mapj = Fittings2.map(e1j[j], e2_j, -1, scv[0], 0);
                        if (mapj == null) {
                            Constants.LOGGER.severe(String.format("MAP %d can not be formed", j + 1));
                            throw new InfeasibleMAPException(String.format("MAP %d can not be formed", j + 1));
                            // TODO handle error codes
                        }
                    }
                }
            }
            if (mapj == null) {
                mapj = Fittings2.mapExp(e1j[j]);
            }
            if (mapj == null || !mapj.isFeasible()) {
                Constants.LOGGER.severe("Trying to compose with null MAP");
                throw new InfeasibleMAPException(String.format("MAP %d is can not be formed", j + 1));
            }
            if (j == 1) {
                kpc = mapj;
            } else {
                kpc = KPC.kroneckerCompose(kpc, mapj);
            }
        }
        Constants.LOGGER.info("Successfully composed map with traits:");
        Constants.LOGGER.info(Arrays.toString(e1j) + ", " + Arrays.toString(e3j));
        return kpc;
    }

    /**
     * Evaluates the fit of a MAP to a trace by computing the objective functions of
     * the MAP now that it has been composed from the fit characteristics.
     *
     * @param map the MAP to be evaluated
     * @return vector with 2 items [ac fit value, bc fit value].
     */
    public static double[] evaluate(Trace t, MAP map, boolean evalBC) {
        double e0_sq = t.moments[0] * t.moments[0];
        double trace_scv = (t.moments[1] - e0_sq) / e0_sq;
        DMatrixRMaj map_acf = new DMatrixRMaj(map.getAcf(t.acLags));
        DMatrixRMaj trace_acf = new DMatrixRMaj(t.ac);
        DMatrixRMaj diff_acf = new DMatrixRMaj();
        CommonOps_DDRM.subtract(trace_acf, map_acf, diff_acf);
        double acfEval = NormOps_DDRM.normP1(diff_acf) / NormOps_DDRM.normP2(trace_acf) +
                (map.getSCV()-trace_scv) * (map.getSCV()-trace_scv) / (trace_scv * trace_scv);
        if (evalBC) {
            DMatrixRMaj trace_bc = new DMatrixRMaj(t.bc);
            double[] bcVals = new double[t.bcLags.length];
            DMatrixRMaj diff_bc = new DMatrixRMaj();
            for (int i = 0; i < t.bcLags.length; i++) {
                double val = map.getJoint(t.bcLags[i], new int[]{1, 1, 1});
                bcVals[i] = val;
            }
            DMatrixRMaj map_bc = new DMatrixRMaj(bcVals);
            CommonOps_DDRM.subtract(trace_bc, map_bc, diff_bc);
            double bcEval = NormOps_DDRM.normP2(diff_bc) / NormOps_DDRM.normP2(trace_bc);
            return new double[] {acfEval, bcEval};
        } else {
            return new double[] {acfEval, 0};
        }
    }
}