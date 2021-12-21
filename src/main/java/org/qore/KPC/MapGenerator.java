package org.qore.KPC;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class MapGenerator {

    /**
     * Generates a Random valid MAP.
     *
     * @param nStates number of states of the MAP to generate.
     * @return Random MAP with desired number of states.
     */
    public static MAP randomMAP(int nStates) {
        DMatrixRMaj D0 = new DMatrixRMaj(nStates, nStates);
        DMatrixRMaj D1 = new DMatrixRMaj(nStates, nStates);
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                D0.set(i, j, Math.random());
                D1.set(i, j, Math.random());
            }
        }
        MAP m = new MAP(D0, D1);
        m.normalize();
        return m;
    }

    /**
     * Generates a random Semi Markov Process with the desired number of states.
     * Uses KPC and diagonal MAPs to ensure valid D0.
     *
     * @param nMaps number of maps to use in composition. Results in 2^nMaps state process.
     * @return Random valid SMP.
     */
    public static SMP randomSMP(int nMaps) {
        boolean goon = true;
        DMatrixRMaj D0;
        DMatrixRMaj P;
        SMP m = null;
        while (goon) {
            D0 = randomMAP(2).D0;
            P = randomP(2);
            SMP sm = new SMP(D0, P);
            double[] moments = sm.getMoments();
            double gamma = sm.getGamma();
            if (gamma > 0 && gamma < 1) {
                m = Fittings2.sm(moments[0], moments[1], moments[2], gamma);
                if (!m.isFeasible()) {
                    goon = false;
                }
            }
        }
        for (int j = 0; j < nMaps-1; j++) {
            MAP diag;
            do {
                MAP rnd = randomMAP(2);
                double[] moments = rnd.getMoments();
                double scv = (moments[1] - moments[0] * moments[0]) / (moments[0] * moments[0]);
                double g2 = rnd.getGamma();
                diag = Fittings2.feasblock(moments[0], moments[2], scv, g2);
            } while (diag == null || !diag.isFeasible());
            m = KPC.kroneckerCompose(m, diag);
        }
        return m;
    }

    /**
     * Generates a random row-wise Markovian probability matrix.
     *
     * @param size number of states of the MC.
     * @return Valid MC.
     */
    public static DMatrixRMaj randomP(int size) {
        double p = Math.random();
        double q = Math.random();
        return new DMatrixRMaj(new double[][]{{p, 1-p}, {1-q, q}});
    }

    /**
     * Generates sample of interarrival times from MAP.
     *
     * @param nSamples Number of iats to generate.
     */
    public static Sample sample(MAP m , int nSamples) {
        int numStates = m.numStates;
        double[][] pdf = new double[numStates][2*numStates];
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                double diag = Math.abs(m.D0.get(i, i));
                if (i == j) {
                    pdf[i][j] = 0;
                } else {
                    pdf[i][j] = m.D0.get(i, j) / diag;
                }
                pdf[i][numStates+j] = m.D1.get(i,j) / diag;
            }
        }
        Sample sample = new Sample(nSamples);
        List<Integer>[] visits = new ArrayList[nSamples];
        for (int i = 0; i < nSamples; i++) {
            visits[i] = new ArrayList<>();
        }
        int currState = chooseState(m.pi.getData());
        for (int i = 0; i < nSamples; i++) {
            boolean arrived = false;
            sample.initials[i] = currState;
            visits[i].add(currState);
            while (!arrived) {
                int destState = chooseState(pdf[currState]);
                if (destState >= numStates) {
                    arrived = true;
                    destState -= numStates;
                    sample.absorbings[i] = currState;
                } else {
                    visits[i].add(destState);
                }
                currState = destState;
            }
        }
        ExponentialDistribution[] expD = new ExponentialDistribution[numStates];
        for (int j = 0; j < numStates; j++) {
            expD[j] = new ExponentialDistribution(-1.0 / m.D0.get(j,j));
        }
        for (int i = 0; i < nSamples; i++) {
            for (int state : visits[i]) {
                sample.iats[i] += expD[state].sample();
            }
        }
        return sample;
    }

    /**
     * Samples IATs from a Semi-Markov distribution
     *
     * @param sm semi-markov MAP to sample from
     * @param nSamples number of iats to generate
     * @return Sample object with iats, initial, and end states.
     */
    public static Sample sample(SMP sm, int nSamples) {
        DMatrixRMaj e = MatrixOps.getOnes(sm.numStates);
        DMatrixRMaj tmp = new DMatrixRMaj(sm.numStates, sm.numStates);
        DMatrixRMaj D1 = new DMatrixRMaj(sm.numStates, sm.numStates);
        CommonOps_DDRM.mult(e, sm.pi, tmp);
        CommonOps_DDRM.mult(-1, sm.D0, tmp, D1);
        MAP map = new MAP(sm.D0, D1);
        return sample(map, nSamples);
    }

    /**
     * Randomly choose a state from a vector of probabilities
     *
     * @param probs Probability vector to choose next state.
     * @throws IllegalArgumentException if input does not sum to 1
     * @return Next state chosen.
     */
    protected static int chooseState(double[] probs) {
        double sum = 0;
        for (double d : probs) {
            sum += d;
        }
        if ((sum - 1) > Constants.ZERO) {
            throw new IllegalArgumentException("Probability vector does not sum to 1");
        }
        double r = Math.random();
        sum = 0;
        int state = 0;
        while (state < probs.length) {
            if (r > sum && r < sum + probs[state]) {
                return state;
            }
            sum += probs[state];
            state++;
        }

        throw new IllegalArgumentException("Something went wrong during sample");
    }
}
