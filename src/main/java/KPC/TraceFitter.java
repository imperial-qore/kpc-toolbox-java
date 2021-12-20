package org.qore.KPC;

import de.xypron.jcobyla.Cobyla;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.*;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.RandomMatrices_DDRM;

import java.util.*;
import java.util.logging.Level;

/**
 * Class that handles the logic of fitting a MAP to a Trace
 * including order selection, AC fitting, BC fitting, and display.
 */
public class TraceFitter {
    Trace t;
    FittingOptions options;

    /**
     * Generates new fitting instance with predetermined options object
     * and Trace to be fit with a MAP.
     *
     * @param t Trace object to be fit.
     * @param options Options employed when fitting the Trace.
     */
    public TraceFitter(Trace t, FittingOptions options) {
        this.t = t;
        this.options = options;
        Constants.LOGGER.setLevel(Level.SEVERE);
    }

    /**
     * Fits a MAP to the Trace.
     *
     * @return List of best MAPs found to fit the trace.
     */
    public List<MAP> fit() {
        if (options.numMAPs < 1) {
            if (options.possibleOrders > 1) {
                int[] states = new int[options.possibleOrders];
                for (int i = 0; i < options.possibleOrders; i++) {
                    states[i] = i;
                }
                options.numMAPs = orderSelection(states);
            } else {
                options.numMAPs = orderSelection();
            }
            System.out.printf("Order determined automatically. Fitting MAP with %d states%n", (int) Math.pow(2, options.numMAPs));
        } else {
            System.out.printf("Order determined manually. Fitting MAP with %d states%n", (int) Math.pow(2, options.numMAPs));
        }
        List<ACFitResult> acResults = new ArrayList<>();
        for (int i = 0; i < options.maxRunsAC; i++) {
            ACFitResult ac = fitAC();
            acResults.add(ac);
            System.out.printf("Adding AC fit %d of %d with fit value %f%n", i+1, options.maxRunsAC, ac.value);
        }
        System.out.printf("%nKeeping top %d AC Fits%n%n", options.maxResAC);
        Collections.sort(acResults);
        acResults = acResults.subList(0, options.maxResAC);
        List<MAP> kpcMAPs = new ArrayList<>();
        if (!options.onlyAc) {
            System.out.println("Starting BC fitting process...");
            int i = 0;
            for (ACFitResult ac : acResults) {
                System.out.printf("Adding BC on top of AC fit value %f with value...\t", ac.value);
                BCFitResult bc = fitBC(ac.scv, ac.gamma);
                if (bc == null) {
                    System.out.println("Null");
                    Constants.LOGGER.severe("Could not fit MAP to the AC Fit " +  i);
                } else {
                    System.out.println(bc.value);
                    MAP result;
                    if (options.allowSM) {
                        result = KPC.composeSMMAP(bc.e1, bc.e3, ac.scv, ac.gamma, options.numMAPs, false);
                    } else {
                        result = KPC.composeMAP(bc.e1, bc.e3, ac.scv, ac.gamma, options.numMAPs, false);
                    }

                    if (!result.isFeasible()) {
                        System.out.println("result is not feasible ");
                    }
                    result.scale(t.moments[0], options.allowSM);
                    kpcMAPs.add(result);
                }
                i++;
            }
        } else {
            System.out.println("Skipping BC fitting process...");
            for (ACFitResult ac : acResults) {
                double[] e1 = new double[options.numMAPs];
                double[] e3 = new double[options.numMAPs];
                for (int i = 0; i < options.numMAPs; i++) {
                    e1[i] = 1;
                    e3[i] = 1.51 * (1 + ac.scv[i]) * (1 + ac.scv[i]);
                }
                if (options.allowSM) {
                    kpcMAPs.add(KPC.composeSMMAP(e1, e3, ac.scv, ac.gamma, options.numMAPs, true));
                } else {
                    kpcMAPs.add(KPC.composeMAP(e1, e3, ac.scv, ac.gamma, options.numMAPs, true));
                }

            }
        }
        return kpcMAPs;
    }

    /**
     * Fits the SCV and Gamma of a trace to a MAP with
     * similar characteristics.
     *
     * @return ACFitResult object comprising of best found values.
     */
    public ACFitResult fitAC() {
        int J = options.numMAPs;
        ObjectiveFunction f = new ObjectiveFunction(new ACObjFunc(t, J));
        SimpleBounds bounds = new SimpleBounds(createLB(J), createUB(J));
        double[] x = new double[2*J];
        for (int i = 0; i < J; i++) {
            x[i] = Math.random() + 1;
            x[i+J] = Math.random();
        }
        InitialGuess x0 = new InitialGuess(x);
        RandomGenerator r = new MersenneTwister();
        MaxEval maxEvals;
        if (options.maxEvalsAC == 0) {
            maxEvals = MaxEval.unlimited();
        } else {
            maxEvals = new MaxEval(options.maxEvalsAC);
        }
        double[] sigma = new double[2*J];
        Arrays.fill(sigma, options.sigma);
        CMAESOptimizer optimizer = new CMAESOptimizer(
                options.maxIterAC, 1e-5, true, 10, 50, r, false, null);
        PointValuePair pv = optimizer.optimize(maxEvals, f, bounds, x0, GoalType.MINIMIZE,
                new CMAESOptimizer.Sigma(sigma),
                new CMAESOptimizer.PopulationSize(10*J));
        double[] scv = Arrays.copyOfRange(pv.getPoint(), 0, J);
        double[] gamma = Arrays.copyOfRange(pv.getPoint(), J, 2*J);
        return new ACFitResult(scv, gamma, pv.getValue());
    }

    /**
     * Fits the bicorrelations for a given AC fit. Runs multiple times and chooses
     * the best result.
     *
     * @param scv The SCV value from the AC fit
     * @param gamma The Gamma value from the AC fit
     * @return BCFitResult object comprising of best found values
     */
    public BCFitResult fitBC(double[] scv, double[] gamma) {
        int J = options.numMAPs;
        double e1val = Math.pow(t.moments[0], 1.0 / options.numMAPs);
        double[] e2 = new double[J];
        double[] x = new double[2*J];
        Arrays.fill(x, 0, J, e1val);
        double r = Math.random() * t.moments[0];
        for (int j = 0; j < J; j++) {
            e2[j] = (1 + scv[j]) * x[j] * x[j];
            x[J + j] = (1.5 + r) * e2[j] * e2[j] / x[j];
        }
        DMatrixRMaj x0 = new DMatrixRMaj(x);
        DMatrixRMaj xNext = new DMatrixRMaj(x.length);

        double bestVal = Double.MAX_VALUE;
        BCFitResult bestFit = null;
        for (int i = 0; i < options.maxRunsBC; i++) {
            BCObjFunc func = new BCObjFunc(t, options.numMAPs, scv, gamma, options.allowSM, options.verbose);
            int m = 4 * (J-1) + 2 + ((scv[0] > 1) ? 2 : 0);
            Cobyla.findMinimum(func, 2 * J, m, x,
                    options.rhobeg, options.rhoend, options.iprint, options.maxIterBC);
            // save value with reference to AC vals
            double[] e1j = Arrays.copyOfRange(x, 0, J);
            double[] e3j = Arrays.copyOfRange(x, J, x.length);
            double[] e2j = new double[J];
            for (int j = 0; j < J; j++) {
                e2j[j] = (1 + scv[j])*e1j[j]*e1j[j];
            }
            double value = func.value(e1j, e3j);
            if (value < bestVal) {
                bestFit = new BCFitResult(e1j, e2j, e3j, value);
                bestVal = value;
            }
            // recreate x0 and E2
            DMatrixRMaj rnd = new DMatrixRMaj(J, 1);
            CommonOps_DDRM.fill(rnd, 0.25);
            RandomMatrices_DDRM.addUniform(rnd, 0, 1.75, new Random());
            DMatrixRMaj one = new DMatrixRMaj(J,1);
            CommonOps_DDRM.fill(one, 1);
            DMatrixRMaj multiplier = new DMatrixRMaj(2*J, 1);
            CommonOps_DDRM.concatRows(rnd, one, multiplier);
            CommonOps_DDRM.elementMult(x0, multiplier, xNext);
            x = xNext.data; // set e1 for next itr
            r = Math.random() + t.moments[0];
            for (int j = 0; j < J; j++) { // set next starting parameters
                e2[j] = (1 + scv[j]) * x[j] * x[j];
                x[J + j] = (1.5 + r) * e2[j] * e2[j] / x[j];
            }
        }
        if (bestVal == Double.MAX_VALUE) {
            return null;
        }
        return bestFit;
    }

    /**
     * Chooses the optimal order for the MAP fit
     *
     * @return Optimal number of states
     */
    public int orderSelection() {
        return orderSelection(new int[] {1, 2, 3, 4, 5, 6});
    }

    /**
     * Chooses the optimal order for the MAP fit
     *
     * @param states Potential orders of MAP fit to choose from.
     * @return Optimal number of states
     */
    public int orderSelection(int[] states) {
        if (options.osCriterion == OSCriterion.BIC) {
            return findBIC(states);
        } else {
            System.out.println("Other comaprison methods not implemented yet");
            return 4;
        }
    }

    /**
     * Performs order selection using Bayseian Information Criterion (BIC).
     *
     * @param states Potential orders of MAP fit to choose from.
     * @return Optimal number of states according to minimal BIC.
     */
    public int findBIC(int[] states) {
        int nlags = t.acFull.length;
        int nLagsEnd = nlags;
        int orderMax = (int) Math.pow(2, states[states.length-1]);
        for (int i = 0; i < nlags-1; i++) {
            if (t.acFull[i+1] < Constants.CONSTRAINT_TOL) {
                nLagsEnd = i - 1 + orderMax;
                break;
            }
        }
        final int NLAGSMAX = 10000;
        int[] SAlags;
        if (nLagsEnd > NLAGSMAX) {
            int[] SAlagsT = Equations.logspacei(1, nLagsEnd - orderMax, NLAGSMAX);
            double temp = NLAGSMAX - SAlagsT.length;
            double step = Math.pow((nLagsEnd - orderMax), (1/temp));
            int i = 0;
            while (i < SAlagsT.length && (SAlagsT[i] - Math.round(Math.pow(step, i))) < temp) {
                i++;
            }
            int newLength = SAlagsT[i] + SAlagsT.length - i;
            SAlags = new int[newLength];
            int jj = i;
            for (int j = 0; j < newLength; j++) {
                if (j < SAlagsT[i]) {
                    SAlags[j] = j+1;
                } else {
                    SAlags[j] = SAlagsT[jj++]+1;
                }
            }
        } else {
            SAlags = new int[nLagsEnd - orderMax];
            for (int i = 0; i < nLagsEnd - orderMax; i++) {
                SAlags[i] = i;
            }
        }
        int nSamples = SAlags.length;
        double[][] allX = new double[orderMax][nSamples];
        double[] y = t.acAtLags(SAlags, nSamples);
        for (int i = 1; i <= orderMax; i++) {
            Equations.elementAdd(SAlags, 1);
            allX[i-1] = t.acAtLags(SAlags, nSamples);
//            int[] temp = Arrays.copyOfRange(SAlags, SAlags.length-50, SAlags.length);
//            System.out.println(Arrays.toString(temp));
        }
        int bestOrder = -1;
        double bestVal = Double.MAX_VALUE;

        OLSMultipleLinearRegression model;
        for (int nMaps : states) {
            int order = (int) Math.pow(2, nMaps);
            double[][] x = Arrays.copyOfRange(allX, 0, order);
            x = Equations.transpose(x);
            model = new OLSMultipleLinearRegression();
            model.newSampleData(y, x);
            double rss = model.calculateResidualSumOfSquares();
            double bic = nSamples*Math.log(rss) - nSamples*Math.log(nSamples) + order*Math.log(nSamples);

            if (bic < bestVal) {
                bestVal = bic;
                bestOrder = nMaps;
            }
        }

        if (bestVal > .05) {
            Constants.LOGGER.severe("The RSS indicates that the data may not be fit well by a MAP");
            Constants.LOGGER.severe("Consider using another mdoel");
        }
        return bestOrder != -1 ? bestOrder : 3;
    }

    /**
     * Prints to console to closeness of fit of some MAPS
     * to the trace.
     *
     * @param maps The MAPs to be compared to the trace.
     */
    public void displayResults(List<MAP> maps) {
        int bestMAP = -1;
        double bestfit = Double.MAX_VALUE;
        for (int i = 0; i < maps.size(); i++) {
            MAP map = maps.get(i);
            double[] result = KPC.evaluate(t, map, true);
            double val = result[0] + result[1];
            String s = String.format("MAP %d: ACF: %f, BCF: %f, Total: %f", i+1, result[0], result[1], val);
            if (!map.isFeasible()) {
                s += "\t(Semi-Markov)";
            }
            if (val < bestfit) {
                bestfit = val;
                bestMAP = i;
                s += "\t**** best ****";
            }
            System.out.println(s);
        }
        MAP best = maps.get(bestMAP);
        System.out.println("\nFitting process completed. Best MAP shown below\n");
        System.out.println("Moments Comparison\n");
        System.out.println("Moment\t\tOriginal Trace\t\t\t\tFit MAP\t\t\t\t\t\tDifference");
        double[] moments = best.getMoments(new int[] {1, 2, 3});
        for (int i = 0; i < 3; i++) {
            System.out.println(i+1 + "\t\t" + t.moments[i] + "\t\t" + moments[i] + "\t\t" + (t.moments[i]- moments[i]));
        }
        System.out.println("\nAutocorrelation Comparison\n");
        System.out.println("Lag\t\t\tOriginal Trace\t\t\t\tFit MAP\t\t\t\t\t\tDifference");
        double[] ac = best.getAcf(new int[] {1,2,3,4,5});
        for (int i = 0; i < 5; i++) {
            System.out.println(i+1 + "\t\t" + t.ac[i] + "\t\t\t" + ac[i] + "\t\t\t" + (t.ac[i]- ac[i]));
        }
    }



    /**
     * Creates the lower bounds of the data array to be used
     * by the Apache Optimizer.
     *
     * @param J Number of MAPS
     * @return list of lower bounds ordered in correspondence with input
     */
    private double[] createLB(int J) {
        // scv(0) > 0.5;    scv(j) > 1;     gamma(j) > 0
        double[] lb = new double[2*J];
        lb[0] = 0.5 - Constants.CONSTRAINT_TOL;
        lb[J] = 0;
        for (int i = 1; i < J; i++) {
            lb[i] = 1 + Constants.CONSTRAINT_TOL;
            lb[J+i] = Constants.CONSTRAINT_TOL;
        }
        return lb;
    }

    /**
     * Creates the upper bounds of the data array to be used
     * by the Apache Optimizer.
     *
     * @param J Number of MAPS
     * @return list of upper bounds ordered in correspondence with input
     */
    private double[] createUB(int J) {
        // gamma(j) < 1
        double[] ub = new double[2*J];
        for (int i = 0; i < J; i++) {
            ub[i] = Double.MAX_VALUE;
            ub[J+i] = 1 + Constants.CONSTRAINT_TOL;
        }
        return ub;
    }

    /**
     * Old AC fitting method implementing JCOBYLA as optimizer, removed due to
     * poor performance. Kept in case needed in future
     */
    @Deprecated
    public ACFitResult oldAC () {
        int J = options.numMAPs;
        System.out.println("Fitting " + J + " MAPs");
        OldACObjFunc func = new OldACObjFunc(t, J);
        double[] x = new double[2*J];
        for (int i = 0; i < J; i++) {
            x[i] = Math.random() + 1;
        }
        for (int i = J; i < 2*J; i++) {
            x[i] = Math.random();
        }
        Cobyla.findMinimum(func, 2*J, 3*J, x,
                options.rhobeg, options.rhoend, options.iprint, options.maxIterAC);
        double[] scv = Arrays.copyOfRange(x, 0, J);
        double[] gamma = Arrays.copyOfRange(x, J, x.length);
        ACFitResult acFit = new ACFitResult(scv, gamma, func.value(scv, gamma));
        System.out.println(acFit.value);
        return acFit;
    }
}