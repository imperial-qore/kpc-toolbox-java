package org.qore.KPC;

import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;

import static org.junit.Assert.assertEquals;

public class Test16 {
    Trace t;
    FittingOptions options;
    TraceFitter fit;

    @Before
    public void setUp() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/bcaug.csv";
        t = new Trace(path);
        options = new FittingOptions(4);
        fit = new TraceFitter(t, options);
        Constants.LOGGER.setLevel(Level.WARNING);
    }

    @Test
    public void testBestFit() {
        double[] e1_e3 = {2.591145363010744e-01, 2.096890904703476e-01, 2.397661366398542e-01, 2.412484536001585e-01,
                2.501740216612849e-01, 1.849545786232304e-01, 8.074768209178093e-02, 8.319772551933895e-02};
        double[] e1 = Arrays.copyOfRange(e1_e3, 0, 4);
        double[] e3 = Arrays.copyOfRange(e1_e3, 4, 8);
        double[] scv_g2 = {1.737872678096437e+00, 1.897418799086275e+00, 1.062984267448193e+00, 1.081043681188437e+00,
                9.547093751230068e-01, 3.796387006758699e-01, 9.991639804644602e-01, 9.999639317576479e-01};
        double[] scv = Arrays.copyOfRange(scv_g2, 0, 4);
        double[] g2 = Arrays.copyOfRange(scv_g2, 4, 8);
        double[] e2 = {1.838217104855952e-01, 1.273980983691944e-01, 1.185964275461960e-01, 1.211184411352844e-01};

        ACObjFunc fA = new ACObjFunc(t, 4);
        double rA = fA.value(scv_g2);
        BCObjFunc fB = new BCObjFunc(t, 4, scv, g2, false, false);
        double rB = fB.value(e1, e3);
    }

    @Test
    public void testLargeDEC() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/DEC-PKT-1-UDP.csv";
        t = new Trace(path);
        options = new FittingOptions(4);
        options.maxEvalsAC = 10000;
        options.maxIterAC = 550;
        options.sigma = 0.35;
        fit = new TraceFitter(t, options);
        ACFitResult result = fit.fitAC();
//        System.out.println("scv: " + Arrays.toString(result.scv));
//        System.out.println("gamma: " + Arrays.toString(result.gamma));
//        System.out.println("value: " + result.value);
    }

    @Test
    public void testGoodAC() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/DEC-PKT-1-UDP.csv";
        t = new Trace(path);
        options = new FittingOptions(4);
        options.maxEvalsAC = 1000;
        fit = new TraceFitter(t, options);
        ACObjFunc f = new ACObjFunc(t, 4);
        double[] x = {1.016711847968116e+00,
                1.183681770337907e+00,
                1.152056626165421e+00,
                1.138838493558732e+00,
                9.999850583519337e-01,
                6.698188348670987e-01,
                9.516043163749496e-01,
                2.484090690622519e-01};
        double val = f.value(x);
        assertEquals(val, 2.942869, 1e-5);

    }

    /*
    @Test
    public void findGoodParameters() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/DEC-PKT-1-UDP.csv";
        t = new Trace(path);
        options = new FittingOptions(4);
        double bestVal = Double.MAX_VALUE;
        int bestI = 0;
        int bestJ = 0;
        double bestK = 0;
        for (int i = 10000; i <= 20000; i += 1000) {
            for (int j = 400; j <= 600; j += 50) {
                for (double k = .3; k < .65; k += .05) {
                    options.maxEvalsAC = i;
                    options.maxIterAC = j;
                    options.sigma = k;
                    fit = new TraceFitting(t, options);
                    ACFitResult result = fit.fitAC();
                    System.out.printf("value = %f, i = %d, j = %d, k = %f%n", result.value, i, j, k);
                    if (result.value < bestVal) {
                        bestI = i;
                        bestJ = j;
                        bestK = k;
                        bestVal = result.value;
                    }
                }
            }
        }
        // best so far: 10,000, 550, .35
        // best i = 15000, best j = 550, best k = 0.450000
        System.out.printf("Best value = %f, best i = %d, best j = %d, best k = %f%n", bestVal, bestI, bestJ, bestK);
    }
    */

    @Test
    public void testLargeEvaluation() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/bcaug.csv";
        t = new Trace(path);
        options = new FittingOptions(4);
//        fit = new TraceFitting(t, options);
        double[] scv = {1.0796000166400115, 1.7366150862756324, 1.8437250215884489, 1.1223816802900712};
        double[] g2 = {0.21621851494800032, 0.9578401944876733, 0.5677839493332579, 0.9999279678640932};
        double[] e1 = {0.22337416518644793, 0.2576790179800146, 0.24077694071886813, 0.2586441743871845};
        double[] e3 = {0.12018050952086921, 0.1911958300666558, 0.16942471534913178, 0.11672722596964862};
        BCObjFunc fB = new BCObjFunc(t, 4, scv, g2, true,false);
        double bcval = fB.value(e1, e3);
//        System.out.println("obj func value: " + bcval);
        MAP m = KPC.composeSMMAP(e1, e3, scv, g2, 4, false);
        double[] vals = KPC.evaluate(t, m, true);
//        System.out.println("evaluation: " + Arrays.toString(vals));
    }

    @Test
    public void testLargeEvaluation2() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/bcaug.csv";
        t = new Trace(path);
//        options = new FittingOptions(4);
//        fit = new TraceFitting(t, options);
        double[] scv = {1.1861809890753094, 1.7776835604147927, 1.1223255573326945, 1.5404566651509126};
        double[] g2 = {0.1933519204626638, 0.9537714918223732, 0.9999120419560055, 0.4440024098932075};
        double[] x = {1.1861809890753094, 1.7776835604147927, 1.1223255573326945, 1.5404566651509126,
                        0.1933519204626638, 0.9537714918223732, 0.9999120419560055, 0.4440024098932075};
        double[] e1 = {0.09373893634139145, 0.20088504776821417, 0.3008832638899168, 0.5209770343231517};
        double[] e3 = {0.010588443321990305, 0.08802270790788456, 0.1824146361421601, 1.368945322715493};
        ACObjFunc fA = new ACObjFunc(t, 4);
        double acval = fA.value(x);
        BCObjFunc fB = new BCObjFunc(t, 4, scv, g2, true,false);
        double bcval = fB.value(e1, e3);
//        System.out.println("obj func ac: " + acval + "\t bc: " + bcval);
        MAP m = KPC.composeSMMAP(e1, e3, scv, g2, 4, false);
        double[] vals = KPC.evaluate(t, m, true);
//        System.out.println("evaluation: " + Arrays.toString(vals));
        double[] s_vals = KPC.evaluate(t, m, true);
//        System.out.println("scaled evaluation: " + Arrays.toString(s_vals));
    }

    @Test
    public void testLargeEvaluation3() {
        // TODO these composed to infeasible MAP but are feasible on MATLAB
        double[] scv = {0.9855495252274138, 1.9294583038292972, 1.1127600486344094, 1.7533553833480633};
        double[] gamma = {0.9339467621427353, 0.45900857401091527, 0.9999343013628589, 0.9600850254344416};
        double[] e1 = {0.32599621982529686, 0.4671832963107999, 0.3177891868290314, 0.4167936540935424};
        double[] e3 = {0.20364052206757546, 1.5259282891505421, 0.21696951737071593, 0.9518194737847541};
        MAP m = KPC.composeSMMAP(e1, e3, scv, gamma, 4, true);
    }
}
