package org.qore.KPC;

import de.xypron.jcobyla.Calcfc;
import de.xypron.jcobyla.Cobyla;
import de.xypron.jcobyla.CobylaExitStatus;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class TraceFitter4Test {
    Trace t;
    FittingOptions options;
    TraceFitter fit;
    double[] scv = {1.616776856193191e+00, 1.990734816761814e+00};
    double[] gamma ={9.791016635338956e-01, 6.318114925844853e-01};
    double[] e2 ={8.224067899492617e-03, 9.399351780498133e-03};
    double[] e1 ={0.056060891348814, 0.056060891348814};
    double[] e3 ={0.001812515516495, 0.002367577293816};

    /* Good fit values
    double[] e1 ={0.059635641301638, 0.052700423274181};
    double[] e3 ={0.002830570564294, 0.002202716859225};
     */


    @Before
    public void setUp() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/bcaug.csv";
        t = new Trace(path);
        options = new FittingOptions();
        options.numMAPs = 2;
        fit = new TraceFitter(t, options);
        Constants.LOGGER.setLevel(Level.WARNING);
    }

    /**
     * Tests used when writing files to capture errors
     */
    @Test
    public void testFitACExecutes() {
        fit.options.numMAPs = 2;
        fit.options.iprint = 2;
        fit.fitAC();
    }

    @Test
    public void testFitBCExecutes() {
        ACFitResult prev = fit.fitAC();
        fit.fitBC(prev.scv, prev.gamma);
    }


    @Test
    public void testACFValue4() {
        ACObjFunc f = new ACObjFunc(t, options.numMAPs);
        double[] x = {scv[0], scv[1], gamma[0], gamma[1]};
        double v = f.value(x);
        double real = 6.163844741297546e+00;
        assertEquals(v, real, 1e-6);
    }

    @Test
    public void testACFValue8() {
        double[] ac_data = {1.108041114528882e+00, 1.700711825145974e+00, 1.900031326673984e+00,
                9.999362300262029e-01, 9.613605133281835e-01, 4.385767777898044e-01};
        ACObjFunc f = new ACObjFunc(t, 3);
        double v = f.value(ac_data);
        assertEquals(v, 1.7559334593523, 1e-5);
    }

    /**
     * Tests that a given call to compose with the parameters from an
     * ac fit result returns the correct MAP
     * values obtained from matlab run
     * Does not check infeasible branch
     */

    /**
     * Tests that the BC Fit objective function returns the correct value
     * on MAP(4)s. Values obtained form MATLAB test run
     */
    @Test
    public void testBCValue4() throws IOException {
        setUp();
        BCObjFunc f = new BCObjFunc(t, options.numMAPs, scv, gamma, false, options.verbose);
        double v = f.value(e1, e3);
        System.out.println("value = " + v);
        double answer = 0.154476443957429;
        assertEquals(v, answer, 1e-5);
    }

    @Test
    public void testBCConstraints4() {
        BCObjFunc f = new BCObjFunc(t, options.numMAPs, scv, gamma, false, options.verbose);
        double[] cons = new double[8];
        double[] x = new double[4];
        x[0] = e1[0];
        x[1] = e1[1];
        x[2] = e3[0];
        x[3] = e3[1];
        f.compute(4, 8, x, cons);
        double[] real = new double[] {
                0.001938420790417,
                0.000002823635781,
                0.003113704671423,
                0.000003688341369,
                1.644248566792492,
                -0.144248566792492,
                -0.000022135752310,
                0.028395993367214 };
        assertArrayEquals(cons, real, 1e-6);
    }

    @Test
    public void test02FindMinimum() {
        System.out.format("%nOutput from test problem 2 (2D unit circle calculation)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[0] * x[0] - x[1] * x[1];
                return x[0] * x[1];
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.findMinimum(calcfc, 2, 1, x, options.rhobeg, options.rhoend, 0, 3000);
        System.out.println(Arrays.toString(x));
        double v = calcfc.compute(2, 1, x, new double[1]);
        System.out.println("val: " + v);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), -Math.sqrt(0.5) }, x, 1.0e-5);
    }
}