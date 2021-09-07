package KPC;

import junit.framework.TestCase;
import org.ejml.data.DMatrixRMaj;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;

import static org.junit.Assert.assertArrayEquals;

public class TraceFitter2Test extends TestCase {
    Trace t;
    FittingOptions options;
    TraceFitter fit;
    double[] scv = {1.310198213431212e+00};
    double[] gamma = {9.760638375345613e-01};
    double[] ac_data = {1.310198213431212e+00, 9.760638375345613e-01};
    double[] e2 = {2.281861275059570e-05};
    double[] e1 = {3.142823538823539e-03};
    double[] e3 = {2.486262256309365e-07};

    @Before
    public void setUp() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/bcaug.csv";
        t = new Trace(path);
        options = new FittingOptions();
        options.numMAPs = 1;
        fit = new TraceFitter(t, options);
        Constants.LOGGER.setLevel(Level.WARNING);
    }

    @Test
    public void testFitACExecutes() {
        fit.fitAC();
    }
    @Test
    public void testFitBCExecutes() {
        ACFitResult prev = fit.fitAC();
        fit.fitBC(prev.scv, prev.gamma);
    }

    @Test
    public void testACFValue2() {
//        OldACObjFunc f = new OldACObjFunc(t, options.numMAPs);
        ACObjFunc f = new ACObjFunc(t, options.numMAPs);
        double v = f.value(ac_data);
        double real = 6.849301658340621e+00;
        assertEquals(v, real, 1e-6);

    }

    /**
     * Tests that a given call to compose with the parameters from an
     * ac fit result returns the correct MAP
     * values obtained from matlab run
     * Does not check infeasible branch
     */


    /**
     * Tests that the BC Fit objective function returns the correct value
     * on MAP(2)s. Values obtained form MATLAB test run
     */
    @Test
    public void testBCValue2() {
        BCObjFunc f = new BCObjFunc(t, options.numMAPs, scv, gamma, false, options.verbose);
        double v = f.value(e1, e3);
        System.out.println("value = " + v);
        double answer = 4.247673974007398e-01;
        assertEquals(v, answer, 1e-6);
    }

    /**
     * Tests that the constraints are generated correctly for a single
     * execution of the optimization problem with scv and gamma values
     * received from an acfit result.
     */
    @Test
    public void testBCConstraints2() {
        BCObjFunc f = new BCObjFunc(t, options.numMAPs, scv, gamma, false, options.verbose);
        double[] cons = new double[4];
        double[] x = new double[2];
        x[0] = e1[0];
        x[1] = e3[0];
        f.compute(2, 4, x, cons);
        double[] real = new double[] {3.063933059455678e-06,
                1.128657966360576e-10,
                1.876331541783999e+00,
                -3.763315417839987e-01};
        System.out.println(Arrays.toString(cons));
        assertArrayEquals(cons, real, 1e-6);
    }

    @Test
    public void testEvaluate2() {
        DMatrixRMaj d0 = new DMatrixRMaj(new double[][]
                {{-1.066650758528454e+00, 8.440228187684076e-01},
                 {4.554557043385951e-01, -1.163239251432586e+00}});
        DMatrixRMaj d1 = new DMatrixRMaj(new double[][]
                {{2.079251269846806e-02, 2.018354270615780e-01},
                 {4.597609231526781e-01, 2.480226239413128e-01}});
        MAP m = new MAP(d0, d1);

        double[] truth = new double[] {1.312927737670788e+01, 1.695537945656511e+08};
        double[] results = KPC.evaluate(t, m, true);
        assertEquals(results[0],truth[0], 0.00001);
        assertEquals(results[1]/100000000,truth[1]/100000000, 0.00001);
    }

    @Test
    public void testBIC() {
        System.out.println(t.acLags[0]);
        fit.findBIC(new int[] {1,2,3,4,5,6,7});
    }
}