package KPC;

import junit.framework.TestCase;
import org.ejml.data.DMatrixRMaj;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

public class MAP2Test {
    MAP map;
    DMatrixRMaj real_d0;
    DMatrixRMaj real_d1;
    DMatrixRMaj real_p;
    DMatrixRMaj real_q;
    double tol = 1e-6;
    Logger logger;

    @Before
    public void initialize() {
        double[][] d0 = {{-1.075103950256064, 0.248102680895525},
                         {0.534970435097334, -1.945304577679050}};
        double[][] d1 = {{0.507576981795755, 0.319424287564784},
                         {0.846333041249930, 0.564001101331786}};
        real_d0 = new DMatrixRMaj(d0);
        real_d1 = new DMatrixRMaj(d1);
        double[][] q = {{-0.5675269684603097,0.5675269684603097},
                        {1.381303476347264, -1.381303476347264}};
        double[][] p = {{0.6113153420397894,0.3886846579602106},
                        {0.6031799283393360,0.3968200716606640}};
        real_p = new DMatrixRMaj(p);
        real_q = new DMatrixRMaj(q);

        map = new MAP(real_d0, real_d1);
    }

    @Before
    public void setUpLogger() {
        logger = Logger.getLogger("KPC");
        logger.setUseParentHandlers(false);
        logger.setLevel(Level.WARNING);
    }

    @Test
    public void MAPConstructor() {

        assertArrayEquals(map.D0.data, real_d0.data, tol);
        assertArrayEquals(map.D1.data, real_d1.data, tol);
        assertArrayEquals(map.Q.data, real_q.data, tol);
        assertArrayEquals(map.P.data, real_p.data, tol);
        Assert.assertEquals(map.numStates, 2);
    }

    @Test
    public void MAP_isFeasible() {
        Assert.assertTrue(map.isFeasible(tol));
    }

    @Test
    public void infeasible_d0() {
        double[][] d0 = {{-1, 1}, {1, 2}};
        double[][] d1 = {{0, 0}, {1, 0}};
        MAP badMap = new MAP(new DMatrixRMaj(d0), new DMatrixRMaj(d1));
        Assert.assertFalse(badMap.feasibleD0());
        // check logs for correct reasoning
    }

    @Test
    public void infeasible_d1() {
        double[][] d0 = {{-1, 1}, {1, -2}};
        double[][] d1 = {{0, 0}, {1, -1}};
        MAP badMap = new MAP(new DMatrixRMaj(d0), new DMatrixRMaj(d1));
        Assert.assertFalse(badMap.feasibleD1());
        // check logs for correct reasoning

    }
    @Test
    public void non_stochastic_P() {
        double[][] d0 = {{-1, 1}, {2, -2}};
        double[][] d1 = {{0, 0}, {1, 0}};
        MAP badMap = new MAP(new DMatrixRMaj(d0), new DMatrixRMaj(d1));
        Assert.assertFalse(badMap.stochasticP(tol));
        // check logs for correct reasoning

    }
    @Test
    public void non_stochastic_Q() {
        double[][] d0 = {{-1, 1}, {1, -2}};
        double[][] d1 = {{0, 0}, {2, 0}};
        MAP badMap = new MAP(new DMatrixRMaj(d0), new DMatrixRMaj(d1));
        Assert.assertFalse(badMap.stochasticQ(tol));
        // check logs for correct reasoning
    }

    @Test
    public void canGetLambda() {
        double real_lam = 0.996876046319468;
        Assert.assertEquals(map.getLambda(), real_lam, tol);
    }

    @Test
    public void canGetMean() {
        double real_mean = 1.003133743349603;
        Assert.assertEquals(map.getMean(), real_mean, tol);
    }

    @Test
    public void canGetVariance() {
        double real_var = 1.066421432939353;
        Assert.assertEquals(map.getVariance(), real_var, tol);
    }

    @Test
    public void canGetSCV() {
        double real_scv = 1.059768937917609;
        Assert.assertEquals(map.getSCV(), real_scv, tol);
    }

    @Test
    public void moments2() {
        double[] moments = map.getMoments(new int[]{1, 2, 3});
        double[] real_moments = {1.003133743349603, 2.072698739985940, 6.515820704802493};
        assertArrayEquals(moments, real_moments, tol);
    }
    @Test
    public void lags2() {
        double[] lags = map.getAcf(new int[]{1, 2, 3, 4, 5});
        double[] real_lags = {0.0002294108739174502, 0.000001866352366783731,
                0.00000001518354858390565, 0.0000000001235243649141085,
                0.000000000001004866143098112};
        assertArrayEquals(lags, real_lags, tol);

    }
    @Test
    public void joint2() {
        double real_a = 1.009432734664093;
        double real_b = 1.009430721865400;
        double bic_a = map.getJoint(new int[]{1,2,3}, new int[]{1,1,1});
        double bic_b = map.getJoint(new int[]{1,10,100}, new int[]{1,1,1});
        Assert.assertEquals(bic_a, real_a, tol);
        Assert.assertEquals(bic_b, real_b, tol);

    }

    @Test
    public void bc2() {
        double[] real = {1.009432734664093, 1.009430721865400};
        double[] mapBC = map.getBCs(new int[][] {{1,2,3}, {1,10,100}});
        assertArrayEquals(real, mapBC, tol);
    }

    @Test
    public void test_ctmc_easy() {
        double[][] d0 = new double[][] {{-1, 1}, {1, -2}};
        double[][] d1 = new double[][] {{0, 0}, {1, 0}};

        MAP m2 = new MAP(new DMatrixRMaj(d0), new DMatrixRMaj(d1));
        System.out.println("pi: " + m2.ctmc(false));
    }

    @Test
    public void test_ctmc_hard() {
        initialize();
        DMatrixRMaj pi = map.ctmc(false);
        double[] real_pi = {0.7087858669427001, 0.2912141330572999};
        assertArrayEquals(real_pi, pi.data, tol);
    }

    @Test
    public void test_dtmc() {
        initialize();
        double[] real_pi = {0.6081272954705266, 0.3918727045294734};
        assertArrayEquals(real_pi, map.dtmc().data, tol);

    }

    @Test
    public void test_matrixPower() {
        initialize();
        DMatrixRMaj pow = MatrixOps.matrixPower(map.P, 100000);
        DMatrixRMaj result = new DMatrixRMaj(new double[] {6.081272953144892e-01,     3.918727046894990e-01,
                6.081272953144891e-01,     3.918727046894989e-01});
        assertArrayEquals(pow.data, result.data, tol);
    }

    @Test
    public void testNormalize() {
        MAP m = new MAP(new DMatrixRMaj(new double[][] {{0.12178954126852182, 0.8792341452085077},
                                                        {0.5501946937270578, 0.15827470035372126}}),
                        new DMatrixRMaj(new double[][] {{0.330546520956415, 0.36025793647074156},
                                                        {0.24376603735923152, 0.7927048900743877}}));
        Assert.assertFalse(m.isFeasible());
        m.normalize();
        m.P.print();
        Assert.assertTrue(m.isFeasible());
    }



}