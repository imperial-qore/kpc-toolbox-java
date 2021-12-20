package org.qore.KPC;

import org.ejml.data.DMatrixRMaj;
import org.junit.Test;

import java.util.Set;

import static org.junit.Assert.*;

public class Fittings2Test {
    double tol = 1e-4;

    /**
     *  MAP2 that should not return any errors, values taken from MATLAB code
     *  Tests the hyper-exponential fitting, ensuring arithmetic is correct
     */
    @Test
    public void testMap2Hyper() {
        double e1 = 0.246448715775199;
        double e2 = 0.149441138009822;
        double e3 = -1;
        double g2 = 0.680709874224242;
        double scv=(e2-e1*e1)/(e1*e1);

        MAP m = Fittings2.map(e1, e2, e3, scv, g2);
        double[] real_d0 = new double[] {-3.296067645038661, 0, 0, -1.139813454069052e+03};
        double[] real_d1 = new double[] {3.097971217888393, 1.980964271502682e-01,
                2.954274375030541e+02, 8.443860165659984e+02};
        assertArrayEquals(m.D0.data, real_d0, tol);
        assertArrayEquals(m.D1.data, real_d1, tol);
    }

    /**
     *  MAP2 that should not return any errors, values taken from MATLAB run.
     *  Tests the hypo-exponential fitting, ensuring arithmetic is correct
     */
    @Test
    public void testMap2Hypo() {
        double e1 = 0.236771812825797;
        double e2 = 0.097218124518012;
        double e3 = 0.055602196080333;
        double g2 = 0;
        double scv = (e2 - e1 * e1) / (e1 * e1);

        MAP m = Fittings2.map(e1, e2, e3, scv, g2);
        double[] real_d0 = new double[] {-6.646677372312877, 4.846548077526329, 0, -6.646966881989579};
        double[] real_d1 = new double[] {1.416515309615199, 0.383613985171349,
                                        5.230474487645011, 1.416492394344568};
        assertArrayEquals(m.D0.data, real_d0, tol);
        assertArrayEquals(m.D1.data, real_d1, tol);
    }

    // TODO: tests for MAP2s with values that don't pass.

    /**
     * General case MAP2 fitting based on traits. Values taken from successful MATLAB run
     */
    @Test
    public void testFeasblockBad() {
        double e1 = 1;
        double e3 = 11.460804888834867;
        double e2 = 1.146080488883487e+01;
        double g2 = 3.396097343698996e-01;

        MAP m = Fittings2.feasblock(e1, e2, e3, g2);

        double[] real_d0 = new double[] {8.511008891136352e-01, 0, 0, -3.149657053634182e-01};
        double[] real_d1 = new double[] {-5.700715179937511e-01,-2.810293711198841e-01,
                                         1.040001429146599e-01, 2.109655624487584e-01};
        assertArrayEquals(m.D0.data, real_d0, tol);
        assertArrayEquals(m.D1.data, real_d1, tol);
//        System.out.println(m.D0.toString());
//        System.out.println(m.D1.toString());

        assertFalse(m.isFeasible(tol));
    }

    @Test
    public void testFeasblockGood() {
        double e1 = 1;
        double e3 = 1.146080488883487e+01;
        double e2 = 2.754984031255612e+00;
        double g2 = 3.396097343698996e-01;

        MAP m = Fittings2.feasblock(e1, e2, e3, g2);

        double[] real_d0 = new double[] {-7.210640882748437e-01, 0, 0, -4.138537690048018e+01};
        double[] real_d1 = new double[] {5.858840913169368e-01, 1.351799969579069e-01,
                1.957186516684821e+01, 2.181351173363197e+01};
        assertArrayEquals(m.D0.data, real_d0, tol);
        assertArrayEquals(m.D1.data, real_d1, tol);

        assertTrue(m.isFeasible(tol));

    }

    @Test
    public void testFeasblock2() {
        double e1 = 5.209770343231517e-01;
        double e3 = 1.368945322715493e+00;
        double e2 = 0.6895233052594170;
        double g2 = 4.440024098932075e-01;

        MAP m = Fittings2.feasblock(e1, e2, e3, g2);

        double[] real_d0 = new double[] {-1.532072833279684e+00, 0, 0, -1.733557019809538e+04};
        double[] real_d1 = new double[] {1.333457383446309e+00, 1.986154498333752e-01,
                                            5.000317237890302e+03, 1.233525296020508e+04};
        assertArrayEquals(m.D0.data, real_d0, tol);
        assertArrayEquals(m.D1.data, real_d1, tol);
        assertTrue(m.isFeasible(tol));
    }


    @Test
    public void testMmpp() {
        double e1 = 3.142823538823539e-03;
        double e2 = 2.281860358311090e-05;
        double e3 = 2.488184018487869e-07;
        double scv = 1.310197285298235e+00;
        double g2 = 9.760638751446140e-01;
        double[] d0 = {-2.751404162444633e+02,     9.370166376635052e-01,
                7.198498528323075e+02,    -3.482645020462810e+04};
        double[] d1 = {2.742033996067998e+02,                         0,
                0,     3.410660035179579e+04};

        MAP m = Fittings2.mmpp(e1, e2, e3, scv, g2);

        assertArrayEquals(m.D0.data, d0, tol);
        assertArrayEquals(m.D1.data, d1, tol);

    }

    @Test
    public void testSM2() {
        double[] e = {7.913074756311242e-01,1.221457299114743e+00,2.788650030325246e+00};
        DMatrixRMaj d0 = Fittings2.d0Fit(e[0], e[1], e[2]);
        double[] realD0 = new double[] {-1.340945569595714e+00, 3.475293876570082e-01, 0, -2.210061741986475e+00};
        assertArrayEquals(d0.data, realD0, tol);
        double g2 = 3.730323196053044e-01;
        MAP m = Fittings2.sm(e[0], e[1], e[2], g2);
        double[] realD1 = new double[] {1.014273361394829e+00, -2.085717945612298e-02, 1.143627865194707e+00, 1.066433876791769e+00};
        double[] realP = new double[] {8.904965349564904e-01, 1.095034650435096e-01, 5.174642153511860e-01, 4.825357846488140e-01};
        assertArrayEquals(m.D0.data, realD0, tol);
        assertArrayEquals(m.D1.data, realD1, tol);
        assertArrayEquals(m.P.data, realP, tol);
    }

    @Test
    public void testErlang() {
        int k = 4;
        double r = Math.random();
        MAP m = Fittings2.erlang(r, k);
        double mu = k/r;
        double[] realD0 = {-mu,mu,0,0, 0,-mu,mu,0,0,0,-mu, mu,0,0,0,-mu};
        double[] realD1 = {0,0,0,0,0,0,0,0,0,0,0,0,mu,0,0,0};
        assertArrayEquals(m.D0.data, realD0, 1e-6);
        assertArrayEquals(m.D1.data, realD1, 1e-6);
    }

    @Test
    public void testMap2exp() {
        double r = Math.random();
        MAP m = Fittings2.mapExp(r);
        double mu = 1/r;
        assertEquals(m.D0.data[0], -mu, 1e-6);
        assertEquals(m.D1.data[0], mu, 1e-6);
        assertEquals(m.D0.getNumElements(), 1);
        assertEquals(m.D1.getNumElements(), 1);
    }
}