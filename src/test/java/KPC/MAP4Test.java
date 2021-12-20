package org.qore.KPC;

import org.junit.Test;
import static org.junit.Assert.*;

public class MAP4Test {

    @Test
    public void testSMP() {
        double[] scv = {0.7331352969761109, 1.0996892434159826, 1.7801877459402649, 2.1121281262981246};
        double[] g2 = {0.8724712592108149, 0.9999384515996077, 0.9581275338000508, 0.6173050401158501};
        double[] e1 = {0.25448872001266254, 0.32364512668396117, 0.2589819039012036, 0.2271961072548022};
        double[] e3 = {0.07049627834106463, 0.22443464740922772, 0.21793351302257924, 0.16755825838324634};
        SMP m = KPC.composeSMMAP(e1, e3, scv, g2, 4, false);
        assertFalse(m.isFeasible());
        assertTrue(m.isValid());
    }

    @Test
    public void testSM2() {
        double[] e1 = {1.7645535017942162, 1.4342666214708872, 0.8474552688202979, 1.2481517517994611};
        double[] e3 = {31.95766599413806, 17.46042811386062, 3.8105615679320453, 14.434946209600717};
        double[] scv = {0.9736431861502529, 1.98824779754105, 1.031981262788616, 1.1923028553595123};
        double[] g2 = {0.3978471167952484, 0.017505387585007828, 0.1720056358096722, -0.07683280006688202};
        SMP m = KPC.composeSMMAP(e1, e3, scv, g2, 4, true);
        assertFalse(m.isFeasible());
        assertTrue(m.isValid());
    }

}
