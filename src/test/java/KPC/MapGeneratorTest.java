package KPC;

import junit.framework.TestCase;
import org.junit.Before;

import java.util.Arrays;
import java.util.logging.Level;

import static org.junit.Assert.assertArrayEquals;

public class MapGeneratorTest extends TestCase {

    @Before
    public void setUp() {
        Constants.LOGGER.setLevel(Level.SEVERE);
    }

    public void testChooseNext() {
        double[] pdf = {.2, .2, .3, .3};
        double[] counts = new double[4];
        int n = 100000;
        for (int i = 0; i < n; i++) {
            int a = MapGenerator.chooseState(pdf);
            assertTrue(a >= 0 && a < 4);
            counts[a] += 1.0/n;
        }
        assertArrayEquals(pdf, counts, .05);
        System.out.println(Arrays.toString(counts));
    }

    public void testRandomMAP() {
        Constants.LOGGER.setLevel(Level.INFO);
        int maxMaps = 6;
        for (int i = 2; i <= Math.pow(2, maxMaps); i *= 2) {
            MAP m = MapGenerator.randomMAP(i);
            assertTrue(m.isFeasible());
        }
    }

    public void testRandomSMP() {
        Constants.LOGGER.setLevel(Level.INFO);
        int maxMaps = 6;
        for (int i = 2; i <= maxMaps; i++) {
            SMP sm = MapGenerator.randomSMP(i);
            assertFalse(sm.isFeasible());
            assertTrue(sm.isValid());

        }

    }

    public void testMAPSample() {
        int n = 1000000;
        int iter = 1;
        int maxSize = 64;
        for (int nStates = 2; nStates <= maxSize; nStates *= 2) {
            for (int j = 0; j < iter; j++) {
                MAP m = MapGenerator.randomMAP(nStates);
                Sample s = MapGenerator.sample(m, n);
                double[] e = new double[3];
                for (int i = 0; i < n; i++) {
                    double v = s.iats[i];
                    e[0] += v / n;
                    e[1] += v * v / n;
                    e[2] += v * v * v / n;
                }
                double[] moments = m.getMoments();
                for (int i = 0; i < 3; i++) {
                    assertTrue(Math.abs((moments[i] - e[i]) / moments[i]) < .05);
                }
            }
        }
    }

    public void testSMPSample() {
        int n = 100000;
        int iter = 1;
        int maxMaps = 6;
        for (int nMaps = 1; nMaps <= maxMaps; nMaps++) {
            for (int j = 0; j < iter; j++) {
                SMP sm = MapGenerator.randomSMP(nMaps);
                Sample s = MapGenerator.sample(sm, n);
                Trace t = new Trace(s.getIATs());
                double[] traceAC = t.getAc();
                double[] mapAC = sm.getAcf(t.getAcLags());
                double[] e = new double[3];
                for (int i = 0; i < n; i++) {
                    double v = s.iats[i];
                    e[0] += v / n;
                    e[1] += v * v / n;
                    e[2] += v * v * v / n;
                }
                double[] moments = sm.getMoments();
                for (int i = 0; i < 3; i++) {
                    assertTrue((Math.abs((moments[i] - e[i]) / e[i]) < .05) || Math.abs(moments[i] - e[i]) < Constants.ZERO);
                }
                for (int i = 0; i < mapAC.length; i++) {
                    double diff = Math.abs(traceAC[i] - mapAC[i]);

                    if (diff > .05) {
                        System.out.println("Lag " + i+1 + "off by " + diff);
                    }
                }
            }
        }
    }
}




// --------------------------------------------------------------------------------------------------



//        MAP m = new MAP(new DMatrixRMaj(new double[][] {{-1.7942550992641235, 0.2895331696644361},
//                                                        {0.0011393773153844933, -0.511960658474298}}),
//                        new DMatrixRMaj(new double[][] {{0.8141749596909594, 0.6905469699087279},
//                                                        {0.5079283174504905, 0.0028929637084229576}}));