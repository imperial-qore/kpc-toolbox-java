package KPC;

import junit.framework.TestCase;
import org.ejml.data.DMatrixRMaj;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;


public class TraceTest extends TestCase {
    Trace t;

    // All hard coded numbers come from BC Aug network traffic data
    String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/bcaug.csv";
    double tol = .00001;

    @Before
    public void setUp() throws IOException {
        t = new Trace(path);
    }
    @Test
    public void testGetData() {
        try {
            DMatrixRMaj temp = t.getData(path);
            assertEquals(temp.getNumElements(), 999999);
            double d1 = .000168000000000000;
            double d2 = .002668000000000000;

            assertEquals(d1, temp.get(0), tol);
            assertEquals(d2, temp.get(1), tol);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void testMoments() {
        double[] real = new double[] {.003142823538823539, .00004171801994357199, .000002010425529820079};
        try {
            DMatrixRMaj data = t.getData(path);
            t.generateMoments(data);
            assertArrayEquals(real, t.getMoments(3), tol);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testACF() {
        double[] real_ac = new double[] {
                .2000484323902797,
                .1892741490928188,
                .1389506214955594,
                .1421333993032179,
                .1171267609012947,
                .1236777287082127,
                .1121206711095049,
                .1005101934892867,
                .1001901288451582,
                .09879628650406407
        };
        try {
            DMatrixRMaj data = t.getData(path);
            t.generateAC(data);
            for (int i = 0; i < 10; i++) {
                assertEquals(real_ac[i], t.ac[i], 1e-7);
            }
            System.out.println(Arrays.toString(t.ac));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testBC() {
        double[] real_bc = new double[] {
                1.2184472845989310e-07,
                7.7265562312809010e-08,
                5.8437705882148630e-08,
                5.4518075525906270e-08,
                5.5591955856097480e-08,
        };
        try {
            DMatrixRMaj data = t.getData(path);
            t.generateBC(data);
            for (int i = 0; i < 5; i++) {
                assertEquals(real_bc[i], t.bc[i], tol);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testDEC() throws IOException {
        path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/DEC-PKT-1-UDP.csv";
        t = new Trace(path);
    }


}