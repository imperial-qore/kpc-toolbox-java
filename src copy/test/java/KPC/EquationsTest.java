package KPC;

import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.assertArrayEquals;


public class EquationsTest extends TestCase {

    @Test
    public void testCumSumInteger() {
        int[] input = {1, 2, 3, 4, 5, 6};
        int[] real = {1, 3, 6, 10, 15, 21};
        int[] output = Equations.cumSum(input);
        assertArrayEquals(output, real);
    }

    @Test
    public void testCumSumDouble() {
        double[] input = {.1, .2, .3, .4, .5, .6};
        double[] real = {.1, .3, .6, 1.0, 1.5, 2.1};
        double[] output = Equations.cumSum(input);
        assertArrayEquals(output, real, 1e-10);
    }

    @Test
    public void testLogspace() {
        double[] real = {1.000000000000000e+00,
                1.668100537200059e+00,
                2.782559402207124e+00,
                4.641588833612778e+00,
                7.742636826811269e+00,
                1.291549665014884e+01,
                2.154434690031884e+01,
                3.593813663804627e+01,
                5.994842503189409e+01,
                1.000000000000000e+02};
        double[] test = Equations.logspace(1, 100, 10);
        assertArrayEquals(real, test, 1e-6);
    }

    @Test
    public void testLogspacei() {
        int[] real = new int[]{1,2,3,5,8,13,22,36,60,100};
        int[] test = Equations.logspacei(1, 100,10);
        assertArrayEquals(real, test);

    }
}