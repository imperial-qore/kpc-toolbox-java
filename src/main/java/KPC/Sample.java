package KPC;

import java.util.Arrays;

/**
 * Class holds the necessary inforamtion of a sample from a
 * Markov process, including the IATs, initial states, and
 * absorbing states.
 */
public class Sample {

    double[] iats;
    int[] initials;
    int[] absorbings;

    /**
     * Constructs empty object to hold data.
     * @param size Number of samples to hold.
     */
    public Sample(int size) {
        iats = new double[size];
        initials = new int[size];
        absorbings = new int[size];
    }

    public void print() {
        if (iats.length < 20) {
            System.out.println(Arrays.toString(iats));
        } else {
            System.out.println("Printing first 20 elements");
            System.out.println(Arrays.toString(Arrays.copyOfRange(iats, 0, 20)));
        }
    }

    public double[] getIATs() {
        return iats;
    }
}
