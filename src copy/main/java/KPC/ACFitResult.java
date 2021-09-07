package KPC;

import java.util.Comparator;

/**
 * Class encapsulates the relevant data of from an AC optimization run,
 * including the best parameter values found and their OF value
 */
public class ACFitResult implements Comparable<ACFitResult> {
    double value;
    double[] scv;
    double[] gamma;

    /**
     * Generates a new instance of a the optimally fit values
     *
     * @param scv SCVs of the composing MAPs
     * @param gamma Gamma(2)s of the composing MAPs
     * @param value Objective function closeness of fit
     */
    public ACFitResult(double[] scv, double[] gamma, double value)  {
        this.value = value;
        this.scv = scv;
        this.gamma = gamma;
    }

    @Override
    public int compareTo(ACFitResult o) {
        return Double.compare(value, o.value);
    }
}

/**
 * Compares two ACFitResults by comparing their OF values.
 */
class ACFitComparator implements Comparator<ACFitResult> {
    @Override
    public int compare(ACFitResult o1, ACFitResult o2) {
        return Double.compare(o1.value, o2.value);
    }
}
