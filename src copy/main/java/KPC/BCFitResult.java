package KPC;

import java.util.Comparator;

/**
 * Class encapsulates the relevant data of from an BC optimization run,
 * including the best parameter values found and their OF value
 */
public class BCFitResult {
    double value;
    double[] e1;
    double[] e2;
    double[] e3;

    /**
     * Generates a new instance of a BC fit result
     *
     * @param e1 Fit e1 values of composing MAP(2)s
     * @param e2 Fit e2 values (implied from e1, and SCV, included for ease of access)
     * @param e3 Fit e3 values of composing MAP(2)s
     * @param value BC objective function closeness of fit
     */
    public BCFitResult(double[] e1, double[] e2, double[] e3, double value) {
        this.value = value;
        this.e1 = e1;
        this.e2 = e2;
        this.e3 = e3;
    }
}

/**
 * Compares two BCFitResults by comparing their OF values.
 */
class BCFitComparator implements Comparator<BCFitResult> {
    @Override
    public int compare(BCFitResult o1, BCFitResult o2) {
        return Double.compare(o1.value, o2.value);
    }
}
