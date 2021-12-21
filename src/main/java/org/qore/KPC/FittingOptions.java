package org.qore.KPC;

/**
 * Enum of the order selection criterion to use (Only BIC implemented).
 */
enum OSCriterion {BIC, AIC, MES}

/**
 * Encapsulating options class for fitting process.
 */
public class FittingOptions {

    /**
     * Only fit AC values (not BC) and generate moments automatically.
     */
    boolean onlyAc = true;

    /**
     * Allow semi-markov models in fitting process.
     */
    boolean allowSM = false;

    /**
     * Employ verbose fitting
     */
    boolean verbose = true;

    /**
     * Maximum number of MAPs to consider during order selection.
     */
    int possibleOrders = 6;

    /**
     * log2 of number of states of final MAP, or alternatively the
     * number of MAP(2)s to use in KPC. Default 0 means use order
     * selection to determine automatically.
     */
    int numMAPs = 0;

    /**
     * The number of AC lags to fit (logarithmically spaced)
     */
    int nLags = 500;

    /**
     * Number of BC lags values to fit (logarithmically spaced)
     */
    int nBCLags = 5;

    /**
     * Number of runs of AC fitting.
     */
    int maxRunsAC = 30;

    /**
     * Max number of optimization evaluations for a single AC fitting.
     * If zero, is set to unlimited.
     */
    int maxEvalsAC = 3000;

    /**
     * Max number of optimization iterations for a single AC fitting.
     */
    int maxIterAC = 300;

    /**
     * Number of AC fitting results considered in BC fitting.
     */
    int maxResAC = 10;

    /**
     * Number of BC fitting runs for each AC fit.
     */
    int maxRunsBC = 5;

    /**
     * Max number of optimization iterations for a single BC fitting.
     */
    int maxIterBC = 30;

    /**
     * AC fit hyper-parameter for jump numStates.
     */
    double sigma = 0.5;

    /**
     * BC fit hyper-parameter.
     */
    double rhobeg = .05;
    /**
     * BC fit hyper-parameter.
     */
    double rhoend = 1.0e-8;

    /**
     * level of print for BC fitting.
     */
    int iprint = 0;

    /**
     * Order selection comparison criterion used.
     */
    OSCriterion osCriterion = OSCriterion.BIC;

    /**
     * Constructor using all default arguments
     */
    public FittingOptions() {}

    /**
     * Constructor specifying number of MAP(2)s to use in KPC.
     * @param numMAPs
     */
    public FittingOptions(int numMAPs) {
        this.numMAPs = numMAPs;
    }

    public void setPossibleOrders(int n) {
        possibleOrders = n;
    }
    public void setOnlyAc(boolean b) {
        onlyAc = b;
    }
    public void setAllowSM(boolean b) {
        allowSM = b;
    }
    public void setVerbose(boolean b) {
        verbose = b;
    }
    public void setNumMAPs(int n) {
        numMAPs = n;
    }
    public void setMaxRunsAC(int n) {
        maxRunsAC = n;
    }
    public void setMaxEvalsAC(int n) {
        maxEvalsAC = n;
    }
    public void setMaxIterAC(int n) {
        maxIterAC = n;
    }
    public void setMaxResAC(int n) {
        maxResAC = n;
    }
    public void setMaxRunsBC(int n) {
        maxRunsBC = n;
    }
    public void setMaxIterBC(int n) {
        maxIterBC = n;
    }
    public void setSigma(double d) {
        sigma = d;
    }
}
