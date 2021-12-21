package org.qore.KPC;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Scanner;

/**
 * Class for representing an Inter-Arrival Time Trace to be fit.
 */
public class Trace {
    double[] moments;
    /**
     * AC values for all lags specified in acLags
     */
    double[] ac;
    /**
     * BC values for all lags specified in bcLags
     */
    double[] bc;
    /**
     * AC values for all lags specified in acLags
     */
    double [] acFull; // index-0 = lag-0 = 1;

    int[] acLags;
    int[] bcLagValues;
    int[][] bcLags;

    int length = 999999; // max length of file
    int nLags = 500; // 500 logarithmically spaced correlations
    int nBCValues = 5; // 5^2 different combos of bicorrelations
    int maxMoments = 10;

    /**
     * Creates new trace object by reading in values and calculating characteristics.
     * Does not keep data stored.
     *
     * @param path Path to the trace data.
     * @throws IOException If can not find or read from file specified.
     */
    public Trace(String path) throws IOException {
        DMatrixRMaj data = getData(path);
        init(data);
    }

    /**
     * Construct a trace object from an array of IATs.
     * @param data Vector of doubles representing IATs.
     */
    public Trace(double[] data) {
        length = data.length;
        init(new DMatrixRMaj(data));
    }

    /**
     * Trace constructor allows user to specify the characterstics to fit.
     * @param path Path to trace file.
     * @param nLags Number of AC lags to capture.
     * @param nBCValues Number of BC lags valeus to capture.
     * @param maxMoments Number of moments to capture.
     * @throws IOException Can not find or read file from path.
     */
    public Trace(String path, int nLags, int nBCValues, int maxMoments) throws IOException {
        this.nLags = nLags;
        this.nBCValues = nBCValues;
        this.maxMoments = maxMoments;
        DMatrixRMaj data = getData(path);
        init(data);
    }

    /**
     * Calculate characteristics of a trace from the given data.
     * @param data Vector of doubles representing IATs.
     */
    private void init(DMatrixRMaj data) {
        int nMinSupportAC = 10; // minimum number of pts to estimate gamma
        acLags = Equations.logspacei(1, (int) Math.ceil((double) length/nMinSupportAC), nLags);
        bcLagValues = Equations.logspacei(1, acLags[acLags.length-1], nBCValues);
        getBCLags();

        moments = new double[maxMoments];
        ac = new double[acLags.length];
        acFull = new double[(int) Math.ceil( (double)length / nMinSupportAC)];
        bc = new double[nBCValues*nBCValues];
        generateMoments(data);
        generateAC(data);
        generateBC(data);
    }

    /**
     * Generates a vector of trace data by reading until end of file. Assumes
     * each line has one double value representing IAT.
     *
     * @param path Path to trace data file.
     * @return Vector of data points in DMatrix object.
     * @throws IOException If can not find or read from file specified.
     */
    public DMatrixRMaj getData(String path) throws IOException {
        double[] data = new double[length];
        try (FileInputStream input = new FileInputStream(path); Scanner reader = new Scanner(input)) {
            int i = 0;
            while (reader.hasNext() && i < length) {
                data[i++] = reader.nextDouble();
            }
            length = i;
            data = Arrays.copyOfRange(data, 0, i);
        }
        return new DMatrixRMaj(data);
    }

    /** Generates the moments of the Trace.
     *
     * @param data Data of IATS read from Trace file.
     */
    public void generateMoments(DMatrixRMaj data) {
        DMatrixRMaj result = data.copy();
        for (int i = 1; i <= maxMoments; i++) {
            if (i != 1) {
                CommonOps_DDRM.elementPower(data, i, result);
            }
            moments[i-1] = CommonOps_DDRM.elementSum(result) / result.getNumElements();
        }
    }

    /** Generates the autocorrelations of the Trace.
     *
     * @param data Data of IATS read from Trace file.
     */
    public void generateAC(DMatrixRMaj data) {
        double mean;
        if (moments[0] != 0) {
            mean = moments[0];
        } else {
            mean = CommonOps_DDRM.elementSum(data) / data.getNumElements();
        }
        int n = data.getNumElements();
        int l = Equations.nextPowerOfTwo(2*n + 1);

        DMatrixRMaj xAdj = new DMatrixRMaj();
        CommonOps_DDRM.subtract(data, mean, xAdj);
        double[] paddedData = Arrays.copyOf(xAdj.data, l);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] z = fft.transform(paddedData, TransformType.FORWARD);
        Complex[] norms = new Complex[l];
        for (int i = 0; i < l; i++) {
            norms[i] = z[i].multiply(z[i].conjugate());
        }
        Complex[] acov = fft.transform(norms, TransformType.INVERSE);
        double[] acorr = new double[n];
        double var = 1;
        for (int j = 0; j < n; j++) {
            if (j == 0) {
                var = acov[j].getReal()/n;
            }
            acorr[j] = acov[j].getReal() / (n*var);
        }
        acFull = acorr;
        for (int i = 0; i < acLags.length; i++) {
            ac[i] = acFull[acLags[i]];
        }
    }

    /** Generates the bicorrelations of the Trace.
     *
     * @param data Data of IATS read from Trace file.
     */
    public void generateBC(DMatrixRMaj data) {
        // Done assuming all orders = 1 (bicorrelations)
        // Can be abstracted to include other orders
        for (int i = 0; i < bcLags.length; i++) {
            int[] lags = Equations.cumSum(bcLags[i]);
            for (int j = 0; j < lags.length; j++) {
                lags[j]--;
            }
            int elementCount = data.getNumElements() - lags[2];
            DMatrixRMaj previous = CommonOps_DDRM.extract(data, 0, elementCount, 0, 1);
            for (int j = 1; j < lags.length; j++) {
                DMatrixRMaj dest = CommonOps_DDRM.extract(data, lags[j], elementCount + lags[j], 0, 1);
                // add in power operation here
                CommonOps_DDRM.elementMult(previous, dest);
            }
            bc[i] = CommonOps_DDRM.elementSum(previous) / elementCount;
        }
    }

    /**
     * Generates the bicorrelation trios form the determined
     * bicorrelation values of interest.
     */
    public void getBCLags() {
        assert bcLagValues != null;
        int nRows = bcLagValues.length*bcLagValues.length;
        bcLags = new int[nRows][3];
        int row = 0;
        for (int i : bcLagValues) {
            for (int j : bcLagValues) {
                bcLags[row++] = new int[]{1, i, j};
            }
        }
    }

    /**
     * Gets the AC values at the lags desired.
     *
     * @param lags List of lags of AC values desired.
     * @param n Number of lags desired.
     * @return
     */
    public double[] acAtLags(int[] lags, int n) {
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            try {
                int lag = lags[i];
                result[i] = acFull[lags[i]];
            } catch (Exception e) {
                int[] temp = Arrays.copyOfRange(lags, lags.length-50, lags.length);
                System.out.println(Arrays.toString(temp));
            }
        }
        return result;
    }

    /**
     * @return Autocorrelation lags used in trace.
     */
    public int[] getAcLags() {
        return acLags;
    }

    /**
     * @return Bicorrelation lags used in trace.
     */
    public int[][] getBcLags() {
        return bcLags;
    }

    /**
     * @return all AC values.
     */
    public double[] getAcFull() {
        return acFull;
    }

    /**
     * @return subset of AC values of trace (log spaced).
     */
    public double[] getAc() {
        return ac;
    }

    /**
     * @return BC values of trace.
     */
    public double[] getBc() {
        return bc;
    }

    /**
     * @param n Number of moments.
     * @return First n moments of a trace.
     */
    public double[] getMoments(int n) {
        if (n > maxMoments) {
            throw new IllegalArgumentException(String.format("Only %d moments available", maxMoments));
        }
        double[] result = new double[n];
        System.arraycopy(moments, 0, result, 0, n);
        return result;
    }
}