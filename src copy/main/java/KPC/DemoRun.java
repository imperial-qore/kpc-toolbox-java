package KPC;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;

/**
 * Executable class that can runs a demo of the MAP fitter on two datasets.
 */
public class DemoRun {

    /**
     * Runs a fit on the Morristown Bellcore Aug 89
     * internet traffic dataset and PKT-UDP dataset.
     * @param args unused
     */
    public static void main(String[] args) {
        Constants.LOGGER.setLevel(Level.SEVERE);
        Constants.LOGGER.addHandler(new ConsoleHandler());
        try {
            BCAug();
//            DEC_UDP();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Demo on BC-Aug89 dataset.
     * @throws IOException if trace data file can not be found.
     */
    public static void BCAug() throws IOException {
        Path path = Paths.get("data/bcaug.csv");
        Trace t = new Trace(path.toAbsolutePath().toString());
        FittingOptions opt = new FittingOptions();
        TraceFitter fitter = new TraceFitter(t, opt);
        List<MAP> maps = fitter.fit();
        fitter.displayResults(maps);
    }

    /**
     * Demo of DEC-PKT-UDP trace dataset.
     * @throws IOException if trace data file can not be found.
     */
    public static void DEC_UDP() throws IOException {
        Path path = Paths.get("data/DEC-PKT-1-UDP.csv");
        FittingOptions opt = new FittingOptions();
        opt.maxEvalsAC = 10000;
        opt.maxIterAC = 550;
        opt.sigma = 0.35;
        Trace t = new Trace(path.toAbsolutePath().toString());
        TraceFitter fitter = new TraceFitter(t, opt);
        List<MAP> maps = fitter.fit();
        fitter.displayResults(maps);
    }
}