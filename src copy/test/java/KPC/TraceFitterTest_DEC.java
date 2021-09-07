package KPC;

import junit.framework.TestCase;
import org.junit.Before;

import java.io.IOException;
import java.util.logging.Level;

public class TraceFitterTest_DEC extends TestCase {
    Trace t;
    FittingOptions options;
    TraceFitter fit;

    @Before
    public void setUp() throws IOException {
        String path = "/Users/kodiakconrad/Imperial/IndividualProject/KPC_Toolbox/data/DEC-PKT-1-UDP.csv";
        t = new Trace(path);
        options = new FittingOptions(5);
        fit = new TraceFitter(t, options);
        Constants.LOGGER.setLevel(Level.SEVERE);
    }

    public void testFitBC() {
        double[] e1 = {0.37453321375040155, 0.32453321375040156, 0.32453321375040156, 0.32453321375040156, 0.32453321375040156};
        double[] e3 = {0.12685329200996345, 0.3350104705379787, 0.31081047277027724, 0.23121392689444184, 0.33009627817338033};
        double[] scv = {0.5723884706878561, 1.5552780836214652, 1.4612560779775272, 1.122833705287266, 1.536467450159754};
        double[] gamma = {0.4979850978887004, 0.8059122014932878, 0.9999922905518491, 0.20441436683906442, 0.6120482871081027};
        BCObjFunc func = new BCObjFunc(t, options.numMAPs, scv, gamma, false, true);
        func.value(e1, e3);
        System.out.println("obj func value: " + func.value(e1, e3));
    }
}