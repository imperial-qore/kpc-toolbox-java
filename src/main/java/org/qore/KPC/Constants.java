package org.qore.KPC;

import java.util.logging.Level;
import java.util.logging.Logger;

public class Constants {
    protected static final Logger LOGGER = Logger.getLogger("KPC");
    protected static final double CONSTRAINT_TOL = 1.0e-10;
    protected static final double FEASIBLE_TOL = 1.0e-5;
    protected static final double ZERO = 1.0e-10;
    protected static final int MAX_MAP_SIZE = 128;
    public static void setLogger(Level l) {
        LOGGER.setLevel(l);
    }
}
