package org.qore.KPC;

import org.ejml.data.DMatrixRMaj;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class KPCTest {
    MAP a;
    MAP b;
    MAP c;
    double tol = .0000000001;

    @Before
    public void initialize() {
        DMatrixRMaj a_d0 = new DMatrixRMaj(new double[][] {{-1.397946868742791,0.280807532966513},
                {0.583702076963942, -1.344518468875679}});

        DMatrixRMaj a_d1 = new DMatrixRMaj(new double[][] {{0.670410898360615, 0.446728437415663},
                {0.595968746857703, 0.164847645054035}});

        DMatrixRMaj b_d0 = new DMatrixRMaj(new double[][] {{-1.994029211808941, 0.484244876766691},
                {0.392305839523815, -1.337624894612294}});

        DMatrixRMaj b_d1 = new DMatrixRMaj(new double[][] {{0.795791459640391, 0.713992875401859},
                {0.867062176050182, 0.078256879038297}});

        DMatrixRMaj c_d0 = new DMatrixRMaj(new double[][] {{-1.193312546793983, 0.761322847922017},
                {0.594321492404827, -2.045190433526966}});

        DMatrixRMaj c_d1 = new DMatrixRMaj(new double[][] {{0.141645466043449, 0.290344232828517},
                {0.635347792243281, 0.815521148878858}});

        a = new MAP(a_d0, a_d1);
        b = new MAP(b_d0, b_d1);
        c = new MAP(c_d0, c_d1);
    }

    @Test
    public void correctMAP() {
        MAP d = KPC.kroneckerCompose(a, b);
        DMatrixRMaj real_d0 = new DMatrixRMaj(new double[][]
                {{-2.787546892829965, 0.676948609180735, 0.559938423631229, -0.135979609196528},
                        {0.548422719951830, -1.869928532975663, -0.110162434965039, 0.375615146690670},
                        {1.163918992459650, -0.282654740327866, -2.681009102754735, 0.651076180271244},
                        {-0.228989733335134, 0.780774429183870, 0.527462446687548, -1.798461375234114}});
        DMatrixRMaj real_d1 = new DMatrixRMaj(new double[][]
                {{0.533507267365220, 0.478668605021239, 0.355502675273882, 0.318960921554189},
                        {0.581287932380312, 0.052464264578963, 0.387341331049122, 0.034959573289805},
                        {0.474266838961946, 0.425517439218574, 0.131184348075831, 0.117700044095355},
                        {0.516741958508340, 0.046638654133449, 0.142933157837299, 0.012900462218742}});
        assertArrayEquals(real_d0.data, d.D0.data, tol);
    }


}
