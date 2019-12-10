package contactrees.test.operators;

import org.junit.Test;

public class SimulationTestConversionOperators extends TestSingleOperator {
  
    @Test
    public void testAddRemoveConversion() throws Exception {
        testStationarity("AddRemoveConversion");
    }

    @Test
    public void testAddRemoveConversionGibbs() throws Exception {
        testStationarity("AddRemoveConversionGibbs");
    }
    
    @Test
    public void testConvertedEdgeSlide() throws Exception {
        testStationarity("ConvertedEdgeSlide");
    }
    
    @Test
    public void testResampleConversionMoves() throws Exception {
        testStationarity("ResampleConversionMoves");
    }
    
    @Test
    public void testGibbsSampleMovesPerConversion() throws Exception {
        testStationarity("GibbsSampleMovesPerConversion");
    }
    
}
