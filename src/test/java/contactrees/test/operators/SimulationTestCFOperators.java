package contactrees.test.operators;

import org.junit.Test;

public class SimulationTestCFOperators extends TestSingleOperator {
	
    @Test
    public void testWilsonBalding() throws Exception {
        testStationarity("WilsonBalding");
    }
    
    @Test
    public void testNarrowExchange() throws Exception {
        testStationarity("NarrowExchange");
    }
    
    @Test
    public void testWideExchange() throws Exception {
        testStationarity("WideExchange");
    }
    
    @Test
    public void testUniform() throws Exception {
        testStationarity("Uniform");
    }

    @Test
    public void testCFConversionSwap() throws Exception {
        testStationarity("CFConversionSwap");
    }
}
