package contactrees.test.operators;

import org.junit.Test;

public class ScaleOperatorsSimulationTest extends SingleOperatorTest {

    @Test
    public void testFullTreeScaler() throws Exception {
        testStationarity("Scale");
    }

    @Test
    public void testRootOnly() throws Exception {
        testStationarity("ScaleRoot");
    }

}
