package contactrees.test.operators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import beast.base.inference.OperatorSchedule;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import contactrees.model.ConversionPrior;
import contactrees.operators.CFSubtreeSlide;
import contactrees.test.ContactreesTest;

public class CFSubtreeSlideTest extends ContactreesTest {

    private static final double SCALE_FACTOR = 0.8;

    private CFSubtreeSlide getOperator(boolean optimise) {
        RealScalar<NonNegativeReal> cRate = new RealScalarParam<>(1.0, NonNegativeReal.INSTANCE);
        RealScalar<UnitInterval> pMove = new RealScalarParam<>(0.4, UnitInterval.INSTANCE);
        ConversionPrior prior = new ConversionPrior();
        prior.initByName("network", acg, "conversionRate", cRate);

        CFSubtreeSlide operator = new CFSubtreeSlide();
        operator.initByName("acg", acg, "conversionRate", cRate, "pMove", pMove,
                            "blockSet", blockSet, "conversionPrior", prior,
                            "scaleFactor", SCALE_FACTOR, "optimise", optimise,
                            "weight", 1.0);

        // calcDelta() consults the schedule, so the operator needs one to be tunable.
        OperatorSchedule schedule = new OperatorSchedule();
        schedule.initByName("autoOptimizeDelay", 0);
        operator.setOperatorSchedule(schedule);

        return operator;
    }

    @Test
    public void testCoercableParameterStartsAtScaleFactor() {
        assertEquals(SCALE_FACTOR, getOperator(true).getCoercableParameterValue(), EPS);
    }

    @Test
    public void testOptimiseAdjustsScaleFactor() {
        CFSubtreeSlide operator = getOperator(true);
        operator.optimize(Math.log(0.9));
        assertTrue(operator.getCoercableParameterValue() != SCALE_FACTOR);
    }

    @Test
    public void testOptimiseCanBeDisabled() {
        CFSubtreeSlide operator = getOperator(false);
        operator.optimize(Math.log(0.9));
        assertEquals(SCALE_FACTOR, operator.getCoercableParameterValue(), EPS);
    }

    @Test
    public void testCoercableParameterStaysInUnitInterval() {
        CFSubtreeSlide operator = getOperator(true);
        operator.setCoercableParameterValue(17.0);
        assertTrue(operator.getCoercableParameterValue() > 0.0);
        assertTrue(operator.getCoercableParameterValue() < 1.0);
    }
}
