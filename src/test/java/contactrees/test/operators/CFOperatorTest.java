package contactrees.test.operators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import contactrees.Conversion;
import contactrees.model.ConversionPrior;
import contactrees.operators.CFOperator;
import contactrees.test.ContactreesTest;

public class CFOperatorTest extends ContactreesTest {
	
    double P_MOVE = 0.4;
    double C_RATE = 1.0;
            
	@Test
	public void testCollapseExpandSymmetry() {
		RealScalar<NonNegativeReal> cRate = new RealScalarParam<>(C_RATE, NonNegativeReal.INSTANCE);
		RealScalar<UnitInterval> pMove = new RealScalarParam<>(P_MOVE, UnitInterval.INSTANCE);
		ConversionPrior prior = new ConversionPrior();
		prior.initByName("network", acg, "conversionRate", cRate);

		CFOperator cfOp = new CollapseExpandOperator();
		cfOp.initByName("acg", acg, "conversionRate", cRate, "pMove", pMove,
						"blockSet", blockSet, "conversionPrior", prior,
						"weight", 1.0);
		double logHGF = cfOp.proposal();
//		System.out.println(logHGF);
		assertEquals(logHGF, 0.0, EPS);
		
	}
	
    /**
     * A conversion that degenerates into a loop while collapsing must never be visible
     * to the block-move Gibbs sampling in its invalid, half-moved state.
     */
    @Test
    public void testCollapseConversionsKeepsBorrowingStateValid() {
        RealScalar<NonNegativeReal> cRate = new RealScalarParam<>(C_RATE, NonNegativeReal.INSTANCE);
        RealScalar<UnitInterval> pMove = new RealScalarParam<>(P_MOVE, UnitInterval.INSTANCE);
        ConversionPrior prior = new ConversionPrior();
        prior.initByName("network", acg, "conversionRate", cRate);

        Conversion loopCandidate = new Conversion(node1, node2, 0.5, acg, 3);
        acg.addConversion(loopCandidate);
        blockSet.getBlocks().get(0).addMove(loopCandidate);

        CollapseValidityCheckingOperator cfOp = new CollapseValidityCheckingOperator();
        cfOp.initByName("acg", acg, "conversionRate", cRate, "pMove", pMove,
                        "blockSet", blockSet, "conversionPrior", prior,
                        "gibbsSampleBlockMoves", true,
                        "weight", 1.0);

        double logHGF = cfOp.proposal();
        assertTrue(Double.isFinite(logHGF));
        assertTrue(cfOp.checkedBorrowingState);
    }

	public class CollapseExpandOperator extends CFOperator {
		@Override
		public double proposal() {
			double logHGF = 0.0;
			
			int nConvs = acg.getConvCount();
            int nBlocks = blockSet.getBlockCount();
            int nMoves= blockSet.countMoves();
			
            logHGF -= expandConversions(node1, node3, 2.0);
            logHGF += collapseConversions(node1, node2, 1.0);
            assert acg.getConvCount() == nConvs;
            assert blockSet.getBlockCount() == nBlocks;
            assert blockSet.countMoves() == nMoves;
            assertEquals(logHGF, 0.0, EPS);
            
            
            logHGF -= expandConversions(node1, node3, 2.0);
            logHGF += collapseConversions(node1, node2, 1.0);
            assert acg.getConvCount() == nConvs;
            assert blockSet.getBlockCount() == nBlocks;
            assert blockSet.countMoves() == nMoves;
            assertEquals(logHGF, 0.0, EPS);
            
            logHGF -= expandConversions(node1, node3, 2.0);
            logHGF += collapseConversions(node1, node2, 1.0);
            assert acg.getConvCount() == nConvs;
            assert blockSet.getBlockCount() == nBlocks;
            assert blockSet.countMoves() == nMoves;
            assertEquals(logHGF, 0.0, EPS);
            
			return logHGF;
		}
	}

    /**
     * Collapses a conversion and records whether the borrowing probability was evaluated,
     * checking that every conversion is valid at that point.
     */
    public class CollapseValidityCheckingOperator extends CFOperator {
        boolean checkedBorrowingState = false;

        @Override
        public double proposal() {
            return collapseConversions(node1, node2, 0.2);
        }

        @Override
        public double getBorrowingsProbGibbs(Conversion conv, boolean mtreesChanged) {
            checkedBorrowingState = true;
            assert conv.isValid() : "Borrowing probability evaluated on invalid conversion " + conv.getID();
            for (Conversion other : acg.getConversions()) {
                assert other.isValid() : "Invalid conversion present during borrowing evaluation: " + other.getID();
            }
            return 0.0;
        }
    }
}
