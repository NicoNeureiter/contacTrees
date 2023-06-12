package contactrees.test.operators;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import contactrees.model.ConversionPrior;
import contactrees.operators.CFOperator;
import contactrees.test.ContactreesTest;

public class CFOperatorTest extends ContactreesTest {

    double P_MOVE = 0.4;
    double C_RATE = 1.0;

	@Test
	public void testCollapseExpandSymmetry() {
		RealParameter cRate = new RealParameter(Double.toString(C_RATE));
		RealParameter pMove = new RealParameter(Double.toString(P_MOVE));
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
}
