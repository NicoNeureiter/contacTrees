package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import contactrees.model.ConversionPrior;
import contactrees.util.Util;

public class ConversionPriorTest extends ContactreesTest {

	@Test
	public void TwoTaxaTest() {
		double gamma = 0.5;
		int N = 10;

		String newick = "(1:1.0,2:1.0)3:0.0;";
		setupFromNewick(newick);
		ConversionPrior prior = new ConversionPrior();
		prior.initByName("network", acg, "conversionRate", "0.5");

		for (int i=0; i<N; i++) {
		    // TODO change network

		    int n = acg.getConvCount();
			assert prior.calculateLogP() == poissonLogPDF(2*gamma, n) + 0.5*n + Util.logFactorial(n);

		}
	}

	@Test
	public void invalidConversionTest() {
		ConversionPrior prior = new ConversionPrior();
		prior.initByName("network", acg, "conversionRate", "0.5");
		assertTrue(Double.isFinite(prior.calculateLogP()));

		// Move tip 2 above conv1 (height 0.5), but keep it below its parent (height 1.0)
		node2.setHeight(0.7);
		assertFalse(conv1.isValid());
		assertEquals(Double.NEGATIVE_INFINITY, prior.calculateLogP(), 0.0);
	}

	/** Bounds on the conversion count renormalise the prior by P(lower <= N <= upper). */
	@Test
	public void conversionCountBoundsTest() {
		int lower = 1, upper = 3;
		double mean = 0.5 * acg.getClonalFramePairedLength();
		double pInBounds = 0.0;
		for (int k = lower; k <= upper; k++)
			pInBounds += Math.exp(poissonLogPDF(mean, k));

		ConversionPrior unbounded = new ConversionPrior();
		unbounded.initByName("network", acg, "conversionRate", "0.5");
		ConversionPrior bounded = new ConversionPrior();
		bounded.initByName("network", acg, "conversionRate", "0.5",
				"lowerConvCountBound", lower, "upperConvCountBound", upper);
		assertEquals(unbounded.calculateLogP() - Math.log(pInBounds), bounded.calculateLogP(), EPS);
	}

	protected double poissonLogPDF(double mean, int value) {
		return value * Math.log(mean) - mean - Util.logFactorial(value);
	}
}
