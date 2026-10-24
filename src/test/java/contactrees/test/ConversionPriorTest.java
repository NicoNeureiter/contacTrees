package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

import contactrees.Conversion;
import contactrees.ConversionGraph;
import beast.base.evolution.tree.Node;
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
		prior.initByName("network", acg, "conversionRates", "0.5");

		for (int i=0; i<N; i++) {
		    // TODO change network

		    int n = acg.getConvCount();
			assert prior.calculateLogP() == poissonLogPDF(2*gamma, n) + 0.5*n + Util.logFactorial(n);

		}
	}
	
	@Test
	public void MultipleEpochsTest() {
		double gamma_0 = 0.25; // the conversion rate in the earliest epoch
		double gamma_1 = 0.5; // the conversion rate in the latest epoch
		// Setup the initial tree topology
		String newick = "((1:1.0,2:1.0)4:2.0,3:3.0)5:0.0;";
		ConversionGraph acg = getACGFromNewick(newick);
		// Assign nodes
		Node root = acg.getRoot();
        Node node4 = root.getRight();
        Node node3 = root.getLeft();
        Node node2 = node3.getRight();
        Node node1 = node3.getLeft();
		ConversionPrior prior = new ConversionPrior();
		prior.initByName("network", acg, "linearContactGrowth", "true", "conversionRateChangeTimes", "2.0", "conversionRates", "0.25 0.5", "reverseChangeTimes", "true");

		double[] times = prior.getTimes();
		double[] rates = prior.getRates();
		double[] clonalFrameLengths = prior.getClonalFrameLengths();
		double[] expectedConversions = prior.getExpectedConversions();
		// The internal times array should consist of change times from lowest to highest,
		// and all other internal arrays should be reversed 
		// so they proceed from the lowest epoch to the highest:
		assertArrayEquals(new double[] {2.0}, times, 0.0);
		assertArrayEquals(new double[] {0.5, 0.25}, rates, 0.0);
		assertArrayEquals(new double[] {5.0, 2.0}, clonalFrameLengths, 0.0);
		assertArrayEquals(new double[] {2.5, 0.5}, expectedConversions, 0.0);
		// Get prior with no conversion edges:
		assertEquals(-(5.0)*gamma_1 - (2.0)*gamma_0, prior.calculateLogP(), 0.0);
		// Add one conversion in the lower epoch and check the prior again:
        Conversion conv1 = new Conversion(node3, node4, 1.5, acg, 1);
        acg.addConversion(conv1);
        assertEquals(Math.log(gamma_1) - (5.0)*gamma_1 - (2.0)*gamma_0, prior.calculateLogP(), 0.0);
        // Add another conversion in the higher epoch and check the prior again:
        Conversion conv2 = new Conversion(node4, node3, 2.5, acg, 2);
        acg.addConversion(conv2);
        assertEquals(Math.log(gamma_0) + Math.log(gamma_1) - (5.0)*gamma_1 - (2.0)*gamma_0, prior.calculateLogP(), 0.0);
	}

	@Test
	public void invalidConversionTest() {
		ConversionPrior prior = new ConversionPrior();
		prior.initByName("network", acg, "conversionRates", "0.5");
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
		unbounded.initByName("network", acg, "conversionRates", "0.5");
		ConversionPrior bounded = new ConversionPrior();
		bounded.initByName("network", acg, "conversionRates", "0.5",
				"lowerConvCountBound", lower, "upperConvCountBound", upper);
		assertEquals(unbounded.calculateLogP() - Math.log(pInBounds), bounded.calculateLogP(), EPS);
	}

	protected double poissonLogPDF(double mean, int value) {
		return value * Math.log(mean) - mean - Util.logFactorial(value);
	}
}
