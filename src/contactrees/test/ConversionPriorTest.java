package contactrees.test;

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

	protected double poissonLogPDF(double mean, int value) {
		return value * Math.log(mean) - mean - Util.logFactorial(value);
	}
}
