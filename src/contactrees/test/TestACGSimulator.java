package contactrees.test;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import contactrees.ACGWithBlocksReader;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.ConversionGraph;
import contactrees.model.ACGSimulator;
import beast.app.BeastMCMC;
import beast.util.Randomizer;

public class TestACGSimulator {

	static ArrayList<ACGWithBlocksReader> samples;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		int nBlocks = 10;

		samples = new ArrayList<>();
		
		// Run the simulation
		String fbase = "simulateACGs2taxon_fixedCF";
		String xmlPath = String.format("examples/ACGsimulations/%s.xml", fbase);
		String workingDir = "simulations/";
		
		System.setProperty("beast.useWindow", "true"); // Trick beast into not stopping application
		String[] args = {"-overwrite", "-prefix", workingDir, xmlPath};
		BeastMCMC.main(args);
		
		// Read the prior samplesSimulator
		String treesPath = workingDir + fbase + ".trees";
		BufferedReader reader = new BufferedReader(new FileReader(treesPath));
		
		String line = "";
		while (!line.startsWith("0\ttree"))
			line = reader.readLine();
		
		while (!line.startsWith("End")) {
			samples.add(parseSampleLine(line, nBlocks));
			line = reader.readLine();
		}
		
		reader.close();
	}
	
	public static ACGWithBlocksReader parseSampleLine(String line, int nBlocks) {
		String newick = line
				.split("&R\\]", 2)[1]
				.replace(";"," ")
				.trim();
		
		assert newick != "";
		assert newick != null;
		
		return ACGWithBlocksReader.newFromNewick(nBlocks, newick);		
	}

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void test() {
		int nSamples = samples.size(); 
		
		// Statistics to collect
		double[] convCounts = new double[nSamples];
		double[] convCountsTarget = new double[nSamples];
		double[] moveCounts = new double[nSamples];
		double[] moveCountsTarget = new double[nSamples];
		
		// Evaluate statistics per sample
		for (int i=0; i<nSamples; i++) {
			ConversionGraph acg = samples.get(i);
			BlockSet blockSet = samples.get(i).blockSet;
			double iProb = i / (double) nSamples;
			
			// Collect empirical and target distribution for conversion counts
			convCounts[i] = acg.getConvCount();
			PoissonDistribution distr = new PoissonDistribution(meanConvCount(acg, 0.5));
			convCountsTarget[i] = (double) distr.inverseCumulativeProbability(iProb);
			
			// Collect empirical and target distribution for block move counts
			moveCounts[i] = 0;
			for (Block b : blockSet.getBlocks()) moveCounts[i] += b.size();
			int nTries = acg.getConvCount()*blockSet.getBlockCount();
			
			BinomialDistribution distrMoves = new BinomialDistribution(nTries, 0.4);
			moveCountsTarget[i] = (double) distrMoves.inverseCumulativeProbability(iProb);
			
//			System.out.print(moveCounts[i]);
//			System.out.print(" -- ");
//			System.out.println(moveCountsTarget[i]);
//			System.out.println(nTries);
		}
		

		KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
		double pValueKSTest = ksTest.kolmogorovSmirnovTest(convCountsTarget, convCounts);
		System.out.println(pValueKSTest);

		ksTest = new KolmogorovSmirnovTest();
		double pValueKSTestMoves = ksTest.kolmogorovSmirnovTest(moveCountsTarget, moveCounts);
		System.out.println(pValueKSTestMoves);
		
		
	}
	
	private double meanConvCount(ConversionGraph acg, double conversionRate) {
		return conversionRate * acg.getClonalFramePairedLength();
	}
	
}
