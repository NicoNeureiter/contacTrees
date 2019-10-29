package contactrees.test;

import static org.junit.Assert.*;

import java.awt.font.NumericShaper;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
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

public class TestACGOperators {

	static ArrayList<ACGWithBlocksReader> samplesSimulator;
	static ArrayList<ACGWithBlocksReader> samplesMCMC;
	
	static int N_SAMPLES = 3000;
	static int N_BLOCKS = 10;
	static double CONV_RATE = 0.02;
	static double MOVE_PROB = 0.0;
	static double POP_SIZE = 12.0;
	static int BURNIN_SAMPLES = 10;
	static int LOG_INTERVAL = 10000;
	
	// Simulator settings
//	static String FBASE_SIM = "simulateACGs2taxon";
//	static String FBASE_SIM = "simulateACGs2taxon_fixedCF";
	static String FBASE_SIM = "simulateACGs5taxon";
//	static String FBASE_SIM = "simulateACGs5taxon_fixedCF";
	// MCMC settings 
//	static String FBASE_MCMC = "contactrees2taxon";
	static String FBASE_MCMC = "contactrees5taxon";
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		int chainLength = LOG_INTERVAL*N_SAMPLES;
		String xmlParams = String.format("nSims=%d,convRate=%f,moveProb=%f,popSize=%f,chainLength=%d,logInterval=%d", 
										 N_SAMPLES, CONV_RATE, MOVE_PROB, POP_SIZE, chainLength, LOG_INTERVAL);
		System.out.println(xmlParams);
		String workingDir = "simulations/";
		System.setProperty("beast.useWindow", "true"); // Trick beast into not stopping application


		// Run the simulation
		samplesSimulator = new ArrayList<>();
		String xmlPath = String.format("examples/ACGsimulations/%s.xml", FBASE_SIM);
		String[] simArgs = {"-overwrite", "-D", xmlParams, "-prefix", workingDir, xmlPath};
		BeastMCMC.main(simArgs);
		// Read the prior samplesSimulator
		String treesPath = workingDir + FBASE_SIM + ".trees";
		samplesSimulator = readSamplesFromNexus(treesPath);

		// Run the simulation
		samplesMCMC = new ArrayList<>();
		xmlPath = String.format("examples/ACGsimulations/%s.xml", FBASE_MCMC);
		String[] mcmcArgs = {"-overwrite", "-D", xmlParams, "-prefix", workingDir, xmlPath};
		BeastMCMC.main(mcmcArgs);
		// Read the prior samplesSimulator
		treesPath = workingDir + FBASE_MCMC + ".trees";
		samplesMCMC = readSamplesFromNexus(treesPath);

	}
	
	protected static ArrayList<ACGWithBlocksReader> readSamplesFromNexus(String path) throws IOException {
		ArrayList<ACGWithBlocksReader> samples = new ArrayList<>();
		BufferedReader reader = new BufferedReader(new FileReader(path));
		
		String line = "";
		while (!line.startsWith("0\ttree"))
			line = reader.readLine();
		
		while (!line.startsWith("End")) {
			samples.add(parseSampleLine(line));
			line = reader.readLine();
		}
		
		return samples;
	}
	
	public static ACGWithBlocksReader parseSampleLine(String line) {
		String newick = line
				.split("&R\\]", 2)[1]
				.replace(";"," ")
				.trim();
		
		assert newick != "";
		assert newick != null;
		
		ACGWithBlocksReader reader = new ACGWithBlocksReader(N_BLOCKS);
		reader.fromExtendedNewick(newick, false, 1);
		
		return reader;
	}

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void test() {
		int nSamples = samplesSimulator.size(); 
		
		// Statistics to collect
		double[] convCounts = new double[nSamples];
		double[] convCountsTarget = new double[nSamples];
		double[] moveCounts = new double[nSamples];
		double[] moveCountsTarget = new double[nSamples];
		double[] rootHeights = new double[nSamples];
		double[] rootHeightsTarget = new double[nSamples];
		List<Double> convHeights = new ArrayList<Double>();
		List<Double> convHeightsTarget = new ArrayList<Double>();
		
		// Evaluate statistics per sample
		for (int i=0; i<nSamples; i++) {
			collectStatistics(i, samplesMCMC, rootHeights, convCounts);
			collectStatistics(i, samplesSimulator, rootHeightsTarget, convCountsTarget);

//			double iProb = i / (double) nSamples;
//
//			if(CONV_RATE > 0.0) {
//				PoissonDistribution distr = new PoissonDistribution(meanConvCount(acg, CONV_RATE));
//				convCountsTarget[i] = (double) distr.inverseCumulativeProbability(iProb);
//			}

////			// Collect empirical and target distribution for block move counts
////			moveCounts[i] = 0;
////			for (Block b : blockSet.getBlocks()) moveCounts[i] += b.size();
////			int nTries = acg.getConvCount()*blockSet.getBlockCount();
////			if (MOVE_PROB > 0.0) {
////				BinomialDistribution distrMoves = new BinomialDistribution(nTries, MOVE_PROB);
////				moveCountsTarget[i] = (double) distrMoves.inverseCumulativeProbability(iProb);
////			}
//			
//			// Collect empirical and target distribution for root heights
//			rootHeights[i] = acg.getRoot().getHeight();
//			ExponentialDistribution distrHeights = new ExponentialDistribution(POP_SIZE);
//			rootHeightsTarget[i] = distrHeights.inverseCumulativeProbability(iProb);
		}
		

//		KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
//		double pValueKSTest = ksTest.kolmogorovSmirnovTest(convCountsTarget, convCounts);
//
////		ksTest = new KolmogorovSmirnovTest();
////		double pValueKSTestMoves = ksTest.kolmogorovSmirnovTest(moveCountsTarget, moveCounts);
////		System.out.println(pValueKSTestMoves);
//
//		ksTest = new KolmogorovSmirnovTest();
//		double pValueKSTestHeights = ksTest.kolmogorovSmirnovTest(rootHeightsTarget, rootHeights);
//		
//		System.out.println();
//		System.out.println("KS-test conv.-count:" + pValueKSTest);
//		System.out.println("KS-test root-height:" + pValueKSTestHeights);

		System.out.println();
		System.out.println("Mean root-height:" + mean(rootHeights));
		System.out.println("Stdev root-height:" + std(rootHeights));
		System.out.println("Stdev of mean:" + stdOfMean(rootHeights));
		System.out.println();
		System.out.println("Mean target root-height:" + mean(rootHeightsTarget));
		System.out.println("Stdev target root-height:" + std(rootHeightsTarget));
		System.out.println("Stdev of mean:" + stdOfMean(rootHeightsTarget));

		System.out.println();
		System.out.println("Mean conv.-count:" + mean(convCounts));
		System.out.println("Stdev of mean:" + stdOfMean(convCounts));
		System.out.println();
		System.out.println("Mean target conv.-count:" + mean(convCountsTarget));
		System.out.println("Stdev of mean:" + stdOfMean(convCountsTarget));
		System.out.println();
		
	}
	
	protected void collectStatistics(int i, List<ACGWithBlocksReader> samples, double[] heightsArray, double[] convCountArray) {
		heightsArray[i] = samples.get(i).acg.getRoot().getHeight();
		convCountArray[i] = samples.get(i).acg.getConvCount();
	}
	
	protected double mean(double[] xs) {
		double sum = 0.0;
		for (double x : xs) sum += x;
		return sum / xs.length;
	}
	
	protected double std(double[] xs) {
		double mu = mean(xs);
		double sum = 0.;
		for (double x : xs) sum += (x-mu)*(x-mu);
		return Math.sqrt(sum / xs.length);
	}
	
	protected double stdOfMean(double[] xs) {
		double mu = mean(xs);
		double sum = 0.;
		for (double x : xs) sum += (x-mu)*(x-mu);
		return Math.sqrt(sum / xs.length / xs.length);
	}
	
	protected double meanConvCount(ConversionGraph acg, double conversionRate) {
		return conversionRate * acg.getClonalFramePairedLength();
	}
	
	protected static String blockNames(int nBlocks) {
		String names = "{";
		for (int i=0; i<nBlocks; i++) {
			names += "block" + i;
			if (i < nBlocks -1)
				names += ",";
		}
		return names + "}";
	}
	
	protected static String taxonNames(int nTaxa) {
		String names = "";
		char taxName = 'A';

		for (int i=0; i<nTaxa; i++) {
			names += taxName;
			
			taxName += 1;
			if (i < nTaxa -1)
				names += ",";
		}
		
		return names;
	}
}
