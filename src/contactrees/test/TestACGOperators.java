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
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.hamcrest.core.IsEqual;
import org.hamcrest.core.IsNot;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import contactrees.ACGWithBlocksReader;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.model.ACGSimulator;
import beast.app.BeastMCMC;
import beast.util.Randomizer;

public class TestACGOperators {

	static ArrayList<ACGWithBlocksReader> samplesSimulator;
	static ArrayList<ACGWithBlocksReader> samplesMCMC;
	
	static int N_SAMPLES = 6000;
	static int N_SAMPLES_SIMU = N_SAMPLES;
	static int N_SAMPLES_MCMC = N_SAMPLES;
//
	static int N_BLOCKS = 10;
	static double CONV_RATE = 0.05;
	static double P_MOVE = 0.15;
	static double POP_SIZE = 50.0;
	static int BURNIN_SAMPLES = 100;
	static int LOG_INTERVAL = 150000;
	
	static boolean SAMPLE_HEIGHT = false;
	static double HEIGHT = 40.;

	static int N_TAXA = 10;
	static boolean FIXED_CF = false;
	
	static String FBASE_SIM = String.format("simulateACGs%dtaxon%s", N_TAXA, FIXED_CF ? "_fixedCF" : "");
//	static String FBASE_SIM = String.format("simulateTrees%dtaxon%s", N_TAXA, FIXED_CF ? "_fixedCF" : "");
	static String FBASE_MCMC = String.format("contactrees%dtaxon%s", N_TAXA, FIXED_CF ? "_fixedCF" : "");
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		int chainLength = LOG_INTERVAL*(N_SAMPLES_MCMC + BURNIN_SAMPLES);
		String xmlParams = String.format("nSims=%d,convRate=%f,pMove=%f,popSize=%f,chainLength=%d,logInterval=%d,height=%f,sampleHeight=%b", 
				N_SAMPLES_SIMU, CONV_RATE, P_MOVE, POP_SIZE, chainLength, LOG_INTERVAL, HEIGHT, SAMPLE_HEIGHT);
		System.out.println(xmlParams);
		String workingDir = "simulations/";
		System.setProperty("beast.useWindow", "true"); // Trick beast into not stopping application
		System.setProperty("beast.debug", "true");

//		String seed = "1440939580";
		String seed = String.format("%d", Randomizer.nextInt());
		System.out.println("Seed: " + seed);

		
		// Run the simulation
		samplesSimulator = new ArrayList<>();
		String xmlPath = String.format("examples/ACGsimulations/%s.xml", FBASE_SIM);
		System.out.println(xmlPath);
		String[] simArgs = {"-overwrite", "-D", xmlParams, "-seed", seed, "-prefix", workingDir, xmlPath};
		BeastMCMC.main(simArgs);
		
		// Read the prior samplesSimulator
		samplesSimulator = readSamplesFromNexus(workingDir + FBASE_SIM + ".trees");
		
		// Run the MCMC
		samplesMCMC = new ArrayList<>();
		xmlPath = String.format("examples/ACGsimulations/%s.xml", FBASE_MCMC);
		System.out.println(xmlPath);
		String[] mcmcArgs = {"-overwrite", "-D", xmlParams, "-seed", seed, "-prefix", workingDir, xmlPath};
		BeastMCMC.main(mcmcArgs);
		
		// Read the prior samplesMCMC and drop the burn-in
		ArrayList<ACGWithBlocksReader> tmp = readSamplesFromNexus(workingDir + FBASE_MCMC + ".trees");
		samplesMCMC = new ArrayList<ACGWithBlocksReader>();
		for (int i=0; i<N_SAMPLES_MCMC; i++) {
			samplesMCMC.add(tmp.get(BURNIN_SAMPLES+i));
		}

	}

	@Test
	public void test() {
		assert samplesMCMC.size() == N_SAMPLES_MCMC;
		assert samplesSimulator.size() == N_SAMPLES_SIMU;
		
		// Statistics to collect
		double[] convCounts = new double[N_SAMPLES_MCMC];
		double[] convCountsTarget = new double[N_SAMPLES_SIMU];
		double[] moveCounts = new double[N_SAMPLES_MCMC];
		double[] moveCountsTarget = new double[N_SAMPLES_SIMU];
		double[] rootHeights = new double[N_SAMPLES_MCMC];
		double[] rootHeightsTarget = new double[N_SAMPLES_SIMU];
		double[] meanConvHeights = new double[N_SAMPLES_MCMC];
		double[] meanConvHeightsTarget = new double[N_SAMPLES_SIMU];
//		List<Double> convHeights = new ArrayList<Double>();
//		List<Double> convHeightsTarget = new ArrayList<Double>();

		
		// Evaluate statistics per sample
		for (int i=0; i<N_SAMPLES_MCMC; i++)
			collectStatistics(i, samplesMCMC, rootHeights, convCounts, 
							  moveCounts, meanConvHeights);
		for (int i=0; i<N_SAMPLES_SIMU; i++) {
			collectStatistics(i, samplesSimulator, rootHeightsTarget, convCountsTarget, 
							  moveCountsTarget, meanConvHeightsTarget);
		}
		
		printStats(rootHeights, rootHeightsTarget, "rootHeight");
		printStats(convCounts, convCountsTarget, "convCounts");
		printStats(meanConvHeights, meanConvHeightsTarget, "meanConvHeights");
		printStats(moveCounts, moveCountsTarget, "moveCounts");

	}
	
	protected static double ksTest(double[] samples1, double[] samples2) {
		KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
		return ksTest.kolmogorovSmirnovTest(samples1, samples2);
	}
	
	protected void collectStatistics(int i, List<ACGWithBlocksReader> samples, double[] heightsArray, double[] convCountArray, 
									 double[] moveCountsArray, double[] meanConvHeightsArray) {
		ACGWithBlocksReader sample = samples.get(i);
		heightsArray[i] = sample.acg.getRoot().getHeight();
		convCountArray[i] = sample.acg.getConvCount();
		moveCountsArray[i] = sample.blockSet.countMoves();
		meanConvHeightsArray[i] = meanConvHeight(sample.acg);
//		for (Conversion conv : sample.acg.getConversions()) {
////			System.out.println(conv.getHeight());
//			convHeightsList.add(conv.getHeight());
//		}
	}
	
	protected void printStats(double[] xMCMC, double[] xSimu, String statName) {
		double ksStatistic = ksTest(xSimu, xMCMC);

		System.out.println();
		System.out.println(String.format("KS-test %s: %.4f", statName, ksStatistic));
		System.out.println(String.format("MCMC %s: %.2f  \t   %.3f\t+/- %.2f", statName, 
									     (new Median()).evaluate(xMCMC), mean(xMCMC), stdOfMean(xMCMC)));
		System.out.println(String.format("Simu %s: %.2f  \t   %.3f\t+/- %.2f", statName, 
										 (new Median()).evaluate(xSimu), mean(xSimu), stdOfMean(xSimu)));
	}
	
	protected double[] meanArr(double[] xs, double[] ys) {
		assert xs.length == ys.length;
		double[] zs = new double[xs.length];
		for (int i=0; i<xs.length; i++) {
			zs[i] = (xs[i] + ys[i])/2.;
		}
		return zs;
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
	
	protected double meanConvHeight(ConversionGraph acg) {
		double sum = 0.0;
		for (Conversion conv : acg.getConversions()) {
			sum += conv.getHeight();
		}
		return sum / acg.getConvCount();
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
	
	protected double[] listToArray(List<Double> list) {
		double[] array = new double[list.size()]; 
		int i = 0;
		for (double value : list) {
			array[i] = value;
			i++;
		}
		return array;
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
	
}
