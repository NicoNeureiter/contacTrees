package contactrees.test.operators;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.junit.Before;
import org.junit.Test;

import contactrees.ACGWithBlocksReader;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.util.Util;
import beast.app.BeastMCMC;
import beast.util.Randomizer;


public abstract class TestSingleOperator {

	static ArrayList<ACGWithBlocksReader> samplesSimulator;
	static ArrayList<ACGWithBlocksReader> samplesMCMC;
	
	static final int N_SAMPLES = 500;
	static final int N_SAMPLES_SIMU = 5000;
    static final int LOG_INTERVAL = 10000;

    static final int N_TAXA = 10;
	static final int N_BLOCKS = 10;
	static final double CONV_RATE = 0.025;
	static final double P_MOVE = 0.15;
	static final double POP_SIZE = 50.0;
 
	static final double KS_THRESHOLD = 0.001;
	
	static String FBASE_SIM = String.format("simulateACGs%dtaxon", N_TAXA);
	static String FBASE_MCMC;  // Will be assigned in setUp, based on subclass specific operator
	

	public void generateSamples(String operator) throws Exception {
	    FBASE_MCMC = String.format("%dtaxa_%s", N_TAXA, operator);
	    
		int chainLength = LOG_INTERVAL*N_SAMPLES;
		String xmlParams = String.format("nSims=%d,convRate=%f,pMove=%f,popSize=%f,chainLength=%d,logInterval=%d", 
		        N_SAMPLES_SIMU, CONV_RATE, P_MOVE, POP_SIZE, chainLength, LOG_INTERVAL);
		System.out.println(xmlParams);
		String workingDir = "simulations/";
		System.setProperty("beast.useWindow", "true"); // Trick beast into not stopping application

		String seed = String.format("%d", Randomizer.nextInt());
		seed = "1575254328";
		System.out.println("Seed: " + seed);

		
		// Run the simulation and read samples from file
		samplesSimulator = new ArrayList<>();
		String xmlPath = String.format("examples/ACGsimulations/%s.xml", FBASE_SIM);
		System.out.println(xmlPath);
		String[] simArgs = {"-overwrite", "-D", xmlParams, "-seed", seed, "-prefix", workingDir, xmlPath};
		BeastMCMC.main(simArgs);
		samplesSimulator = readSamplesFromNexus(workingDir + FBASE_SIM + ".trees");
		
		// Run the MCMC and read samples from file
		samplesMCMC = new ArrayList<>();
		xmlPath = String.format("examples/operatorTests/%s.xml", FBASE_MCMC);
		String[] mcmcArgs = {"-overwrite", "-D", xmlParams, "-seed", seed, "-prefix", workingDir, xmlPath};
		System.out.println(xmlPath);
		
		BeastMCMC.main(mcmcArgs);
		samplesMCMC = readSamplesFromNexus(workingDir + FBASE_MCMC + ".trees");
		samplesMCMC.remove(0);
	}

	public void testStationarity(String operatorName) throws Exception {
	    generateSamples(operatorName);
	    
		assert samplesMCMC.size() == N_SAMPLES : samplesMCMC.size();
		assert samplesSimulator.size() == N_SAMPLES_SIMU;
		
		
		// Statistics to collect
		double[] convCounts = new double[N_SAMPLES];
		double[] convCountsTarget = new double[N_SAMPLES_SIMU];
		double[] moveCounts = new double[N_SAMPLES];
		double[] moveCountsTarget = new double[N_SAMPLES_SIMU];
		double[] rootHeights = new double[N_SAMPLES];
		double[] rootHeightsTarget = new double[N_SAMPLES_SIMU];
		double[] meanConvHeights = new double[N_SAMPLES];
		double[] meanConvHeightsTarget = new double[N_SAMPLES_SIMU];
		
		// Evaluate statistics per sample
		for (int i=0; i<N_SAMPLES; i++)
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
		heightsArray[i] = sample.getRoot().getHeight();
		convCountArray[i] = sample.getConvCount();
		moveCountsArray[i] = sample.blockSet.countMoves();
		meanConvHeightsArray[i] = meanConvHeight(sample);
	}
	
	protected void printStats(double[] xMCMC, double[] xSimu, String statName) {
		
		double ksStatistic = ksTest(xSimu, xMCMC);
		
		System.out.println();
		System.out.println(String.format("KS-test %s: %.4f", statName, ksStatistic));
		System.out.println(String.format("MCMC %s: %.2f  \t   %.3f\t+/- %.2f", statName, 
									     (new Median()).evaluate(xMCMC), mean(xMCMC), stdOfMean(xMCMC)));
		System.out.println(String.format("Simu %s: %.2f  \t   %.3f\t+/- %.2f", statName, 
										 (new Median()).evaluate(xSimu), mean(xSimu), stdOfMean(xSimu)));
		

		// We test whether a bootstrap sample of the KS-Statistic averages out to 0.5 (as it should)-
		assertTrue(ksStatistic > KS_THRESHOLD);
//        int nFolds = 500;
//        double ksStatSum = 0.0;
//        for (int iFold=0; iFold<nFolds; iFold++) {
//            double[] xMCMCFold = Util.sampleSubset(xMCMC, 200);
//            double[] xSimuFold = Util.sampleSubset(xSimu, 200);
//            ksStatSum += ksTest(xSimuFold, xMCMCFold);
//        }
//        double ksStatMean = ksStatSum / nFolds;
//        System.out.println(ksStatMean);
//        assertTrue(ksStatMean >= 0.1);
        
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
		
		reader.close();
		return samples;
	}
	
	public static ACGWithBlocksReader parseSampleLine(String line) {
		String newick = line
				.split("&R\\]", 2)[1]
				.replace(";"," ")
				.trim();
		
		assert newick != "";
		assert newick != null;
		
		ACGWithBlocksReader reader = ACGWithBlocksReader.newFromNewick(N_BLOCKS, newick);
		//reader.fromExtendedNewick(newick, false, 1);
		
		return reader;
	}
}
