package contactrees.test.operators;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.junit.Before;
import org.junit.Test;

import contactrees.ACGWithBlocks;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.util.Util;
import beast.app.BeastMCMC;
import beast.util.Randomizer;


public abstract class TestSingleOperator {

	static ArrayList<ACGWithBlocks> samplesSimulator;
	static ArrayList<ACGWithBlocks> samplesMCMC;

    static final int N_TAXA = 10;
    static final int N_BLOCKS = 10;
	static final double KS_THRESHOLD = 0.001;
	
	static String FBASE_SIM = String.format("simulateACGs%dtaxon", N_TAXA);
	static String FBASE_MCMC;  // Will be assigned in setUp, based on subclass specific operator
	
	public void generateSamples(String operator) throws Exception {
	    FBASE_MCMC = String.format("%dtaxa_%s", N_TAXA, operator);
	    
		String workingDir = "simulations/";
		System.setProperty("beast.useWindow", "true"); // Trick beast into not stopping application

		String seed = String.format("%d", Randomizer.nextInt());
//		seed = "1575254328";
		System.out.println("Seed: " + seed);

		
		// Run the simulation and read samples from file
		samplesSimulator = new ArrayList<>();
		String xmlPath = String.format("examples/ACGsimulations/%s.xml", FBASE_SIM);
		String paramsPath = "examples/operatorTests/parameters.json";
		System.out.println(xmlPath);
		String[] simArgs = {"-overwrite", "-DF", paramsPath, "-seed", seed, "-prefix", workingDir, xmlPath};
		System.out.println(simArgs);
		BeastMCMC.main(simArgs);
		samplesSimulator = readSamplesFromNexus(workingDir + FBASE_SIM + ".trees");
		
		// Run the MCMC and read samples from file
		samplesMCMC = new ArrayList<>();
		xmlPath = String.format("examples/operatorTests/%s.xml", FBASE_MCMC);
		String[] mcmcArgs = {"-overwrite", "-DF", paramsPath, "-seed", seed, "-prefix", workingDir, xmlPath};
		System.out.println(xmlPath);
		
		BeastMCMC.main(mcmcArgs);
		samplesMCMC = readSamplesFromNexus(workingDir + FBASE_MCMC + ".trees");
        samplesMCMC.remove(0);
        samplesMCMC.remove(1);
	}

	public void testStationarity(String operatorName) throws Exception {
	    generateSamples(operatorName);

        int n_samples_simu = samplesSimulator.size();
	    int n_samples_mcmc = samplesMCMC.size();
		
		// Statistics to collect
        double[] convCounts = new double[n_samples_mcmc];
        double[] convCountsTarget = new double[n_samples_simu];
		double[] moveCounts = new double[n_samples_mcmc];
		double[] moveCountsTarget = new double[n_samples_simu];
		double[] rootHeights = new double[n_samples_mcmc];
		double[] rootHeightsTarget = new double[n_samples_simu];
		double[] meanConvHeights = new double[n_samples_mcmc];
		double[] meanConvHeightsTarget = new double[n_samples_simu];
        double[] obsConvRate = new double[n_samples_mcmc];
        double[] obsConvRateTarget= new double[n_samples_simu];
        double[] movesPerConv = new double[n_samples_mcmc];
        double[] movesPerConvTarget= new double[n_samples_simu];
        
		// Evaluate statistics per sample
		for (int i=0; i<n_samples_mcmc; i++)
			collectStatistics(i, samplesMCMC, rootHeights, convCounts, 
							  moveCounts, meanConvHeights, obsConvRate,
							  movesPerConv);
		for (int i=0; i<n_samples_simu; i++) {
			collectStatistics(i, samplesSimulator, rootHeightsTarget, convCountsTarget, 
							  moveCountsTarget, meanConvHeightsTarget, obsConvRateTarget, 
							  movesPerConvTarget);
		}
		
		printStats(rootHeights, rootHeightsTarget, "rootHeight");
        printStats(obsConvRate, obsConvRateTarget, "obsConvRate");
		printStats(convCounts, convCountsTarget, "convCounts");
		printStats(meanConvHeights, meanConvHeightsTarget, "meanConvHeights");
		printStats(moveCounts, moveCountsTarget, "moveCounts");
        printStats(movesPerConv, movesPerConvTarget, "movesPerConv");
	}
	
	protected static double ksTest(double[] samples1, double[] samples2) {
		KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
		return ksTest.kolmogorovSmirnovTest(samples1, samples2);
	}
	
	protected void collectStatistics(int i, List<ACGWithBlocks> samples, double[] heightsArray, double[] convCountArray, 
									 double[] moveCountsArray, double[] meanConvHeightsArray, double[] obsConvRate, double[] movesPerConv) {
		ACGWithBlocks sample = samples.get(i);
		heightsArray[i] = sample.getRoot().getHeight();
		convCountArray[i] = sample.getConvCount();
		moveCountsArray[i] = sample.blockSet.countMoves();
		meanConvHeightsArray[i] = meanConvHeight(sample);
		obsConvRate[i] = (double) sample.getConvCount() / ((double) sample.getClonalFramePairedLength());
		movesPerConv[i] = (double) sample.blockSet.countMoves() / ((double) sample.getConvCount());

	}
	
	protected void printStats(double[] xMCMC, double[] xSimu, String statName) {
	    xMCMC = filterNaNs(xMCMC);
	    xSimu = filterNaNs(xSimu);
	    
	    double ksStatistic = ksTest(xSimu, xMCMC);
		
		System.out.println();
		System.out.println(String.format("KS-test %s: %.4f", statName, ksStatistic));
		System.out.println(String.format("MCMC %s: %.2f  \t   %.3f\t+/- %.2f", statName, 
									     (new Median()).evaluate(xMCMC), mean(xMCMC), std(xMCMC)));
		System.out.println(String.format("Simu %s: %.2f  \t   %.3f\t+/- %.2f", statName, 
										 (new Median()).evaluate(xSimu), mean(xSimu), std(xSimu)));

		// We test whether a bootstrap sample of the KS-Statistic averages out to 0.5 (as it should)-
		assertTrue(ksStatistic > KS_THRESHOLD);
	}
	
	double[] filterNaNs(double[] x) {
	    List<Double> xList = asList(x);
	    xList.removeIf(v -> Double.isNaN(v));
        return asArray(xList); 
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
	
	protected double[] asArray(Collection<Double> list) {
		double[] array = new double[list.size()]; 
		int i = 0;
		for (double value : list) {
			array[i] = value;
			i++;
		}
		return array;
	}
	
	protected ArrayList<Double> asList(double[] array) {
	    ArrayList<Double> list = new ArrayList<>();
	    for (Double x : array) {
	        list.add(x);
	    }
	    return list;
	}
	
	protected static ArrayList<ACGWithBlocks> readSamplesFromNexus(String path) throws IOException {
		ArrayList<ACGWithBlocks> samples = new ArrayList<>();
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
	
	public static ACGWithBlocks parseSampleLine(String line) {
		String newick = line
				.split("&R\\]", 2)[1]
				.replace(";"," ")
				.trim();
		
		assert newick != "";
		assert newick != null;
		
		ACGWithBlocks reader = ACGWithBlocks.newFromNewick(N_BLOCKS, newick);
		//reader.fromExtendedNewick(newick, false, 1);
		
		return reader;
	}
}
