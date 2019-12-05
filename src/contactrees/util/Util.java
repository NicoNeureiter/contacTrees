/**
 * 
 */
package contactrees.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * 
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class Util {

	/**
	 * 
	 */
	public Util() {}

	static public ArrayList<Integer> deepCopyIntegers(Iterable<Integer> original){
	    ArrayList<Integer> copy = new ArrayList<>();
	    for (Integer x : original)
	    	copy.add(new Integer(x));

	    return copy;
	}
	
	static public <T> T sampleFrom(Collection<T> population) {
		int z = Randomizer.nextInt(population.size());
		for (T candidate : population) {
			if (z == 0)
				return candidate;
			z--;
		}
		assert population.size() == 0;
		throw new RuntimeException("Can not sample from empty set.");
	}
	
	static public void sortByHeight(List<Node> nodes, boolean reverse) {
		if (reverse) {
	        Collections.sort(nodes, (Node n1, Node n2) -> {
	            if (n1.getHeight()<n2.getHeight())
	                return -1;
	            if (n1.getHeight()>n2.getHeight())
	                return 1;
	            return 0;
	        });
		} else {
	        Collections.sort(nodes, (Node n1, Node n2) -> {
	            if (n1.getHeight()<n2.getHeight())
	                return 1;
	            if (n1.getHeight()>n2.getHeight())
	                return -1;
	            return 0;
	        });
		}
		
	}
	
	static public double logFactorial(int n) {
		double res = 0;
		for (int i=1; i<=n; i++) {
			res += Math.log(i);
		}
		return res;
	}
	
	static public double sum(double[] values) {
		double total = 0;
		for (double x : values) total += x; 
		return total;
	}
	
	static public double sum(Iterable<Double> values) {
		double total = 0;
		for (double x : values) total += x; 
		return total;
	}
	
	static public int sampleCategorical(double[] weights) {
		double weightSum = sum(weights);
		double u = weightSum * Randomizer.nextDouble();
		
		int i = 0;
		for (double w : weights) {
			if (u < w)
				return i;
			
			u -= w;
			i++;
		}
	
		throw new RuntimeException("Can not sample from an empty collection.");
	}
	
	static public int sampleCategorical(Iterable<Double> weights) {
		double weightSum = sum(weights);
		double u = weightSum * Randomizer.nextDouble();
		
		int i = 0;
		for (double w : weights) {
			if (u < w)
				return i;
			
			u -= w;
			i++;
		}
	
		throw new RuntimeException("Can not sample from an empty collection.");
	}
    
    static public double[] sampleSubset(double[] samples, int subsetSize) {
        assert subsetSize <= samples.length; 
        
        double[] rndSubset = new double[subsetSize];
        int i = 0;
        
        for (int sampleIdx : Randomizer.shuffled(samples.length)) {
            rndSubset[i] = samples[sampleIdx];
            i++;
            if (i >= subsetSize)
                break;
        }
        
        return rndSubset;       
    }
    
    static public <T> ArrayList<T> sampleSubset(List<T> samples, int subsetSize) {
        assert subsetSize <= samples.size(); 
        
        ArrayList<T> rndSubset = new ArrayList<>(subsetSize);
        int i = 0;
        
        for (int sampleIdx : Randomizer.shuffled(samples.size())) {
            rndSubset.set(i, samples.get(sampleIdx));
            i++;
            if (i >= subsetSize)
                break;
        }
        
        return rndSubset;       
    }
}
