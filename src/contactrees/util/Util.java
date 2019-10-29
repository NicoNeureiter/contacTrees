/**
 * 
 */
package contactrees.util;

import java.util.ArrayList;
import java.util.Collection;

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
	
}
