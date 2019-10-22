/**
 * 
 */
package contactrees.util;

import java.util.ArrayList;

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

}
