/**
 * 
 */
package contactrees;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import javax.naming.directory.InvalidAttributesException;
import javax.xml.crypto.dsig.keyinfo.RetrievalMethod;

import beast.util.Randomizer;

/**
 * A container class for conversion edges.
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class ConversionList implements Iterable<Conversion> {

	HashMap<Integer, Conversion> convs;
	ConversionGraph acg;
	
	public ConversionList(ConversionGraph acg) {
		this.convs = new HashMap<>();
		this.acg = acg;
	}
	
	/**
	 * Find a free hash-code/conversion ID.
	 * @return The chosen conversion ID.
	 */
	private int getFreeKey() {
		int key = Randomizer.nextInt(Integer.MAX_VALUE);
		while (convs.containsKey(key)) {
			key = Randomizer.nextInt(Integer.MAX_VALUE);
		}
		return key;
	}
	
	/**
	 * Create a new conversion, put it in the list and return it.
	 * @return The new conversion.
	 */
	public Conversion addNewConversion() {
		startEditing();
		int key = getFreeKey();
		Conversion conv = new Conversion(key);
		convs.put(key, conv);
		return conv;
	}
	
	/**
	 * Add a conversion to the list and set the hash code as the conversion ID.
	 * @param Conversion to be added.
	 */
	public void add(Conversion conv) {
		assert conv.id == 0;
		startEditing();
		int key = getFreeKey();
		conv.setID(key);
		convs.put(key, conv);
	}
	
	/**
	 * Obtain the conversion at the given key.
	 * @param key
	 * @return The requested conversion.
	 */
	public Conversion get(int key) {
		return convs.get(key);
	}

	/**
	 * Remove the conversion at the given key from the list.
	 * @param The key of the conversion to be removed. 
	 */
	public void remove(int key) {
		startEditing();
		convs.remove(key);
	}
	
	/**
	 * Remove a conversion from the list.
	 * @param The conversion to be removed. 
	 */
	public void remove(Conversion conv) {
		remove(conv.getID());
	}

	/**
	 * Obtain the number of conversions in the list.
	 * @return Number of conversions.
	 */
	public int size() {
		return convs.size();
	}
	
	/**
	 * Obtain Iterator object to itarate over the conversions in the list. 
	 */
	@Override
	public Iterator<Conversion> iterator() {
		return convs.values().iterator();
	}
	
	/**
	 * Obtain the conversions in the list as a collection (not indexed).
	 * @return Conversions 
	 */
	public Collection<Conversion> getConversions() {
		return convs.values();
	}
	
	/**
	 * Obtain IDs (= hash map keys) for all conversions in the list.
	 * @return
	 */
	public Set<Integer> getKeys() {
		return convs.keySet();
	}
	
	/**
	 * Remove all conversions from the list.
	 */
	public void clear() {
		startEditing();
		convs.clear();
	}
	
	/**
	 * Obtain the number of conversions in the list.
	 * @return Number of conversions.
	 */
	public int getConvCount() {
		return convs.size();
	}
	
	/**
	 * Copy all conversions into a new conversion list (with same IDs).
	 * Note: The conversion attributes are not copied (e.g. they reference the same Node objects).
	 * @return Copied conversion list.
	 */
	public ConversionList copy() {
        ConversionList convListCopy = new ConversionList(acg);
		
        for (Conversion conv : convs.values()) {
        	Conversion convCopy = conv.getCopy();
            convListCopy.convs.put(convCopy.getID(), convCopy);
        }
		
		return convListCopy;
	}
	
	/**
	 * Choose a random conversion from the list (uniformly).
	 * @return Random conversion.
	 * @throws InvalidAttributesException
	 */
	public Conversion getRandomConversion() throws InvalidAttributesException {
		int z = Randomizer.nextInt(getConvCount());
		for (Conversion conv : convs.values()) {
			if (z == 0)
				return conv;
			z--;
		}
		assert convs.size() == 0;
		throw new InvalidAttributesException("No conversions to choose from.");
	}


    /**
     * Mark ACG statenode as dirty if available.
     */
    public void startEditing() {
        if (acg != null)
            acg.startEditing(null);
    }
	
}
