/**
 *
 */
package contactrees;

import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;

import javax.naming.directory.InvalidAttributesException;

import beast.base.util.Randomizer;
import contactrees.util.Util;

/**
 * A container class for conversion edges.
 *
 * @author Nico Neureiter
 */
public class ConversionList implements Iterable<Conversion> {

    HashMap<Integer, Conversion> convs;
    ConversionGraph acg;
    Conversion _lastAdded;
    Deque<Conversion> trashCan;
    static int MAX_TRASH_CAN_SIZE = 100;

    public ConversionList(ConversionGraph acg) {
        this.convs = new HashMap<>();
        this.acg = acg;
        this.trashCan = new LinkedList<>();
    }

    /**
     * Find a free hash-code/conversion ID.
     * @return The chosen conversion ID.
     */
    private Integer getFreeKey() {
        int key = Randomizer.nextInt(10000);
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

        Conversion conv = getNewOrRecycledConversion(getFreeKey());
        add(conv);

        return conv;
    }

    /**
     * If there are old conversions in the trash-can, take one of
     * them, clear any previous information and return it. Otherwise
     * create a new conversion.
     * @param key The conversion ID for the new/recycled conversion.
     * @return conv The new/recycled conversion.
     */
    public Conversion getNewOrRecycledConversion(Integer key) {
        if (trashCan.isEmpty())
            return new Conversion(key);
        else {
            Conversion conv = trashCan.pop();
            conv.clear();
            conv.setID(key);
            return conv;
        }
    }

    /**
     * Add a conversion to the list and set the hash code as the conversion ID.
     * @param Conversion to be added.
     */
    public void add(Conversion conv) {
        startEditing();

        if (conv.id == 0)
            conv.setID(getFreeKey());

        if (convs.containsKey(conv.id))
            throw new RuntimeException("Conversion " + conv.id + " is already in the ConversionList.");

        convs.put(conv.id, conv);
        _lastAdded = conv;
    }


    /**
     * Create a duplicate of the speicified conversion and add it to the list
     * @param Conversion to be duplicated
     * @return The duplicate
     */
    public Conversion duplicateConversion(Conversion conv) {
        assert convs.containsKey(conv.id);
        startEditing();

        // Copy the original conversion
        Conversion newConv= getNewOrRecycledConversion(conv.id);
        conv.copyTo(newConv);

        // Assign a new ID and add the copied conversion to the list
        newConv.setID(getFreeKey());
        add(newConv);

        return newConv;
    }

	/**
	 * Obtain the conversion at the given key.
	 * @param key
	 * @return The requested conversion
	 */
	public Conversion get(Integer key) {
		return convs.get(key);
	}

	/**
	 * Remove the conversion at the given key from the list.
	 * @param The key of the conversion to be removed.
	 */
	public void remove(Integer key) {
		startEditing();
		if (trashCan.size() < MAX_TRASH_CAN_SIZE)
		    trashCan.add(convs.get(key));
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
	 * Obtain Iterator object to iterate over the conversions in the list.
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
		for (Conversion conv : convs.values()) {
		    if (trashCan.size() >= MAX_TRASH_CAN_SIZE)
		        break;
		    trashCan.add(conv);
		}
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
        copyTo(convListCopy);
        return convListCopy;
    }

    public void copyTo(ConversionList other) {
        for (Conversion conv : convs.values()) {
            Conversion convCopy = getNewOrRecycledConversion(conv.getID());
            conv.copyTo(convCopy);
            other.convs.put(convCopy.getID(), convCopy);
        }
    }

	/**
	 * Choose a random conversion from the list (uniformly).
	 * @return Random conversion.
	 * @throws InvalidAttributesException
	 */
	public Conversion getRandomConversion() {
		return Util.sampleFrom(convs.values());
	}

    /**
     * Mark ACG statenode as dirty if available.
     */
    public void startEditing() {
        if (acg != null)
            acg.startEditing(null);
    }

    /**
     * @return An array containing the conversions.
     */
    public Conversion[] asArray() {
      Conversion[] convArray = new Conversion[size()];
      int i = 0;
      for (Conversion conv : this)
    	  convArray[i++] = conv;

      return convArray;
    }

    /**
     * Return the conversions as an array, sorted by height in ascending order
     */
    public Conversion[] asSortedArray() {
    	Conversion[] convArray = asArray();
    	Arrays.sort(convArray,
        		(c1, c2) -> {
        			if (c1.height < c2.height) return -1;
        			if (c1.height > c2.height) return 1;
        			return 0;
		});
    	return convArray;
    }
}
