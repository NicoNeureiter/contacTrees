/**
 *
 */
package contactrees;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;

import javax.naming.directory.InvalidAttributesException;

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

    /** Dense, run-local slot index (see {@link Conversion#getSlot()}) for O(1) access by callers
     *  that already hold a Conversion object; never serialized. */
    private final ArrayList<Conversion> bySlot = new ArrayList<>();
    private final Deque<Integer> freeSlots = new ArrayDeque<>();

    /** Monotonically increasing, never reused for the life of this ConversionGraph -- closes a
     *  silent id-rebind hazard the old random-then-escalate scheme had (it only checked
     *  collisions against LIVE keys, so a Block holding a stale id could, after enough churn, be
     *  silently rebound to an unrelated conversion instead of failing loudly). */
    protected int nextConversionId = 1;

    public ConversionList(ConversionGraph acg) {
        this.convs = new HashMap<>();
        this.acg = acg;
        this.trashCan = new LinkedList<>();
    }

    /**
     * Find a free conversion ID.
     * @return The chosen conversion ID.
     */
    private Integer getFreeKey() {
        return nextConversionId++;
    }

    /**
     * Claim a slot for a conversion that doesn't have one yet: reuse a freed slot if one is
     * available, otherwise extend the dense array.
     */
    private int allocateSlot(Conversion conv) {
        int slot;
        if (!freeSlots.isEmpty()) {
            slot = freeSlots.pop();
            bySlot.set(slot, conv);
        } else {
            slot = bySlot.size();
            bySlot.add(conv);
        }
        conv.setSlot(slot);
        return slot;
    }

    /** Release a conversion's slot back to the free list. No-op if it has none. */
    private void releaseSlot(Conversion conv) {
        int slot = conv.getSlot();
        if (slot < 0)
            return;
        bySlot.set(slot, null);
        freeSlots.push(slot);
        conv.setSlot(-1);
    }

    /**
     * Obtain the conversion at the given slot -- O(1) dense array access, for callers that
     * already hold a Conversion object (and hence its slot) rather than a bare id.
     * @param slot
     * @return The requested conversion.
     */
    public Conversion getBySlot(int slot) {
        return bySlot.get(slot);
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

        if (conv.id == 0) {
            conv.setID(getFreeKey());
        } else if (conv.id >= nextConversionId) {
            // A caller assigned an explicit id (Newick parsing, test fixtures, ...), bypassing
            // getFreeKey() entirely -- push the counter past it so a later auto-generated id can
            // never collide with it.
            nextConversionId = conv.id + 1;
        }

        if (convs.containsKey(conv.id))
            throw new RuntimeException("Conversion " + conv.id + " is already in the ConversionList.");

        convs.put(conv.id, conv);
        if (conv.getSlot() < 0)
            allocateSlot(conv);
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

        // Assign a new ID and slot (copyTo copied the original's, which add() would otherwise
        // keep as-is) and add the copied conversion to the list
        newConv.setID(getFreeKey());
        newConv.setSlot(-1);
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
		Conversion conv = convs.get(key);
		if (conv != null)
		    releaseSlot(conv);
		if (trashCan.size() < MAX_TRASH_CAN_SIZE)
		    trashCan.add(conv);
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
		bySlot.clear();
		freeSlots.clear();
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
        other.nextConversionId = Math.max(other.nextConversionId, this.nextConversionId);

        for (Conversion conv : convs.values()) {
            Conversion convCopy = getNewOrRecycledConversion(conv.getID());
            conv.copyTo(convCopy);
            other.convs.put(convCopy.getID(), convCopy);

            int slot = convCopy.getSlot();
            while (other.bySlot.size() <= slot)
                other.bySlot.add(null);
            other.bySlot.set(slot, convCopy);
        }

        // Derive the free-slot list from occupancy rather than mirroring this list's free-list
        // order -- HashMap iteration order isn't meaningful, and this is simpler and just as
        // correct. (other.bySlot/freeSlots must already be empty here: copyTo() is only ever
        // called right after other.clear().)
        other.freeSlots.clear();
        for (int i = 0; i < other.bySlot.size(); i++) {
            if (other.bySlot.get(i) == null)
                other.freeSlots.push(i);
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
            acg.startEditingConversion();
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
