/**
 * 
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Node;

import beast.core.Description;
import beast.core.StateNode;
import contactrees.util.Util;


/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("A block which is affected by a marginal tree and defined through a set of conversions wihtin a network.")
public class Block extends StateNode { 

	// The indices (in the ConversionList) of the conversions affecting this block. 
	List<Integer> convIds, convIdsStored;

	@Override
	public void initAndValidate() {
		convIds = new ArrayList<Integer>();
        convIdsStored = new ArrayList<>();
	}

	/**
	 * Add a conversion over which this block moves in its history.
	 * @param The new conversion.
	 */
	public void addMove(Conversion conv) {
		convIds.add(conv.getID());
	}
	
	
	/**
	 * Remove a conversion move from this blocks history.
	 * @param The conversion to be removed.
	 */
	public void removeMove(Conversion conv) {
		convIds.remove(conv.getID());
	}

	/**
	 * Evaluate whether this block was affected by the given conversion.
	 * @param The conversion in question.
	 * @return "true" iff the block moved over the conversion.
	 */
	public boolean isAffected(Conversion conv) {
		// TODO resolve overshadowed conversions somewhere
		return convIds.contains(conv.getID()); 
	}
	
	/**
	 * Number of conversions over which this block moves.
	 * @return Number of conversions.
	 */
	public int size() {
		return convIds.size();
	}

	
	/*
	 * Loggable implementations 
	 */
	
	/**
	 * init() and close() can stay empty.
	 */
	@Override
	public void init(PrintStream out) {}
	@Override
	public void close(PrintStream out) {}
	
	@Override
	public int getDimension() {
		return 1;
	}

	/*
	 * StateNode implementations
	 */
	@Override
	public Block copy() {
		Block copy = new Block();
		
		copy.setID(getID());
		copy.convIds = Util.deepCopyIntegers(convIds);
		copy.convIdsStored = Util.deepCopyIntegers(convIdsStored);
		
		return copy;
	}

	@Override
	public void assignTo(StateNode other) {
		final Block block = (Block) other;
		block.setID(getID());
		block.convIds = Util.deepCopyIntegers(convIds);
		block.convIdsStored = Util.deepCopyIntegers(convIdsStored);
	}

	@Override
	public void assignFrom(StateNode other) {
		final Block block = (Block) other;
		setID(block.getID());
		convIds = Util.deepCopyIntegers(block.convIds);
		convIdsStored = Util.deepCopyIntegers(block.convIdsStored);
	}

	@Override
	protected void store() {
		convIdsStored = Util.deepCopyIntegers(convIds);
	}

	@Override
	public void restore() {
		List<Integer> tmp = convIds;
		convIds = convIdsStored;
		convIdsStored = tmp;
		
//		ConversionList convList = acg.getConversions();
//		List<Integer> toStore = convIdsStored;
//		
//		// Remove all active conversion that are not stored
//		for (int id : convIds) {
//			if (!convIdsStored.contains(id)) {
//				toStore.add(id);
//				removeConversion(convList.get(id));
//			}
//		}
//		
//		// Add all stored conversions that are not active
//		for (int id : convIdsStored) {
//			if (!convIds.contains(id)) {
//				addConversion(convList.get(id));
//			}
//		}
//		
//		// Add the removed conversions to the stored ones
//		for (int id : toStore) {
//			convIdsStored.add(id);
//		}
	}

	@Override
	public void assignFromFragile(StateNode other) {
		assignFrom(other);
	}
	
	@Override
	public double getArrayValue(int dim) {
		return convIds.size();
	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
		setSomethingIsDirty(isDirty);
	}

	@Override
	public void fromXML(Node node) {
		// TODO Auto-generated method stub	
	}

	@Override
	public int scale(double scale) {
		throw new UnsupportedOperationException("StateNode method scale() not applicable to Block."); 
	}	
	
}
