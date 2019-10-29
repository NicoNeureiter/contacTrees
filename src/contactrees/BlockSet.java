package contactrees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import contactrees.util.parsers.ExtendedNewickBaseVisitor;
import contactrees.util.parsers.ExtendedNewickLexer;
import contactrees.util.parsers.ExtendedNewickParser;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("A container for multiple Block objects.")
public class BlockSet extends CalculationNode implements Iterable<Block> {

	final public Input<ArrayList<Block>> blocksInput = new Input<>(
			"block",
			"A containter for multiple Block objects.",
			new ArrayList<>());
	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The network on which the conversion moves are defined.",
			Input.Validate.REQUIRED);
	
	protected ArrayList<Block> blocks;
	protected ConversionGraph acg;
	
	@Override
	public void initAndValidate() {
		blocks = blocksInput.get();	
		acg = networkInput.get();
	}
	
	public ArrayList<Block> getBlocks(){
		return blocks;
	}
	
	public int getBlockCount() {
		return blocks.size();
	}
	
	@Override
    protected boolean requiresRecalculation() {
		// TODO Remove this
		return true;
    }
	
	public void removeConversion(Conversion conv) {
		for (Block block : blocks) {
			if (block.convIds.contains(conv.id)) {
				block.convIds.remove(conv.id);
			}
		}
	}
	
	public void addBlockMove(Conversion conv, Block block) {
		block.convIds.add(conv.id);
	}
	
	public void removeBlockMove(Conversion conv, Block block) {
		block.convIds.remove(conv.id);
	}
	
	public void addBlockMove(Conversion conv, int blockID) {
		addBlockMove(conv, blocks.get(blockID));
	}
	
	public void removeBlockMove(Conversion conv, int blockID) {
		removeBlockMove(conv, blocks.get(blockID));
	}

	/**
	 * Obtain the blocks currently affected by the given conversion.
	 * @return List of affected blocks.
	 */

	public List<Integer> getAffectedBlocks(Conversion conv){
		List<Integer> affectedBlocks = new ArrayList<>();
		
		for (int i = 0; i<blocks.size(); i++) {
			if (blocks.get(i).convIds.contains(conv.id))
				affectedBlocks.add(i);
		}
		
		return affectedBlocks;
	}

	/**
	 * Obtain a hash-map, mapping each conversion to the currently affected blocks.
	 * @return Map from conversion to affected blocks.
	 */
	public HashMap<Conversion, List<Integer>> getAffectedBlocks(){
		HashMap<Conversion, List<Integer>> affectedBlocksByConv = new HashMap<>();
		
		for (Conversion conv : acg.convs) {
			affectedBlocksByConv.put(conv, getAffectedBlocks(conv));
		}

		return affectedBlocksByConv;
	}
	
	/**
     * Count the number of blocks affected by any conversion.
     * 
     * @return Number of blocks
     */
    public int countAffectedBlocks() {
    	int affectedCount = 0;
    	for (Block block : blocks) {
    		if (!block.convIds.isEmpty())
    			affectedCount++;
    	}
    	
    	return affectedCount;
    }

    /**
     * Obtain the set of conversion edges which do not affect any blocks.
     * 
     * @return The edges without an affected blocks.
     */
    public HashSet<Conversion> getUselessConversions() {
        HashSet<Conversion> uselessConvs = new HashSet<>(acg.convs.getConversions());
        
        for (Block block : blocks) {
	        for (int convId : block.convIds) {
	            uselessConvs.remove(acg.convs.get(convId));
	        }
        }
        
        return uselessConvs;
    }

    /**
     * Count the number of conversion edges which do not affect any block.
     * 
     * @return The number of edges without an affected block.
     */
    public int countUselessConversion() {
        return getUselessConversions().size();
    }

	@Override
	public Iterator<Block> iterator() {
		return blocks.iterator();
	}

	
	/*
	 * TESTING INTERFACE
	 */
	static public BlockSet getBlockSet(ConversionGraph acg, ArrayList<Block> blocks) {
		BlockSet bs = new BlockSet();
		bs.initAndValidate();
		bs.acg = acg;
		bs.blocks = blocks;
		return bs;
	}

	static public BlockSet getBlockSet(ConversionGraph acg) {
		return getBlockSet(acg, new ArrayList<>());
	}

	
}


