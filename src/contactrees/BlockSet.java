package contactrees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("A container for multiple Block objects.")
public class BlockSet extends CalculationNode {

	final public Input<List<Block>> blocksInput = new Input<>(
			"block",
			"A containter for multiple Block objects.",
			new ArrayList<>());
	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The network on which the conversion moves are defined.",
			Input.Validate.REQUIRED);
	
	List<Block> blocks;
	ConversionGraph acg;
	int blockCount;
	
	@Override
	public void initAndValidate() {
		blocks = blocksInput.get();	
		acg = networkInput.get();
		blockCount = blocks.size(); 
	}
	
	/**
	 * Obtain the block currently affected by the given conversion.
	 * @return List of affected blocks.
	 */
	public List<Block> getAffectedBlocks(Conversion conv){
		List<Block> affectedBlocks = new ArrayList<>();
		
		for (Block block : blocks) {
			if (block.convs.contains(conv)) {
				affectedBlocks.add(block);
			}
		}
		return affectedBlocks;
	}
	
	/**
	 * Obtain a hash-map, mapping each conversion to the currently affected blocks.
	 * @return Map from conversion to affected blocks.
	 */
	public HashMap<Conversion, List<Block>> getAffectedBlocks(){
		HashMap<Conversion, List<Block>> affectedBlocksByConv = new HashMap<>();
		
		for (Conversion conv : acg.convs) {
			affectedBlocksByConv.put(conv, getAffectedBlocks(conv));
		}
		
		return affectedBlocksByConv;
	}
	
	/**
	 * Count the number of blocks affected by the given conversion.
	 * @return Number of affected blocks.
	 */
	public int countAffectedBlocks(Conversion conv){
		return getAffectedBlocks(conv).size();
	}
	
	/**
     * Count the number of blocks affected by any conversion.
     * 
     * @return Number of blocks
     */
    public int countAffectedBlocks() {
    	int affectedCount = 0;
    	for (Block block : blocks) {
    		if (!block.convs.isEmpty())
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
        HashSet<Conversion> uselessConvs = new HashSet<>(acg.convs);
        
        for (Block block : blocks) {
	        for (Conversion conv : block.convs) {
	            uselessConvs.remove(conv);
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

}
