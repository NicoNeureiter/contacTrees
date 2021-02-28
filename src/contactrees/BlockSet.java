package contactrees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;

/**
 * @author Nico Neureiter
 */
@Description("A container for multiple Block objects.")
public class BlockSet extends CalculationNode implements Iterable<Block> {

	final public Input<ArrayList<Block>> blocksInput = new Input<>(
			"block",
			"A containter for multiple Block objects.",
			new ArrayList<>());
	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The network on which the conversion moves are defined.");
	public Input<Boolean> deferNetworkSpecificationInput = new Input<>(
	        "deferNetworkSpecification",
	        "Set this flag to allow deferring the specification of the network input (usually to be set in ACGWithBlocks).",
	        Input.Validate.XOR, networkInput);

	protected ArrayList<Block> blocks;
	protected ConversionGraph acg;

	@Override
	public void initAndValidate() {
		blocks = blocksInput.get();
		for (Block block : blocks)
		    block.initAndValidate();

        if (networkInput.get() != null)
            acg = networkInput.get();
        else
            assert deferNetworkSpecificationInput.get() == true;
	}

	public void setNetwork(ConversionGraph acg) {
	    this.acg = acg;
	}

	public ArrayList<Block> getBlocks(){
		return blocks;
	}

	public int getBlockCount() {
		return blocks.size();
	}

	@Override
    protected boolean requiresRecalculation() {
        return true;
//        for (Block block : this) {
//            if (block.somethingIsDirty())
//                return true;
//        }
//        return false;
    }

	public void removeConversion(Conversion conv) {
		for (Block block : blocks) {
			if (block.isAffected(conv)) {
				removeBlockMove(conv, block);
			}
		}

	}

	public void addBlockMove(Conversion conv, Block block) {
		block.addMove(conv);
	}

	public void removeBlockMove(Conversion conv, Block block) {
		block.removeMove(conv);
	}

	public void addBlockMove(Conversion conv, int blockID) {
		addBlockMove(conv, blocks.get(blockID));
	}

    public void addBlockMove(Conversion conv, String blockName) {
        addBlockMove(conv, getBlockByName(blockName));
    }

    public Block getBlockByName(String blockName) {
        for (Block block : blocks)
            if (block.getID().equals(blockName))
                return block;
        System.out.println("blockName: " + blockName);
        throw new RuntimeException("Block not found: " + blockName);
    }

	public void removeBlockMove(Conversion conv, int blockID) {
		removeBlockMove(conv, blocks.get(blockID));
	}

	/**
	 * Obtain the IDs of the blocks currently affected by the given conversion.
	 * @return List of affected blocks.
	 */

	public List<Integer> getAffectedBlockIDs(Conversion conv){
		List<Integer> affectedBlocks = new ArrayList<>();

		for (int i = 0; i<blocks.size(); i++) {
			if (blocks.get(i).isAffected(conv))
				affectedBlocks.add(i);
		}

		return affectedBlocks;
	}

    /**
     * Obtain the blocks currently affected by the given conversion.
     * @return List of affected blocks.
     */
    public List<Block> getAffectedBlocks(Conversion conv){
        List<Block> affectedBlocks = new ArrayList<>();
        for (Block block : this) {
            if (block.isAffected(conv))
                affectedBlocks.add(block);
        }
        return affectedBlocks;
    }

    /**
     * Obtain the IDs of the blocks currently affected by the given conversion.
     * @return List of affected blocks.
     */

    public List<String> getAffectedBlockNames(Conversion conv){
        List<String> affectedBlocks = new ArrayList<>();

        for (Block block : this) {
            if (block.isAffected(conv))
                affectedBlocks.add(block.getID());
        }

        return affectedBlocks;
    }

    /**
     * Obtain a hash-map, mapping each conversion to the IDs of the currently affected blocks.
     * @return Map from conversion to affected block IDs.
     */
    public HashMap<Conversion, List<Integer>> getAffectedBlockIDs(){
        HashMap<Conversion, List<Integer>> affectedBlocksByConv = new HashMap<>();

        assert acg != null;
        for (Conversion conv : acg.getConversions()) {
            affectedBlocksByConv.put(conv, getAffectedBlockIDs(conv));
        }

        return affectedBlocksByConv;
    }

    /**
     * Obtain a hash-map, mapping each conversion to the currently affected blocks.
     * @return Map from conversion to affected blocks.
     */
    public HashMap<Conversion, List<Block>> getAffectedBlocks(){
        HashMap<Conversion, List<Block>> affectedBlocksByConv = new HashMap<>();

        for (Conversion conv : acg.getConversions()) {
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
            if (block.countMoves() > 0)
                affectedCount++;
        }

        return affectedCount;
    }

    /**
     * Count the number of blocks affected by the given conversion.
     *
     * @conv the conversion edge.
     * @return Number of blocks
     */
    public int countAffectedBlocks(Conversion conv) {
        int affectedCount = 0;
        for (Block block : blocks) {
            if (block.isAffected(conv)) {
                affectedCount++;
            }
        }
        return affectedCount;
    }

	/**
     * Count the number of times a block is affected by a conversion edge.
     *
     * @return Number of moves
     */
    public int countMoves() {
        int moveCount = 0;

        for (Block block : blocks)
            moveCount += block.countMoves();

        return moveCount;
    }

    /**
     * Obtain the set of conversion edges which do not affect any blocks.
     *
     * @return The edges without an affected blocks.
     */
    public HashSet<Conversion> getUselessConversions() {
        HashSet<Conversion> uselessConvs = new HashSet<>(acg.getConversions().getConversions());

        for (Block block : blocks) {
	        for (int convId : block.getConversionIDs()) {
	            uselessConvs.remove(acg.getConversions().get(convId));
	        }
        }

        return uselessConvs;
    }

    /**
     * List the names of all blocks.
     *
     * @return list of block names
     */
    public String[] getBlockNames() {
        String[] blockNames = new String[size()];
        for (int i=0; i<size(); i++) {
            blockNames[i] = blocks.get(i).getID();
        }
        return blockNames;
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

	public int size() {
	    return blocks.size();
	}

	public Block get(int i) {
	    return blocks.get(i);
	}

	/*
	 * TESTING INTERFACE
	 */
	static public BlockSet create(ConversionGraph acg, ArrayList<Block> blocks) {
		BlockSet bs = new BlockSet();
		bs.initByName("network", acg);
//		bs.initAndValidate();
//		bs.acg = acg;
		bs.blocks = blocks;
		return bs;
	}

	static public BlockSet create(ConversionGraph acg) {
		return create(acg, new ArrayList<>());
	}

//    public BlockSet copy() {
//        return copy(acg);
//    }

    public BlockSet copy(ConversionGraph newACG) {
        ArrayList<Block> otherBlocks = new ArrayList<>();
        for (Block b : blocks) {
            otherBlocks.add(b.copy());
        }

        BlockSet other = new BlockSet();
        other.initByName("block", otherBlocks,
                         "network", newACG);

        return other;
    }

}
