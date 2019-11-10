/**
 * 
 */
package contactrees.model;

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.ConversionGraph;

/**
 * 
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class ConversionMovePrior extends Distribution {

	final public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The conversion graph containing the conversion edges.",
			Input.Validate.REQUIRED);

	final public Input<BlockSet> blockSetInput = new Input<>(
			"blockSet",
			"The blockSet containing all blocks that evolve along the clonal frame and "
			+ "move over conversion edges in the network.",
			Input.Validate.REQUIRED);
	
	final public Input<RealParameter> pMoveInput = new Input<>(
			"pMove",
			"The rate at which a pair of lineages will get in contact and form a conversion.",
			Input.Validate.REQUIRED);

	ConversionGraph acg;
	BlockSet blockSet;
	double pMove;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		acg = networkInput.get();
		blockSet = blockSetInput.get();
		pMove = pMoveInput.get().getValue();
	}
	
	@Override
	public double calculateLogP() {
		logP = 0.0;
		int moveCount = blockSet.countMoves();
		int n = acg.getConvCount() * blockSet.getBlockCount();
		
		for (Block b : blockSet) {
		    for (int cID : b.getConversionIDs()) {
		        assert acg.getConversions().getKeys().contains(cID) : cID;
		    }
		}
		
        if (pMove == 0.) {
        	assert moveCount == 0;
        	return 0;
        }
		
		logP = moveCount * Math.log(pMove) + (n-moveCount) * Math.log(1-pMove);
		
		return logP;
	}
	
	@Override
    protected boolean requiresRecalculation() {
		// For now we use the safe version (always recalculate)
		return true;
		// TODO: Use the version below when sure that dirty logic is fine in ACG and BlockSet.
		// return acg.somethingIsDirty() || blockSet.somethingIsDirty();
    }
	
	@Override
	public List<String> getArguments() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public List<String> getConditions() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public void sample(State state, Random random) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

}
