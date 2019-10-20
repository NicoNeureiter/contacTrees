package contactrees;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.alignment.FilteredAlignment;


public class MarginalTree extends Tree {

	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The conversion graph (network) from which the marginal tree is sampled.",
			Input.Validate.REQUIRED);
	public Input<Block> movesInput = new Input<>(
			"moves",
			"The moves (conversion edges) this marginal tree follows along the conversion graph.",
			Input.Validate.REQUIRED);

	public ConversionGraph network;
	public Block moves;
	
	public void initAndValidate() {
		network = networkInput.get();
		moves = movesInput.get();
		
		update();
	}
	
	public void update() {
		// Build marginal tree by following <moves> in the <network>.
	}
	
}


