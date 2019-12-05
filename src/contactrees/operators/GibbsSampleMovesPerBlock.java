/**
 * 
 */
package contactrees.operators;

import beast.core.Input;
import beast.evolution.likelihood.TreeLikelihood;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalTree;

/**
 * 
 * 
 * @author Nico Neureiter
 */
public class GibbsSampleMovesPerBlock extends GibbsBlockMovesOperator {

    final public Input<Block> blockInput = new Input<>(
            "block",
            "Block of sites which are either inherited or passed via a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<MarginalTree> marginalTreeInput = new Input<>(
            "marginalTree",
            "MarginalTree defined by the moves of the given block.",
            Input.Validate.REQUIRED);

    public Input<TreeLikelihood> treeLHInput = new Input<>(
            "treeLikelihood",
            "BEASTObject computing the tree likelihood.",
            Input.Validate.REQUIRED);

    protected Block block;
    protected TreeLikelihood treeLH;
    protected MarginalTree marginalTree;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        block= blockInput.get();
        marginalTree = marginalTreeInput.get();
        treeLH = treeLHInput.get();
    }
    
    @Override
    public double proposal() {
        for (Conversion conv : acg.getConversions()) {
            sampleBlockMove(block, conv, marginalTree, treeLH);
        }
        
        return Double.POSITIVE_INFINITY;
    }
   
    
}
