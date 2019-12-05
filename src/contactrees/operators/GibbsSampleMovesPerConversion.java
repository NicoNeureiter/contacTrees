/**
 * 
 */
package contactrees.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.evolution.likelihood.TreeLikelihood;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.MarginalTree;

/**
 * 
 * 
 * @author Nico Neureiter
 */
public class GibbsSampleMovesPerConversion extends GibbsBlockMovesOperator {

    final public Input<BlockSet> blockSetInput = new Input<>(
            "blockSet",
            "Block set containing all blocks affected by the the given network.",
            Input.Validate.REQUIRED);

    public Input<List<TreeLikelihood>> treeLHsInput = new Input<>(
            "treeLikelihood",
            "BEASTObject computing the tree likelihood.",
            new ArrayList<TreeLikelihood>());

    protected BlockSet blockSet;
    protected List<TreeLikelihood> treeLHs;
    protected List<MarginalTree> marginalTrees;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        blockSet = blockSetInput.get();
        treeLHs = treeLHsInput.get();
    }
    
    @Override
    public double proposal() {
        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;
        
        // Sample conversion to change
        Conversion conv = acg.getConversions().getRandomConversion();
        
        // Iterate over all blocks and resample the move over the chosen conversion
        for (TreeLikelihood treeLH : treeLHs) {
            
            // Select marginalTree and block corresponding to treeLH
            MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
            Block block = marginalTree.block;
            assert blockSet.getBlocks().contains(block);
            
            sampleBlockMove(block, conv, marginalTree, treeLH);
        }
        
        return Double.POSITIVE_INFINITY;
    }
    
}
