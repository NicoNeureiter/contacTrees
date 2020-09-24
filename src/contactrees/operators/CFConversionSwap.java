package contactrees.operators;

import contactrees.Block;
import contactrees.Conversion;
import contactrees.util.Util;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;

/**
 * Swap attachment points of a subtree and a conversion, keeping heights constant.
 * This can preserve the topology of all marginal trees, 
 * but changes marginal node heights (depends on the sampled block moves). 
 *
 * @author Nico Neureiter
 */
@Description("CF/conversion swap operator.")
public class CFConversionSwap extends CFOperator {

    Input<Boolean> resampleMovesInput = new Input<>(
            "resampleMoves",
            "Whether to sample new blocks to be moved over the conversion (set to ´true´) or simply invert the previos moves (set to ´false´)",
            false);
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public double proposal() {
        double logHGF = 0.0;

        // Determine whether we can apply this operator:
        if (acg.getLeafNodeCount()<3 || acg.getConvCount()==0)
            return Double.NEGATIVE_INFINITY;

        // Acquire list of conversions compatible with swap:
        Conversion conv = acg.getConversions().getRandomConversion();

        // Cancel if conv is a cherry-edge
        if (conv.getNode1().getParent() == conv.getNode2().getParent())
            return Double.NEGATIVE_INFINITY;
        
        Node child = conv.getNode1();
        Node parent = child.getParent();
        Node sibling = Util.getOtherChild(parent, child);
        Node newSibling = conv.getNode2();
        
        // Cancel if a speciation happens below parent and above conv.height 
        if (sibling.getHeight() > conv.getHeight())
            return Double.NEGATIVE_INFINITY;

        // Cancel if speciation happens above conv and below parent.height
        if (newSibling.getParent().getHeight() < parent.getHeight())
            return Double.NEGATIVE_INFINITY;

        // No nodes in the way  =>  We can simply:
        //    * Re-attach parent above new sibling
        //    * Re-attach conv above old sibling
        
        // Re-attach parent above new sibling
        Node grandParent = parent.getParent();
        Node newGrandParent = newSibling.getParent();
        
        grandParent.removeChild(parent);
        grandParent.addChild(sibling);
        newGrandParent.removeChild(newSibling);
        newGrandParent.addChild(parent);
        parent.removeChild(sibling);
        parent.addChild(newSibling);
        
        assert Util.getSibling(child) == newSibling;
        assert parent.getParent() == newGrandParent;
        assert sibling.getParent() == grandParent;
        assert sibling != child;
        
        // Re-attach conv above old sibling
        conv.setNode2(sibling);
        
        // Move other conversions from parent to sibling
        for (Conversion c : acg.getConversions()) {
            if (c.getNode1() == parent)
                c.setNode1(sibling);
            if (c.getNode2() == parent)
                c.setNode2(sibling);
        }
        
        // Update conversion which were previously attached to "newSibling", but are now above "parent".
        for (Conversion c : acg.getConversions()) {
            if (c == conv)
                continue;
            if (c.getNode1() == newSibling && c.getHeight() > parent.getHeight())
                c.setNode1(parent);
            if (c.getNode2() == newSibling && c.getHeight() > parent.getHeight())
                c.setNode2(parent);
        }
        
        if (resampleMovesInput.get()) {
            logHGF += drawAffectedBlocksGibbs(conv, true);
        } else {
            invertAffectedBlocks(conv);
        }
        
        assert !acg.isInvalid() : "CFCS proposed invalid state.";

        return logHGF;
    }
    
    /**
     * Invert all block-moves of the given conversion (i.e. remove the present ones and add all other block moves).
     * @param conv
     */
    protected void invertAffectedBlocks(Conversion conv) {
        for (Block block : blockSet) {
            if (block.isAffected(conv)) {
                block.removeMove(conv);
            } else {
                block.addMove(conv);
            }
        }
    }
}
