package contactrees.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.MarginalTree;

/**
 * Uniform operator for clonal frame nodes. This operator is capable of
 * shifting an internal CF node past conversions, getting better acceptance
 * rates than the standard uniform operator when a large number of conversions
 * is present.
 *
 * @author Nico Neureiter
 */
@Description("Uniform operator for clonal frame nodes.")
public class CFUniform extends CFOperator {

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Root height proposal parameter.", 0.8);

    @Override
    public double doProposal() {

        double logHGF = 0.0;

        // Select internal non-root node at random.
        Node node = acg.getNode(acg.getLeafNodeCount()
                + Randomizer.nextInt(acg.getInternalNodeCount()));
        // Choice of height is symmetric -> no effect on HGF

        Node leftChild = node.getLeft();
        Node rightChild = node.getRight();

        double oldHeight = node.getHeight();
        double maxChildHeight = Math.max(leftChild.getHeight(), rightChild.getHeight());

        // Choose new height:
        double newHeight;
        if (node.isRoot()) {
            double fMin = Math.min(scaleFactorInput.get(), 1.0/scaleFactorInput.get());
            double fMax = 1.0/fMin;

            double f = Randomizer.uniform(fMin, fMax);
            newHeight = node.getHeight() * f;
            logHGF += Math.log(1.0/f);

            if (newHeight < maxChildHeight)
                return Double.NEGATIVE_INFINITY;
        } else {
            Node parent = node.getParent();
            newHeight = Randomizer.uniform(maxChildHeight, parent.getHeight());
            // Choice of height is symmetric -> no effect on HGF
        }

        // If no conversion is attached to the branches around the node within the range of the
        // move (and none is added or removed), the move is a pure height change: the order of
        // events along every branch stays the same, so the marginal trees keep their topology.
        boolean pureHeightChange = !hasConversionNear(node,
                Math.min(oldHeight, newHeight), Math.max(oldHeight, newHeight));
        int convCountBefore = acg.getConvCount();

        if (newHeight>oldHeight) {
            logHGF -= expandConversions(leftChild, rightChild, newHeight);
        } else {
            logHGF += collapseConversions(leftChild, rightChild, newHeight);
        }

        if (logHGF > Double.NEGATIVE_INFINITY) {
            assert !acg.isInvalid() : "CFUniform proposed invalid state.";

            if (pureHeightChange && acg.getConvCount() == convCountBefore) {
                for (MarginalTree tree : acg.getMarginalTrees()) {
                    tree.updateCreatorHeight(node);
                    tree.markUpToDate();
                }
            }
            // Otherwise the marginal trees resweep themselves when they are next used.
        }

        return logHGF;
    }

    /**
     * @return true iff a conversion attached to the branch above the node or above one of its
     * children lies within [lower, upper].
     */
    private boolean hasConversionNear(Node node, double lower, double upper) {
        Node left = node.getLeft();
        Node right = node.getRight();
        for (Conversion conv : acg.getConversions()) {
            double h = conv.getHeight();
            if (h < lower || h > upper)
                continue;
            Node n1 = conv.getNode1();
            Node n2 = conv.getNode2();
            if (n1 == node || n1 == left || n1 == right || n2 == node || n2 == left || n2 == right)
                return true;
        }
        return false;
    }

}
