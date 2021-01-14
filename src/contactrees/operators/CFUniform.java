package contactrees.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

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
    public double proposal() {

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

        if (newHeight>oldHeight) {
            logHGF -= expandConversions(leftChild, rightChild, newHeight);
        } else {
            logHGF += collapseConversions(leftChild, rightChild, newHeight);
        }

        if (logHGF > Double.NEGATIVE_INFINITY)
            assert !acg.isInvalid() : "CFUniform proposed invalid state.";

        return logHGF;
    }

}
