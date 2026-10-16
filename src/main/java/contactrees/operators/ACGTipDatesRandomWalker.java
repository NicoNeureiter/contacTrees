package contactrees.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.operator.TipDatesRandomWalker;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.ConversionGraph;

/**
 * @author Nico Neureiter
 */
@Description("Tip date random walk on a conversion graph. Conversion edges attached to the moved tip "
        + "keep their relative height within their range of valid heights.")
public class ACGTipDatesRandomWalker extends TipDatesRandomWalker {

    protected ConversionGraph acg;

    @Override
    public void initAndValidate() {
        if (!(treeInput.get() instanceof ConversionGraph))
            throw new IllegalArgumentException("ACGTipDatesRandomWalker requires a ConversionGraph as tree.");
        acg = (ConversionGraph) treeInput.get();
        super.initAndValidate();
    }

    @Override
    public double proposal() {
        // Randomly select leaf node
        int i = Randomizer.nextInt(taxonIndices.length);
        Node node = acg.getNode(taxonIndices[i]);

        double oldHeight = node.getHeight();
        double newHeight = oldHeight;
        if (useGaussian) {
            newHeight += Randomizer.nextGaussian() * windowSize;
        } else {
            newHeight += Randomizer.nextDouble() * 2 * windowSize - windowSize;
        }

        double parentHeight = node.getParent().getHeight();
        if (newHeight > parentHeight) {
            if (reflectValue) {
                // Same reflection bounds as in TipDatesRandomWalker
                newHeight = reflectValue(newHeight, Math.min(0.0, oldHeight), parentHeight);
            } else {
                return Double.NEGATIVE_INFINITY;
            }
        }
        if (newHeight == oldHeight)
            return Double.NEGATIVE_INFINITY;

        return moveTip(node, newHeight);
    }

    /**
     * Move the given tip to a new height and map every conversion attached to its
     * edge linearly from its old valid range [tMin, tMax] to the new one. The upper
     * bound tMax does not depend on the tip height, only tMin does.
     *
     * @param node tip to move
     * @param newHeight new height of the tip
     * @return log Jacobian of the conversion height transformation (log Hastings ratio)
     */
    public double moveTip(Node node, double newHeight) {
        double oldHeight = node.getHeight();
        double logHR = 0.0;

        // Compute all new conversion heights before changing the state
        List<Conversion> attachedConvs = new ArrayList<>();
        List<Double> newConvHeights = new ArrayList<>();
        for (Conversion conv : acg.getConversions()) {
            Node other;
            if (conv.getNode1() == node)
                other = conv.getNode2();
            else if (conv.getNode2() == node)
                other = conv.getNode1();
            else
                continue;

            double tMax = Math.min(node.getParent().getHeight(), other.getParent().getHeight());
            double tMinOld = Math.max(oldHeight, other.getHeight());
            double tMinNew = Math.max(newHeight, other.getHeight());

            // No valid position for this conversion at the new tip height
            if (tMinNew >= tMax)
                return Double.NEGATIVE_INFINITY;

            double relHeight = (conv.getHeight() - tMinOld) / (tMax - tMinOld);
            attachedConvs.add(conv);
            newConvHeights.add(tMinNew + relHeight * (tMax - tMinNew));
            logHR += Math.log((tMax - tMinNew) / (tMax - tMinOld));
        }

        acg.startEditing(this);
        node.setHeight(newHeight);
        for (int j = 0; j < attachedConvs.size(); j++)
            attachedConvs.get(j).setHeight(newConvHeights.get(j));

        return logHR;
    }

}
