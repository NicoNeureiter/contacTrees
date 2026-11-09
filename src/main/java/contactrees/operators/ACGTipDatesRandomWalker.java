package contactrees.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.operator.TipDatesRandomWalker;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.MarginalTree;

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

        double[] oldConvHeights = new double[attachedConvs.size()];
        for (int j = 0; j < attachedConvs.size(); j++)
            oldConvHeights[j] = attachedConvs.get(j).getHeight();

        acg.beginProposal();
        acg.startEditing(this);
        node.setHeight(newHeight);
        for (int j = 0; j < attachedConvs.size(); j++)
            attachedConvs.get(j).setHeight(newConvHeights.get(j));

        updateMarginalTrees(node, attachedConvs, oldConvHeights);

        return logHR;
    }

    /**
     * Update the marginal trees after the tip move. The attached conversions stay on their
     * branches, so a tree keeps its topology unless one of them changed its order with another
     * conversion of the tree's block on a shared branch; then only heights change (the tip and
     * the attached conversions the block follows). Otherwise the tree is reswept.
     */
    protected void updateMarginalTrees(Node tip, List<Conversion> attached, double[] oldHeights) {
        List<Conversion> followed = new ArrayList<>();
        for (MarginalTree tree : acg.getMarginalTrees()) {
            followed.clear();
            for (Conversion conv : attached)
                if (tree.block.isAffected(conv))
                    followed.add(conv);

            if (keepsTopology(tree, followed, attached, oldHeights)) {
                tree.updateCreatorHeight(tip);
                for (Conversion conv : followed)
                    tree.updateCreatorHeight(conv);
                tree.markUpToDate();
            } else {
                tree.ensureUpToDate();
            }
        }
    }

    /**
     * @return true iff every attached conversion the tree's block follows lies strictly inside its
     * branches and kept its order relative to the block's other conversions on those branches.
     */
    private boolean keepsTopology(MarginalTree tree, List<Conversion> followed, List<Conversion> attached, double[] oldHeights) {
        for (Conversion conv : followed) {
            if (!isStrictlyInsideBranches(conv))
                return false;
            double oldConv = oldHeights[attached.indexOf(conv)];
            Node n1 = conv.getNode1(), n2 = conv.getNode2();
            for (int id : tree.block.getConversionIDs()) {
                if (id == conv.getID())
                    continue;
                Conversion other = acg.getConversions().get(id);
                Node o1 = other.getNode1(), o2 = other.getNode2();
                if (o1 != n1 && o1 != n2 && o2 != n1 && o2 != n2)
                    continue;
                int j = attached.indexOf(other);
                double oldOther = (j >= 0) ? oldHeights[j] : other.getHeight();
                int before = Double.compare(oldConv, oldOther);
                int after = Double.compare(conv.getHeight(), other.getHeight());
                if (before == 0 || before != after)
                    return false;
            }
        }
        return true;
    }

    /**
     * @return true iff the conversion lies strictly between the bottom and the top of both its
     * branches (no tie with a CF event on them).
     */
    private static boolean isStrictlyInsideBranches(Conversion conv) {
        double h = conv.getHeight();
        Node n1 = conv.getNode1(), n2 = conv.getNode2();
        return h > n1.getHeight() && h > n2.getHeight()
                && (n1.isRoot() || h < n1.getParent().getHeight())
                && (n2.isRoot() || h < n2.getParent().getHeight());
    }

}
