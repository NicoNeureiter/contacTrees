package contactrees.operators;

import java.util.HashSet;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.MarginalTree;
import contactrees.util.Util;

/**
 * @author Nico Neureiter
 */
@Description("Cause recombinant edge to hop between clonal frame edges.")
public class ConvertedEdgeHopGibbs extends ConversionCreationOperator {

    public Input<Boolean> sourceOnlyInput = new Input<>("sourceOnly",
                                                        "Only move the source end (node2) of the conversion edge.",
                                                        false);

    public Input<Integer> nClosestRelativesInput = new Input<>("nClosestRelatives",
                                                        "The maximum number of closest relatives (lowest tMRCA to the oringinal attachment point) to consider for the hop.",
                                                        -1);

    public Input<Double> maxHopDistanceInput = new Input<>("maxHopDistance",
            "The maximum distance from the original attachment point to the new one, measured in time to their MRCA.",
            Double.MAX_VALUE);

    public ConvertedEdgeHopGibbs() { }

    @Override
    public double proposal() {
        double logHGF = 0.;

        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;

        for (TreeLikelihood tLH : treeLHsInput.get()) {
            tLH.store();
            ((MarginalTree) tLH.treeInput.get()).store();
        }

        // Select recombination at random
        Conversion conv = chooseConversion();
        double height = conv.getHeight();

        // Compute back probability
        logHGF += getBorrowingsProbGibbs(conv, true);

        // Choose whether to move departure or arrival point
        boolean moveDeparture = Randomizer.nextBoolean();
        if (sourceOnlyInput.get()) {
            moveDeparture = false;
        }
        Node nodeToMove = moveDeparture ? conv.getNode1() : conv.getNode2();
        Node nodeToStay = moveDeparture ? conv.getNode2() : conv.getNode1();

        // Find list of CF edges alive at pointHeight
        HashSet<Node> candidates;
        int nClosest = nClosestRelativesInput.get();
        if (nClosest > 0) {
            candidates = Util.getClosestRelatives(nodeToMove, height, nClosest+1);
            // We add 1, so that we can remove the other end of the conversion from the candidates.
            // TODO find a more elegant way to do this

        } else {
            // Non-positive numbers are interpreted as no constraint on the candidate lineages
            candidates = acg.getLineagesAtHeight(conv.getHeight());
        }

        candidates.remove(conv.getNode1());
        candidates.remove(conv.getNode2());

        if (candidates.isEmpty())
            return Double.NEGATIVE_INFINITY;


        Node newNode = Util.sampleFrom(candidates);
        logHGF -= Math.log(1. / candidates.size());

        // Select new attachment point:
        if (moveDeparture)
            conv.setNode1(newNode);
        else
            conv.setNode2(newNode);

        assert !acg.isInvalid() : "ConvertedEdgeHopGibbs produced invalid state.";

        // Resample affected blocks
        logHGF -= drawBorrowingsGibbs(conv, true);

        if (nClosest > 0) {
            HashSet<Node> candidatesBack = Util.getClosestRelatives(newNode, height, nClosest+1);

            // Directly reject if there is no back-move (nodeToMove is no neighbor of newNode)
            if (!candidatesBack.contains(nodeToMove))
                return Double.NEGATIVE_INFINITY;

            candidatesBack.remove(newNode);
            candidatesBack.remove(nodeToStay);

            logHGF += Math.log(1. / candidatesBack.size());
        } else {
            logHGF += Math.log(1. / candidates.size());
        }

        for (TreeLikelihood tLH : treeLHsInput.get()) {
            tLH.restore();
            ((MarginalTree) tLH.treeInput.get()).restore();
        }

        return logHGF;
    }


}
