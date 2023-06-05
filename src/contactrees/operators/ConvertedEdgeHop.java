package contactrees.operators;

import java.util.HashSet;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.util.Util;

/**
 * @author Nico Neureiter
 */
@Description("Cause recombinant edge to hop between clonal frame edges.")
public class ConvertedEdgeHop extends ACGOperator {

    public ConvertedEdgeHop() { }

    @Override
    public double proposal() {

        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;

        // Select recombination at random
        Conversion conv = chooseConversion();

        // Find list of CF edges alive at pointHeight
        HashSet<Node> activeLineages = acg.getLineagesAtHeight(conv.getHeight());
        activeLineages.remove(conv.getNode1());
        activeLineages.remove(conv.getNode2());

        if (activeLineages.isEmpty())
            return Double.NEGATIVE_INFINITY;

        Node newNode = Util.sampleFrom(activeLineages);

        // Choose whether to move departure or arrival point
        boolean moveDeparture = conv.getNode2().isRoot() || Randomizer.nextBoolean();
        // Select new attachment point:
        if (moveDeparture)
            conv.setNode1(newNode);
        else
            conv.setNode2(newNode);

        assert !acg.isInvalid() : "CEHContemp produced invalid state.";

        return 0.0;
    }

}
