package contactrees.operators;

import contactrees.Conversion;
import contactrees.util.Util;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.HashSet;

/**
 * @author Nico Neureiter
 */
@Description("Cause recombinant edge to hop between clonal frame edges.")
public class ConvertedEdgeHopGibbs extends ConversionCreationOperator {

    public ConvertedEdgeHopGibbs() { }
    
    @Override
    public double proposal() {
        double logHGF = 0.;
        
        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;
        
        // Select recombination at random
        Conversion conv = chooseConversion();
        
        // Compute back probability
        logHGF += getAffectedBlocksProbGibbs(conv);

        // Choose whether to move departure or arrival point
        boolean moveDeparture = conv.getNode2().isRoot() || Randomizer.nextBoolean();

        double height = conv.getHeight();
        Node convNode = moveDeparture ? conv.getNode1() : conv.getNode2();

        // Find list of CF edges alive at pointHeight
        HashSet<Node> activeLineages = acg.getLineagesAtHeight(height);
        activeLineages.remove(conv.getNode1());
        activeLineages.remove(conv.getNode2());

        if (activeLineages.isEmpty())
            return Double.NEGATIVE_INFINITY;
        
        Node newNode = Util.sampleFrom(activeLineages);
        
        // Select new attachment point:
        if (moveDeparture)
            conv.setNode1(newNode);
        else
            conv.setNode2(newNode);

        assert !acg.isInvalid() : "CEHContemp produced invalid state.";
        
        // Resample affected blocks
        logHGF -= drawAffectedBlocksGibbs(conv);

        return logHGF;
    }
    
}
