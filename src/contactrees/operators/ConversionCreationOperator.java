package contactrees.operators;

import java.util.List;
import java.util.Set;

import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.model.ACGDistribution;
import contactrees.model.ConversionPrior;

/**
 * Abstract class of ACG operators that add new converted edges
 * and their affected sites to an existing ConversionGraph.
 *
 * @author Nico Neureiter
 */
public abstract class ConversionCreationOperator extends BorrowingOperator {

    public Input<RealParameter> conversionRateInput = new Input<>(
            "conversionRate",
            "Rate at which conversions happen along pairs of edges on the clonal frame.",
            Input.Validate.REQUIRED);

    public Input<ACGDistribution> networkPriorInput = new Input<>(
            "networkPrior",
            "The network prior (or conversion prior) defines how new edges should be sampled.");

    public Input<ConversionPrior> conversionPriorInput = new Input<>(
            "conversionPrior",
            "The conversion prior (or network prior) defines how new edges should be sampled.",
            Input.Validate.XOR, networkPriorInput);

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    /**
     * Add and return a new conversion edge.
     * @return The new conversion.
     */
    protected Conversion addNewConversion() {
        return acg.addNewConversion();
    }

    /**
     * Remove the specified conversion
     * @param The conversion to be removed
     */
    protected void removeConversion(Conversion conv) {
        acg.removeConversion(conv);
        blockSet.removeConversion(conv);
    }

    /**
     * Attach chosen recombination to the clonal frame.  Note that only the
     * attachment points (nodes and heights) are set, the affected region of
     * the alignment is not modified.
     *
     * @param conv conversion
     * @return log probability density of chosen attachment.
     */
    public double attachEdge(Conversion conv) {
        if (networkPriorInput.get() != null)
            return networkPriorInput.get().attachEdge(conv);
        else
            return conversionPriorInput.get().attachEdge(conv);
    }

    /**
     * Retrieve probability density for both attachment points of the given
     * recombinant edge.
     *
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeAttachmentProb(Conversion conv) {
        if (networkPriorInput.get() != null)
            return networkPriorInput.get().getEdgeAttachmentProb(conv);
        else
            return conversionPriorInput.get().getEdgeAttachmentProb(conv);
    }

    /**
     * Take a recombination with an existing departure point and determine
     * the arrival point by allowing it to coalesce with the clonal frame.
     *
     * @param conv recombination to modify
     * @return log probability density of coalescent point chosen.
     */
    public double coalesceEdge(Conversion conv) {
        double height = conv.getHeight();
        double logP = 0.0;

        // Find the other lineages at the same height as node2.
        Set<Node> activeLineages = acg.getLineagesAtHeight(height);
        activeLineages.remove(conv.getNode1());

        // Sample a second node uniformly at random
        int choice = Randomizer.nextInt(activeLineages.size());
        int i = 0;
        for (Node node : activeLineages) {
            if (i == choice) {
                conv.setNode2(node);
                break;
            }
            i++;
        }

        // The only random choice was the lineage
        logP -= Math.log(activeLineages.size());

        return logP;
    }

    @Override
    public List<StateNode> listStateNodes() {
    	List<StateNode> list = super.listStateNodes(); 
        list.remove(conversionRateInput.get());
    	return list;
    }
}
