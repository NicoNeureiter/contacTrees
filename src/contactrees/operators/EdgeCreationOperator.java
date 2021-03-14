/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package contactrees.operators;

import java.util.List;
import java.util.Set;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.model.ACGDistribution;

/**
 * Abstract class of ACG operators that add new conversion
 * edges to an existing ConversionGraph.
 *
 * @author Nico Neureiter
 */
public abstract class EdgeCreationOperator extends ACGOperator {

	final public Input<BlockSet> blockSetInput = new Input<>(
			"blockSet",
			"Block of site which are either inherited or passed via a conversion edge.",
			Input.Validate.REQUIRED);

	public Input<RealParameter> conversionRateInput = new Input<>(
            "conversionRate",
            "Rate at which conversions happen along pairs of edges on the clonal frame.",
            Input.Validate.REQUIRED);

	public Input<ACGDistribution> networkPriorInput = new Input<>(
	        "networkPrior",
	        "The network prior defines how new edges should be sampled.",
	        Input.Validate.REQUIRED);

	protected BlockSet blockSet;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        blockSet = blockSetInput.get();
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
        return networkPriorInput.get().attachEdge(conv);
    }

    /**
     * Retrieve probability density for both attachment points of the given
     * recombinant edge.
     *
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeAttachmentProb(Conversion conv) {
        return networkPriorInput.get().getEdgeAttachmentProb(conv);
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

    /**
     * Include all block StateNodes in the list of affected state nodes.
     * This does not work automatically, since blocks are only indirect
     * inputs through the blockSet.
     */
    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = super.listStateNodes();
        for (Block block : blockSet.getBlocks()) {
            list.add(block);
        }
        return list;
    }

}
