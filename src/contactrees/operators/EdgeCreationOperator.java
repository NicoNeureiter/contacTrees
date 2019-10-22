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

import contactrees.BlockSet;
import contactrees.Conversion;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import java.util.Set;

/**
 * Abstract class of ACG operators that add new conversion
 * edges to an existing ConversionGraph.
 *
 * @author Nico Neureiter <nico.neureiter@gmail.com>
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
    	Conversion conv = acg.addNewConversion();
    	blockSet.addConversion(conv);
    	return conv;
    }
    
    /**
     * Remove the specified conversion
     * @param The conversion to be removed
     */
    protected void removeConversion(Conversion conv) {
    	acg.deleteConversion(conv);
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
        
        double logP = 0.0;
        
        // Select departure point
        double u = Randomizer.nextDouble()*acg.getClonalFrameLength();
        logP += Math.log(1.0/acg.getClonalFrameLength());
        
        for (Node node : acg.getNodesAsArray()) {
            if (node.isRoot())
                continue;
            
            if (u<node.getLength()) {
                conv.setHeight(node.getHeight() + u);
                conv.setNode1(node);
//                System.out.print(conv.getHeight());
//                System.out.print("  \t");
//                System.out.println(acg.getRoot().getHeight());
                break;
            } else
                u -= node.getLength();
        }
        
        // Select arrival point
        logP += coalesceEdge(conv);
        
        return logP;
    }
    
    /**
     * Retrieve probability density for both attachment points of the given
     * recombinant edge.
     * 
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeAttachmentProb(Conversion conv) {
        double logP = 0.0;
        
        logP += Math.log(1.0/acg.getClonalFrameLength());
        logP += getEdgeCoalescenceProb(conv);
        
        return logP;
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
        
    	// Find the other lineages at the same height as node1. 
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
     * Get probability density for the arrival time of the given recombinant
     * edge under ClonalOrigin's coalescent model.
     * 
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeCoalescenceProb(Conversion conv) {
    	int nLineages = acg.countLineagesAtHeight(conv.getHeight());
    	return - Math.log(nLineages - 1);

        
//        List<CFEventList.Event> events = acg.getCFEvents();
//        
//        // Find event immediately below departure point
//        int startIdx = 0;
//        while (events.get(startIdx+1).getHeight()<conv.getHeight1())
//            startIdx += 1;
//        
//        // Compute probability of edge length and arrival
//        for (int i=startIdx; i<events.size() && events.get(i).getHeight()<conv.getHeight2(); i++) {           
//            double t1 = Math.max(conv.getHeight1(), events.get(i).getHeight());
//            double t2 = conv.getHeight2();
//            if (i<events.size()-1)
//                t2 = Math.min(t2, events.get(i+1).getHeight());
//        
//            double intervalArea = popFunc.getIntegral(t1, t2);
//            logP += -intervalArea*events.get(i).getLineageCount();
//        }
//        
//        // Probability of single coalescence event
//        logP += Math.log(1.0/popFunc.getPopSize(conv.getHeight2()));
//        
//        return logP;
    }
    
}
