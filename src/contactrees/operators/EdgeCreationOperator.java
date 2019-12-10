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

import contactrees.Block;
import contactrees.BlockSet;
import contactrees.CFEventList;
import contactrees.CFEventList.Event;
import contactrees.Conversion;
import contactrees.util.Util;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.HashSet;
import java.util.List;
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
//    	blockSet.addConversion(conv);
    	return conv;
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
    	CFEventList cfEventList = acg.getCFEventList();
    	List<Event> cfEvents = cfEventList.getCFEvents();
        double logP = 0.0;

        // Choose event interval
    	double[] intervalVolumes = cfEventList.getIntervalVolumes();
    	int iEvent = Util.sampleCategorical(intervalVolumes);
        Event event = cfEvents.get(iEvent);
        
    	// Choose height within interval
        double height = Randomizer.uniform(event.getHeight(), cfEvents.get(iEvent+1).getHeight());
    	conv.setHeight(height);
    	
    	// Choose source lineage (given the height)
    	Set<Node> activeLineages = acg.getLineagesAtHeight(height);
    	Node node1 = Util.sampleFrom(activeLineages);
    	conv.setNode1(node1);
    	assert node1.getHeight() < height;
        
        // Choose destination lineage (given the height and node1)
        activeLineages.remove(node1);
        Node node2 = Util.sampleFrom(activeLineages);
        conv.setNode2(node2);
        
        assert conv.getNode1() != null;
        assert conv.getNode2() != null;
        assert conv.getNode1() != conv.getNode2();
        assert !conv.getNode1().isRoot();
        assert !conv.getNode2().isRoot();
        assert conv.isValid();
        
        return Math.log(1.0/acg.getClonalFramePairedLength());
        
//        
//        // Select departure point
//        double u = Randomizer.nextDouble()*acg.getClonalFramePairedLength();
//        logP += Math.log(1.0/acg.getClonalFramePairedLength());
//        
//        List<Event> eventList = acg.getCFEvents(); 
//        for (int i=0; i<eventList.size()-1; i++) {
//        	Event start = eventList.get(i);
//        	Event end = eventList.get(i+1);
//        	int k = start.getLineageCount();
//        	
//        	double intervalLength = end.getHeight() - start.getHeight(); 
//        	double pairedLength = intervalLength * k * (k-1); 
//        	
//        	if (u < pairedLength) {
//        		System.out.println(start.getNode().getNr());
//        		System.out.println(start.getNode());
//        		System.out.println(intervalLength);
//        		System.out.println(start.getNode().getLength());
//        		System.out.println(u + " < "+ pairedLength);
//        		// Pick height uniformly at random
//        		double height = Randomizer.uniform(0, intervalLength);
//        		conv.setHeight(height);
//        		
//        		// Pick pair of lineages uniformly from active interval
//        		HashSet<Node> lineages = acg.getLineagesAtHeight(height);
//        		Node node1 = Util.sampleFrom(lineages);
//        		conv.setNode1(node1);
//        		assert node1.getParent().getHeight() > height;
//        		
//        		lineages.remove(node1);
//        		Node node2 = Util.sampleFrom(lineages);
//        		conv.setNode2(node2);
//        		assert node2.getParent().getHeight() > height : "Parent too young: "
//        				+ node2.getParent().getHeight() + " < " + height;
//        		
//        		
//        		break;
//        	}
//        	u -= pairedLength;
//        }
//        
//        
//        return logP;
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
        
        logP += Math.log(1.0/acg.getClonalFramePairedLength());
//        logP += getEdgeCoalescenceProb(conv);
        
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
     * Get probability density for the arrival time of the given recombinant
     * edge under ClonalOrigin's coalescent model.
     * 
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeCoalescenceProb(Conversion conv) {
    	int nLineages = acg.countLineagesAtHeight(conv.getHeight());
    	return - Math.log(nLineages - 1);
    }
    
    /**
     * Include all block StateNodes in the list of affected state nodes.
     * This does not work automatically, since blocks are only indirect 
     * inputs through the blockSet.
     */
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = super.listStateNodes();
        for (Block block : blockSet.getBlocks()) {
            list.add(block);
        }
        return list;
    }
    
}
