package contactrees.operators;

import contactrees.Conversion;
import contactrees.util.Util;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public abstract class CFOperator extends ConversionCreationOperator {

    /**
     * Take conversions which connect to edge above srcNode at times greater than
     * destTime and attach them instead to the lineage above destNode.
     *
     * Assumes topology has not yet been altered.
     *
     * @param srcNode   source node for move
     * @param destNode  dest node for move
     * @param destTime  new time of attachment of edge above srcNode to edge
     *                  above destNode
     * @return log probability of the collapsed attachments.
     */
    protected double collapseConversions(Node srcNode, Node destNode, double destTime) {
        double logP = 0.0;

        boolean reverseRootMove = srcNode.getParent().isRoot();
        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getSibling(srcNode);
        double maxChildHeight = getMaxRootChildHeight();
        
        // Conversions which degenerate to point back to their source node are removed
        List<Conversion> toRemove = new LinkedList<>();

        // Collapse non-root conversions

        Node node = destNode;
        while (!node.isRoot() && node.getHeight() < srcNodeP.getHeight()) {

            double lowerBound = Math.max(destTime, node.getHeight());
            double upperBound = Math.min(node.getParent().getHeight(),
                    srcNodeP.getHeight());

            for (Conversion conv : acg.getConversions()) {
                if (conv.getHeight() > lowerBound && conv.getHeight() < upperBound) {

                    if (conv.getNode1() == srcNode)
                        conv.setNode1(node);

                    if (conv.getNode2() == srcNode)
                        conv.setNode2(node);


                    if (conv.getNode1() == conv.getNode2())
                    	// Conversion degenerated to a loop -> remove and adjust logP later on
                    	toRemove.add(conv);
                    else {
                        if ((conv.getNode1() == node) || (conv.getNode2() == node))
                        	// No loop, but an edge could have been moved -> adjust logP
                            logP += Math.log(0.5);                    	
                    }               	
                }
            }
            
            node = node.getParent();
        }
        // Remove degenerated conversions
        double L = 2.0 * (srcNode.getParent().getHeight() - destTime);
        double Nexp = L * conversionRateInput.get().getValue();
        if (Nexp == 0.) { 
        	assert toRemove.isEmpty();
        } else {
        	// We remove "toRemove.size()" edges
        	// In reverse move these are added according to prior:
	        logP += -Nexp + toRemove.size() * Math.log(Nexp); // Factorial cancels
        }
        
        // individual conversion states
        for (Conversion conv : toRemove) {
        	logP += Math.log(1.0/L) + getAffectedBlocksProb(conv);

            removeConversion(conv);
        }
        
        // Apply topology modifications.
        disconnectEdge(srcNode);
        connectEdge(srcNode, destNode, destTime);

        if (reverseRootMove && destTime < maxChildHeight) {
            acg.setRoot(srcNodeS);
        }

        return logP;
    }

    /**
     * Take length of new edge above srcNode that is greater than the
     * original height of srcNode.parent and shifts a random fraction of
     * conversion attachments to it from the lineage above destNode.
     *
     * In the case that destNode was the root, the conversions starting
     * above destNode are drawn from the prior.
     *
     * Assumes topology has not yet been altered.
     *
     * @param srcNode source node for the move
     * @param destNode dest node for the move
     * @param destTime new time drawn for srcNode.P.
     * @return log probability of new conversion configuration.
     */
    protected double expandConversions(Node srcNode, Node destNode, double destTime) {
        double logP = 0.0;

        double volatileHeight = acg.getRoot().getHeight();
        boolean forwardRootMove = destTime > volatileHeight;
        double parentHeight = srcNode.getParent().getHeight();
        Node sibling = getSibling(srcNode);

        Node node = srcNode.getParent();
        // TODO make more efficient by iterating over conversions per interval.
        // TODO stop at destTime
        while (!node.isRoot()) {
            for (Conversion conv : acg.getConversions()) {
            	if (conv.getHeight() < destTime) {
	                if (conv.getNode1() == node) {
	                    if (Randomizer.nextBoolean())  // TODO Would it make sense to parameterize this 50/50 choice?
	                        conv.setNode1(srcNode);
	                    logP += Math.log(0.5);
	                }
	
	                if (conv.getNode2() == node) {
	                    if (Randomizer.nextBoolean())
	                        conv.setNode2(srcNode);
	                    logP += Math.log(0.5);
	                }
            	}
            }

            node = node.getParent();
        }

        // Apply topology modifications.
        disconnectEdge(srcNode);
        connectEdge(srcNode, destNode, destTime);

        if (forwardRootMove) {
        	acg.setRoot(srcNode.getParent());
        }

        // Randomly add edges which would be pruned in a corresponding collapse move
        double L = 2.0*(destTime - parentHeight);
        double Nexp = L * conversionRateInput.get().getValue();
        int N = 0;
        if (Nexp > 0.) {
	        // Choose number of new conversions according to Poisson distribution
	        N = (int)Randomizer.nextPoisson(Nexp);
	        logP += -Nexp + N*Math.log(Nexp); // Factorial cancels
        }

        // Randomly place conversions between the new source and an ancestor of the old source. 
        for (int i=0; i<N; i++) {
            Conversion conv = addNewConversion();

            double convHeight = Randomizer.uniform(destTime, parentHeight);
            boolean convDirection = Randomizer.nextBoolean();
            logP += Math.log(1.0/L);

            conv.setHeight(convHeight);
            		
            Node other = getAncestorAtHeight(sibling, convHeight);
            if (convDirection) {
            	conv.setNode1(srcNode);
            	conv.setNode2(other);
            } else {
                conv.setNode1(other);
                conv.setNode2(srcNode);
            }

            logP += drawAffectedBlocks(conv);
        }
        
        return logP;
    }
	
	protected static Node getAncestorAtHeight(Node node, final double height) {
		assert height > node.getHeight();
		
		while (node.getParent().getHeight() < height) {
			node = node.getParent();
		}
		
		return node;
	}
		
}
