package contactrees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hamcrest.core.IsInstanceOf;

import contactrees.MarginalNode;
import contactrees.CFEventList;
import contactrees.CFEventList.Event;
import contactrees.Conversion;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;


/**
 * A calculation node extracting the marginal tree for a block from the ACG.
 * The resulting marginal trees are the basis for the likelihood computation 
 * of the corresponding block.  
 * 
 * @author Nico Neureiter
 */

public class MarginalTree extends Tree {

	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The conversion graph (network) from which the marginal tree is sampled.",
			Input.Validate.REQUIRED);
	public Input<Block> blockInput = new Input<>(
			"block",
			"The block object cinitArraysontaining the moves this marginal tree follows along the conversion graph.",
			Input.Validate.REQUIRED);

    public Input<RealParameter> ratesCFInput = new Input<>(
            "ratesCF",
            "The rates associated with nodes in the clonal frame.");
    public Input<RealParameter> ratesMarginalTreeInput = new Input<>(
            "ratesMarginalTree",
            "The rates associated with nodes in the marginal tree, computed from ratesCF.");
	
	public ConversionGraph acg;
	public Block block;
	boolean outdated;
	
	public void initAndValidate() {
	    
		acg = networkInput.get();
		block = blockInput.get();
		
		// Initialize to clonal frame of acg
		assignFrom(acg);
		recalculate();
        super.initAndValidate();
		    
		outdated = true;
		
//		nodeTypeInput.set(MarginalNode.class.getName());
		assert nodeTypeInput.get().equals(MarginalNode.class.getName());
	}
	
	@Override
	public boolean requiresRecalculation() {
	    
		outdated = true;
		
		// Recalculate right away if needed
		if (outdated) {
			recalculate();
			outdated = false;
			return true;
		} else {
            return false;
		}
	}
	
	public void recalculate() {
	    startEditing(null);
	    
	    Node startRoot = acg.getRoot();
	    
        List<Event> cfEvents = acg.getCFEvents();
        Map<Node, MarginalNode> activeCFlineages = new HashMap<>();
        ArrayList<Conversion> convs = getBlockConversions();
        for (Conversion c : convs) {
            assert c != null;
        }
        convs.sort((c1, c2) -> {
            if (c1.height < c2.height)
                return -1;

            if (c1.height > c2.height)
                return 1;

            return 0;
        });
        
        int iConv = 0;
        int nextNonLeafNr = acg.getLeafNodeCount();
        for (int iEvent = 0; iEvent < cfEvents.size(); iEvent++) {
            Event event = cfEvents.get(iEvent);
            Node node = event.getNode();
            
            int nActive = activeCFlineages.size();
            
            // Process the current CF-event
            switch (event.getType()) {
                case SAMPLE:
                    MarginalNode marginalLeaf = new MarginalNode(acg, node.getNr(), event.getHeight());
                    marginalLeaf.setID(node.getID());
                    marginalLeaf.cfNodeNr = node.getNr();
                    activeCFlineages.put(node, marginalLeaf);
                    
                    assert activeCFlineages.size() == nActive + 1;
                    
                    break;

                case COALESCENCE:
                    Node left = node.getLeft();
                    Node right = node.getRight();
                    
                    if (activeCFlineages.containsKey(left) && activeCFlineages.containsKey(right)) {
                        MarginalNode marginalLeft = activeCFlineages.get(left);
                        MarginalNode marginalRight = activeCFlineages.get(right);
                        
                        // Create a new marginal node at the coalescence event
                        MarginalNode marginalNode = new MarginalNode(acg, nextNonLeafNr++, event.getHeight(), 
                                                                     marginalLeft, marginalRight);
                        marginalNode.cfNodeNr = node.getNr();

                        // Remove the old and add the new marginal node to the active lineages.
                        activeCFlineages.remove(left);
                        activeCFlineages.remove(right);
                        activeCFlineages.put(node, marginalNode);
                        
                        assert activeCFlineages.size() == nActive - 1;
                        
                    } else {
                        // Only one side is active -> no coalescence in marginal tree (i.e. no marginal node)
                            
                        if (activeCFlineages.containsKey(left)) {
                            MarginalNode marginalLeft = activeCFlineages.get(left);
                            activeCFlineages.remove(left);
                            activeCFlineages.put(node, marginalLeft);
                            break;
                        }

                        if (activeCFlineages.containsKey(right)) {
                            MarginalNode marginalRight = activeCFlineages.get(right);
                            activeCFlineages.remove(right);
                            activeCFlineages.put(node, marginalRight);
                            break;
                        }
                        
                        assert activeCFlineages.size() == nActive;
                        
                    }
                    break;
            }

            // Process all conversion below the next CF-event 
            while (iConv < convs.size() &&
                    (event.node.isRoot() || convs.get(iConv).height < cfEvents.get(iEvent + 1).getHeight())) {
                
                nActive = activeCFlineages.size();
                
                Conversion conv = convs.get(iConv++);
                Node node1 = conv.getNode1();
                Node node2 = conv.getNode2();
                
                if (activeCFlineages.containsKey(node1) && activeCFlineages.containsKey(node2)) {
                    // Both lineages at the conversion are active --> coalescence in the marginal tree
                    
                    // Pop active lineages of node1 and node2
                    MarginalNode left = activeCFlineages.get(node2);
                    MarginalNode right = activeCFlineages.get(node1);
                    activeCFlineages.remove(node1);
                    activeCFlineages.remove(node2);
                    
                    // Create a MarginalNode at the point of the conversion and add it as a new lineage
                    MarginalNode convNode = new MarginalNode(acg, nextNonLeafNr++, conv.height, left, right);
                    convNode.convID = conv.getID();
                    activeCFlineages.put(node2, convNode);
                
                    assert activeCFlineages.size() == nActive - 1;
                    
                } else {
                    // node1 or node2 already moved away (overshadowed by another conversion)
                    
                    if (activeCFlineages.containsKey(node1)) {
                        // node1 passes conversion, but node2 branched away --> CF lineage of node2 is continued by node1
                        MarginalNode margNode1 = activeCFlineages.get(node1);
                        activeCFlineages.remove(node1);
                        activeCFlineages.put(node2, margNode1);
                    }    
                    // else: node1 already branched away --> conversion has no effect
                    
                    assert activeCFlineages.size() == nActive;
                }
            }
        }

        // A single active CF lineage (the root) should remain:
        assert activeCFlineages.size() == 1;
        assert activeCFlineages.containsKey(acg.getRoot());
        assert m_nodes.length == acg.getNodeCount();
        root = activeCFlineages.get(acg.getRoot());
        setRoot(root);
        initArrays();
        
        assert root.isRoot();
        assert root == getRoot();
        
        int rootCount = 0;
        for (Node node : getNodesAsArray())
            if (node.isRoot())
                rootCount += 1;
        assert rootCount == 1;
        
        for (Node node : getNodesAsArray()) {
            if (node.isRoot())
                assert node == getRoot();
        }

//        acg.setEverythingDirty(true); // TODO Get around this!
//        setSomethingIsDirty(true);
        
        
    }
	
	@Override
    public boolean somethingIsDirty() {
	    // Whould be dirty if either ACG or BlockSet is dirty
        boolean ALWAYS_DIRTY = true;
        
        if (ALWAYS_DIRTY) 
            return true;
        
        return (acg.somethingIsDirty() || block.somethingIsDirty());
	}
	
	@Override
    public void store() {
//	    super.store();
    }
	
	@Override
    public void restore() {
//        for (Node node : getNodesAsArray())
//            assert node.getClass() == MarginalNode.class;
//
//        super.restore();
//
//        for (Node node : getNodesAsArray())
//            assert node.getClass() == MarginalNode.class;
        
        outdated = true;
	}
	
    /**
     * Obtain the list of conversions which affect this block, i.e. the ones defining this marginal tree.
     * 
     * @return List of relevant conversions. 
     */
    public ArrayList<Conversion> getBlockConversions() {
        ArrayList<Conversion> blockConvs = new ArrayList<>();
        ConversionList convList = acg.getConversions();
        
        for (int cID : block.getConversionIDs()) {
            Conversion c = convList.get(cID);
            assert c != null;
            blockConvs.add(c);
        }
        
        return blockConvs;
    }

    /**
     * The marginal tree in Newick format.
     */
    @Override
    public String toString() {
        return root.toString();
    }	
	
}

