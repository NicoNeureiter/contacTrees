package contactrees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import contactrees.MarginalNode;
import contactrees.CFEventList.Event;
import contactrees.Conversion;
import beast.core.Input;
import beast.core.StateNode;
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
    
    boolean nextTimeDirty;
    
//    Map<String, MarginalNode> oldNodes, newNodes;
    
    public void initAndValidate() {
        
        acg = networkInput.get();
        block = blockInput.get();
        
//        oldNodes = new HashMap<>();
//        newNodes = new HashMap<>();
        
        // Initialize to clonal frame of acg
        String beastID = ID;
        assignFrom(acg);
        setID(beastID);
        
        recalculate();
        super.initAndValidate();
            
        outdated = true;
        
        nextTimeDirty = true;
        
        // assert nodeTypeInput.get().equals(MarginalNode.class.getName()); 
        //   \--> Not needed, since we override newNode()
    }
    
    @Override
    public boolean requiresRecalculation() {
//        recalculate();
//        return true;

////      System.out.println(outdated);
////        outdated = true;
        outdated = outdated || checkOutdated();
        
        // Recalculate right away if needed
        if (outdated) {
            recalculate();
            outdated = false;
            return true;
        } else {
            System.out.println("Not outdated");
            return false;
        }
    }
    
    public boolean checkOutdated() {
        return (acg.somethingIsDirty() || block.somethingIsDirty());
//      // Check whether a clonal frame edge is dirty
//      if (!outdated) {
//          for (Node cfNode : acg.getNodesAsArray()) {
//              if (cfNode.isDirty() > Tree.IS_CLEAN) {
//                  outdated = true;
//                  break;
//              }
//          }
//      }
//      
//      // Check whether a conversion is dirty
//      if (!outdated) {
//          // Check whether recalculation is necessary
//          for (Node mNodeAsNode : getNodesAsArray()) {
//              MarginalNode mNode = (MarginalNode) mNodeAsNode;
//              if (mNode.isConversionNode()) {
//                  try {
//                      Conversion conv = mNode.getConversion();
//                        if (conv.isDirty() > Tree.IS_CLEAN) {
////                            System.out.println("A conversion is dirty");
//                            outdated = true;
//                            break;
//                        }   
//                  } catch(RuntimeException e) {
//                      outdated = true;
//                      break;
//                  }
//              } else {
//                  Node cfNode = mNode.getCFNode();
//                  if (cfNode.isDirty() >= Tree.IS_FILTHY) {
////                        System.out.println("ACG is dirty");
//                      outdated = true;
//                      break;
//                  }   
//              }
//          }
//      }
//      
//      // Check whether a block is dirty
//      if (!outdated) {
//          if (block.somethingIsDirty()) {
////                System.out.println("Block is dirty");
//              outdated = true;
//          }
//        }
    }
    
    public void recalculate() {
        startEditing(null);
        
//      System.out.println("                                 [[[ " + block.index + " ]]]");

//        newNodes = new HashMap<>();
//        oldNodes = new HashMap<>();
        
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
        int nLeafs = acg.getLeafNodeCount();
        int nextNonLeafNr = nLeafs;
        
        for (int iEvent = 0; iEvent < cfEvents.size(); iEvent++) {
            Event event = cfEvents.get(iEvent);
            Node node = event.getNode();
            
            int nActive = activeCFlineages.size();
            
            // Process the current CF-event
            switch (event.getType()) {
                case SAMPLE:
                    MarginalNode marginalLeaf = registerLeafNode(node);
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
                        MarginalNode marginalNode = registerNode(node, marginalLeft, marginalRight, nextNonLeafNr++); 

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
                    MarginalNode convNode = registerNode(conv, left, right, nextNonLeafNr++);
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

        // A single active CF lineage (the root) should remain:;
        root = activeCFlineages.get(acg.getRoot());
        setRootOnly(root);

//        assert activeCFlineages.size() == 1;
//        assert activeCFlineages.containsKey(acg.getRoot());
//        assert m_nodes.length == acg.getNodeCount()
//        assert m_nodes[getNodeCount() - 1] == root;
//        assert root.isRoot();
//        assert root == getRoot();
//        int rootCount = 0;
//        for (Node node : getNodesAsArray()) {
//            assert node != null;
//            if (node.isRoot())
//                rootCount += 1;
//        }
//        assert rootCount == 1;

//        if (nextTimeDirty) {
//            setEverythingDirty(true);
//            nextTimeDirty = false;
////            System.out.println("This time dirty!");
//        }
        
//        oldNodes = newNodes;        
    }
    
    public MarginalNode registerLeafNode(Node node) {
        MarginalNode marginalNode = (MarginalNode) m_nodes[node.getNr()];
        marginalNode.setHeight(node.getHeight());
        marginalNode.setID(node.getID());
        marginalNode.cfNodeNr = node.getNr();
        marginalNode.convID = -1;
        marginalNode.makeDirty(Tree.IS_FILTHY);
        
        return marginalNode;
    }
    
    public MarginalNode registerNode(Node node, MarginalNode left, MarginalNode right, int nodeNr) {
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];
        marginalNode.update(node.getHeight(), left, right);
//        marginalNode.setID(id);
        marginalNode.cfNodeNr = node.getNr();
        marginalNode.convID = -1;
        marginalNode.makeDirty(Tree.IS_FILTHY);
        
//        newNodes.put(id, marginalNode);
        return marginalNode;
    }
    
    public MarginalNode registerNode(Conversion conv, MarginalNode left, MarginalNode right, int nodeNr) {
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];
        marginalNode.update(conv.getHeight(), left, right);
//        marginalNode.setID(id);
        marginalNode.cfNodeNr = -1;
        marginalNode.convID = conv.getID();
        marginalNode.makeDirty(Tree.IS_FILTHY);
        
//        newNodes.put(id, marginalNode);
        return marginalNode;
    }

    @Override
    public boolean somethingIsDirty() {
        return true;
//        // Whould be dirty if either ACG or BlockSet is dirty
//        return (acg.somethingIsDirty() || block.somethingIsDirty());
    }
    
    @Override
    public void store() {}
    
    @Override
    public void restore() {
//        for (Node node : getNodesAsArray())
//            assert node.getClass() == MarginalNode.class;
//
//        super.restore();
//
//        for (Node node : getNodesAsArray())
//            assert node.getClass() == MarginalNode.class;
        
//      System.out.println("################################### RESTORE ###################################");
//        super.restore();
        postCache = null;
//        oldNodes = new HashMap<>();
//        newNodes = new HashMap<>();
        block.setSomethingIsDirty(true);
        acg.setSomethingIsDirty(true);
        nextTimeDirty = true;
        outdated = true;
    }
    

    protected void initArraysSlim() { 
        // initialise tree-as-array representation + its stored variant
        m_nodes = new MarginalNode[nodeCount];
        listMarginalNodes((MarginalNode) root, m_nodes);
        postCache = null;
    }
    
    /**
     * convert tree to array representation *
     */
    void listMarginalNodes(final MarginalNode node, final Node[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);
        
        for (final Node child : node.getChildren()) {
            listMarginalNodes((MarginalNode) child, nodes);
        }
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
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Tree tree = (Tree) other;
        
        final MarginalNode[] nodes = new MarginalNode[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = newNode();
        }
        
        setID(tree.getID());

        root = (MarginalNode) nodes[tree.getRoot().getNr()];
        root.assignFrom(nodes, tree.getRoot());
        root.setParent(null);
        
        nodeCount = tree.getNodeCount();
        internalNodeCount = tree.getInternalNodeCount();
        leafNodeCount = tree.getLeafNodeCount();
        
        initArraysSlim();
    }
    
    @Override
    protected MarginalNode newNode() {
        return new MarginalNode();
    }

    /**
     * The marginal tree in Newick format.
     */
    @Override
    public String toString() {
        return root.toString();
    }   
    
}

