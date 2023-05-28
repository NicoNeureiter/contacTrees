package contactrees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import contactrees.CFEventList.Event;


/**
 * A calculation node extracting the marginal tree for a block from the ACG.
 * The resulting marginal trees are the basis for the likelihood computation
 * of the corresponding block.
 *
 * @author Nico Neureiter
 */

public class MarginalTreeSlow extends Tree {

    public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph (network) from which the marginal tree is sampled.",
            Input.Validate.REQUIRED);
    public Input<Block> blockInput = new Input<>(
            "block",
            "The block object cinitArraysontaining the moves this marginal tree follows along the conversion graph.",
            Input.Validate.REQUIRED);
    public Input<BranchRateModel.Base> branchRateModelInput = new Input<>(
            "branchRateModel",
            "A model describing the rates on the branches of the clonal frame tree.");

    public ConversionGraph acg;
    public Block block;
    protected BranchRateModel.Base branchRateModel;
    private boolean hasBranchRates;

    protected boolean outdated;
    protected boolean manuallyUpdated = false;
    String lastBlockState;
    protected boolean changed = true;

    boolean customdebug = false;

    public void setManuallyUpdated() {
        manuallyUpdated = true;
    }

    @Override
    public void initAndValidate() {

        acg = networkInput.get();
        block = blockInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
            hasBranchRates = true;
        } else {
            branchRateModel = new StrictClockModel();
            hasBranchRates = false;
        }

        // Initialize to clonal frame of acg
        String beastID = ID;
        assignFrom(acg);
        setID(beastID);

        super.initAndValidate();
        recalculate();

        outdated = true;
        manuallyUpdated = false;
    }

    @Override
    public boolean requiresRecalculation() {
        outdated = checkOutdated();
        changed = outdated;

        if (outdated) {
            recalculate();

            if (manuallyUpdated) {
                manuallyUpdated = false;
                return false;
            } else
                return true;


        } else {
            if (customdebug) {
                String oldNewick = root.toNewick();
                recalculate();
                String newNewick = root.toNewick();
                assert oldNewick == newNewick;
            }
            return false;
        }
    }

    @Override
    public boolean somethingIsDirty() {
        return this.changed;
    }

    public boolean checkOutdated() {
        if (acg.somethingIsDirty())
            return true;

        if (block.isAheadOfMarginalTree())
            return true;

        if (branchRateModel.isDirtyCalculation())
            return true;

        return outdated;
    }


    public void recalculate() {
        startEditing(null);
        if (customdebug) System.out.print("*");

        List<Event> cfEvents = acg.getCFEvents();
        Map<Node, MarginalNodeSlow> activeCFlineages = new HashMap<>();
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
                    MarginalNodeSlow marginalLeaf = registerLeafNode(node);
                    activeCFlineages.put(node, marginalLeaf);

                    assert activeCFlineages.size() == nActive + 1;
                    break;

                case COALESCENCE:
                    Node left = node.getChild(0);
                    Node right = node.getChild(1);

                    if (activeCFlineages.containsKey(left) && activeCFlineages.containsKey(right)) {
                        // Create a new marginal node at the coalescence event
                        MarginalNodeSlow marginalNodeSlow = registerNode(node, activeCFlineages, nextNonLeafNr++);

                        // Remove the old and add the new marginal node to the active lineages.
                        activeCFlineages.remove(left);
                        activeCFlineages.remove(right);
                        activeCFlineages.put(node, marginalNodeSlow);

                        assert activeCFlineages.size() == nActive - 1;

                    } else {
                        // Only one side is active -> no coalescence in marginal tree (i.e. no marginal node)

                        if (activeCFlineages.containsKey(left)) {
                            MarginalNodeSlow marginalLeft = activeCFlineages.get(left);
                            updateTimeLength(node.getHeight(), left, marginalLeft);

                            activeCFlineages.remove(left);
                            activeCFlineages.put(node, marginalLeft);
                            break;
                        }

                        if (activeCFlineages.containsKey(right)) {
                            MarginalNodeSlow marginalRight = activeCFlineages.get(right);
                            updateTimeLength(node.getHeight(), right, marginalRight);

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

                    // Create a MarginalNodeSlow at the point of the conversion and add it as a new lineage
                    MarginalNodeSlow convNode = registerNode(conv, activeCFlineages, nextNonLeafNr++);

                    // Remove child lineages and add new one
                    activeCFlineages.remove(node1);
                    activeCFlineages.remove(node2);
                    activeCFlineages.put(node2, convNode);

                    assert activeCFlineages.size() == nActive - 1;

                } else {
                    // node1 or node2 already moved away (overshadowed by another conversion)

                    if (activeCFlineages.containsKey(node1)) {
                        // node1 passes conversion, but node2 branched away --> CF lineage of node2 is continued by node1
                        MarginalNodeSlow margNode1 = activeCFlineages.get(node1);
                        updateTimeLength(conv.getHeight(), node1, margNode1);

                        activeCFlineages.remove(node1);
                        activeCFlineages.put(node2, margNode1);
                    }
                    // else: node1 already branched away --> conversion has no effect

                    assert activeCFlineages.size() == nActive;
                }
            }
        }

        // A single active CF lineage (the root) should remain:;
        MarginalNodeSlow newRoot = activeCFlineages.get(acg.getRoot());
        setRootOnly(newRoot);

        if (hasBranchRates) {
            rollOutTime(newRoot, newRoot.getHeight());
        }

        outdated = false;
//        block.updatedMarginalTree();
//        lastBlockState = block.toString();
    }

    public MarginalNodeSlow registerLeafNode(Node node) {
        MarginalNodeSlow marginalNode = (MarginalNodeSlow) m_nodes[node.getNr()];
        marginalNode.setHeight(node.getHeight());
        marginalNode.lastEventHeight = node.getHeight();
        marginalNode.setID(node.getID());
        marginalNode.timeLength = 0;

        marginalNode.makeDirty(Tree.IS_FILTHY);

        return marginalNode;
    }

    public MarginalNodeSlow registerNode(Node node, Map<Node, MarginalNodeSlow> activeCFlineages, int nodeNr) {
        Node left = node.getChild(0);
        Node right = node.getChild(1);
        double height = node.getHeight();
        MarginalNodeSlow marginalLeft = activeCFlineages.get(left);
        MarginalNodeSlow marginalRight = activeCFlineages.get(right);

        MarginalNodeSlow marginalNode = (MarginalNodeSlow) m_nodes[nodeNr];
        marginalNode.update(height, marginalLeft, marginalRight);
        marginalNode.timeLength = 0;
        marginalNode.lastEventHeight = height;
        marginalNode.makeDirty(Tree.IS_FILTHY);

        // Update time-length of both children
        if (hasBranchRates) {
            updateTimeLength(height, left, marginalLeft);
            updateTimeLength(height, right, marginalRight);
        }

        return marginalNode;
    }

    public MarginalNodeSlow registerNode(Conversion conv, Map<Node, MarginalNodeSlow> activeCFlineages, int nodeNr) {
        Node left = conv.getNode1();
        Node right = conv.getNode2();
        double height = conv.getHeight();
        MarginalNodeSlow marginalLeft = activeCFlineages.get(left);
        MarginalNodeSlow marginalRight = activeCFlineages.get(right);

        MarginalNodeSlow marginalNode = (MarginalNodeSlow) m_nodes[nodeNr];
        marginalNode.update(height, marginalLeft, marginalRight);
        marginalNode.timeLength = 0;
        marginalNode.lastEventHeight = height;
        marginalNode.makeDirty(Tree.IS_FILTHY);

        // Update time-length of both children
        if (hasBranchRates) {
            updateTimeLength(height, left, marginalLeft);
            updateTimeLength(height, right, marginalRight);
        }

        return marginalNode;
    }

    void updateTimeLength(double parentHeight, Node child, MarginalNodeSlow marginalChild) {
        marginalChild.timeLength += (parentHeight - marginalChild.lastEventHeight) * branchRateModel.getRateForBranch(child);
        marginalChild.lastEventHeight = parentHeight;
    }

    /**
     * Assign new heights to the nodes in the tree by
     * rolling out the branch time-lengths from the
     * root down (starting at rootHeight).
     *
     * @param root
     * @param rootHeight
     */
    void rollOutTime(MarginalNodeSlow root, double rootHeight) {
        root.setHeight(rootHeight);

        for (Node childAsNode : root.getChildren()) {
            MarginalNodeSlow child = (MarginalNodeSlow) childAsNode;
            rollOutTime(child, rootHeight - child.timeLength);
        }
    }

    @Override
    public void store() {}

    @Override
    public void restore() {
        postCache = null;

        manuallyUpdated = false;
        if (customdebug) System.out.print("r");
        recalculate();
    }

    protected void initArraysSlim() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new MarginalNodeSlow[nodeCount];
        listMarginalNodes((MarginalNodeSlow) root, m_nodes);
        postCache = null;
    }

    /**
     * convert tree to array representation *
     */
    void listMarginalNodes(final MarginalNodeSlow node, final Node[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);

        for (final Node child : node.getChildren()) {
            listMarginalNodes((MarginalNodeSlow) child, nodes);
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

        final MarginalNodeSlow[] nodes = new MarginalNodeSlow[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = newNode();
        }

        setID(tree.getID());

        root = (MarginalNodeSlow) nodes[tree.getRoot().getNr()];
        root.assignFrom(nodes, tree.getRoot());
        root.setParent(null);

        nodeCount = tree.getNodeCount();
        internalNodeCount = tree.getInternalNodeCount();
        leafNodeCount = tree.getLeafNodeCount();

        initArraysSlim();
    }

    @Override
    protected MarginalNodeSlow newNode() {
        return new MarginalNodeSlow();
    }

    /**
     * The marginal tree in Newick format.
     */
    @Override
    public String toString() {
        return root.toString();
    }


    /**
     * Marginal trees are currently implemented as lean calculation
     * nodes, hence the state can be inconsistent at intermediate stages.
     */
    @Override
    public int getChecksum() {
        return 0;
    }

}

