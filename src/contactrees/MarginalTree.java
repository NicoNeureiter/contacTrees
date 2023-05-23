package contactrees;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contactrees.CFEventList.Event;
import contactrees.model.likelihood.CTreeLikelihood;


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
    public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel", "", new StrictClockModel());
    public Input<ArrayList<String>> frozenTaxaInput = new Input<>(
            "frozenTaxa",
            "Taxa for which the last branch should have a fixed branch rate of 0.",
            new ArrayList<>());

    public ConversionGraph acg;
    public Block block;
    protected BranchRateModel.Base branchRateModel;
    protected boolean hasBranchRates;
    protected ArrayList<String> frozenTaxa;

    protected boolean outdated;
    protected boolean manuallyUpdated = false;
    protected String lastBlockState;
    protected boolean changed = true;
    protected boolean[] isFree;
    protected int nextNonLeafNr;
    ArrayList<MarginalNode> activeCFlineages = new ArrayList<>();

    protected boolean customdebug = false;

    public void setManuallyUpdated() {
        manuallyUpdated = true;
    }

    public void setBranchRateModel(BranchRateModel.Base brm) {
        branchRateModelInput.setValue(brm, this);
        branchRateModel = brm;
        hasBranchRates = (branchRateModel != null);
    }

    private CTreeLikelihood getLikelihood() {
        for (BEASTInterface beastObject : getOutputs())
            if (beastObject instanceof CTreeLikelihood)
                return (CTreeLikelihood) beastObject;
        throw new RuntimeException("No CTreeLikelihood found in the outputs of the marginal tree.");
    }

    @Override
    public void initAndValidate() {

        acg = networkInput.get();
        block = blockInput.get();
        frozenTaxa = frozenTaxaInput.get();

        branchRateModel = branchRateModelInput.get();
//        hasBranchRates = !(branchRateModel instanceof StrictClockModel);
        hasBranchRates = (branchRateModel != null);

        // Initialize to clonal frame of acg
        String beastID = ID;
        assignFrom(acg);
        setID(beastID);

        super.initAndValidate();
        activeCFlineages = new ArrayList<>();
        for (Node node : acg.getNodesAsArray())
            activeCFlineages.add(null);

        recalculate();
        makeOutdated();
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

    public void makeOutdated() {
        outdated = true;
        setEverythingDirty(true);
        manuallyUpdated = false;
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


    public int getNextFreeIndex() {
        // Skip all indices that have been used
        while (!isFree[nextNonLeafNr])
            nextNonLeafNr += 1;

        // Once we return the index, it will be marked as used up
        isFree[nextNonLeafNr] = false;

        // Return the index
        return nextNonLeafNr;
    }

    public void recalculate() {
        startEditing(null);
        if (customdebug) System.out.print("*");

        List<Event> cfEvents = acg.getCFEvents();
        activeCFlineages.replaceAll((oldNode) -> {
            return null;
        });
        List<Conversion> convs = getBlockConversions();
        convs.sort((c1, c2) -> {
            if (c1.height < c2.height)
                return -1;

            if (c1.height > c2.height)
                return 1;

            return 0;
        });

        int iConv = 0;
        int nLeafs = acg.getLeafNodeCount();

        nextNonLeafNr = nLeafs;

        for (int iEvent = 0; iEvent < cfEvents.size(); iEvent++) {
            Event event = cfEvents.get(iEvent);
            Node node = event.getNode();

            // Process the current CF-event
            switch (event.getType()) {
                case SAMPLE:
                    MarginalNode marginalLeaf = registerLeafNode(node);
                    activeCFlineages.set(node.getNr(), marginalLeaf);
                    break;

                case COALESCENCE:
                    Node left = node.getChild(0);
                    Node right = node.getChild(1);

                    if (activeCFlineages.get(left.getNr()) != null && activeCFlineages.get(right.getNr()) != null) {
                        // Create a new marginal node at the coalescence event
                        MarginalNode marginalNode = registerNode(node, activeCFlineages, nextNonLeafNr++);

                        // Remove the old and add the new marginal node to the active lineages.
                        activeCFlineages.set(left.getNr(), null);
                        activeCFlineages.set(right.getNr(), null);
                        activeCFlineages.set(node.getNr(), marginalNode);
                    } else {
                        // Only one side is active -> no coalescence in marginal tree (i.e. no marginal node)

                        if (activeCFlineages.get(left.getNr()) != null) {
                            MarginalNode marginalLeft = activeCFlineages.get(left.getNr());
                            updateTimeLength(node.getHeight(), left, marginalLeft);

                            activeCFlineages.set(left.getNr(), null);
                            activeCFlineages.set(node.getNr(), marginalLeft);
                            break;
                        }

                        if (activeCFlineages.get(right.getNr()) != null) {
                            MarginalNode marginalRight = activeCFlineages.get(right.getNr());
                            updateTimeLength(node.getHeight(), right, marginalRight);

                            activeCFlineages.set(right.getNr(), null);
                            activeCFlineages.set(node.getNr(), marginalRight);
                            break;
                        }
                    }
                    break;
            }

            // Process all conversion below the next CF-event
            while (iConv < convs.size() &&
                    (event.node.isRoot() || convs.get(iConv).height < cfEvents.get(iEvent + 1).getHeight())) {

                Conversion conv = convs.get(iConv++);
                Node node1 = conv.getNode1();
                Node node2 = conv.getNode2();

                if ((activeCFlineages.get(node1.getNr()) != null) && (activeCFlineages.get(node2.getNr()) != null)) {
                    // Both lineages at the conversion are active --> coalescence in the marginal tree

                    // Create a MarginalNode at the point of the conversion and add it as a new lineage
                    MarginalNode convNode = registerNode(conv, activeCFlineages, nextNonLeafNr++);

                    // Remove child lineages and add new one
                    activeCFlineages.set(node1.getNr(), null);
                    activeCFlineages.set(node2.getNr(), null);
                    activeCFlineages.set(node2.getNr(), convNode);
                } else {
                    // node1 or node2 already moved away (overshadowed by another conversion)

                    if (activeCFlineages.get(node1.getNr()) != null) {
                        // node1 passes conversion, but node2 branched away --> CF lineage of node2 is continued by node1
                        MarginalNode margNode1 = activeCFlineages.get(node1.getNr());
                        updateTimeLength(conv.getHeight(), node1, margNode1);

                        activeCFlineages.set(node1.getNr(), null);
                        activeCFlineages.set(node2.getNr(), margNode1);
                    }
                    // else: node1 already branched away --> conversion has no effect
                }
            }
        }

        // A single active CF lineage (the root) should remain:;
        MarginalNode newRoot = activeCFlineages.get(acg.getRoot().getNr());
        setRootOnly(newRoot);

        if (hasBranchRates) {
            rollOutTime(newRoot, newRoot.getHeight());
        }

        outdated = false;
//        block.updatedMarginalTree();
//        lastBlockState = block.toString();

    }

    public MarginalNode registerLeafNode(Node node) {
        MarginalNode marginalNode = (MarginalNode) m_nodes[node.getNr()];
//        if (marginalNode.getHeight() != node.getHeight()) marginalNode.makeDirty(Tree.IS_FILTHY);
//        else if (marginalNode.getLength() != node.getLength()) marginalNode.makeDirty(Tree.IS_FILTHY);
        if ((!marginalNode.equalsNode(node))) {
            marginalNode.setHeight(node.getHeight());
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }
        marginalNode.lastEventHeight = node.getHeight();
        marginalNode.setID(node.getID());
        marginalNode.timeLength = 0;

        return marginalNode;
    }

    public MarginalNode registerNode(Node node, ArrayList<MarginalNode> activeCFlineages, int nodeNr) {
        Node left = node.getChild(0);
        Node right = node.getChild(1);
        double height = node.getHeight();
        MarginalNode newLeft = activeCFlineages.get(left.getNr());
        MarginalNode newRight = activeCFlineages.get(right.getNr());

        // Take the old marginal node from the given index
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];
        Node oldLeft = marginalNode.getLeft();
        Node oldRight = marginalNode.getRight();

        // Mark the node as filthy if one of the children changed
        if (oldLeft.getNr() != newLeft.getNr())
            marginalNode.makeDirty(Tree.IS_FILTHY);

        if (oldRight.getNr() != newRight.getNr())
            marginalNode.makeDirty(Tree.IS_FILTHY);

        // Mark the node and children as dirty if the height changed
        if (marginalNode.getHeight() != node.getHeight()) {
            marginalNode.makeDirty(Tree.IS_FILTHY);
            List<Node> children = marginalNode.getChildren();
            if (children.size() > 0) children.get(0).makeDirty(Tree.IS_FILTHY);
            if (children.size() > 1) children.get(1).makeDirty(Tree.IS_FILTHY);
        }

        // Update the marginal node and meta information
        marginalNode.update(height, newLeft, newRight);
        marginalNode.timeLength = 0;
        marginalNode.lastEventHeight = height;

        // Update time-length of both children
        if (hasBranchRates) {
            updateTimeLength(height, left, newLeft);
            updateTimeLength(height, right, newRight);
        }

        return marginalNode;
    }

    public MarginalNode registerNode(Conversion conv, ArrayList<MarginalNode> activeCFlineages, int nodeNr) {
        Node left = conv.getNode1();
        Node right = conv.getNode2();
        double height = conv.getHeight();
        MarginalNode marginalLeft = activeCFlineages.get(left.getNr());
        MarginalNode marginalRight = activeCFlineages.get(right.getNr());

        // Take the old marginal node from the given index
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];

        // Mark the node as filthy if one of the children changed
        if (marginalNode.getLeft().getNr() != marginalLeft.getNr())
            marginalNode.makeDirty(Tree.IS_FILTHY);
        if (marginalNode.getRight().getNr() != marginalRight.getNr())
            marginalNode.makeDirty(Tree.IS_FILTHY);

        // Mark the node and children as dirty if the height changed
        if (marginalNode.getHeight() != conv.getHeight()) {
            marginalNode.makeDirty(Tree.IS_FILTHY);
            List<Node> children = marginalNode.getChildren();
            if (children.size() > 0) children.get(0).makeDirty(Tree.IS_FILTHY);
            if (children.size() > 1) children.get(1).makeDirty(Tree.IS_FILTHY);
        }

        marginalNode.update(height, marginalLeft, marginalRight);
        marginalNode.timeLength = 0;
        marginalNode.lastEventHeight = height;

        // Update time-length of both children
        if (hasBranchRates) {
            updateTimeLength(height, left, marginalLeft);
            updateTimeLength(height, right, marginalRight);
        }

//        marginalNode.makeDirty(Tree.IS_FILTHY);

        return marginalNode;
    }

    void updateTimeLength(double parentHeight, Node child, MarginalNode marginalChild) {
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
    void rollOutTime(MarginalNode root, double rootHeight) {
        root.setHeight(rootHeight);

        List<Node> children = root.getChildren();
        if (children.size() > 0) {
            MarginalNode leftChild = (MarginalNode) children.get(0);
            rollOutTime(leftChild, rootHeight - leftChild.timeLength);
        }
        if (children.size() > 1) {
            MarginalNode rightChild = (MarginalNode) children.get(1);
            rollOutTime(rightChild, rootHeight - rightChild.timeLength);
        }

////        for (Node childAsNode : root.getChildren()) {
//        for (int i=0; i< root.getChildren().size(); i++) {
//            MarginalNode child = (MarginalNode) root.getChildren().get(i);
//            rollOutTime(child, rootHeight - child.timeLength);
//        }
    }

    @Override
    protected void accept() {
        super.accept();
        setEverythingDirty(false);
    }

    @Override
    public void store() {}

    @Override
    public void restore() {
        postCache = null;

        manuallyUpdated = false;
        if (customdebug) System.out.print("r");
        recalculate();
//        setEverythingDirty(false);
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

    ArrayList<Conversion> _blockConvs = new ArrayList<>();

    /**
     * Obtain the list of conversions which affect this block, i.e. the ones defining this marginal tree.
     *
     * @return List of relevant conversions.
     */
    public List<Conversion> getBlockConversions() {
        _blockConvs.clear();
        ConversionList convList = acg.getConversions();
        List<Integer> cIDs = block.getConversionIDs();

        for (int i=0; i < cIDs.size(); i++) {
            Conversion c = convList.get(cIDs.get(i));
            assert c != null;
            _blockConvs.add(c);
        }

        return _blockConvs;
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


    /**
     * Marginal trees are currently implemented as lean calculation
     * nodes, hence the state can be inconsistent at intermediate stages.
     */
    @Override
    public int getChecksum() {
        return 0;
    }

}

