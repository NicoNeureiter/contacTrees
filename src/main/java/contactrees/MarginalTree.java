package contactrees;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import contactrees.CFEventList.Event;
import contactrees.CFEventList.EventType;


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
    public Input<BranchRateModel> branchRateModelInput = new Input<>(
            "branchRateModel",
            "The branch rate model for rolling out the branch lengths of the marginal tree (if not given, a strict clock is assumed).",
            new StrictClockModel());
    public Input<ArrayList<String>> frozenTaxaInput = new Input<>(
            "frozenTaxa",
            "Taxa for which the last branch should have a fixed branch rate of 0.",
            new ArrayList<>());
    public Input<Integer> validationLevelInput = new Input<>("validationLevel",
            "0 = off; 1 = check this marginal tree against a fresh rebuild after every proposal.", 0);

    public ConversionGraph acg;
    public Block block;
    protected BranchRateModel branchRateModel;
    protected boolean hasBranchRates;
    protected ArrayList<String> frozenTaxa;
    protected int validationLevel;

    protected boolean outdated;
    protected long seenCfVersion = -1;
    protected long seenConvVersion = -1;
    protected long seenBlockVersion = -1;
    protected boolean manuallyUpdated = false;
    protected String lastBlockState;
    protected boolean changed = true;
    protected boolean[] isFree;
    protected int nextNonLeafNr;
    ArrayList<MarginalNode> activeCFlineages = new ArrayList<>();

    public void setManuallyUpdated() {
        manuallyUpdated = true;
    }

    boolean isInBeauti() {
        return ProgramStatus.name.equals("BEAUti");
    }

    // Diagnostics for single-threaded benchmark runs; not synchronized.
    public static long nRecalculate = 0;
    public static long nFilthyNodesMarked = 0;
    public static long nRestore = 0;
    public static long nOutdatedChecks = 0;
    public static long nOutdatedTrue = 0;

    public static void resetCounters() {
        nRecalculate = 0;
        nFilthyNodesMarked = 0;
        nRestore = 0;
        nOutdatedChecks = 0;
        nOutdatedTrue = 0;
    }

    public static String counterSummary() {
        return "MarginalTree counters: recalculate=" + nRecalculate
                + ", filthyNodesMarked=" + nFilthyNodesMarked
                + ", restore=" + nRestore
                + ", outdatedChecks=" + nOutdatedChecks
                + ", outdatedTrue=" + nOutdatedTrue;
    }

    public void setBranchRateModel(Base brm) {
        branchRateModelInput.setValue(brm, this);
        branchRateModel = brm;
        hasBranchRates = (branchRateModel != null);
    }

    @Override
    public void initAndValidate() {

        acg = networkInput.get();
        if (!isInBeauti())
            acg.registerMarginalTree(this);
        block = blockInput.get();
        frozenTaxa = frozenTaxaInput.get();
        validationLevel = validationLevelInput.get();

        branchRateModel = branchRateModelInput.get();
//        hasBranchRates = !(branchRateModel instanceof StrictClockModel);
        hasBranchRates = (branchRateModel != null);

        // Initialize to clonal frame of acg
        String beastID = ID;
        assignFrom(acg);
        setID(beastID);

        super.initAndValidate();
        activeCFlineages = new ArrayList<>();
        for (int i=0; i<acg.getNodeCount(); i++)
            activeCFlineages.add(null);

        recalculate();
        makeOutdated();
    }

    @Override
    public boolean requiresRecalculation() {
        outdated = checkOutdated();
        // OR in, don't overwrite: a direct recalculate() call earlier this proposal (the
        // Gibbs-sampling code in BorrowingOperator calls it directly, bypassing this method)
        // may already have set changed=true to signal a real change to TreeLikelihood, which
        // reads only somethingIsDirty()/changed, never re-derives dirtiness itself. Overwriting
        // that with a fresh (possibly false) checkOutdated() result here would hide the change.
        changed = changed || outdated;

        if (outdated) {
            recalculate();

            if (manuallyUpdated) {
                manuallyUpdated = false;
                return false;
            } else
                return true;


        } else {
            return false;
        }
    }

    public boolean isValidating() {
        return validationLevel > 0;
    }

    /**
     * Verify that this marginal tree matches a rebuild from the current ACG and block state.
     * A mismatch means the tree was stale, i.e. the dirty-tracking logic failed to invalidate it.
     *
     * Calls {@link #recalculate()}, which mutates the tree into the correct state as a side
     * effect; that self-healing is intentional, not a bug.
     *
     * @param phase description of the calling context, used in the failure message
     */
    public void assertConsistent(String phase) {
        int n = getNodeCount();
        String oldNewick = getRoot().toSortedNewick(new int[1], false);
        double[] oldHeights = new double[n];
        double[] oldTimeLengths = new double[n];
        for (int i = 0; i < n; i++) {
            Node node = getNode(i);
            oldHeights[i] = node.getHeight();
            oldTimeLengths[i] = ((MarginalNode) node).timeLength;
        }

        recalculate();

        String newNewick = getRoot().toSortedNewick(new int[1], false);
        if (!oldNewick.equals(newNewick)) {
            throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): stale tree detected"
                    + " for marginal tree " + getID() + " (block " + block.getID() + "):"
                    + " sorted Newick mismatch."
                    + " before=" + truncate(oldNewick) + " after=" + truncate(newNewick));
        }

        for (int i = 0; i < n; i++) {
            Node node = getNode(i);
            double newHeight = node.getHeight();
            if (Math.abs(oldHeights[i] - newHeight) > 1e-12) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): stale tree detected"
                        + " for marginal tree " + getID() + " (block " + block.getID() + "):"
                        + " height mismatch at node " + i
                        + " before=" + oldHeights[i] + " after=" + newHeight);
            }

            double newTimeLength = ((MarginalNode) node).timeLength;
            if (Math.abs(oldTimeLengths[i] - newTimeLength) > 1e-12) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): stale tree detected"
                        + " for marginal tree " + getID() + " (block " + block.getID() + "):"
                        + " timeLength mismatch at node " + i
                        + " before=" + oldTimeLengths[i] + " after=" + newTimeLength);
            }
        }
    }

    private static String truncate(String s) {
        return s.length() > 500 ? s.substring(0, 500) + "..." : s;
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
        nOutdatedChecks++;

        boolean result;
        if (outdated) {
            result = true;
        } else if (branchRateModel instanceof CalculationNode brm && brm.somethingIsDirty()) {
            result = true;
        } else if (acg.getCFVersion() != seenCfVersion) {
            result = true;
        } else if (block.getMoveSetVersion() != seenBlockVersion) {
            result = true;
        } else if (acg.getConvVersion() != seenConvVersion) {
            // Some conversion changed somewhere in the ACG; check precisely whether it was
            // one of THIS block's own conversions before concluding we're outdated.
            result = false;
            for (Conversion conv : getBlockConversions()) {
                if (conv.getVersion() > seenConvVersion) {
                    result = true;
                    break;
                }
            }
            if (!result) {
                // Nothing relevant changed -- re-baseline so we don't re-scan on every
                // subsequent check until convVersion moves again.
                seenConvVersion = acg.getConvVersion();
            }
        } else {
            result = false;
        }

        if (result)
            nOutdatedTrue++;

        return result;
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
        nRecalculate++;
        startEditing(null);

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
                    (event.node.isRoot() || isBelowNextCFEvent(convs.get(iConv), cfEvents.get(iEvent + 1)))) {

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
        seenCfVersion = acg.getCFVersion();
        seenConvVersion = acg.getConvVersion();
        seenBlockVersion = block.getMoveSetVersion();

        // Some operators (the Gibbs-sampling borrowing code in BorrowingOperator) call
        // recalculate() directly, bypassing requiresRecalculation(). somethingIsDirty() is
        // the ONLY way TreeLikelihood learns this tree changed (its own requiresRecalculation()
        // reads treeInput.get().somethingIsDirty(), never re-derives it), so a direct call must
        // still flip it -- otherwise the change is invisible to TreeLikelihood and it silently
        // keeps a stale cached logP. Always safe: the next requiresRecalculation() call
        // overwrites this with a fresh checkOutdated() result regardless.
        changed = true;
//        block.updatedMarginalTree();
//        lastBlockState = block.toString();

    }

    /**
     * Decide whether a conversion should still be processed before the next CF event.
     *
     * A conversion exactly at the height of the next CF event only belongs below that event
     * if the event is the coalescence of one of the conversion's own lineages.
     */
    private boolean isBelowNextCFEvent(Conversion conv, Event nextEvent) {
        double convHeight = conv.getHeight();
        double nextHeight = nextEvent.getHeight();

        if (convHeight < nextHeight)
            return true;

        if (convHeight > nextHeight)
            return false;

        if (nextEvent.getType() != EventType.COALESCENCE)
            return false;

        Node nextNode = nextEvent.getNode();
        return conv.getNode1().getParent() == nextNode || conv.getNode2().getParent() == nextNode;
    }

    public MarginalNode registerLeafNode(Node node) {
        MarginalNode marginalNode = (MarginalNode) m_nodes[node.getNr()];
//        if (marginalNode.getHeight() != node.getHeight()) marginalNode.makeDirty(Tree.IS_FILTHY);
//        else if (marginalNode.getLength() != node.getLength()) marginalNode.makeDirty(Tree.IS_FILTHY);
        if ((!marginalNode.equalsNode(node))) {
            marginalNode.setHeight(node.getHeight());
            nFilthyNodesMarked++;
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

        if (nodeNr >= m_nodes.length) {
            throw binaryMergeCapacityError(nodeNr, newLeft, newRight,
                "coalescence at node " + node.getNr() + " height=" + height);
        }

        // Take the old marginal node from the given index
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];
        Node oldLeft = marginalNode.getLeft();
        Node oldRight = marginalNode.getRight();

        // Mark the node as filthy if one of the children changed
        if (oldLeft.getNr() != newLeft.getNr()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }

        if (oldRight.getNr() != newRight.getNr()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }

        // Mark the node and children as dirty if the height changed
        if (marginalNode.getHeight() != node.getHeight()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
            List<Node> children = marginalNode.getChildren();
            if (children.size() > 0) {
                nFilthyNodesMarked++;
                children.get(0).makeDirty(Tree.IS_FILTHY);
            }
            if (children.size() > 1) {
                nFilthyNodesMarked++;
                children.get(1).makeDirty(Tree.IS_FILTHY);
            }
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

        if (nodeNr >= m_nodes.length) {
            throw binaryMergeCapacityError(nodeNr, marginalLeft, marginalRight,
                "conversion " + conv.getID() + " at height=" + height
                    + " between nodes " + left.getNr() + " and " + right.getNr());
        }

        // Take the old marginal node from the given index
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];

        // Mark the node as filthy if one of the children changed
        if (marginalNode.getLeft().getNr() != marginalLeft.getNr()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }
        if (marginalNode.getRight().getNr() != marginalRight.getNr()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }

        // Mark the node and children as dirty if the height changed
        if (marginalNode.getHeight() != conv.getHeight()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
            List<Node> children = marginalNode.getChildren();
            if (children.size() > 0) {
                nFilthyNodesMarked++;
                children.get(0).makeDirty(Tree.IS_FILTHY);
            }
            if (children.size() > 1) {
                nFilthyNodesMarked++;
                children.get(1).makeDirty(Tree.IS_FILTHY);
            }
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

    /**
     * Build a descriptive error for a merge that would write past the end of the node array,
     * instead of letting it surface as an ArrayIndexOutOfBoundsException. Callers check the
     * capacity themselves so this (and the context string) is only ever built on the failure
     * path, not on every registerNode() call.
     */
    private IllegalStateException binaryMergeCapacityError(int nodeNr, MarginalNode left, MarginalNode right, String context) {
        String leftId = left == null ? "null" : Integer.toString(left.getNr());
        String rightId = right == null ? "null" : Integer.toString(right.getNr());
        return new IllegalStateException("MarginalTree exceeded node capacity while registering " + context
                + "; nodeNr=" + nodeNr + ", capacity=" + m_nodes.length
                + ", leftMarginal=" + leftId + ", rightMarginal=" + rightId
                + ", sameLineage=" + (left != null && left == right));
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
        // Reset the sticky "changed" signal (see requiresRecalculation()) now that this
        // whole proposal cycle -- and any TreeLikelihood read of it -- is done.
        changed = false;
    }

    @Override
    public void store() {}

    @Override
    public void restore() {
        nRestore++;
        refresh();
    }

    /**
     * Rebuild this marginal tree against the current ACG and block state.
     *
     * This is NOT an undo: use it from operators that mutate the ACG mid-proposal and need the
     * derived tree brought up to date.
     */
    public void refresh() {
        postCache = null;

        manuallyUpdated = false;
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
            validateConversionRef(c);
            _blockConvs.add(c);
        }

        return _blockConvs;
    }

    /**
     * Reject conversions that a block refers to but that are missing or degenerate.
     */
    private void validateConversionRef(Conversion conversion) {
        if (conversion == null) {
            throw new IllegalStateException(
                    "MarginalTree encountered a missing conversion referenced by block " + block.getID());
        }

        if (!conversion.isValid()) {
            throw new IllegalStateException("MarginalTree encountered invalid conversion " + conversion.getID()
                    + " referenced by block " + block.getID());
        }
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

