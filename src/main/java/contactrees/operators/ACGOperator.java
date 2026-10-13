package contactrees.operators;

import javax.naming.directory.InvalidAttributesException;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.evolution.tree.Node;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.MarginalTree;

/**
 * Abstract class of operators which act on the ConversionGraph state.
 *
 * @author Nico Neureiter
 */
public abstract class ACGOperator extends Operator {

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Ancestral conversion graph.",
            Input.Validate.REQUIRED);

	protected ConversionGraph acg;

    @Override
    public void initAndValidate() {
        acg = acgInput.get();
    }

    @Override
    public final double proposal() {
        double logHGF = doProposal();
        if (logHGF > Double.NEGATIVE_INFINITY && operatorUpdatesMarginalTree()) {
            for (MarginalTree marginalTree : acg.getMarginalTrees()) {
                if (marginalTree.isValidating())
                    marginalTree.assertConsistent(getID());
            }
        }
        return logHGF;
    }

    /** Perform the actual proposal. See {@link #proposal()}. */
    protected abstract double doProposal();

    /**
     * @return true if this operator's doProposal() already leaves every marginal tree it could
     *         have affected fully, FINALLY recalculated (synced) by the time it returns.
     *         {@link MarginalTree#assertConsistent} compares the tree's CURRENT,
     *         not-yet-recalculated display against a fresh rebuild, so calling it from
     *         {@link #proposal()} is only valid when this returns true -- otherwise the
     *         operator's proposal would report a false "stale tree".
     *         <p>
     *         Default false, and -- as of this writing -- no operator overrides it, including the
     *         Gibbs-sampling ones in {@link BorrowingOperator}: {@code sampleBlockMove()}'s own
     *         choreography (evaluate old state -&gt; flip -&gt; recalculate for the flipped state
     *         -&gt; maybe flip back with NO third recalculate) deliberately leaves a marginal tree
     *         displaying the flipped topology whenever the block's move reverts, relying on
     *         BEAST's lazy {@code requiresRecalculation()} to self-heal on the next dirty-check --
     *         this is correct, load-bearing, and extensively tested behaviour (see
     *         {@code MarginalTreeGibbsCacheTest}), not a bug, but it means even the Gibbs
     *         operators cannot honestly claim "already synced" at the point {@link #proposal()}
     *         would ask. Every non-Gibbs operator relies on the same lazy self-healing even more
     *         directly (a CF-height change alone would otherwise force recalculating every
     *         marginal tree, not just the ones actually consulted before the next accept/reject).
     *         This hook stays real, general infrastructure for a hypothetical future operator
     *         that genuinely never needs that self-healing -- flip it on for one only after
     *         checking it doesn't go through {@code sampleBlockMove}'s revert path.
     */
    protected boolean operatorUpdatesMarginalTree() {
        return false;
    }

    /**
     * Report that a block's move set changed at an existing, unmoved conversion's height (a
     * plain borrowing flip, not a re-attachment or a newly created conversion) to the block's
     * own marginal tree, if it has one registered on this ACG.
     */
    protected void noteBlockMoveChanged(Block block, Conversion conv) {
        MarginalTree marginalTree = acg.getMarginalTree(block);
        if (marginalTree != null)
            marginalTree.noteChangedAtHeight(conv.getHeight());
    }

    /**
     * Return sister of node.
     *
     * @param node to return sister of
     * @return sister of node
     */
    protected Node getSibling(Node node) {
        Node parent = node.getParent();
        if (parent.getLeft() == node)
            return parent.getRight();
        else
            return parent.getLeft();
    }

    /**
     * Disconnect edge <node, node.parent> from its sister and
     * grandparent (if the grandparent exists), forming a new edge
     * <sister, grandparent>. All conversions originally above node.parent
     * are re-attached to sister.
     *
     * Conversions on edge above node are not modified.
     *
     * @param node base of edge to detach.
     */
    protected void disconnectEdge(Node node) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to disconnectEdge().");

        Node parent = node.getParent();
        Node sister = getSibling(node);

        // Remove parent and replace edges [grandparent -> parent -> sister] by [grandparent -> sister]
        if (parent.isRoot()) {
            parent.removeChild(sister);
            sister.setParent(null);
        } else {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            parent.setParent(null);
            parent.removeChild(sister);
            grandParent.addChild(sister);
        }

        // Move all conversions that attached to removed "parent" branch to the "sister" branch
        for (Conversion conv : acg.getConversions()){
            if (conv.getNode1() == parent)
                conv.setNode1(sister);

            if (conv.getNode2() == parent)
                conv.setNode2(sister);
        }
    }

    /**
     * Connect edge <node, node.parent> above destEdgeBase, forming new
     * edge <destEdgeBase, node.parent> and <node.parent, destEdgeBase.parent>.
     * All conversions above destEdgeBase that are older than destTime
     * are transferred to the new edge above node.parent.
     *
     * Conversions on edge above node are not modified.
     *
     * @param node base of edge to attach
     * @param destEdgeBase base of edge to be bisected
     * @param destTime height at which bisection will occur
     */
    protected void connectEdge(Node node, Node destEdgeBase, double destTime) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to connectEdge().");

        // Old parent of "node" is moved to the destination
        Node parent = node.getParent();

        if (destEdgeBase.isRoot()) {
            parent.addChild(destEdgeBase);
        } else {
            Node grandParent = destEdgeBase.getParent();
            grandParent.removeChild(destEdgeBase);
            grandParent.addChild(parent);
            parent.addChild(destEdgeBase);
        }

        parent.setHeight(destTime);

        // Some conversion may have to be moved from "destinationBase" to "parent"
        for (Conversion conv : acg.getConversions()) {
            if (conv.getNode1() == destEdgeBase && conv.getHeight() > destTime)
                conv.setNode1(parent);

            if (conv.getNode2() == destEdgeBase && conv.getHeight() > destTime)
                conv.setNode2(parent);
        }
    }

    /**
     * @return conversion selected uniformly at random
     * @throws InvalidAttributesException
     */
    protected Conversion chooseConversion() {
		return acg.getConversions().getRandomConversion();
    }

    /**
     * @return maximum height of children of root
     */
    protected double getMaxRootChildHeight() {
        double max = 0.0;
        for (Node child : acg.getRoot().getChildren())
            max = Math.max(max, child.getHeight());

        return max;
    }

}
