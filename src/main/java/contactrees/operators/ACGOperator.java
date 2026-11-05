package contactrees.operators;

import java.util.HashSet;
import java.util.Set;

import javax.naming.directory.InvalidAttributesException;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.evolution.tree.Node;
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

    /**
     * Wraps {@link #doProposal()}: marginal trees modified by a proposal that fails directly
     * (returns -infinity) are rolled back here, since BEAST does not restore calculation nodes
     * in that case.
     */
    @Override
    public final double proposal() {
        acg.beginProposal();
        double logHR;
        try {
            logHR = doProposal();
        } catch (RuntimeException | Error e) {
            acg.rollbackMarginalTrees();
            throw e;
        }
        if (logHR == Double.NEGATIVE_INFINITY)
            acg.rollbackMarginalTrees();
        return logHR;
    }

    /**
     * The actual proposal of the operator.
     *
     * @return log Hastings ratio, or -infinity if the proposal is rejected directly.
     */
    protected abstract double doProposal();

    /*
     * Operator-driven marginal tree updates.
     *
     * After an operator finished all its edits, it may tell the marginal trees directly which of
     * them are affected. Unaffected trees are marked up to date (O(1)), affected trees are
     * updated (resweep). Trees that are not handled here update themselves lazily when they are
     * next used. Declaring a tree unaffected when it is not leads to a stale marginal tree, so
     * these helpers must only be called once all edits of the proposal are done.
     */

    /**
     * @return the marginal trees whose blocks move over any of the given conversions.
     */
    protected Set<MarginalTree> treesBorrowing(Conversion... convs) {
        Set<MarginalTree> trees = new HashSet<>();
        for (MarginalTree tree : acg.getMarginalTrees()) {
            for (Conversion conv : convs) {
                if (tree.block.isAffected(conv)) {
                    trees.add(tree);
                    break;
                }
            }
        }
        return trees;
    }

    /**
     * Update the given (affected) marginal trees and mark all others as unaffected. Only for
     * proposals that change nothing but conversions: a tree whose block follows none of the
     * changed conversions keeps its structure and its branch rates.
     */
    protected void updateMarginalTrees(Set<MarginalTree> affected) {
        for (MarginalTree tree : acg.getMarginalTrees()) {
            if (affected.contains(tree))
                tree.ensureUpToDate();
            else
                tree.markUnaffected();
        }
    }

    /**
     * Mark all marginal trees as unaffected, for edits that cannot affect any of them (e.g. adding
     * a conversion that no block follows).
     */
    protected void markAllMarginalTreesUpToDate() {
        for (MarginalTree tree : acg.getMarginalTrees())
            tree.markUnaffected();
    }

    /**
     * Update the marginal trees after a conversion changed its height but kept its CF branches
     * (node1 and node2). In a tree that follows the conversion, the topology is unchanged if the
     * conversion did not pass another conversion of the block on one of its branches; then only
     * its height is updated. Otherwise the tree is reswept. Trees that do not follow the
     * conversion are unaffected.
     */
    protected void updateMarginalTreesAfterConversionHeightChange(Conversion conv, double oldHeight) {
        double lower = Math.min(oldHeight, conv.getHeight());
        double upper = Math.max(oldHeight, conv.getHeight());
        boolean insideBranches = isStrictlyInsideBranches(conv);
        for (MarginalTree tree : acg.getMarginalTrees()) {
            if (!tree.block.isAffected(conv)) {
                tree.markUnaffected();
            } else if (insideBranches && !passesBlockConversion(tree, conv, lower, upper)) {
                tree.updateCreatorHeight(conv);
                tree.markUpToDate(lower);
            } else {
                tree.ensureUpToDate();
            }
        }
    }

    /**
     * @return true iff the conversion lies strictly between the bottom and the top of both its
     * branches (no tie with a CF event on them).
     */
    protected static boolean isStrictlyInsideBranches(Conversion conv) {
        double h = conv.getHeight();
        Node n1 = conv.getNode1(), n2 = conv.getNode2();
        return h > n1.getHeight() && h > n2.getHeight()
                && (n1.isRoot() || h < n1.getParent().getHeight())
                && (n2.isRoot() || h < n2.getParent().getHeight());
    }

    /**
     * @return true iff another conversion followed by the tree's block is attached to one of the
     * branches of conv within [lower, upper].
     */
    protected boolean passesBlockConversion(MarginalTree tree, Conversion conv, double lower, double upper) {
        Node n1 = conv.getNode1(), n2 = conv.getNode2();
        for (int id : tree.block.getConversionIDs()) {
            if (id == conv.getID())
                continue;
            Conversion other = acg.getConversions().get(id);
            double h = other.getHeight();
            if (h < lower || h > upper)
                continue;
            Node o1 = other.getNode1(), o2 = other.getNode2();
            if (o1 == n1 || o1 == n2 || o2 == n1 || o2 == n2)
                return true;
        }
        return false;
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
