package contactrees;

import java.util.List;
import java.util.TreeMap;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * Adds a field to the Node class specifying whether a given
 * node in a marginal tree corresponds to a node in the clonal frame.
 *
 * @author Nico Neureiter
 */
public class MarginalNode extends Node {

    /**
     * The length of this branch in time, i.e. scaled by the branch-rate.
     * This is used in the creation of the marginal tree, since scaling
     * branch lengths is easier than scaling heights.
     */
    public double timeLength = 0.0;

    /**
     *
     */
    public double lastEventHeight = 0.0;

    /**
     * @return unmodifiable list of children of this node
     */
    @Override
    public List<Node> getChildren() {
        return children;
    }

    private ConversionGraph acg;

    public MarginalNode() {
        super();
        makeDirty(Tree.IS_FILTHY);
    }

    public MarginalNode(ConversionGraph acg) {
        super();
        this.acg = acg;
        makeDirty(Tree.IS_FILTHY);
    }

    public MarginalNode(ConversionGraph acg, int nr, double height) {
        this(acg);
        setNr(nr);
        this.height = height;
    }

    public MarginalNode(ConversionGraph acg, int nr, double height, MarginalNode child1, MarginalNode child2) {
        this(acg, nr, height);

        children.add(child1);
        children.add(child2);
        child1.parent = this;
        child2.parent = this;
    }

    protected void setTree(MarginalTree tree) {
        m_tree = tree;
    }

    @Override
    public void setParent(final Node parent) {
        this.parent = parent;
    }

    @Override
    protected void setParent(final Node parent, final boolean inOperator) {
        this.parent = parent;
    }

    public void update(double height, MarginalNode child1, MarginalNode child2) {
        this.height = height;
        children.set(0, child1);
        children.set(1, child2);
        child1.parent = this;
        child2.parent = this;
        makeDirty(Tree.IS_FILTHY);
    }

    @Override
    public void setHeight(final double height) {
        this.height = height;
        makeDirty(Tree.IS_DIRTY);
    }

    /**
     * @return (deep) copy of node
     */
    @Override
    public Node copy() {
        final MarginalNode node = new MarginalNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.lengthMetaDataString = lengthMetaDataString;
        node.metaData = new TreeMap<>(metaData);
        node.lengthMetaData = new TreeMap<>(lengthMetaData);
        node.setParent(null);
        node.setID(getID());

        for (final Node child : getChildren()) {
            node.addChild(child.copy());
        }

        node.acg = acg;

        return node;
    }

    public void setTimeLength(double timeLength) {
        this.timeLength = timeLength;
    }

    public void addToTimeLength(double timeLength) {
        this.timeLength += timeLength;
    }

    public void addToTimeLength(double length, double branchRate) {
        this.timeLength += length * branchRate;
    }

}
