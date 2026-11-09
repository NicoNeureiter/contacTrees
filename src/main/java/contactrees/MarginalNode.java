package contactrees;

import java.util.List;
import java.util.TreeMap;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

/**
 * A node of a {@link MarginalTree}. Marginal nodes have a fixed node number for their whole
 * life; their parent/children links and heights are rewritten by the marginal tree without
 * any StateNode bookkeeping (the marginal tree tracks dirtiness and snapshots itself).
 *
 * @author Nico Neureiter
 */
public class MarginalNode extends Node {

    /**
     * @return unmodifiable list of children of this node
     */
    @Override
    public List<Node> getChildren() {
        return children;
    }

    public MarginalNode() {
        super();
        makeDirty(Tree.IS_FILTHY);
    }

    protected void setTree(MarginalTree tree) {
        m_tree = tree;
    }

    @Override
    public void setParent(final Node parent) {
        this.parent = parent;
    }

    @Override
    public void setParent(final Node parent, final boolean inOperator) {
        this.parent = parent;
    }

    public void update(double height, MarginalNode child1, MarginalNode child2) {
        this.height = height;
        children.set(0, child1);
        children.set(1, child2);
        child1.parent = this;
        child2.parent = this;
    }

    @Override
    public void setHeight(final double height) {
        this.height = height;
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

        return node;
    }

}
