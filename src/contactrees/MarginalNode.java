package contactrees;

import java.util.TreeMap;

import javax.naming.directory.InvalidAttributeValueException;

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
     * If >= 0, denotes the node Nr of the CF node that
     * this node corresponds to.  If <0, this node does
     * not correspond to any CF node.
     */
    public int cfNodeNr = -1;
    
    /**
     * If >= 0, denotes the ID of the conversion
     * that this node corresponds to. If <0, this
     * node comes from a CF nodes (not a conversion).
     */
    public int convID = -1;
    
    private ConversionGraph acg;
    
    public MarginalNode() {
        super(); 
//        makeDirty(Tree.IS_FILTHY);
    }
    
    public MarginalNode(ConversionGraph acg) {
        super();
        this.acg = acg;
//        makeDirty(Tree.IS_FILTHY);
    }
    
    public MarginalNode(ConversionGraph acg, int nr, double height) {
        this(acg);
        setNr(nr);
        setHeight(height);
    }
    
    public MarginalNode(ConversionGraph acg, int nr, double height, MarginalNode child1, MarginalNode child2) {
        this(acg, nr, height);
        addChild(child1);
        addChild(child2);
    }
    
    @Override
    public int isDirty() {
        return Tree.IS_FILTHY;
//        int acgNodeDirty = Tree.IS_CLEAN;
//        if (isConversionNode())
//            acgNodeDirty = getConversion().isDirty();
//        else
//            acgNodeDirty = getCFNode().isDirty();
//        
//        return Math.max(super.isDirty(), acgNodeDirty);
    }
    
    public boolean isConversionNode() {
        if (cfNodeNr >= 0) {
            assert convID == -1;
            return false;
        } else {
            assert convID >= 0;
            return true;
        }
    }

    public Conversion getConversion() {
        assert isConversionNode();
        Conversion conv = acg.getConversions().get(convID);
        if (conv == null) {
            throw new RuntimeException("Conversion ID is obsolete (probably conversion was removed from ACG).");
        }
        return conv;
    }

    public Node getCFNode() {
        assert !isConversionNode();
        Node node = acg.getNode(cfNodeNr);
        assert node != null;
        return node;
    }

    public int getCFNodeNr() {
        return cfNodeNr;
    }
    
    public int getConvID() {
        return convID;
    }
    

    /**
     * @return (deep) copy of node
     */
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
        node.cfNodeNr = cfNodeNr;
        node.convID = convID;
        
        return node;
    }
}
