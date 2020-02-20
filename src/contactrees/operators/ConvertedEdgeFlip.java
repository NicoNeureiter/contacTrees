package contactrees.operators;

import contactrees.Conversion;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Nico Neureiter
 */
@Description("Operator which reverses the nodes that an edge corresponds to "
             + "leaving everything else unchanged.")
public class ConvertedEdgeFlip extends ACGOperator {

    @Override
    public double proposal() {

        if (acg.getConvCount()==0)
            return Double.NEGATIVE_INFINITY;
        
        Conversion conv = chooseConversion();
        
        Node node1 = conv.getNode1();
        Node node2 = conv.getNode2();

        conv.setNode1(node2);
        conv.setNode2(node1);

        assert !acg.isInvalid() : "CEF produced invalid state.";
        
        return 0.0;
    }
    
}
