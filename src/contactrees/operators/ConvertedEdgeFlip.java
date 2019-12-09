package contactrees.operators;

import contactrees.Conversion;
import beast.core.Description;
import beast.evolution.tree.Node;

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
        
        Conversion recomb = chooseConversion();
        
        Node node1 = recomb.getNode1();
        Node node2 = recomb.getNode2();
        
        recomb.setNode1(node2);
        recomb.setNode2(node1);

        assert !acg.isInvalid() : "CEF produced invalid state.";
        
        return 0.0;
    }
    
}
