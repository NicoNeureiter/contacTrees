package contactrees.operators;

import contactrees.Conversion;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Nico Neureiter
 */
@Description("Operator which moves contact edges about on clonal frame.")
public class ConvertedEdgeSlide extends ACGOperator {

    public Input<Double> apertureSizeInput = new Input<>("apertureSize",
            "Window size as a fraction of the clonal frame tree height."
                    + "Default is 0.1.", 0.1);

    public ConvertedEdgeSlide() { }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }
    
    @Override
    public double proposal() {
        
        double logHR = 0.0;
        
        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;

        // Select edge at random:
        Conversion conv = chooseConversion();
        
        // Get current (old) attachment height
        double oldHeight;
        oldHeight = conv.getHeight();
        
        // Choose window:
        double w = apertureSizeInput.get()*acg.getRoot().getHeight();

        // Set new height
        double newHeight = oldHeight + (Randomizer.nextDouble() - 0.5)*w;
        
        // Check for boundary violation
        if (newHeight>acg.getRoot().getHeight())
            return Double.NEGATIVE_INFINITY;
        
        // Get node below current (old) attachment point
        Node newNode1, newNode2;
        newNode1 = conv.getNode1();
        newNode2 = conv.getNode2();

        // Choose nodes below new attachment point
        if (newHeight<oldHeight) {
            // Fix node1          
            while (newHeight < newNode1.getHeight()) {
                if (newNode1.isLeaf())
                    return Double.NEGATIVE_INFINITY;
                
                if (Randomizer.nextBoolean())
                    newNode1 = newNode1.getLeft();
                else
                    newNode1 = newNode1.getRight();                
                logHR += -Math.log(0.5);
            }

            // Fix node2
            while (newHeight < newNode2.getHeight()) {
                if (newNode2.isLeaf())
                    return Double.NEGATIVE_INFINITY;
                
                if (Randomizer.nextBoolean())
                    newNode2 = newNode2.getLeft();
                else
                    newNode2 = newNode2.getRight();                
                logHR += -Math.log(0.5);
            }
        } else {
            while (!newNode1.isRoot() && newHeight > newNode1.getParent().getHeight()) {
                newNode1 = newNode1.getParent();
                logHR += Math.log(0.5);
            }
            while (!newNode2.isRoot() && newHeight > newNode2.getParent().getHeight()) {
                newNode2 = newNode2.getParent();
                logHR += Math.log(0.5);
            }
        }
        
        if (newNode1 == newNode2) {
            return Double.NEGATIVE_INFINITY;
        }
        
        // Write changes back to recombination object
        conv.setHeight(newHeight);
        conv.setNode1(newNode1);
        conv.setNode2(newNode2);

        assert !acg.isInvalid() : "CESlide produced invalid state";
        
        return logHR;
    }

}
