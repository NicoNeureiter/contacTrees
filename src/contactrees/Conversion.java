package contactrees;

import contactrees.ConversionGraph;
import beast.evolution.tree.Node;

/**
 * A class representing contact events that are one-edge 
 * modifications of the clonal frame. A conversion edge connects 
 * two lineages at the same time on the clonal frame.     
 * 
 * @author Nico Neureiter
 */

public class Conversion {
	
	// TODO Remove "acg" from Conversion?
	protected ConversionGraph acg;
	
    /**
     * Nodes below branches to which the conversion edge connects.
     */
    protected Node node1, node2;
    
    /**
     * Height at which the conversion edge connects.
     */
    protected double height;
    
    /**
     * ID of the conversion in the conversion list.
     */
    protected Integer id;
    
    /**
     * Keep track whether this conversion changed within this step.
     */
    protected boolean hasStartedEditing = false;
    
    public Conversion(int id) {
    	this.id = id;
    }

    /**
     * Construct new conversion event with specified properties.
     *
     * @param node1
     * @param node2
     * @param height
//     * @param affectedSites
     * @param acg
     */
    public Conversion(Node node1, Node node2, double height, ConversionGraph acg, int id) {
        this.node1 = node1;
        this.node2 = node2;
        this.height = height;
        this.acg = acg;
        this.id = id;
    }

    /**
     * Obtain node below the starting-point of the conversion edge.
     * 
     * @return node
     */
    public Node getNode1() {
        return node1;
    }

    /**
     * Obtain node below the end-point of the conversion edge.
     * 
     * @return node
     */
    public Node getNode2() {
        return node2;
    }
    
    /**
     * Return height of the conversion edge in the clonal frame.
     * 
     * @return height
     */
    public double getHeight() {
        return height;
    }

    /**
     * Return ID of the conversion edge in the conversion list.
     * 
     * @return ID
     */
    public int getID() {
    	return id;
    }
    
    /**
     * Set node below the starting-point of the conversion edge.
     * 
     * @param node1 
     */
    public void setNode1(Node node1) {
        startEditing();
        this.node1 = node1;
    }

    /**
     * Set node below the end-point of the conversion edge.
     * 
     * @param node2 
     */
    public void setNode2(Node node2) {
        startEditing();
        this.node2 = node2;
    }

    /**
     * Set height of the conversion edge in the clonal frame.
     * 
     * @param height
     */
    public void setHeight(double height) {
        startEditing();
        this.height = height;
    }
    
    /**
     * Set ID of the conversion edge in the conversion list.
     * 
     * @param ID
     */
    public void setID(int id) {
    	startEditing();
    	this.id = id;
    }
    
    /**
     * Check validity of recombination specification: whether specified heights
     * belong to edges above specified nodes.
     * 
     * @return true if specification is valid
     */
    public boolean isValid() {
        if ((node1.getHeight() > height) || (node2.getHeight() > height))
            return false;
        
        if (node1.getParent().getHeight() < height)
        	return false;
        if (node2.getParent().getHeight() < height)
        	return false;
        
        if (node1 == node2)
        	return false;
        
        // General sanity checks
        assert !(node1.isRoot() || node2.isRoot());
        assert acg.getConversions().get(id) == this;
        
        return true;
    }
    
    /**
     * Assign conversion graph.
     * @param acg 
     */
    public void setConversionGraph(ConversionGraph acg) {
        this.acg = acg;
    }

    /**
     * Mark ARG statenode as dirty if available.
     */
    public void startEditing() {
        // TODO: hasStartedEditing = true <<<<<<<<<<<<<<<<
        // ...
        if (acg != null)
            acg.startEditing(null);
    }
    
    /**
     * Obtain new recombination with exactly the same
     * field values as this one.
     *
     * @return copy of Conversion object
     */
    public Conversion getCopy() {
        Conversion copy = new Conversion(id);
        
        copy.node1 = node1;
        copy.node2 = node2;
        copy.height = height;
        copy.acg = acg;
        
        copy.newickMetaDataBottom = newickMetaDataBottom;
        copy.newickMetaDataMiddle = newickMetaDataMiddle;
        copy.newickMetaDataTop = newickMetaDataTop;
        
        return copy;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Conversion that = (Conversion) o;

        if (!node1.equals(that.node1)) return false;
        if (!node2.equals(that.node2)) return false;
        if (Double.compare(that.height, height) != 0) return false;
        if (acg != null ? !acg.equals(that.acg) : that.acg != null)
            return false;
        
        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = acg != null ? acg.hashCode() : 0;
        result = 31 * result + node1.hashCode();
        result = 31 * result + node2.hashCode();
        result = 31 * result + id; 				// TODO is this right?
        temp = Double.doubleToLongBits(height);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
                
        return result;
    }

    @Override
    public String toString() {
        return String.format("Conversion edge at height %f from node %d to node %d. ",
                			 height, node1.getNr(), node2.getNr());
    }
    

    /*
    * Bacter-style compatibility methods for tools like ACGAnnotator.
    */

    /**
     * Used by ACGAnnotator to incoroporate additional metadata into
     * the summary ACG.
     */
    public String newickMetaDataBottom, newickMetaDataMiddle, newickMetaDataTop;
    
}
