/**
 * 
 */
package contactrees;

import java.util.Iterator;
import java.util.List;

import contactrees.ConversionGraph;
import beast.evolution.tree.Node;

/**
 * A class representing contact events that are one-edge 
 * modifications of the clonal frame. A conversion edge connects 
 * two lineages at the same time on the clonal frame.     
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 *
 */
public class Conversion {
	
	protected ConversionGraph acg;
	
    /**
     * Nodes below branches to which the conversion edge connects.
     */
    protected Node node1, node2;
    
    /**
     * Height at which the conversion edge connects.
     */
    protected double height;
    
//    /**
//     * A list of indices of the affected sites in the alignment. 
//     */
//    protected List<Integer> affectedSites;

    public Conversion() { }

    /**
     * Construct new conversion event with specified properties.
     *
     * @param node1
     * @param node2
     * @param height
//     * @param affectedSites
     * @param acg
     */
    public Conversion(Node node1, Node node2, double height,
//            List<Integer> affectedSites,
            ConversionGraph acg) {
        this.node1 = node1;
        this.node2 = node2;
        this.height = height;
//        this.affectedSites = affectedSites;
        this.acg = acg;
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
    
//    /**
//     * Return the list of indices of the sites affected by the conversion event.
//     * @return affectedSites
//     */
//    public List<Integer> getAffectedSites() {
//        return affectedSites;
//    }
//
//    /**
//     * Set the list of indices of the sites affected by the conversion event.
//     * 
//     * @param affectedSites
//     */
//    public void setAffectedSites(List<Integer> affectedSites) {
//        startEditing();
//        this.affectedSites = affectedSites;
//    }
//
//    /**
//     * @return total number of sites affected by the conversion event.
//     */
//    public int getSiteCount() {
//        return affectedSites.size();
//    }
    
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
        
        if (node1.isRoot() || node2.isRoot())
            return false;
        
//        for (Integer idx : affectedSites) {
//	        if (idx < 0)
//	            return false;
//        	if (idx >= acg.siteCount)
//        		return false;
//        }
//
//        // affectedSite list must be sorted (and unique elements)
//        // TODO Is this really necessary? 
//        // TODO Use guava or move to util!
//        if (affectedSites.size() > 1) {
//	        Iterator<Integer> it = affectedSites.iterator();
//	        Integer i, j;
//	        i = it.next();
//	        while (it.hasNext()) {
//	        	j = it.next();
//	        	if (i >= j)
//	        		return false;
//	        	i = j;
//	        }
//        }
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
        Conversion copy = new Conversion();
        copy.node1 = node1;
        copy.node2 = node2;
        copy.height = height;
//        copy.affectedSites = affectedSites;
        copy.acg = acg;
        
//        copy.newickMetaDataBottom = newickMetaDataBottom;
//        copy.newickMetaDataMiddle = newickMetaDataMiddle;
//        copy.newickMetaDataTop = newickMetaDataTop;
        
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
        
//        // Check whether affectedSites match (length and all elements)
//        if (that.affectedSites.size() != affectedSites.size()) return false;
//        Iterator<Integer> iter = affectedSites.iterator();
//        Iterator<Integer> iterThat = that.affectedSites.iterator();
//        while (iterThat.hasNext() && iter.hasNext()) {
//        	if (iter.next() != iterThat.next())
//        		return false;
//        }
        
        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = acg != null ? acg.hashCode() : 0;
        result = 31 * result + node1.hashCode();
        result = 31 * result + node2.hashCode();
        temp = Double.doubleToLongBits(height);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        
//        // Implement caching for affectedSites list. 
//        throw new java.lang.UnsupportedOperationException("Not supported yet.");
        
        return result;
    }

    @Override
    public String toString() {
        return String.format("Conversion edge at height %d from node %d to node %d. ",
//                + "Affected sites: %s" ,
                height, node1.getNr(), node2.getNr()
//                , affectedSites.toString()
                );
    }
    

    /*
    * Bacter-style compatibility methods for tools like ACGAnnotator.
    */

    /**
     * Used by ACGAnnotator to incoroporate additional metadata into
     * the summary ACG.
     */
    public String newickMetaDataBottom, newickMetaDataMiddle, newickMetaDataTop;
    
    /**
     * Both heights (start and end of node) are the same in the contactree model. 
     * 
     * @return height
     */
    public double getHeight1() {
        return height;
    }
    public double getHeight2() {
        return height;
    }
}
