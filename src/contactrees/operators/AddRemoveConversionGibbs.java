package contactrees.operators;

import contactrees.Conversion;

import beast.core.Description;
import beast.util.Randomizer;

/**
 * @author Nico Neureiter
 */
@Description("Operator which adds and removes conversions to/from an ACG.")
public class AddRemoveConversionGibbs extends ConversionCreationOperator {
    
    public AddRemoveConversionGibbs() { }
    
    @Override
    public double proposal() {
        double logHGF = 0;

        if (Randomizer.nextBoolean()) {
            
            // Add
            
            logHGF += Math.log(1.0/(acg.getConvCount()+1));
            logHGF -= drawNewConversion();
            
        } else {
            
            // Remove
            
            if (acg.getConvCount()==0)
                return Double.NEGATIVE_INFINITY;
            
            // Select conversion to remove:
            Conversion conv = chooseConversion();

            // Calculate HGF
            logHGF += getConversionProb(conv);
            logHGF -= Math.log(1.0/acg.getConvCount());
            
            // Remove conversion
            removeConversion(conv);
            
        }
        
        assert !acg.isInvalid() : "AddRemoveConv produced invalid state.";
        
        return logHGF;
    }
    
    /**
     * Add new conversion to ACG, returning the probability density of the
     * new edge and converted region location.
     *
     * @return log of proposal density
     */
    public double drawNewConversion() {
        Conversion newConversion = addNewConversion();

        double logP = attachEdge(newConversion);
        logP += drawAffectedBlocksGibbs(newConversion);
        
        // TODO Remove after testing is done (quite expensive)!!!
        double logP2 = getConversionProb(newConversion);
        assert Math.abs(logP - logP2) < 1E-5;

        return logP;
    }
      
    /**
     * Obtain proposal density for the move which results in the current state
     * by adding the conversion conv to a state without that recombination.
     * 
     * @param conv conversion
     * @return log of proposal density
     */
    public double getConversionProb(Conversion conv) {
        double logP; 
        logP = getEdgeAttachmentProb(conv) + getAffectedBlocksProbGibbs(conv);
        
        return logP;
    }

}
