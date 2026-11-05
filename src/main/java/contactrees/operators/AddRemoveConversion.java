package contactrees.operators;

import java.util.Set;

import beast.base.core.Description;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.MarginalTree;

/**
 * @author Nico Neureiter
 */
@Description("Operator which adds and removes conversions to/from an ACG.")
public class AddRemoveConversion extends ConversionCreationOperator {

    public AddRemoveConversion() { }

    @Override
    public double doProposal() {
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

            // Remove conversion (only the trees following it change)
            Set<MarginalTree> affected = treesBorrowing(conv);
            removeConversion(conv);
            updateMarginalTrees(affected);
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
        double logP = attachEdge(newConversion) + drawBorrowings(newConversion);

        // Only the trees following the new conversion change
        updateMarginalTrees(treesBorrowing(newConversion));
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
        return getEdgeAttachmentProb(conv) + getBorrowingsProb(conv);
    }

}
