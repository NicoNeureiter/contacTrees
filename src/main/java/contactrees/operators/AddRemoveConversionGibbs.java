package contactrees.operators;

import beast.base.core.Description;
import beast.base.util.Randomizer;
import contactrees.Conversion;

/**
 * @author Nico Neureiter
 */
@Description("Operator which adds and removes conversions to/from an ACG.")
public class AddRemoveConversionGibbs extends ConversionCreationOperator {

    public AddRemoveConversionGibbs() { }

    boolean activateSanityChecks = false;

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

        double logPBlocks = drawBorrowingsGibbs(newConversion, false);
        if (activateSanityChecks) {
            double err = Math.abs(logPBlocks - getBorrowingsProbGibbs(newConversion, true));
            assert err < 1E-7;
        }
        logP += logPBlocks;

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
        double logP = 0;
        logP = getEdgeAttachmentProb(conv);

        double logPBlocks = getBorrowingsProbGibbs(conv, false);
        if (activateSanityChecks) {
            double err = Math.abs(logPBlocks - getBorrowingsProbGibbs(conv, true));
            assert err < 1E-7;
        }
        logP += logPBlocks;

        return logP;
    }

    /*
     *  FOR TESTING: alternative proposal function with fixed conversion
     */

    public double _proposal(Conversion conv) {
        double logHGF = 0;

        if (Randomizer.nextBoolean()) {
            if (acg.getConvCount() > 0)
                return Double.NEGATIVE_INFINITY;

            // Add
            logHGF += Math.log(1.0/(acg.getConvCount()+1));
            logHGF -= _drawNewConversion(conv);

        } else {
            // Remove
            if (acg.getConvCount()==0)
                return Double.NEGATIVE_INFINITY;

            // Calculate HGF
            logHGF += getConversionProb(conv);
            logHGF -= Math.log(1.0/acg.getConvCount());
            assert acg.getConvCount() == 1;

            // Remove conversion
            removeConversion(conv);
        }

        assert !acg.isInvalid() : "AddRemoveConv produced invalid state.";

        return logHGF;
    }

    public double _drawNewConversion(Conversion newConversion) {
        acg.addConversion(newConversion);
        return drawBorrowingsGibbs(newConversion, true);
    }

}
