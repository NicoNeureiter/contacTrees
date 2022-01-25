package contactrees.operators;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.MarginalTree;

/**
 * Gibbs operator to resample borrowings at a random conversion (contact edge).
 * @author Nico Neureiter
 */
@Description("Gibbs operator to resample borrowings at a random conversion (contact edge).")
public class GibbsSampleMovesPerConversion extends BorrowingOperator {

    public Input<Boolean> mcmcmcInput = new Input<>(
            "mcmcmc",
            "Set this to `true` when using the operator in MCMCMC (otherwise the operator samples from the cold/unheated likelihood).",
            Boolean.FALSE);

    protected BlockSet blockSet;
    protected List<MarginalTree> marginalTrees;

    @Override
    public double proposal() {
        double logHGF = 0;

        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;

        // Sample conversion to change
        Conversion conv = acg.getConversions().getRandomConversion();

        // Iterate over all blocks and resample the borrowings at the chosen contact edge
        logHGF -= drawBorrowingsGibbs(conv, true);

        if (mcmcmcInput.get())
            return logHGF;
        else
            return Double.POSITIVE_INFINITY;
    }

}
