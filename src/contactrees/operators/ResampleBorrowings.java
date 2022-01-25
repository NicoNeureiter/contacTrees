/**
 *
 */
package contactrees.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;

/**
 * An operator re-sampling a certain proportion of the Block moves over a
 * randomly selected conversion edges (according to the prior).
 *
 * @author Nico Neureiter
 */
public class ResampleBorrowings extends ACGOperator {


    final public Input<BlockSet> blockSetInput = new Input<>(
            "blockSet",
            "Block of site which are either inherited or passed via a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> pMoveInput = new Input<>(
            "pMove",
            "Probability for a block to follow a conversion edge.",
            Input.Validate.REQUIRED);

    protected BlockSet blockSet;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        blockSet = blockSetInput.get();
    }

    @Override
    public double proposal() {
        double logHGF = 0.0;

        if (acg.getConvCount() == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        Conversion conv = acg.getConversions().getRandomConversion();

        logHGF += getBorrowingsProb(conv);
        logHGF -= drawBorrowings(conv);

        return logHGF;
    }

    /**
     * Choose set of borrowings for this conversion.
     *
     * @param conv Conversion object whose region is to be set.
     * @return log probability density of chosen attachment.
     */
    public double drawBorrowings(Conversion conv) {
        double pMove = pMoveInput.get().getValue();
        double logP = 0;

        if (pMove == 0.) {
            assert blockSet.getAffectedBlockIDs(conv).isEmpty();
            return 0;
        }

        for (Block block : blockSet.getBlocks())
            if (block.isAffected(conv))
                block.removeMove(conv);

        for (Block block : blockSet.getBlocks()) {
            if (Randomizer.nextDouble() < pMove) {
                block.addMove(conv);
                logP += Math.log(pMove);
            } else {
                logP += Math.log(1 - pMove);
            }
        }

        return logP;
    }

    /**
     * Calculate probability of choosing the borrowings at the
     * given conversion.
     *
     * @param conv conversion borrowings are associated with
     * @return log probability density
     */
    public double getBorrowingsProb(Conversion conv) {
        double pMove = pMoveInput.get().getValue();
        int borrowedCount = blockSet.getAffectedBlockIDs(conv).size();
        int nonborrowedCount = blockSet.getBlockCount() - borrowedCount;

        if (pMove == 0.) {
            assert blockSet.getAffectedBlockIDs(conv).isEmpty();
            return 0;
        }

        return borrowedCount*Math.log(pMove) + nonborrowedCount*Math.log(1-pMove);
    }

}
