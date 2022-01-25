/**
 *
 */
package contactrees.operators;

import beast.core.Input;
import beast.evolution.likelihood.TreeLikelihood;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalTree;

/**
 *
 *
 * @author Nico Neureiter
 */
public class GibbsSampleMovesPerBlock extends BorrowingOperator {

    public Input<TreeLikelihood> treeLHInput = new Input<>(
            "treeLikelihood",
            "BEASTObject computing the tree likelihood.",
            Input.Validate.REQUIRED);

    protected Block block;
    protected TreeLikelihood treeLH;
    protected MarginalTree marginalTree;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        treeLH = treeLHInput.get();
    }

    @Override
    public double proposal() {
        for (Conversion conv : acg.getConversions()) {
            sampleBlockMove(conv, treeLH, false);
        }

        return Double.POSITIVE_INFINITY;
    }


}
