package contactrees.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.TreeLikelihood;
import beast.util.Randomizer;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalTree;
import contactrees.util.Util;

/**
 * Abstract class of ACG operators that add new converted edges
 * and their affected sites to an existing ConversionGraph.
 *
 * @author Nico Neureiter
 */
public abstract class ConversionCreationOperator extends EdgeCreationOperator {

    public Input<RealParameter> pMoveInput = new Input<>(
            "pMove",
            "Probability for a block to follow a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<List<TreeLikelihood>> treeLHsInput = new Input<>(
            "treeLikelihood",
            "BEASTObject computing the tree likelihood.",
            new ArrayList<TreeLikelihood>());

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    /**
     * Choose set of blocks to be affected by this conversion.
     *
     * @param conv Conversion object whose region is to be set.
     * @return log probability density of chosen attachment.
     */
    public double drawAffectedBlocks(Conversion conv) {
        double pMove = pMoveInput.get().getValue();
        double logP = 0;

        if (pMove == 0.) {
            assert blockSet.getAffectedBlockIDs(conv).isEmpty();
            return 0;
        }

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
     * Calculate probability of choosing blocks that are
     * affected by the given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getAffectedBlocksProb(Conversion conv) {
        double pMove = pMoveInput.get().getValue();
        int affectedBlockCount = blockSet.getAffectedBlockIDs(conv).size();
        int unaffectedBlockCount = blockSet.getBlockCount() - affectedBlockCount;

        if (pMove == 0.) {
        	assert blockSet.getAffectedBlockIDs(conv).isEmpty();
        	return 0;
        }

        return affectedBlockCount*Math.log(pMove) + unaffectedBlockCount*Math.log(1-pMove);
    }


    /* ****
     * Gibbs sampling for affected blocks
     * ****/

    /**
     * Choose set of blocks to be affected by the given
     * conversion according to the posterior distribution.
     *
     * @param conv
     * @return
     */
    public double drawAffectedBlocksGibbs(Conversion conv, boolean mtreesChanged) {
        double logP = 0;
        for (TreeLikelihood treeLH : treeLHsInput.get()) {
            logP += sampleBlockMove(conv, treeLH, mtreesChanged);
        }
        return logP;
    }

    /**
     * Sample the move of specific block over a specific conversion.
     * @param conv
     * @param treeLH
     * @param mtreeChanged
     * @return log-prob. difference between state before and after the move.
     */
    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH, boolean mtreeChanged) {
        Block block = getBlock(treeLH);
        double logP;

        // Get posterior for current block move
        double logPosteriorOld = calcLogBlockPosterior(conv, treeLH, mtreeChanged);

        // Compute posterior for flipped block move
        GibbsBlockMovesOperator.flipBlockMove(block, conv);
        double logPosteriorNew = calcLogBlockPosterior(conv, treeLH, true);

        // Revert flip with probability $p_old / (p_old + p_new)$
        double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
        double pRevert = Math.exp(logPRevert);

        // Compute the hastings ratio, to remove effect of block move from acceptance-ratio
        // (We can't just return positive infinity, because of combined moves)
        if (Randomizer.nextDouble() < pRevert) {
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
            logP = logPRevert;
        } else {
            logP = Math.log(1 - pRevert);
        }

        return logP;
    }

    /**
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getAffectedBlocksProbGibbs(Conversion conv, boolean mtreesChanged) {
        double logP = 0;

        for (TreeLikelihood treeLH : treeLHsInput.get()) {
            // Get block corresponding to treeLH
            MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
            Block block = marginalTree.block;

            // Get prior and likelihood for current block move
            double logPosteriorOld = calcLogBlockPosterior(conv, treeLH, mtreesChanged);

            // Compute prior and likelihood for flipped block move
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
            double logPosteriorNew = calcLogBlockPosterior(conv, treeLH, true);

            // Revert flip with probability $p_old / (p_old + p_new)$
            double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
            logP += logPRevert;

            // Flip the block-move back to its original state
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
        }

        return logP;
    }

    /*
     * The following methods define mtreesChanged=true as a default.
     */

    public double drawAffectedBlocksGibbs(Conversion conv) {
        return drawAffectedBlocksGibbs(conv, true);
    }

    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH) {
        return sampleBlockMove(conv, treeLH, true);
    }

    public double getAffectedBlocksProbGibbs(Conversion conv) {
        return getAffectedBlocksProbGibbs(conv, true);
    }

    public double calcLogBlockPosterior(Conversion conv, TreeLikelihood treeLH, boolean mtreeChanged) {
        double pMove = pMoveInput.get().getValue();
        MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
        Block block = marginalTree.block;

        double prior = block.isAffected(conv) ? pMove : (1 - pMove);
        double logLH;
        if (mtreeChanged) {
            marginalTree.recalculate();
            logLH = treeLH.calculateLogP();
        } else {
            logLH= treeLH.getCurrentLogP();
        }
        return Math.log(prior) + logLH;
    }

    public Block getBlock(TreeLikelihood treeLH) {
        MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
        return marginalTree.block;
    }

}
