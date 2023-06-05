package contactrees.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.MarginalTree;
import contactrees.util.Util;

/**
 * Abstract class of ACG operators that change the borrowings
 * of existing or new conversions (contact edges).
 *
 * @author Nico Neureiter
 */
public abstract class BorrowingOperator extends ACGOperator {

    final public Input<BlockSet> blockSetInput = new Input<>(
            "blockSet",
            "Block of site which are either inherited or passed via a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> pMoveInput = new Input<>(
            "pMove",
            "Probability for a block to follow a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<List<TreeLikelihood>> treeLHsInput = new Input<>(
            "treeLikelihood",
            "BEASTObject computing the tree likelihood.",
            new ArrayList<TreeLikelihood>());

    protected BlockSet blockSet;
    protected List<TreeLikelihood> treeLHs;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        blockSet = blockSetInput.get();
        treeLHs = treeLHsInput.get();
    }

    /**
     * Choose set of blocks to be borrowed at this conversion.
     *
     * @param conv Conversion whose borrowings are to be drawn.
     * @return log Probability density of chosen borrowings.
     */
    public double drawBorrowings(Conversion conv) {
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
     * Calculate probability of the borrowings at that given conversion.
     *
     * @param conv Conversion at which borrowings are evaluated.
     * @return log Probability density
     */
    public double getBorrowingsProb(Conversion conv) {
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
     * Gibbs sampling for borrowing
     * ****/


    /**
     * Utility method for calculating the posterior likelihood of borrowing
     * blocks at the given conversion.
     * @param conv Conversion edge to evaluate.
     * @param treeLH The TreeLikelihood defining the block and used for calculating the likelihood.
     * @param mtreeChanged Did the marginal tree change before calling the method (relevant for caching)?
     * @return Log posterior probability of the borrowings on the given conversion.
     */
    public double calcLogBorrowingPosterior(Conversion conv, TreeLikelihood treeLH, boolean mtreeChanged) {
        double pMove = pMoveInput.get().getValue();
        MarginalTree marginalTree = getMarginalTree(treeLH);
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

    /**
     * Choose set of blocks to be borrowed on the given
     * conversion according to the posterior distribution.
     * @param conv Conversion edge to evaluate.
     * @param mtreesChanged Did the marginal trees change before calling the method (relevant for caching)?
     * @return Log probability of the sampled borrowings.
     */
    public double drawBorrowingsGibbs(Conversion conv, boolean mtreesChanged) {
        double logP = 0;
        assert treeLHs.size() > 0;
        for (TreeLikelihood treeLH : treeLHs) {
            logP += sampleBlockMove(conv, treeLH, mtreesChanged);
        }
        return logP;
    }

    /**
     * Sample the borrowing of a specific block over a specific conversion.
     * @param conv
     * @param treeLH
     * @param mtreeChanged
     * @return log-prob. difference between state before and after the move.
     */
    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH, boolean mtreeChanged) {
        Block block = getBlock(treeLH);
        double logP;

        // Store fat calculation nodes, so that we can restore the old state
        // before state.storeCalculationNodes() is called in the MCMC object.
        treeLH.store();

        // Get posterior for current block move
        double logPosteriorOld = calcLogBorrowingPosterior(conv, treeLH, mtreeChanged);

        // Compute posterior for flipped block move
        flipBorrowing(block, conv);
        double logPosteriorNew = calcLogBorrowingPosterior(conv, treeLH, true);

        // Revert flip with probability $p_old / (p_old + p_new)$
        double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
        double pRevert = Math.exp(logPRevert);

        // Compute the hastings ratio, to remove effect of block move from acceptance-ratio
        // (We can't just return positive infinity, because of combined moves)
        if (Randomizer.nextDouble() < pRevert) {
            flipBorrowing(block, conv);
            logP = logPRevert;
        } else {
            logP = Math.log(1 - pRevert);
        }

        // Restore fat calculation nodes, so that state.storeCalculationNodes()
        // (in the MCMC object) stores the correct (old) likelihood.
        treeLH.restore();
//        ((MarginalTree) treeLH.treeInput.get()).restore();

        return logP;
    }

    /**
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getBorrowingsProbGibbs(Conversion conv, boolean mtreesChanged) {
        double logP = 0;

        assert treeLHs.size() > 0;

        for (TreeLikelihood treeLH : treeLHs) {
            // Get block corresponding to treeLH
            MarginalTree marginalTree = getMarginalTree(treeLH);
            Block block = marginalTree.block;

            // Store fat calculation nodes, so that we can restore the old state
            // before state.storeCalculationNodes() is called in the MCMC object.
            treeLH.store();

            // Get prior and likelihood for current block move
            double logPosteriorOld = calcLogBorrowingPosterior(conv, treeLH, mtreesChanged);

            // Compute prior and likelihood for flipped block move
            flipBorrowing(block, conv);
            double logPosteriorNew = calcLogBorrowingPosterior(conv, treeLH, true);

            // Revert flip with probability $p_old / (p_old + p_new)$
            double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
            logP += logPRevert;

            // Flip the block-move back to its original state
            flipBorrowing(block, conv);

            // Store fat calculation nodes, so that we can restore the old state
            // before state.storeCalculationNodes() is called in the MCMC object.
            treeLH.restore();
//            marginalTree.restore();
        }

        return logP;
    }

    /**
     * Extract the marginal tree of the given TreeLikelihood object.
     * @param treeLH
     * @return marginalTree
     */
    MarginalTree getMarginalTree(TreeLikelihood treeLH) {
        return (MarginalTree) treeLH.treeInput.get();
    }

    /**
     * Extract the block of the given TreeLikelihood object.
     * @param treeLH
     * @return block
     */
    public Block getBlock(TreeLikelihood treeLH) {
        return getMarginalTree(treeLH).block;
    }

    /**
     * Add a borrowing of the given block at the given conversion if it
     * is absent, remove the borrowing if it is present.
     * @param block
     * @param conv
     */
    public void flipBorrowing(Block block, Conversion conv) {
        if (block.isAffected(conv)) {
            block.removeMove(conv);
        } else {
            block.addMove(conv);
        }
    }

    /**
     * Include all block StateNodes in the list of affected state nodes.
     * This does not work automatically, since blocks are only indirect
     * inputs through the blockSet.
     */
    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = super.listStateNodes();
        for (Block block : blockSet.getBlocks()) {
            list.add(block);
        }
        return list;
    }
}
