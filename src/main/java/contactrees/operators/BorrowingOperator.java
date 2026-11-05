package contactrees.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.type.RealScalar;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
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

    public Input<RealScalar<UnitInterval>> pMoveInput = new Input<>(
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
        double pMove = pMoveInput.get().get();
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
        double pMove = pMoveInput.get().get();
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
     *
     * Caching protocol: the Gibbs moves evaluate tree likelihoods inside the proposal, before
     * BEAST stores the calculation nodes. Every evaluation therefore runs in its own
     * store()/restore() window of the TreeLikelihood (see evaluateInWindow), so that
     *  - the likelihood caches are back at the pre-proposal state when BEAST stores them, and
     *  - no node is recomputed twice between a store() and a restore(): the likelihood cores
     *    flip between two buffers per node on every recomputation, so a second recomputation
     *    would overwrite the stored buffer.
     * The marginal tree itself keeps the changes (and its dirty flags), so the final state is
     * evaluated once more by the MCMC machinery.
     * ****/

    /**
     * Compute the log likelihood of the current state of the given tree likelihood, leaving its
     * caches unchanged.
     */
    protected double evaluateLogLikelihood(TreeLikelihood treeLH) {
        MarginalTree marginalTree = getMarginalTree(treeLH);
        marginalTree.ensureUpToDate();

        // Tree unchanged since the likelihood was last computed by BEAST: use the cached value
        if (marginalTree.isLikelihoodCacheCurrent())
            return treeLH.getCurrentLogP();

        return evaluateInWindow(treeLH);
    }

    /**
     * Evaluate the tree likelihood inside its own store()/restore() window.
     */
    protected static double evaluateInWindow(TreeLikelihood treeLH) {
        treeLH.store();
        try {
            return treeLH.calculateLogP();
        } finally {
            treeLH.restore();
        }
    }

    /**
     * Utility method for calculating the posterior likelihood of borrowing
     * blocks at the given conversion.
     * @param conv Conversion edge to evaluate.
     * @param treeLH The TreeLikelihood defining the block and used for calculating the likelihood.
     * @return Log posterior probability of the borrowings on the given conversion.
     */
    public double calcLogBorrowingPosterior(Conversion conv, TreeLikelihood treeLH) {
        double pMove = pMoveInput.get().get();
        Block block = getBlock(treeLH);

        double prior = block.isAffected(conv) ? pMove : (1 - pMove);
        return Math.log(prior) + evaluateLogLikelihood(treeLH);
    }

    /**
     * Choose set of blocks to be borrowed on the given
     * conversion according to the posterior distribution.
     * @param conv Conversion edge to evaluate.
     * @return Log probability of the sampled borrowings.
     */
    public double drawBorrowingsGibbs(Conversion conv) {
        double logP = 0;
        assert treeLHs.size() > 0;
        for (TreeLikelihood treeLH : treeLHs) {
            logP += sampleBlockMove(conv, treeLH);
        }
        return logP;
    }

    /**
     * Sample the borrowing of a specific block over a specific conversion.
     * @param conv
     * @param treeLH
     * @return log-prob. difference between state before and after the move.
     */
    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH) {
        Block block = getBlock(treeLH);
        double logP;

        // Get posterior for current block move
        double logPosteriorOld = calcLogBorrowingPosterior(conv, treeLH);
        boolean checkpoint = getMarginalTree(treeLH).isCheckpoint();

        // Compute posterior for flipped block move
        flipBorrowing(block, conv);
        double logPosteriorNew = calcLogBorrowingPosterior(conv, treeLH);

        // Revert flip with probability $p_old / (p_old + p_new)$
        double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
        double pRevert = Math.exp(logPRevert);

        // Compute the hastings ratio, to remove effect of block move from acceptance-ratio
        // (We can't just return positive infinity, because of combined moves)
        if (Randomizer.nextDouble() < pRevert) {
            flipBorrowing(block, conv);
            restoreAfterRevertedFlip(getMarginalTree(treeLH), checkpoint);
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
    public double getBorrowingsProbGibbs(Conversion conv) {
        double logP = 0;

        assert treeLHs.size() > 0;

        for (TreeLikelihood treeLH : treeLHs) {
            Block block = getBlock(treeLH);

            // Get prior and likelihood for current block move
            double logPosteriorOld = calcLogBorrowingPosterior(conv, treeLH);
            boolean checkpoint = getMarginalTree(treeLH).isCheckpoint();

            // Compute prior and likelihood for flipped block move
            flipBorrowing(block, conv);
            double logPosteriorNew = calcLogBorrowingPosterior(conv, treeLH);

            // Revert flip with probability $p_old / (p_old + p_new)$
            double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
            logP += logPRevert;

            // Flip the block-move back to its original state
            flipBorrowing(block, conv);
            restoreAfterRevertedFlip(getMarginalTree(treeLH), checkpoint);
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
    /**
     * Bring the marginal tree back in line with its block after a flip was flipped back. If the
     * tree was unchanged before the flip (a checkpoint), it returns to that state without dirty
     * flags, so its likelihood is not recomputed; otherwise it is reswept.
     */
    protected static void restoreAfterRevertedFlip(MarginalTree tree, boolean checkpoint) {
        if (checkpoint)
            tree.revertToCheckpoint();
        else
            tree.ensureUpToDate();
    }

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
