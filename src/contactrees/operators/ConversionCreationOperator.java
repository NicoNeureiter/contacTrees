package contactrees.operators;

import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.MarginalTree;
import contactrees.util.Util;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.TreeLikelihood;
import beast.util.Randomizer;

/**
 * Abstract class of ACG operators that add new converted edges   
 * and their affected sites to an existing ConversionGraph.
 *
 * @author Nico Neureiter (nico.neureiter@gmail.com)
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

    double pMove;
    protected List<TreeLikelihood> treeLHs;
    
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        pMove = pMoveInput.get().getValue();
        treeLHs = treeLHsInput.get();
    }

    /**
     * Choose set of blocks to be affected by this conversion.
     *
     * @param conv Conversion object whose region is to be set.
     * @return log probability density of chosen attachment.
     */
    public double drawAffectedBlocks(Conversion conv) {
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
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getAffectedBlocksProb(Conversion conv) {
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
     * Choose set of blocks to be affected by this conversion according 
     * to the conditional posterior distribution.
     * 
     * @param conv
     * @return
     */
    public double drawAffectedBlocksGibbs(Conversion conv) {
        double logP = 0;
        for (TreeLikelihood treeLH : treeLHs) {
            logP += sampleBlockMove(conv, treeLH);
        }
        return logP;
    }

    /**
     * Sample the move of specific block over a specific conversion. 
     * @param conv
     * @param treeLH
     * @return log-prob. difference between state before and after the move.
     */
    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH) {
        MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
        Block block = marginalTree.block;
        
        // Get prior and likelihood for current block move
        boolean moveOld = block.isAffected(conv); 
        double priorOld = moveOld ? pMove : (1 - pMove);
        marginalTree.recalculate();
        double logLHOld = treeLH.calculateLogP();
        double logPosteriorOld = Math.log(priorOld) + logLHOld;
        
        // Compute prior and likelihood for flipped block move
        double priorNew = 1 - priorOld;
        GibbsBlockMovesOperator.flipBlockMove(block, conv);
        marginalTree.recalculate();
        double logLHNew = treeLH.calculateLogP();
        double logPosteriorNew = Math.log(priorNew) + logLHNew;
        
        // Revert flip with probability $p_old / (p_old + p_new)$
        double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
        boolean revert = (Math.log(Randomizer.nextDouble()) < logPRevert);
        
        // Compute the hastings ratio, to remove effect of block move from acceptance-ratio
        // (We can't just return positive infinity, because of combined moves)
        if (revert) {
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
            return logPRevert;
        } else {
            double pRevert = Math.exp(logPRevert);
            return Math.log(1 - pRevert);
        }
    }


    /**
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getAffectedBlocksProbGibbs(Conversion conv) {
//        return treeLHs.size() * Math.log(1-pMove);
        double logP = 0;
        for (TreeLikelihood treeLH : treeLHs) {            
            // Get block corresponding to treeLH
            MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
            Block block = marginalTree.block;

            
            // Get prior and likelihood for current block move
            boolean moveOld = block.isAffected(conv); 
            double priorOld = moveOld ? pMove : (1 - pMove);
            marginalTree.recalculate();
            double logLHOld = treeLH.calculateLogP();
            double logPosteriorOld = Math.log(priorOld) + logLHOld;
            
            // Compute prior and likelihood for flipped block move
            double priorNew = 1 - priorOld;
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
            marginalTree.recalculate();
            // marginalTree.outdated = true;
            double logLHNew = treeLH.calculateLogP();
            double logPosteriorNew = Math.log(priorNew) + logLHNew;
            
            // Revert flip with probability $p_old / (p_old + p_new)$
            double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
            logP += logPRevert;
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
                        
            assert block.isAffected(conv) == moveOld;
        }
        
        return logP;
    }
    
}
