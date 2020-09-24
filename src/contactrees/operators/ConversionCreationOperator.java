package contactrees.operators;

import contactrees.Block;
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
     * Calculate probability of choosing region affected by the
     * given conversion.
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
     * Choose set of blocks to be affected by this conversion according 
     * to the conditional posterior distribution.
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
    public double drawAffectedBlocksGibbs(Conversion conv) {
        return drawAffectedBlocksGibbs(conv, true);
    }
    

//    protected double evaluateTreeLH() {
//        double logP = 0.0;
//
//        State state = stateInput.get();
//        state.storeCalculationNodes();
//        state.checkCalculationNodesDirtiness();
//
//        try {
//            logP = treeLG.calculateLogP();
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }
//
//        state.restore();
//        state.store(currentState);
//
//        return logP;
//    };

    /**
     * Sample the move of specific block over a specific conversion. 
     * @param conv
     * @param treeLH
     * @return log-prob. difference between state before and after the move.
     */
    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH, boolean mtreeChanged) {
        double pMove = pMoveInput.get().getValue();
        
        MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
        Block block = marginalTree.block;

        treeLH.store();
//        block.store();
        
        // Get prior and likelihood for current block move
        boolean moveOld = block.isAffected(conv); 
        double priorOld = moveOld ? pMove : (1 - pMove);
        double logLHOld;
        
        if (mtreeChanged) {
            marginalTree.recalculate();
            logLHOld = treeLH.calculateLogP();
        } else {
            logLHOld= treeLH.getCurrentLogP();
        }

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
        double p = 0.0;
        if (revert) {
//            System.out.println("°");
            GibbsBlockMovesOperator.flipBlockMove(block, conv);
            treeLH.restore();
//            marginalTree.setManuallyUpdated();
            p = logPRevert;
        } else {
//            System.out.println("'");
            double pRevert = Math.exp(logPRevert);
            // Not optimal, but we need to restore here, so that beast.core.MCMC doesn't store the updated likelihood.
            // (would lead to wrong restored likelihood, if the operator is rejected)
            treeLH.restore();   
            marginalTree.setManuallyUpdated();
            p = Math.log(1 - pRevert);
        }
        
//        block.setEverythingDirty(true);
//        skipNextStateInitialisation = true;
        
        return p;
    }
    
    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH) {
        return sampleBlockMove(conv, treeLH, true);
    }
    
    boolean skipNextStateInitialisation = false;

    @Override
    public boolean requiresStateInitialisation() {

        // TODO See whether we can use this to avoid unnecesary restore
        // Store ACG and Blocks before!
        
        if (skipNextStateInitialisation) {
            skipNextStateInitialisation = false;
            return false;
        } else {
            return true;  
        }
    }

    /**
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */   
    public double getAffectedBlocksProbGibbs(Conversion conv, boolean mtreesChanged) {
        double pMove = pMoveInput.get().getValue();
        double logP = 0;
        
        for (TreeLikelihood treeLH : treeLHsInput.get()) {            
            // Get block corresponding to treeLH
            MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
            Block block = marginalTree.block;

            treeLH.store();
            
            // Get prior and likelihood for current block move
            boolean moveOld = block.isAffected(conv); 
            double priorOld = moveOld ? pMove : (1 - pMove);
            double logLHOld;
            
            if (mtreesChanged) {
                marginalTree.recalculate();
                logLHOld = treeLH.calculateLogP();
            } else {
                logLHOld= treeLH.getCurrentLogP();
            }

//            marginalTree.recalculate();
//            double logLHOld = treeLH.calculateLogP();
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
            
            treeLH.restore();
//            marginalTree.setManuallyUpdated();
        }
        
        return logP;
    }

    public double getAffectedBlocksProbGibbs(Conversion conv) {
        return getAffectedBlocksProbGibbs(conv, true);
    }
    
}
