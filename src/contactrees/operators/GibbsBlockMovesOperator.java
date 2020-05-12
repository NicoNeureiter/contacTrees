/**
 * 
 */
package contactrees.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.TreeLikelihood;
import beast.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.MarginalTree;
import contactrees.util.Util;


/**
 * Operator that resamples whether a given block is affected by a given conversion. 
 * 
 * @author Nico Neureiter
 */
public abstract class GibbsBlockMovesOperator extends ACGOperator {

    public Input<RealParameter> pMoveInput = new Input<>(
            "pMove",
            "Probability of a block moving over a conversion edge.",
            Input.Validate.REQUIRED);
    

    protected double pMove;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        pMove = pMoveInput.get().getValue();
    }

    public double sampleBlockMove(Conversion conv, TreeLikelihood treeLH) {
        MarginalTree marginalTree = (MarginalTree) treeLH.treeInput.get();
        Block block = marginalTree.block;
        return sampleBlockMove(block, conv, marginalTree, treeLH);
    }
    
    public double sampleBlockMove(Block block, Conversion conv, MarginalTree marginalTree, TreeLikelihood treeLH) {

        // Get prior and likelihood for current block move
        boolean moveOld = block.isAffected(conv); 
        double priorOld = moveOld ? pMove : (1 - pMove);
        marginalTree.recalculate();
        double logLHOld = treeLH.calculateLogP();
        double logPosteriorOld = Math.log(priorOld) + logLHOld;
        
        // Compute prior and likelihood for flipped block move
        double priorNew = 1 - priorOld;
        flipBlockMove(block, conv);
        marginalTree.recalculate();
        double logLHNew = treeLH.calculateLogP();
        double logPosteriorNew = Math.log(priorNew) + logLHNew;
        
        // Revert flip with probability $p_old / (p_old + p_new)$
        double logPRevert = logPosteriorOld - Util.logAddExp(logPosteriorOld, logPosteriorNew);
        boolean revert = (Math.log(Randomizer.nextDouble()) < logPRevert);
        if (revert)
            flipBlockMove(block, conv);
        
        return Double.POSITIVE_INFINITY;
    }
        
    
    static public void flipBlockMove(Block block, Conversion conv) {
        if (block.isAffected(conv)) {
            block.removeMove(conv);
        } else {
            block.addMove(conv);
        }
    }
 
    
    
}
