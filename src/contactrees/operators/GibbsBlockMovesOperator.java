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

    public double sampleBlockMove(Block block, Conversion conv, MarginalTree blockTree, TreeLikelihood treeLH) {
        
        boolean moveOld = block.isAffected(conv); 
        double priorOld = moveOld ? pMove : (1 - pMove);
        double logLHOld = treeLH.getCurrentLogP();
        double logPosteriorOld = Math.log(priorOld) + logLHOld;
        
        // Compute prior and likelihood for flipped block move
        double priorNew = 1 - priorOld;
        // Flip the block move
        flipBlockMove(block, conv);
        // Update the marginal tree and the likelihood
        blockTree.recalculate();
        double logLHNew = treeLH.calculateLogP();
        double logPosteriorNew = Math.log(priorNew) + logLHNew;
        
        // Revert flip with probability $p_old / (p_old + p_new)$
        double logPRevert = logPosteriorOld - logAddExp(logPosteriorOld, logPosteriorNew);
        if (Math.log(Randomizer.nextDouble()) < logPRevert)
            flipBlockMove(block, conv);
        
        return Double.POSITIVE_INFINITY;
    }
    
    static public double logAddExp(double logP1, double logP2) {
        double logPMin, logPMax;
        if (logP1 < logP2) {
            logPMin = logP1;
            logPMax = logP2;
        } else {
            logPMin = logP2;
            logPMax = logP1;
        }
        
        return logPMax + Math.log1p(Math.exp(logPMin - logPMax));
    }
        
    
    static public void flipBlockMove(Block block, Conversion conv) {
        if (block.isAffected(conv)) {
            block.removeMove(conv);
        } else {
            block.addMove(conv);
        }
    }
 
    
    
}
