/**
 * 
 */
package contactrees.operators;

import java.util.ArrayList;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.util.Util;

/**
 * An operator re-sampling a certain proportion of the Block moves over a 
 * randomly selected conversion edges (according to the prior).
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class ResampleConversionMoves extends ACGOperator {


    final public Input<BlockSet> blockSetInput = new Input<>(
            "blockSet",
            "Block of site which are either inherited or passed via a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> pMoveInput = new Input<>(
            "pMove",
            "Probability for a block to follow a conversion edge.",
            Input.Validate.REQUIRED);
    
//    public Input<Double> resampleFractionInput = new Input<>(
//            "resampleFraction",
//            "The fraction of blocks which are randomly resampled to be affected by the conversion edge or not.",
//            0.1); 

    protected BlockSet blockSet;
    protected double pMove, resampleFraction;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        blockSet = blockSetInput.get();
        pMove = pMoveInput.get().getValue();
//        resampleFraction = resampleFractionInput.get();
    }
    
    @Override
    public double proposal() {
        double logHGF = 0.0;
        
        if (acg.getConvCount() == 0) {
            return Double.NEGATIVE_INFINITY;
        }
        
        Conversion conv = acg.getConversions().getRandomConversion();
        
        logHGF += getAffectedBlocksProb(conv);
        logHGF -= drawAffectedBlocks(conv);          
        
        return logHGF;
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
            assert blockSet.getAffectedBlocks(conv).isEmpty();
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
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getAffectedBlocksProb(Conversion conv) {
        int affectedBlockCount = blockSet.getAffectedBlocks(conv).size();
        int unaffectedBlockCount = blockSet.getBlockCount() - affectedBlockCount;
        
        if (pMove == 0.) {
            assert blockSet.getAffectedBlocks(conv).isEmpty();
            return 0;
        }
        
        return affectedBlockCount*Math.log(pMove) + unaffectedBlockCount*Math.log(1-pMove);
    }

}
