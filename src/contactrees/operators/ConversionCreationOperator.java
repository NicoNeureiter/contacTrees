/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package contactrees.operators;

import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import beast.core.Input;
import beast.core.parameter.RealParameter;
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

    double pMove;
    BlockSet blockSet;
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
    	blockSet = blockSetInput.get();
        pMove = pMoveInput.get().getValue();
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
