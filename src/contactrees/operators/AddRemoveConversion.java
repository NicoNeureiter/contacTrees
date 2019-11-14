/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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

import contactrees.Conversion;

import javax.naming.directory.InvalidAttributesException;

import beast.core.Description;
import beast.util.Randomizer;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("Operator which adds and removes conversions to/from an ACG.")
public class AddRemoveConversion extends ConversionCreationOperator {
    
    public AddRemoveConversion() { }
    
    @Override
    public double proposal() {
        double logHGF = 0;

        if (Randomizer.nextBoolean()) {
            
            // Add
            
            logHGF += Math.log(1.0/(acg.getConvCount()+1));
            logHGF -= drawNewConversion();
            
        } else {
            
            // Remove
            
            if (acg.getConvCount()==0)
                return Double.NEGATIVE_INFINITY;
            
            // Select conversion to remove:
            Conversion conv = chooseConversion();

            // Calculate HGF
            logHGF += getConversionProb(conv);
            logHGF -= Math.log(1.0/acg.getConvCount());
            
            // Remove conversion
            removeConversion(conv);
            
        }
        
        assert !acg.isInvalid() : "AddRemoveConv produced invalid state.";
        
        return logHGF;
    }
    
    /**
     * Add new conversion to ACG, returning the probability density of the
     * new edge and converted region location.
     *
     * @return log of proposal density
     */
    public double drawNewConversion() {
        Conversion newConversion = addNewConversion();

        double logP = attachEdge(newConversion) + drawAffectedBlocks(newConversion);

        return logP;
    }
      
    /**
     * Obtain proposal density for the move which results in the current state
     * by adding the conversion conv to a state without that recombination.
     * 
     * @param conv conversion
     * @return log of proposal density
     */
    public double getConversionProb(Conversion conv) {
        return getEdgeAttachmentProb(conv) + getAffectedBlocksProb(conv);
    }

}
