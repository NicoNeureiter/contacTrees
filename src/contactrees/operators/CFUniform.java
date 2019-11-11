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

import contactrees.Conversion; 
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * Uniform operator for clonal frame nodes. This operator is capable of
 * shifting an internal CF node past conversions, getting better acceptance
 * rates than the standard uniform operator when a large number of conversions
 * is present.
 *
 * @author Nico Neureiter (nico.neureiter@gmail.com)
 */
@Description("Uniform operator for clonal frame nodes.")
public class CFUniform extends CFOperator {

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Root height proposal parameter.", 0.8);

    @Override
    public double proposal() {

        double logHGF = 0.0;

        double logHalf = Math.log(0.5);

        // Select internal non-root node at random.
        Node node = acg.getNode(acg.getLeafNodeCount()
                + Randomizer.nextInt(acg.getInternalNodeCount()));
        // Choice of height is symmetric -> no effect on HGF

        Node leftChild = node.getLeft();
        Node rightChild = node.getRight();

        double oldHeight = node.getHeight();
        double maxChildHeight = Math.max(leftChild.getHeight(), rightChild.getHeight());

        // Choose new height:
        double newHeight;
        if (node.isRoot()) {
            double fMin = Math.min(scaleFactorInput.get(), 1.0/scaleFactorInput.get());
            double fMax = 1.0/fMin;

            double f = Randomizer.uniform(fMin, fMax);
            newHeight = node.getHeight() * f;
            logHGF += Math.log(1.0/f);
            
            if (newHeight < maxChildHeight)
                return Double.NEGATIVE_INFINITY;
        } else {
            Node parent = node.getParent();
            newHeight = Randomizer.uniform(maxChildHeight, parent.getHeight());
            // Choice of height is symmetric -> no effect on HGF 
        }

        if (newHeight>oldHeight) {
        	logHGF -= expandConversions(node.getLeft(), node.getRight(), newHeight);
        } else {
        	logHGF += collapseConversions(node.getLeft(), node.getRight(), newHeight);
        }
        
        if (logHGF > Double.NEGATIVE_INFINITY)
            assert !acg.isInvalid() : "CFUniform proposed invalid state.";

        return logHGF;
    }
    
}
