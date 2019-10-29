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
            double delta = node.getHeight() - maxChildHeight;
            newHeight = maxChildHeight + delta*f;
            logHGF += Math.log(1.0/f);
            
            // TODO just re-scale distance to the older child (to avoid the rejection below)?
            assert newHeight > maxChildHeight;
//            if (newHeight < maxChildHeight)
//                return Double.NEGATIVE_INFINITY;

        } else {
            Node parent = node.getParent();
            newHeight = Randomizer.uniform(maxChildHeight, parent.getHeight());
        }

        if (newHeight>oldHeight) {
        	logHGF -= expandConversions(node.getLeft(), node.getRight(), newHeight);
//            for (Conversion conv : acg.getConversions()) {
//                if (conv.getNode1() == node && conv.getHeight() < newHeight) {
//                    conv.setNode1(Randomizer.nextBoolean() ? leftChild : rightChild);
//                    logHGF -= logHalf;
//                }
//
//                if (conv.getNode2() == node && conv.getHeight() < newHeight) {
//                    conv.setNode2(Randomizer.nextBoolean() ? leftChild : rightChild);
//                    logHGF -= logHalf;
//                }
//            }
//
//            node.setHeight(newHeight);
//
//            if (node.isRoot()) {
//                // Draw a number of conversions
//                double L = 2.0*(newHeight-oldHeight);
//                double Nexp = L*conversionRateInput.get().getValue();
//                int N = (int)Randomizer.nextPoisson(Nexp);
//                logHGF -= -Nexp + N*Math.log(Nexp); // N! cancels
//
//                for (int i=0; i<N; i++) {
//
//                    Conversion conv = addNewConversion();
//
//                    double u = L*Randomizer.nextDouble();
//                    logHGF -= Math.log(1.0/L); 
//                    		
//                    if (u < 0.5*L) {
//                        conv.setNode1(leftChild);
//                        conv.setNode2(rightChild);
//                        conv.setHeight(oldHeight + u);
//                    } else {
//                        conv.setNode1(rightChild);
//                        conv.setNode2(leftChild);
//                        conv.setHeight(oldHeight + u - 0.5*L);
//                    }
//
//                    logHGF -= drawAffectedBlocks(conv);
//
//                }
//            }
        } else {
        	logHGF += collapseConversions(node.getLeft(), node.getRight(), newHeight);

//            List<Conversion> toRemove = new ArrayList<>();
//
//            for (Conversion conv : acg.getConversions()) {
//                if ((conv.getNode1() == leftChild || conv.getNode1() == rightChild)
//                        && conv.getHeight() > newHeight) {
//                    if (node.isRoot()) {
//                        toRemove.add(conv);
//                        continue;
//                    } else {
//                        conv.setNode1(node);
//                        logHGF += logHalf;
//                    }
//                }
//
//                if ((conv.getNode2() == leftChild || conv.getNode2() == rightChild)
//                        && conv.getHeight() > newHeight) {
//                    conv.setNode2(node);
//                    logHGF += logHalf;
//                }
//                
//                if (conv.getNode1() == conv.getNode2())
//                	return Double.NEGATIVE_INFINITY;
//            }
//
//            if (node.isRoot()) {
//                double L = 2.0*(oldHeight-newHeight);
//                double Nexp = L*conversionRateInput.get().getValue();
//                logHGF += -Nexp + toRemove.size()*Math.log(Nexp); // N! cancels
//
//                for (Conversion conv : toRemove) {
//                    logHGF += Math.log(1.0/L) + getAffectedBlocksProb(conv);
//                    acg.removeConversion(conv);
//                }
//            }
//
//            node.setHeight(newHeight);
        }

        assert !acg.isInvalid() : "CFUniform proposed invalid state.";

        return logHGF;
    }
    
}
