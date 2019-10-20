/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("Scaling operator for recombination graphs.")
public class ACGScaler extends ACGOperator {

    public Input<List<RealParameter>> parametersInput = new Input<>(
            "parameter", "Parameter to scale with ARG.",
            new ArrayList<>());

    public Input<List<RealParameter>> parametersInverseInput = new Input<>(
            "parameterInverse", "Parameter to scale inversely with ARG.",
            new ArrayList<>());
    
    public Input<Double> scaleParamInput = new Input<>("scaleFactor",
            "Scale factor tuning parameter.  Must be < 1.",
            Input.Validate.REQUIRED);
    
    public Input<Boolean> rootOnlyInput = new Input<>(
            "rootOnly", "Scale root node and connections which attach directly "
                    + "below and above the root only.",
            false);
    
    private double scaleParam;
    private boolean rootOnly;
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        scaleParam = scaleParamInput.get();
        rootOnly = rootOnlyInput.get();
    }
    
    @Override
    public double proposal() {
        
        // Keep track of number of positively scaled elements minus
        // negatively scaled elements.
        int count = 0;

        // Choose scaling factor:
        double f = scaleParam + Randomizer.nextDouble()*(1.0/scaleParam - scaleParam);
        
        // Scale clonal frame:
        if (rootOnly) {
            acg.getRoot().setHeight(acg.getRoot().getHeight()*f);
            count += 1;
        } else {
            for (Node node : acg.getInternalNodes()) {
                node.setHeight(node.getHeight()*f);
                count += 1;
            }
        }       
        
        // Scale conversion edges:
        for (Conversion conv : acg.getConversions()) {
        	Node node1 = conv.getNode1();
        	Node node2 = conv.getNode2();
        	
        	boolean isRootChild = node1.getParent().isRoot(); 
        	
        	if (rootOnly) {
        		if (node2.getParent().isRoot() != isRootChild) {
        			// One node gets rescaled, the other does not -> Illegal move
        			return Double.NEGATIVE_INFINITY;
        		}
        	}
        	
            if (!rootOnly || isRootChild) {
                conv.setHeight(conv.getHeight() * f);
                count += 1;
            }

            if (conv.getHeight() < node1.getHeight()
                    || conv.getHeight() < node2.getHeight())
            	// Conversion below its node -> Illegal move
                return Double.NEGATIVE_INFINITY;
        }
        
        // Check for illegal node heights:
        if (rootOnly) {
            for (Node node : acg.getRoot().getChildren()) {
                if (node.getHeight()>node.getParent().getHeight())
                    return Double.NEGATIVE_INFINITY;
            }
        } else {
            for (Node node : acg.getExternalNodes()) {
                if (node.getHeight()>node.getParent().getHeight())
                    return Double.NEGATIVE_INFINITY;
            }
        }
        
        // Scale parameters
        
        for (RealParameter param : parametersInput.get()) {
            try {
                param.startEditing(null);
                param.scale(f);
            } catch (Exception e) {
                
                // Scale throws Exception if param has been scaled outside its
                // bounds.  Needs to change!
                return Double.NEGATIVE_INFINITY;
            }
            
            count += param.getDimension();
        }
        
        for (RealParameter paramInv : parametersInverseInput.get()) {
            try {
                paramInv.startEditing(null);
                paramInv.scale(1.0/f);
            } catch (Exception e) {
                
                // Scale throws Exception if param has been scaled outside its
                // bounds.  Needs to change!
                return Double.NEGATIVE_INFINITY;
            }
            
            count -= paramInv.getDimension();
        }

        assert !acg.isInvalid() : "ACGScaler produced invalid state.";
        
        // Return log of Hastings ratio:
        return (count-2)*Math.log(f);
    }
    
}
