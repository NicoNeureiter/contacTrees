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

import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;

import javax.naming.directory.InvalidAttributesException;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * Abstract class of operators which act on the ConversionGraph state.
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public abstract class ACGOperator extends Operator {
    
    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Ancestral conversion graph.",
            Input.Validate.REQUIRED);
	
	protected ConversionGraph acg;

    @Override
    public void initAndValidate() {
        acg = acgInput.get();
    }
    
    /**
     * Return sister of node.
     * 
     * @param node to return sister of
     * @return sister of node
     */
    protected Node getSibling(Node node) {
        Node parent = node.getParent();
        if (parent.getLeft() == node)
            return parent.getRight();
        else
            return parent.getLeft();
    }

    /**
     * Disconnect edge <node, node.parent> from its sister and
     * grandparent (if the grandparent exists), forming a new edge
     * <sister, grandparent>. All conversions originally above node.parent
     * are re-attached to sister.
     *
     * Conversions on edge above node are not modified.
     * 
     * @param node base of edge to detach.
     */
    protected void disconnectEdge(Node node) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to disconnectEdge().");

        Node parent = node.getParent();
        Node sister = getSibling(node);

        // Remove parent and replace edges [grandparent -> parent -> sister] by [grandparent -> sister]
        if (parent.isRoot()) {
            parent.removeChild(sister);
            sister.setParent(null);
        } else {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            parent.setParent(null);
            parent.removeChild(sister);
            grandParent.addChild(sister);
        }

        // Move all conversions that attached to removed "parent" branch to the "sister" branch
        for (Conversion conv : acg.getConversions()){
            if (conv.getNode1() == parent)
                conv.setNode1(sister);

            if (conv.getNode2() == parent)
                conv.setNode2(sister);
        }
    }

    /**
     * Connect edge <node, node.parent> above destEdgeBase, forming new
     * edge <destEdgeBase, node.parent> and <node.parent, destEdgeBase.parent>.
     * All conversions above destEdgeBase that are older than destTime
     * are transferred to the new edge above node.parent.
     *
     * Conversions on edge above node are not modified.
     * 
     * @param node base of edge to attach
     * @param destEdgeBase base of edge to be bisected
     * @param destTime height at which bisection will occur
     */
    protected void connectEdge(Node node, Node destEdgeBase, double destTime) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to connectEdge().");

        // Old parent of "node" is moved to the destination 
        Node parent = node.getParent();
        
        if (destEdgeBase.isRoot()) {
            parent.addChild(destEdgeBase);
        } else {
            Node grandParent = destEdgeBase.getParent();
            grandParent.removeChild(destEdgeBase);
            grandParent.addChild(parent);
            parent.addChild(destEdgeBase);
        }

        parent.setHeight(destTime);

        // Some conversion may have to be moved from "destinationBase" to "parent"
        for (Conversion conv : acg.getConversions()) {
            if (conv.getNode1() == destEdgeBase && conv.getHeight() > destTime)
                conv.setNode1(parent);

            if (conv.getNode2() == destEdgeBase && conv.getHeight() > destTime)
                conv.setNode2(parent);
        }
    }

    /**
     * @return conversion selected uniformly at random
     * @throws InvalidAttributesException 
     */
    protected Conversion chooseConversion() throws InvalidAttributesException {
		return acg.getConversions().getRandomConversion();
    }

    /**
     * @return maximum height of children of root
     */
    protected double getMaxRootChildHeight() {
        double max = 0.0;
        for (Node child : acg.getRoot().getChildren())
            max = Math.max(max, child.getHeight());

        return max;
    }
}
