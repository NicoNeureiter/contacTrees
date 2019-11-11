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

package contactrees.test;

import contactrees.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

import java.util.List;

/**
 * Unit test for marginal tree traversal.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MarginalTreeTest extends ContactreesTest {
    
    public MarginalTreeTest() {
    }

    @Test
    public void testNonOverlapping() throws Exception {

        // Test all marginals against truth
        // (I have eyeballed each of these trees and claim that they are correct.)
        String[] correctNewickStrings = {
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
            "(1:2.5,(2:0.5,3:0.5)4:2.0)5:0.5",  // Conv 1
            "((1:1.0,2:1.0)4:0.5,3:1.5)5:1.5",  // Conv 2
            "(1:1.5,(2:0.5,3:0.5)4:1.0)5:1.5",  // Conv 1 & 2
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
        };
        
        
        List<Block> blocks = blockSet.getBlocks();
        blocks.get(1).addMove(conv1);
        blocks.get(2).addMove(conv2);
        blocks.get(3).addMove(conv1);
        blocks.get(3).addMove(conv2);
        
        for (int b=0; b<N_BLOCKS; b++) {         
            MarginalTree marginalTree = new MarginalTree();
            marginalTree.initByName("network", acg, "block", blocks.get(b), "nodetype", MarginalNode.class.getName());

            //System.out.println(marginalTree + ";");
            String newickStr = correctNewickStrings[b];
            Tree correctTree = new TreeParser(newickStr, false, true, false, 1);
            assertTrue(treesEquivalent(marginalTree.getRoot(), correctTree.getRoot(), 1e-15));
        }
    }
    
}
