package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.BitSet;

import org.junit.Test;

import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Tests for {@link MarginalTree#computeAffectedCFBranches}, the conservative
 * over-estimate of which CF branches a borrowing-flip toggle can perturb.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeAffectedBranchesTest extends ContactreesTest {

    private MarginalTree newMarginalTree(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName());
        return marginalTree;
    }

    @Test
    public void testSiblingLeavesGiveSmallSpan() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        Conversion conv = new Conversion(node2_1, node2_2, 0.5, acg2, 101);
        acg2.addConversion(conv);
        block.addMove(conv);

        BitSet affected = marginalTree.computeAffectedCFBranches(conv, marginalTree.getBlockConversions());

        BitSet expected = new BitSet();
        expected.set(node2_1.getNr());
        expected.set(node2_2.getNr());
        assertEquals(expected, affected);
    }

    @Test
    public void testClosureGrowsThroughHigherConversion() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        Conversion base = new Conversion(node2_1, node2_2, 0.5, acg2, 101);
        Conversion higher = new Conversion(node2_1, node2_4, 0.9, acg2, 102);
        acg2.addConversion(base);
        acg2.addConversion(higher);
        block.addMove(base);
        block.addMove(higher);

        BitSet affected = marginalTree.computeAffectedCFBranches(base, marginalTree.getBlockConversions());

        // Seed span from `base` alone is just {leaf1, leaf2}; `higher` shares leaf1 and is
        // above base, so the closure must pull in its own span too.
        assertTrue(affected.get(node2_1.getNr()));
        assertTrue(affected.get(node2_2.getNr()));
        assertTrue(affected.get(node2_6.getNr()));
        assertTrue(affected.get(node2_8.getNr()));
        assertTrue(affected.get(node2_4.getNr()));
        assertTrue(affected.get(node2_7.getNr()));
        // The CF root is the MRCA of `higher`'s endpoints, so it stays excluded.
        assertFalse(affected.get(root2.getNr()));
    }

    @Test
    public void testHeightRestrictionExcludesLowerConversions() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        Conversion base = new Conversion(node2_1, node2_2, 0.5, acg2, 101);
        // Same endpoint (leaf1) as the closure test above, but BELOW `base`'s height: must be
        // ignored even though it would otherwise pull leaf4/node7/node8/root into the span.
        Conversion lower = new Conversion(node2_1, node2_4, 0.3, acg2, 102);
        acg2.addConversion(base);
        acg2.addConversion(lower);
        block.addMove(base);
        block.addMove(lower);

        BitSet affected = marginalTree.computeAffectedCFBranches(base, marginalTree.getBlockConversions());

        BitSet expected = new BitSet();
        expected.set(node2_1.getNr());
        expected.set(node2_2.getNr());
        assertEquals(expected, affected);
    }

}
