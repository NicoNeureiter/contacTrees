package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import contactrees.Block;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Tests for {@link MarginalTree}'s "lowest height changed this proposal" bookkeeping
 * ({@link MarginalTree#noteChangedAtHeight}, {@link MarginalTree#getLowestChangedHeight},
 * {@link MarginalTree#makeIneligibleForIncrementalUpdate}): plumbing for a future
 * height-localized incremental rebuild. Nothing reads it yet to change behaviour.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeChangedHeightTest extends ContactreesTest {

    private MarginalTree newMarginalTree(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName());
        return marginalTree;
    }

    @Test
    public void testNothingReportedIsPositiveInfinity() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        assertEquals(Double.POSITIVE_INFINITY, marginalTree.getLowestChangedHeight(), 0.0);
    }

    @Test
    public void testMinimumKeptAcrossSeveralNotifications() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        marginalTree.noteChangedAtHeight(2.0);
        assertEquals(2.0, marginalTree.getLowestChangedHeight(), EPS);

        // A higher report must not raise the running minimum.
        marginalTree.noteChangedAtHeight(3.0);
        assertEquals(2.0, marginalTree.getLowestChangedHeight(), EPS);

        // A lower report lowers it.
        marginalTree.noteChangedAtHeight(0.5);
        assertEquals(0.5, marginalTree.getLowestChangedHeight(), EPS);
    }

    @Test
    public void testCFHeightChangeMakesTreeIneligible() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);
        // initAndValidate() ends with recalculate(); makeOutdated() -- clear the forced
        // outdated=true baseline first, as MarginalTreeTest.testCheckOutdatedPrecision does,
        // so checkOutdated() below is exercising the CF-version branch, not the stale flag.
        marginalTree.recalculate();

        marginalTree.noteChangedAtHeight(0.5);
        assertEquals(0.5, marginalTree.getLowestChangedHeight(), EPS);

        // A CF change is not height-localizable: it must force ineligibility, and
        // checkOutdated() is where that gets applied. Node.setHeight()'s own startEditing() is
        // gated on the tree having a State attached (Node.startEditing(): "if (m_tree != null
        // && m_tree.getState() != null)"), which this standalone fixture doesn't have, so bump
        // cfVersion directly the way ConversionGraph.startEditing() would from a real CF edit.
        acg2.startEditing(null);
        assertTrue(marginalTree.checkOutdated());
        assertEquals(Double.POSITIVE_INFINITY, marginalTree.getLowestChangedHeight(), 0.0);
    }

    @Test
    public void testConversionHeightChangeMakesTreeIneligible() {
        Block block = blockSet2.getBlocks().get(0);
        block.addMove(conv2_1);
        MarginalTree marginalTree = newMarginalTree(block);
        marginalTree.recalculate();

        marginalTree.noteChangedAtHeight(0.5);
        assertEquals(0.5, marginalTree.getLowestChangedHeight(), EPS);

        // Re-attaching/moving one of THIS block's own conversions is not height-localizable
        // either, even though it only bumps convVersion (not cfVersion).
        conv2_1.setHeight(conv2_1.getHeight() + 0.01);
        assertTrue(marginalTree.checkOutdated());
        assertEquals(Double.POSITIVE_INFINITY, marginalTree.getLowestChangedHeight(), 0.0);
    }

    @Test
    public void testUnrelatedConversionHeightChangeLeavesTreeEligible() {
        Block block = blockSet2.getBlocks().get(0);
        // block does NOT carry conv2_2.
        MarginalTree marginalTree = newMarginalTree(block);
        marginalTree.recalculate();

        marginalTree.noteChangedAtHeight(0.5);

        conv2_2.setHeight(conv2_2.getHeight() + 0.01);
        assertEquals("Editing a conversion this block doesn't carry must not touch this tree",
                0.5, marginalTree.getLowestChangedHeight(), EPS);
    }

    @Test
    public void testMakeOutdatedForcesIneligible() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        marginalTree.noteChangedAtHeight(0.5);
        marginalTree.makeOutdated();
        assertEquals(Double.POSITIVE_INFINITY, marginalTree.getLowestChangedHeight(), 0.0);
    }

    @Test
    public void testClearedAfterAccept() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        marginalTree.noteChangedAtHeight(0.5);
        assertEquals(0.5, marginalTree.getLowestChangedHeight(), EPS);

        marginalTree.accept();
        assertEquals(Double.POSITIVE_INFINITY, marginalTree.getLowestChangedHeight(), 0.0);
    }

    @Test
    public void testClearedAfterRestore() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);
        // initAndValidate() already ran recalculate() once, so this tree is mid-proposal
        // (snapshotTakenThisProposal == true) and restore() takes the full-restore path, not
        // the "untouched this proposal" fast path.
        marginalTree.noteChangedAtHeight(0.5);
        assertEquals(0.5, marginalTree.getLowestChangedHeight(), EPS);

        marginalTree.restore();
        assertEquals(Double.POSITIVE_INFINITY, marginalTree.getLowestChangedHeight(), 0.0);
    }

}
