package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import org.junit.Before;
import org.junit.Test;

import contactrees.Block;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Tests for the Stage C replay-and-compare validator (plan file
 * {@code in-coalre-we-recently-enchanted-feather.md}, "Stage C in detail"): with
 * {@code validationLevel=3}, every recalculate() whose {@link MarginalTree#getLowestChangedHeight()}
 * is finite independently re-derives the occupancy transitions from that height upward -- seeded
 * from the PREVIOUS rebuild's occupancy history -- and throws if the result disagrees with the
 * real rebuild's ground truth.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeReplayTest extends ContactreesTest {

    @Before
    public void resetReplayCounter() {
        MarginalTree.resetCounters();
    }

    private MarginalTree newValidatingMarginalTree(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName(),
                "validationLevel", 3);
        return marginalTree;
    }

    /**
     * Adding conv2_1 (node2_2 -> node2_3, height 0.5) to an empty block produces BOTH transition
     * kinds from the plan's table in a single rebuild: a coalescence at the conversion itself
     * (node2_2 and node2_3 both active leaves), and -- as a direct consequence -- a "lineage
     * moves" pass-through at CF node 6 immediately above it (node2_2's own leaf slot is vacated,
     * so only node2_1 is active there). A correctly reported hStart must validate cleanly against
     * both.
     */
    @Test
    public void testReplayValidatesCoalescenceAndMoveInOneFlip() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newValidatingMarginalTree(block);

        block.addMove(conv2_1);
        marginalTree.noteChangedAtHeight(conv2_1.getHeight());
        long before = MarginalTree.nReplayValidated;
        marginalTree.recalculate();

        assertEquals("An eligible rebuild must actually invoke the validator, not silently skip it",
                before + 1, MarginalTree.nReplayValidated);
    }

    /**
     * Flip, then flip back (a Gibbs-style revert): the second recalculate()'s hStart is the same
     * conversion height, and it must validate cleanly against the first rebuild's occupancy
     * history (recovering the original, pre-flip identities).
     */
    @Test
    public void testReplayValidatesRevert() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newValidatingMarginalTree(block);

        block.addMove(conv2_1);
        marginalTree.noteChangedAtHeight(conv2_1.getHeight());
        marginalTree.recalculate();

        block.removeMove(conv2_1);
        marginalTree.noteChangedAtHeight(conv2_1.getHeight());
        long before = MarginalTree.nReplayValidated;
        marginalTree.recalculate();

        assertEquals(before + 1, MarginalTree.nReplayValidated);
    }

    /**
     * A CF-height change makes the tree ineligible (see MarginalTreeChangedHeightTest); the
     * validator must not run at all in that case (there is no valid hStart to seed from).
     */
    @Test
    public void testReplaySkippedWhenIneligible() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newValidatingMarginalTree(block);
        marginalTree.recalculate();

        marginalTree.noteChangedAtHeight(0.5);
        acg2.startEditing(null);
        marginalTree.checkOutdated();

        long before = MarginalTree.nReplayValidated;
        marginalTree.recalculate();

        assertEquals("An ineligible rebuild must not run the validator",
                before, MarginalTree.nReplayValidated);
    }

    /**
     * Sanity check that the validator is not vacuously passing: if the reported hStart is wrong
     * (too high -- above the conversion's real height, but still within the region whose identity
     * the flip actually changed), the replay seeds from stale pre-flip history and must disagree
     * with the real rebuild's ground truth.
     *
     * conv2_1 (node2_2 -> node2_3, height 0.5) vacates node2_2's own leaf slot from height 0.5
     * onward, turning CF node 6's coalescence (height 1.0) into a "lineage moves" pass-through
     * (only node2_1 active). Reporting hStart=0.8 instead of the true 0.5 makes the replay seed
     * node2_2's branch from the PRE-flip history (still "active", own leaf identity) at the point
     * CF6 is replayed, disagreeing with the real (post-flip) ground truth there.
     */
    @Test
    public void testReplayCatchesWrongHStart() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newValidatingMarginalTree(block);
        marginalTree.recalculate();

        block.addMove(conv2_1);
        marginalTree.noteChangedAtHeight(0.8);

        assertThrows(IllegalStateException.class, marginalTree::recalculate);
    }

}
