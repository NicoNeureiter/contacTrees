package contactrees.test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import contactrees.Block;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Tests for {@link MarginalTree}'s per-CF-node occupancy history and its query
 * ({@link MarginalTree#getOccupantAtHeight}).
 *
 * @author Nico Neureiter
 */
public class MarginalTreeOccupancyHistoryTest extends ContactreesTest {

    private MarginalTree newMarginalTree(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName());
        return marginalTree;
    }

    @Test
    public void testLeafBranchOccupiedFromOwnHeight() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        int leafNr = node2_1.getNr();
        assertEquals(leafNr, marginalTree.getOccupantAtHeight(leafNr, node2_1.getHeight()));
        assertEquals(leafNr, marginalTree.getOccupantAtHeight(leafNr, 0.7));
        // Below the leaf's own height: no occupant yet.
        assertEquals(-1, marginalTree.getOccupantAtHeight(leafNr, node2_1.getHeight() - 0.1));
    }

    @Test
    public void testConversionVacatesOldOccupantAndInstallsNew() {
        Block block = blockSet2.getBlocks().get(0);
        block.addMove(conv2_3);
        MarginalTree marginalTree = newMarginalTree(block);

        int leaf3Nr = node2_3.getNr();
        int leaf4Nr = node2_4.getNr();
        int convNodeNr = marginalTree.getMarginalNr(conv2_3);
        double convHeight = conv2_3.getHeight();

        // leaf4's branch: leaf4's own lineage below the conversion, the new coalescence
        // (leaf3 + leaf4) at and above it.
        assertEquals(leaf4Nr, marginalTree.getOccupantAtHeight(leaf4Nr, convHeight - 0.3));
        assertEquals(convNodeNr, marginalTree.getOccupantAtHeight(leaf4Nr, convHeight));
        assertEquals(convNodeNr, marginalTree.getOccupantAtHeight(leaf4Nr, convHeight + 0.1));

        // leaf3's branch: leaf3's own lineage below the conversion, permanently vacated above
        // it (leaf3 is absorbed into the conversion node, which lives on leaf4's branch).
        assertEquals(leaf3Nr, marginalTree.getOccupantAtHeight(leaf3Nr, convHeight - 0.3));
        assertEquals(-1, marginalTree.getOccupantAtHeight(leaf3Nr, convHeight));
    }

    @Test
    public void testHistoryResetNotAccumulatedAcrossRecalculates() {
        Block block = blockSet2.getBlocks().get(0);
        block.addMove(conv2_3);
        MarginalTree marginalTree = newMarginalTree(block);

        int leaf3Nr = node2_3.getNr();
        int leaf4Nr = node2_4.getNr();
        int lenLeaf3First = marginalTree.getOccupancyHistoryLength(leaf3Nr);
        int lenLeaf4First = marginalTree.getOccupancyHistoryLength(leaf4Nr);

        marginalTree.recalculate();
        assertEquals(lenLeaf3First, marginalTree.getOccupancyHistoryLength(leaf3Nr));
        assertEquals(lenLeaf4First, marginalTree.getOccupancyHistoryLength(leaf4Nr));

        marginalTree.recalculate();
        marginalTree.recalculate();
        assertEquals(lenLeaf3First, marginalTree.getOccupancyHistoryLength(leaf3Nr));
        assertEquals(lenLeaf4First, marginalTree.getOccupancyHistoryLength(leaf4Nr));

        // The query answers are still correct after repeated rebuilds, not just the lengths.
        int convNodeNr = marginalTree.getMarginalNr(conv2_3);
        double convHeight = conv2_3.getHeight();
        assertEquals(convNodeNr, marginalTree.getOccupantAtHeight(leaf4Nr, convHeight));
        assertEquals(-1, marginalTree.getOccupantAtHeight(leaf3Nr, convHeight));
    }

    /**
     * MarginalTree.restore() reverts node fields, cfBranchOccupant, and the stable-numbering
     * mapping arrays -- but the per-CF-node occupancy history log (occHeight/occNr/occPrev/
     * occLastOfNode) has no corresponding stored-field/restore support at all, so a rejected
     * proposal's recalculate() permanently overwrites it with a walk that never actually
     * happened. Caught by Stage C's replay validator empirically (bantu-contactrees.xml,
     * seed 12345): after a reject, a CF node's occupancy history still showed the REJECTED
     * proposal's vacate/install events instead of reverting to what was true before it.
     */
    @Test
    public void testOccupancyHistoryRestoredAfterReject() {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newMarginalTree(block);

        int leaf3Nr = node2_3.getNr();
        double convHeight = conv2_3.getHeight();

        marginalTree.recalculate();
        marginalTree.accept();

        int occupantBefore = marginalTree.getOccupantAtHeight(leaf3Nr, convHeight);
        int historyLenBefore = marginalTree.getOccupancyHistoryLength(leaf3Nr);

        // Proposal: add conv2_3 to the block (vacates leaf3's branch at convHeight), then reject.
        acg2.store();
        block.addMove(conv2_3);
        marginalTree.recalculate();
        assertEquals(-1, marginalTree.getOccupantAtHeight(leaf3Nr, convHeight));

        acg2.restore();
        block.restore();
        marginalTree.restore();

        assertEquals("occupancy history must revert to the pre-proposal state after a reject, "
                + "not keep reflecting the rejected proposal's walk",
                occupantBefore, marginalTree.getOccupantAtHeight(leaf3Nr, convHeight));
        assertEquals("occupancy history length must also revert after a reject",
                historyLenBefore, marginalTree.getOccupancyHistoryLength(leaf3Nr));
    }

}
