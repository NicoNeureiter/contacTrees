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

}
