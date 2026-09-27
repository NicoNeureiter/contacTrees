package contactrees.test;

import static org.junit.Assert.assertThrows;

import java.util.List;

import org.junit.Test;

import contactrees.Block;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Unit tests for {@link MarginalTree#assertConsistent(String)}.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeValidatorTest extends ContactreesTest {

    @Test
    public void testAssertConsistentPassesOnFreshTree() throws Exception {
        List<Block> blocks = blockSet.getBlocks();
        blocks.get(1).addMove(conv1);

        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg, "block", blocks.get(1), "nodetype", MarginalNode.class.getName());

        marginalTree.assertConsistent("test");
    }

    @Test
    public void testAssertConsistentThrowsOnStaleTree() throws Exception {
        List<Block> blocks = blockSet.getBlocks();

        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg, "block", blocks.get(0), "nodetype", MarginalNode.class.getName());

        // Mutate the block's conversion set without going through the dirty-tracking path
        // that would normally trigger a rebuild. blocks.get(1) with conv1 added produces a
        // different topology in MarginalTreeTest.testNonOverlapping, so this is a genuine
        // topology change, not a no-op.
        blocks.get(0).addMove(conv1);

        assertThrows(IllegalStateException.class, () -> marginalTree.assertConsistent("test"));
    }

}
