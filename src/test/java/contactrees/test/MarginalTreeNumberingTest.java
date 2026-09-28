package contactrees.test;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Regression test for the Stage 3 stable marginal-node numbering: unlike the height/event-order
 * numbering it replaces, the marginal Nr assigned to a given CF node or conversion must stay the
 * same across rebuilds for as long as that creator keeps producing a marginal coalescence.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeNumberingTest extends ContactreesTest {

    @Test
    public void testStableNumberingAcrossRandomToggles() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion[] convs = {conv2_1, conv2_2, conv2_3};

        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName());

        Randomizer.setSeed(7);
        for (int step = 0; step < 200; step++) {
            // Snapshot the current CF-node / conversion -> marginal Nr mapping for every
            // creator that currently has a marginal coalescence.
            Map<Integer, Integer> cfNrBefore = new HashMap<>();
            for (int cfNr = 0; cfNr < acg2.getNodeCount(); cfNr++) {
                int nr = marginalTree.getMarginalNr(cfNr);
                if (nr >= 0)
                    cfNrBefore.put(cfNr, nr);
            }
            Map<Conversion, Integer> convNrBefore = new HashMap<>();
            for (Conversion conv : convs) {
                int nr = marginalTree.getMarginalNr(conv);
                if (nr >= 0)
                    convNrBefore.put(conv, nr);
            }

            // Toggle a random conversion's membership in the block.
            Conversion conv = convs[Randomizer.nextInt(convs.length)];
            if (block.isAffected(conv)) {
                block.removeMove(conv);
            } else {
                block.addMove(conv);
            }
            marginalTree.recalculate();

            // Any creator that had AND still has a marginal coalescence must keep its Nr.
            for (Map.Entry<Integer, Integer> e : cfNrBefore.entrySet()) {
                int nrAfter = marginalTree.getMarginalNr(e.getKey());
                if (nrAfter >= 0) {
                    assertEquals("CF node " + e.getKey() + " changed marginal Nr across a rebuild"
                            + " that did not remove its coalescence (step " + step + ")",
                            (int) e.getValue(), nrAfter);
                }
            }
            for (Map.Entry<Conversion, Integer> e : convNrBefore.entrySet()) {
                int nrAfter = marginalTree.getMarginalNr(e.getKey());
                if (nrAfter >= 0) {
                    assertEquals("Conversion " + e.getKey().getID() + " changed marginal Nr across"
                            + " a rebuild that did not remove its coalescence (step " + step + ")",
                            (int) e.getValue(), nrAfter);
                }
            }

            marginalTree.assertMappingConsistent("numbering-test-step-" + step);
        }
    }
}
