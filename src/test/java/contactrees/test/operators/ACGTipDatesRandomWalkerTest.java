package contactrees.test.operators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Conversion;
import contactrees.operators.ACGTipDatesRandomWalker;
import contactrees.test.ContactreesTest;

public class ACGTipDatesRandomWalkerTest extends ContactreesTest {

    ACGTipDatesRandomWalker getOperator() {
        ACGTipDatesRandomWalker op = new ACGTipDatesRandomWalker();
        op.initByName("tree", acg, "windowSize", 0.5, "weight", 1.0);
        return op;
    }

    @Test
    public void testConversionKeepsRelativeHeight() {
        ACGTipDatesRandomWalker op = getOperator();

        // conv1 (node2 -> node3) has valid range [0, 1] before and [0.4, 1] after the move
        double logHR = op.moveTip(node2, 0.4);
        assertEquals(0.4, node2.getHeight(), EPS);
        assertEquals(0.7, conv1.getHeight(), EPS);
        assertEquals(Math.log(0.6), logHR, EPS);

        // conv2 is not attached to node2
        assertEquals(1.5, conv2.getHeight(), EPS);
        assertFalse(acg.isInvalid());
    }

    @Test
    public void testLowerBoundFromOtherNode() {
        ACGTipDatesRandomWalker op = getOperator();

        // conv2 (node3 -> node4) is bounded below by node4 (height 1.0), so moving node3 to 0.3
        // only moves conv1 (valid range [0, 1] -> [0.3, 1])
        double logHR = op.moveTip(node3, 0.3);
        assertEquals(1.5, conv2.getHeight(), EPS);
        assertEquals(0.3 + 0.5 * 0.7, conv1.getHeight(), EPS);
        assertEquals(Math.log(0.7), logHR, EPS);
        assertFalse(acg.isInvalid());
    }

    @Test
    public void testRoundTrip() {
        ACGTipDatesRandomWalker op = getOperator();

        double logHR = op.moveTip(node2, 0.8);
        logHR += op.moveTip(node2, 0.0);
        assertEquals(0.0, node2.getHeight(), EPS);
        assertEquals(0.5, conv1.getHeight(), EPS);
        assertEquals(0.0, logHR, EPS);
    }

    @Test
    public void testRejectEmptyRange() {
        ACGTipDatesRandomWalker op = getOperator();

        // conv1 can not be above the parent of node2 (height 1.0), so node3 can not move above it
        double logHR = op.moveTip(node3, 1.2);
        assertEquals(Double.NEGATIVE_INFINITY, logHR, 0.0);
        assertEquals(0.0, node3.getHeight(), EPS);
        assertEquals(0.5, conv1.getHeight(), EPS);
        assertEquals(1.5, conv2.getHeight(), EPS);
    }

    @Test
    public void testProposalKeepsACGValid() {
        Randomizer.setSeed(1);
        for (int i = 0; i < 1000; i++) {
            setup();
            double logHR = getOperator().proposal();
            if (logHR > Double.NEGATIVE_INFINITY)
                assertFalse(acg.isInvalid());
        }
    }

    /**
     * MCMC on a two-taxon ACG with one conversion between the tips and a target that is uniform on
     * {0 < h < H_MAX, h < t < 1}, where h is the height of tip 1 and t the conversion height.
     * The tip operator keeps the relative position of t fixed, so the chain of h given that relative
     * position has density proportional to (1 - h). Without the Jacobian term, h would be uniform.
     */
    @Test
    public void testStationaryDistribution() {
        final double H_MAX = 0.8;
        final int N = 200000;

        setupFromNewick("(1:1.0,2:1.0)3:0.0;");
        Node tip = acg.getRoot().getLeft();
        Node otherTip = acg.getRoot().getRight();
        Conversion conv = new Conversion(tip, otherTip, 0.3, acg, 1);
        acg.addConversion(conv);
        ACGTipDatesRandomWalker op = getOperator();

        Randomizer.setSeed(1);
        double sumH = 0;
        for (int i = 0; i < N; i++) {
            double oldH = tip.getHeight();
            double oldT = conv.getHeight();
            double newH = oldH + Randomizer.nextDouble() * 0.6 - 0.3;

            boolean accept = false;
            if (newH > 0.0 && newH < H_MAX) {
                double logHR = op.moveTip(tip, newH);
                accept = Math.log(Randomizer.nextDouble()) < logHR;
                if (!accept) {
                    tip.setHeight(oldH);
                    conv.setHeight(oldT);
                }
            }
            assertTrue(conv.isValid());
            sumH += tip.getHeight();
        }

        // E[h] for density proportional to (1 - h) on [0, H_MAX]
        double expected = (H_MAX * H_MAX / 2 - H_MAX * H_MAX * H_MAX / 3) / (H_MAX - H_MAX * H_MAX / 2);
        assertEquals(expected, sumH / N, 0.01);
    }
}
