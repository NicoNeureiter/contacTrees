package contactrees.test.model;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.junit.Test;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.DirectSimulator;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.speciation.YuleModel;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.model.ACGDistribution;
import contactrees.model.ConversionMovePrior;
import contactrees.model.ConversionPrior;
import contactrees.test.ContactreesTest;

/**
 * Direct simulation of the network (clonal frame + conversions) and the conversion moves
 * must not depend on the order in which the distributions appear in the CompoundDistribution.
 */
public class NetworkSamplingTest extends ContactreesTest {

    static final int N_SAMPLES = 3000;
    static final double CONV_RATE = 0.3;
    static final double P_MOVE = 0.4;

    /** Factorized model: tree prior, conversion prior and move prior side by side. */
    @Test
    public void factorizedPriorIsOrderIndependent() {
        for (List<Integer> order : permutations(3)) {
            Model m = new Model(false);
            List<Distribution> distrs = Arrays.asList(m.yule, m.networkPrior, m.movePrior);
            simulateAndCheck(m, reorder(distrs, order), "order " + order);
        }
    }

    /** Joint model: ACGDistribution wrapping the tree prior, next to the move prior. */
    @Test
    public void jointPriorIsOrderIndependent() {
        for (List<Integer> order : permutations(2)) {
            Model m = new Model(true);
            List<Distribution> distrs = Arrays.asList(m.networkPrior, m.movePrior);
            simulateAndCheck(m, reorder(distrs, order), "order " + order);
        }
    }

    /** A distribution listing the network as a condition makes sampling order dependent. */
    @Test
    public void networkConditionIsRejected() {
        Model m = new Model(false);
        NetworkConditioned other = new NetworkConditioned();
        other.initByName("x", m.acg);
        assertThrows(IllegalStateException.class, () -> sample(m, Arrays.asList(m.yule, m.networkPrior, other)));
    }

    /** Same, but the network is reached indirectly through a CalculationNode (the BlockSet). */
    @Test
    public void indirectNetworkConditionIsRejected() {
        Model m = new Model(false);
        NetworkConditioned other = new NetworkConditioned();
        other.initByName("x", m.blockSet);
        assertThrows(IllegalStateException.class, () -> sample(m, Arrays.asList(m.yule, m.networkPrior, other)));
    }

    /* Helpers */

    class Model {
        ConversionGraph acg;
        BlockSet blockSet;
        YuleModel yule;
        Distribution networkPrior;
        ConversionMovePrior movePrior;

        Model(boolean joint) {
            acg = getACGFromNewick("((1:1.0,2:1.0)5:1.0,(3:0.5,4:0.5)6:1.5)7:0.0;");
            acg.setID("acg");
            blockSet = getBlockSet(N_BLOCKS, acg);
            blockSet.setID("blockSet");

            RealScalarParam<PositiveReal> birthRate = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
            birthRate.setID("birthRate");
            RealScalarParam<NonNegativeReal> convRate = new RealScalarParam<>(CONV_RATE, NonNegativeReal.INSTANCE);
            convRate.setID("conversionRate");
            RealScalarParam<UnitInterval> pMove = new RealScalarParam<>(P_MOVE, UnitInterval.INSTANCE);
            pMove.setID("pMove");

            yule = new YuleModel();
            yule.initByName("tree", acg, "birthDiffRate", birthRate, "conditionalOnRoot", false);
            yule.setID("yule");

            if (joint) {
                networkPrior = new ACGDistribution();
                networkPrior.initByName("network", acg, "cfModel", yule, "conversionRate", convRate);
            } else {
                networkPrior = new ConversionPrior();
                networkPrior.initByName("network", acg, "conversionRate", convRate);
            }
            networkPrior.setID("networkPrior");

            movePrior = new ConversionMovePrior();
            movePrior.initByName("network", acg, "blockSet", blockSet, "pMove", pMove);
            movePrior.setID("movePrior");
        }
    }

    void sample(Model m, List<Distribution> distrs) {
        Randomizer.setSeed(1);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", new ArrayList<>(distrs));
        new DirectSimulator().clearSampledFlags(prior);
        prior.sample(null, new Random(1));
    }

    void simulateAndCheck(Model m, List<Distribution> distrs, String label) {
        Randomizer.setSeed(1);
        Random random = new Random(1);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", new ArrayList<>(distrs));
        DirectSimulator sim = new DirectSimulator();

        double prevRootHeight = Double.NaN;
        double sumDiff = 0, sumSqDiff = 0;  // conversion count minus its expectation given the tree
        long nConvTotal = 0, nMovesTotal = 0;

        for (int i = 0; i < N_SAMPLES; i++) {
            sim.clearSampledFlags(prior);
            prior.sample(null, random);
            ConversionGraph acg = m.acg;

            // Clonal frame was resampled
            double rootHeight = acg.getRoot().getHeight();
            assertNotEquals(label, prevRootHeight, rootHeight, 0.0);
            prevRootHeight = rootHeight;

            // Conversions are attached to the current clonal frame
            assertFalse(label, acg.isInvalid());
            for (Conversion conv : acg.getConversions()) {
                assertTrue(label, conv.isValid());
                assertSame(label, acg.getNode(conv.getNode1().getNr()), conv.getNode1());
                assertSame(label, acg.getNode(conv.getNode2().getNr()), conv.getNode2());
            }

            // Moves refer to current conversions only
            for (Block block : m.blockSet) {
                for (int convID : block.getConversionIDs())
                    assertNotNull(label, acg.getConversions().get(convID));
                nMovesTotal += block.countMoves();
            }

            int nConv = acg.getConvCount();
            double diff = nConv - CONV_RATE * acg.getClonalFramePairedLength();
            sumDiff += diff;
            sumSqDiff += diff * diff;
            nConvTotal += nConv;
        }

        // Mean conversion count matches the Poisson expectation given the tree
        double meanDiff = sumDiff / N_SAMPLES;
        double seDiff = Math.sqrt(sumSqDiff / N_SAMPLES / N_SAMPLES);
        assertEquals(label, 0.0, meanDiff, 4 * seDiff);

        // Move frequency matches pMove
        long nTrials = nConvTotal * N_BLOCKS;
        double seMoves = Math.sqrt(P_MOVE * (1 - P_MOVE) / nTrials);
        assertEquals(label, P_MOVE, (double) nMovesTotal / nTrials, 4 * seMoves);
    }

    static <T> List<T> reorder(List<T> items, List<Integer> order) {
        List<T> result = new ArrayList<>();
        for (int i : order)
            result.add(items.get(i));
        return result;
    }

    static List<List<Integer>> permutations(int n) {
        List<List<Integer>> result = new ArrayList<>();
        permute(new ArrayList<>(), n, result);
        return result;
    }

    static void permute(List<Integer> prefix, int n, List<List<Integer>> result) {
        if (prefix.size() == n) {
            result.add(new ArrayList<>(prefix));
            return;
        }
        for (int i = 0; i < n; i++) {
            if (prefix.contains(i))
                continue;
            prefix.add(i);
            permute(prefix, n, result);
            prefix.remove(prefix.size() - 1);
        }
    }

    /** Dummy distribution that is conditioned on some object x. */
    public static class NetworkConditioned extends Distribution {
        public Input<BEASTInterface> xInput = new Input<>("x", "Object to condition on.", Input.Validate.REQUIRED);

        @Override
        public void initAndValidate() {}

        @Override
        public double calculateLogP() {
            return 0.0;
        }

        @Override
        public List<String> getArguments() {
            return new ArrayList<>();
        }

        @Override
        public List<String> getConditions() {
            return Arrays.asList(xInput.get().getID());
        }

        @Override
        public void sample(State state, Random random) {
            if (sampledFlag)
                return;
            sampledFlag = true;
            sampleConditions(state, random);
        }
    }
}
