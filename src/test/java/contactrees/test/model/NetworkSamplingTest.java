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
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.DirectSimulator;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.speciation.YuleModel;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;
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
            Model m = new Model();
            List<Distribution> distrs = Arrays.asList(m.yule, m.networkPrior, m.movePrior);
            simulateAndCheck(m, reorder(distrs, order), "order " + order);
        }
    }

    /** A distribution listing the network as a condition makes sampling order dependent. */
    @Test
    public void networkConditionIsRejected() {
        Model m = new Model();
        NetworkConditioned other = new NetworkConditioned();
        other.initByName("x", m.acg);
        assertThrows(IllegalStateException.class, () -> sample(m, Arrays.asList(m.yule, m.networkPrior, other)));
    }

    /** Same, but the network is reached indirectly through a CalculationNode (the BlockSet). */
    @Test
    public void indirectNetworkConditionIsRejected() {
        Model m = new Model();
        NetworkConditioned other = new NetworkConditioned();
        other.initByName("x", m.blockSet);
        assertThrows(IllegalStateException.class, () -> sample(m, Arrays.asList(m.yule, m.networkPrior, other)));
    }

    /** A tree prior depending on the network only indirectly (through TreeIntervals) needs clonalFramePrior. */
    @Test
    public void explicitIndirectClonalFramePrior() {
        for (List<Integer> order : permutations(3)) {
            Model m = new Model(Option.INDIRECT_EXPLICIT);
            List<Distribution> distrs = Arrays.asList(m.indirectPrior, m.networkPrior, m.movePrior);
            simulateAndCheck(m, reorder(distrs, order), "order " + order);
        }
    }

    /** Without clonalFramePrior the indirect tree prior is not found, and the estimated network can not be simulated. */
    @Test
    public void indirectClonalFramePriorIsNotAutoDetected() {
        Model m = new Model(Option.INDIRECT_IMPLICIT);
        assertThrows(IllegalStateException.class,
                () -> sample(m, Arrays.asList(m.indirectPrior, m.networkPrior, m.movePrior)));
    }

    /** Network with estimate=false and no tree prior: the clonal frame stays fixed. */
    @Test
    public void fixedClonalFrame() {
        for (List<Integer> order : permutations(2)) {
            Model m = new Model(Option.FIXED);
            List<Distribution> distrs = Arrays.asList(m.networkPrior, m.movePrior);
            simulateAndCheck(m, reorder(distrs, order), "order " + order, true);
        }
    }

    /** No tree prior and an estimated network is an error. */
    @Test
    public void missingClonalFramePriorIsRejected() {
        Model m = new Model(Option.NO_TREE_PRIOR);
        assertThrows(IllegalStateException.class, () -> sample(m, Arrays.asList(m.networkPrior, m.movePrior)));
    }

    /** clonalFramePrior must not be a conversion prior. */
    @Test
    public void conversionPriorAsClonalFramePriorIsRejected() {
        Model m = new Model();
        ConversionPrior other = new ConversionPrior();
        RealVectorParam<NonNegativeReal> rates = new RealVectorParam<>(new double[] {CONV_RATE}, NonNegativeReal.INSTANCE);
        // initByName wraps exceptions from initAndValidate() in a RuntimeException
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> other.initByName("network", m.acg, "conversionRates", rates, "clonalFramePrior", m.networkPrior));
        assertTrue(e.getMessage(), e.getMessage().contains("clonalFramePrior"));
    }

    /* Helpers */

    enum Option { DEFAULT, INDIRECT_EXPLICIT, INDIRECT_IMPLICIT, FIXED, NO_TREE_PRIOR }

    class Model {
        ConversionGraph acg;
        BlockSet blockSet;
        YuleModel yule;
        Distribution networkPrior;
        ConversionMovePrior movePrior;
        IndirectTreePrior indirectPrior;

        Model() {
            this(Option.DEFAULT);
        }

        Model(Option option) {
            acg = getACGFromNewick("((1:1.0,2:1.0)5:1.0,(3:0.5,4:0.5)6:1.5)7:0.0;");
            acg.setID("acg");
            blockSet = getBlockSet(N_BLOCKS, acg);
            blockSet.setID("blockSet");

            RealScalarParam<PositiveReal> birthRate = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
            birthRate.setID("birthRate");
            RealVectorParam<NonNegativeReal> convRates = new RealVectorParam<>(new double[] {CONV_RATE}, NonNegativeReal.INSTANCE);
            convRates.setID("conversionRates");
            RealScalarParam<UnitInterval> pMove = new RealScalarParam<>(P_MOVE, UnitInterval.INSTANCE);
            pMove.setID("pMove");

            if (option == Option.DEFAULT) {
                yule = new YuleModel();
                yule.initByName("tree", acg, "birthDiffRate", birthRate, "conditionalOnRoot", false);
                yule.setID("yule");
            }
            if (option == Option.FIXED)
                acg.isEstimatedInput.setValue(false, acg);
            if (option == Option.INDIRECT_EXPLICIT || option == Option.INDIRECT_IMPLICIT) {
                TreeIntervals intervals = new TreeIntervals();
                intervals.initByName("tree", acg);
                intervals.setID("intervals");
                indirectPrior = new IndirectTreePrior();
                indirectPrior.initByName("intervals", intervals);
                indirectPrior.setID("indirectPrior");
            }

            networkPrior = new ConversionPrior();
            if (option == Option.INDIRECT_EXPLICIT)
                networkPrior.initByName("network", acg, "conversionRates", convRates,
                        "clonalFramePrior", indirectPrior);
            else
                networkPrior.initByName("network", acg, "conversionRates", convRates);
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
        simulateAndCheck(m, distrs, label, false);
    }

    void simulateAndCheck(Model m, List<Distribution> distrs, String label, boolean fixedClonalFrame) {
        Randomizer.setSeed(1);
        Random random = new Random(1);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", new ArrayList<>(distrs));
        DirectSimulator sim = new DirectSimulator();

        double prevRootHeight = Double.NaN;
        double initialRootHeight = m.acg.getRoot().getHeight();
        double sumDiff = 0, sumSqDiff = 0;  // conversion count minus its expectation given the tree
        long nConvTotal = 0, nMovesTotal = 0;

        for (int i = 0; i < N_SAMPLES; i++) {
            sim.clearSampledFlags(prior);
            prior.sample(null, random);
            ConversionGraph acg = m.acg;

            double rootHeight = acg.getRoot().getHeight();
            if (fixedClonalFrame)
                assertEquals(label, initialRootHeight, rootHeight, 0.0);
            else  // Clonal frame was resampled
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

    /**
     * Tree prior which depends on the network only through TreeIntervals, like a coalescent.
     * Sampling rescales the internal node heights by a random factor.
     */
    public static class IndirectTreePrior extends Distribution {
        public Input<TreeIntervals> intervalsInput = new Input<>("intervals", "Intervals of the tree.", Input.Validate.REQUIRED);

        double[] originalHeights;

        @Override
        public void initAndValidate() {}

        @Override
        public double calculateLogP() {
            return 0.0;
        }

        @Override
        public List<String> getArguments() {
            return Arrays.asList(intervalsInput.get().getID());
        }

        @Override
        public List<String> getConditions() {
            return new ArrayList<>();
        }

        @Override
        public void sample(State state, Random random) {
            if (sampledFlag)
                return;
            sampledFlag = true;

            beast.base.evolution.tree.Tree tree = intervalsInput.get().treeInput.get();
            Node[] nodes = tree.getNodesAsArray();
            if (originalHeights == null) {
                originalHeights = new double[nodes.length];
                for (Node node : nodes)
                    originalHeights[node.getNr()] = node.getHeight();
            }
            double factor = 0.5 + 1.5 * random.nextDouble();
            for (Node node : nodes)
                if (!node.isLeaf())
                    node.setHeight(originalHeights[node.getNr()] * factor);
        }
    }
}
