package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.spec.domain.PositiveReal;
import beast.base.evolution.alignment.Taxon;
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.spec.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.evolution.tree.coalescent.ConstantPopulation;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.MarginalNode;
import contactrees.MarginalTree;
import contactrees.MarginalTreeSlow;
import contactrees.model.FreezableClock;
import contactrees.model.SimulatedACGWithBlocks;
import contactrees.model.likelihood.CTreeLikelihood;

/**
 * Tests of the creator-mapped marginal trees (true heights, stable node numbers, snapshots)
 * on random conversion graphs.
 */
public class MarginalTreeMappingTest extends ContactreesTest {

    static final String CF_NEWICK = "((((1:0.4,2:0.4)11:0.6,3:1.0)12:0.7,(4:0.8,5:0.8)13:0.9)14:0.8,"
            + "(((6:0.3,7:0.3)15:0.9,8:1.2)16:0.6,(9:0.5,10:0.5)17:1.3)18:1.3)19:0.0;";
    static final int N_RANDOM_BLOCKS = 6;

    ConversionGraph randomAcg;
    BlockSet randomBlockSet;

    /**
     * Simulate conversions and block moves on a fixed 10-taxon clonal frame.
     */
    void simulateRandomACG(long seed, double conversionRate, double moveProb) {
        Randomizer.setSeed(seed);

        ConversionGraph clonalFrame = new ConversionGraph();
        clonalFrame.initAndValidate();
        clonalFrame.fromExtendedNewick(CF_NEWICK);

        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));

        randomAcg = new ConversionGraph();
        ArrayList<Block> blocks = new ArrayList<>();
        for (int i = 0; i < N_RANDOM_BLOCKS; i++) {
            Block block = new Block();
            block.setID("block" + i);
            blocks.add(block);
        }
        randomBlockSet = new BlockSet();
        randomBlockSet.initByName("network", randomAcg, "block", blocks);

        SimulatedACGWithBlocks sim = new SimulatedACGWithBlocks();
        sim.initByName(
                "network", randomAcg,
                "blockSet", randomBlockSet,
                "conversionRate", conversionRate,
                "moveProb", moveProb,
                "populationModel", popFunc,
                "clonalFrame", clonalFrame);
        randomAcg = sim.getACG();
        randomBlockSet = sim.getBlockSet();
    }

    MarginalTree newMarginalTree(ConversionGraph acg, Block block) {
        MarginalTree tree = new MarginalTree();
        tree.initByName("network", acg, "block", block, "nodetype", MarginalNode.class.getName());
        tree.setID("mtree." + block.getID());
        return tree;
    }

    MarginalTreeSlow newReferenceTree(ConversionGraph acg, Block block) {
        MarginalTreeSlow tree = new MarginalTreeSlow();
        tree.initByName("network", acg, "block", block, "nodetype", contactrees.MarginalNodeSlow.class.getName());
        return tree;
    }

    /**
     * Without rolled-out rates, the marginal trees must equal the independent reference
     * implementation (clades and heights) and be internally consistent.
     */
    @Test
    public void testMatchesReferenceOnRandomACGs() {
        int nConvsSeen = 0;
        for (long seed = 1; seed <= 30; seed++) {
            simulateRandomACG(seed, 0.4, 0.5);
            nConvsSeen += randomAcg.getConvCount();
            for (Block block : randomBlockSet.getBlocks()) {
                MarginalTree tree = newMarginalTree(randomAcg, block);
                MarginalTreeSlow reference = newReferenceTree(randomAcg, block);

                tree.checkConsistency("test");
                assertTrue("seed " + seed + ", " + block.getID() + ": " + tree + " vs " + reference,
                        treesEquivalent(tree, reference, 1e-12));
            }
        }
        assertTrue("simulation produced no conversions", nConvsSeen > 0);
    }

    /**
     * With a strict clock (rate != 1) on the likelihood, the likelihood of the true-height
     * marginal tree must equal the likelihood of the rolled-out reference tree.
     */
    @Test
    public void testLikelihoodMatchesRolledOutReference() {
        double rate = 0.37;
        Alignment alignment = getAlignment(10);
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", jc);

        for (long seed = 1; seed <= 10; seed++) {
            simulateRandomACG(seed, 0.4, 0.5);
            for (Block block : randomBlockSet.getBlocks()) {
                MarginalTree tree = newMarginalTree(randomAcg, block);
                StrictClockModel clock = new StrictClockModel();
                clock.initByName("clock.rate", new RealScalarParam<>(rate, PositiveReal.INSTANCE));
                CTreeLikelihood likelihood = new CTreeLikelihood();
                likelihood.initByName("data", alignment, "tree", tree, "siteModel", siteModel,
                        "branchRateModel", clock);

                StrictClockModel refClock = new StrictClockModel();
                refClock.initByName("clock.rate", new RealScalarParam<>(rate, PositiveReal.INSTANCE));
                MarginalTreeSlow reference = new MarginalTreeSlow();
                reference.initByName("network", randomAcg, "block", block,
                        "nodetype", contactrees.MarginalNodeSlow.class.getName(), "branchRateModel", refClock);
                TreeLikelihood refLikelihood = new TreeLikelihood();
                refLikelihood.initByName("data", alignment, "tree", reference, "siteModel", siteModel);

                double logP = likelihood.calculateLogP();
                double refLogP = refLikelihood.calculateLogP();
                assertEquals("seed " + seed + ", " + block.getID(), refLogP, logP, 1e-9 * Math.abs(refLogP));
            }
        }
    }

    /**
     * A plain TreeLikelihood must not apply a branch-dependent clonal frame clock to a marginal
     * tree (it would look up rates by marginal node number). Strict clocks remain allowed.
     */
    @Test
    public void testPlainLikelihoodRejectsRelaxedClock() {
        Alignment alignment = getAlignment(10);
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", jc);
        simulateRandomACG(1, 0.4, 0.5);
        Block block = randomBlockSet.getBlocks().get(0);

        MarginalTree relaxedTree = newMarginalTree(randomAcg, block);
        TreeLikelihood relaxed = new TreeLikelihood();
        relaxed.initByName("data", alignment, "tree", relaxedTree, "siteModel", siteModel,
                "branchRateModel", randomRelaxedClock(randomAcg, 1));
        assertThrows(IllegalArgumentException.class, relaxedTree::requiresRecalculation);

        MarginalTree strictTree = newMarginalTree(randomAcg, block);
        StrictClockModel strict = new StrictClockModel();
        strict.initByName("clock.rate", new RealScalarParam<>(0.5, PositiveReal.INSTANCE));
        TreeLikelihood plainStrict = new TreeLikelihood();
        plainStrict.initByName("data", alignment, "tree", strictTree, "siteModel", siteModel, "branchRateModel", strict);
        strictTree.requiresRecalculation();

        MarginalTree wrappedTree = newMarginalTree(randomAcg, block);
        CTreeLikelihood wrapped = new CTreeLikelihood();
        wrapped.initByName("data", alignment, "tree", wrappedTree, "siteModel", siteModel,
                "branchRateModel", randomRelaxedClock(randomAcg, 1));
        wrappedTree.requiresRecalculation();
        assertTrue(Double.isFinite(wrapped.calculateLogP()));
    }

    /**
     * Random relaxed clock on the clonal frame: one rate per CF branch, drawn in [0.3, 3).
     */
    static UCRelaxedClockModel randomRelaxedClock(ConversionGraph acg, long seed) {
        java.util.Random random = new java.util.Random(seed);
        double[] rates = new double[acg.getNodeCount() - 1];
        for (int i = 0; i < rates.length; i++)
            rates[i] = 0.3 + 2.7 * random.nextDouble();
        UCRelaxedClockModel clock = new UCRelaxedClockModel();
        clock.initByName("rates", new RealVectorParam<>(rates, PositiveReal.INSTANCE),
                "distr", new Uniform(), "tree", acg);
        return clock;
    }

    /**
     * With a relaxed clock (optionally with a frozen taxon), the likelihood of the true-height
     * marginal tree with path-averaged rates must equal the likelihood of the rolled-out reference.
     */
    @Test
    public void testRelaxedClockLikelihoodMatchesRolledOutReference() {
        Alignment alignment = getAlignment(10);
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", jc);

        int nWithConversions = 0;
        for (long seed = 1; seed <= 10; seed++) {
            simulateRandomACG(seed, 0.4, 0.5);
            nWithConversions += randomAcg.getConvCount() > 0 ? 1 : 0;
            for (boolean frozen : new boolean[] {false, true}) {
                UCRelaxedClockModel relaxed = randomRelaxedClock(randomAcg, seed);
                FreezableClock clock = new FreezableClock();
                clock.initByName("clock", relaxed, "frozenTaxa", frozen ? List.of(new Taxon("4")) : List.of());

                for (Block block : randomBlockSet.getBlocks()) {
                    MarginalTree tree = newMarginalTree(randomAcg, block);
                    CTreeLikelihood likelihood = new CTreeLikelihood();
                    likelihood.initByName("data", alignment, "tree", tree, "siteModel", siteModel,
                            "branchRateModel", clock);
                    tree.checkConsistency("init");

                    MarginalTreeSlow reference = new MarginalTreeSlow();
                    reference.initByName("network", randomAcg, "block", block,
                            "nodetype", contactrees.MarginalNodeSlow.class.getName(), "branchRateModel", clock);
                    TreeLikelihood refLikelihood = new TreeLikelihood();
                    refLikelihood.initByName("data", alignment, "tree", reference, "siteModel", siteModel);

                    double logP = likelihood.calculateLogP();
                    double refLogP = refLikelihood.calculateLogP();
                    assertEquals("seed " + seed + ", frozen " + frozen + ", " + block.getID(),
                            refLogP, logP, 1e-9 * Math.abs(refLogP));
                }
            }
        }
        assertTrue(nWithConversions > 0);
    }

    /**
     * Flipping one borrowing must leave the node numbers of all creators that keep firing unchanged.
     */
    @Test
    public void testNodeNumbersStableUnderBorrowingFlip() {
        int nFlipsChecked = 0;
        for (long seed = 1; seed <= 20; seed++) {
            simulateRandomACG(seed, 0.4, 0.5);
            if (randomAcg.getConvCount() == 0)
                continue;

            for (Block block : randomBlockSet.getBlocks()) {
                MarginalTree tree = newMarginalTree(randomAcg, block);
                Map<Long, Integer> before = creatorNumbers(tree);

                Conversion conv = randomAcg.getConversions().getRandomConversion();
                if (block.isAffected(conv))
                    block.removeMove(conv);
                else
                    block.addMove(conv);
                tree.ensureUpToDate();
                tree.checkConsistency("after flip");

                Map<Long, Integer> after = creatorNumbers(tree);
                for (Map.Entry<Long, Integer> e : after.entrySet()) {
                    Integer oldNr = before.get(e.getKey());
                    if (oldNr != null)
                        assertEquals("creator kept firing but changed its node number", oldNr, e.getValue());
                }
                nFlipsChecked++;
            }
        }
        assertTrue(nFlipsChecked > 0);
    }

    static Map<Long, Integer> creatorNumbers(MarginalTree tree) {
        Map<Long, Integer> numbers = new HashMap<>();
        for (int nr = tree.getLeafNodeCount(); nr < tree.getNodeCount(); nr++)
            numbers.put(((long) tree.getCreatorKind(nr) << 32) | (tree.getCreatorKey(nr) & 0xffffffffL), nr);
        return numbers;
    }

    /**
     * After a change, restore() must reinstate the exact pre-change tree (structure, numbering
     * and mapping), as seen by the checksum.
     */
    @Test
    public void testRestoreReinstatesExactTree() {
        int nRestores = 0;
        for (long seed = 1; seed <= 20; seed++) {
            simulateRandomACG(seed, 0.4, 0.5);
            if (randomAcg.getConvCount() == 0)
                continue;

            for (Block block : randomBlockSet.getBlocks()) {
                MarginalTree tree = newMarginalTree(randomAcg, block);
                int checksumBefore = tree.getChecksum();
                String newickBefore = tree.getRoot().toNewick();

                List<Conversion> convs = new ArrayList<>(randomAcg.getConversions().getConversions());
                for (Conversion conv : convs) {
                    if (block.isAffected(conv))
                        block.removeMove(conv);
                    else
                        block.addMove(conv);
                }
                tree.ensureUpToDate();
                tree.checkConsistency("after flips");

                tree.restore();
                assertEquals(checksumBefore, tree.getChecksum());
                assertEquals(newickBefore, tree.getRoot().toNewick());
                assertTrue(!tree.isModifiedInProposal());

                // Undo the flips on the block; the tree then catches up again
                for (Conversion conv : convs) {
                    if (block.isAffected(conv))
                        block.removeMove(conv);
                    else
                        block.addMove(conv);
                }
                tree.ensureUpToDate();
                tree.checkConsistency("after undo");
                tree.restore();
                nRestores++;
            }
        }
        assertTrue(nRestores > 0);
    }

    /**
     * Height-only updates through the creator mapping must agree with a resweep.
     */
    @Test
    public void testCreatorHeightUpdate() {
        simulateRandomACG(3, 0.4, 0.5);
        Block block = randomBlockSet.getBlocks().get(0);
        MarginalTree tree = newMarginalTree(randomAcg, block);

        // Move the root up: never changes the order of events along any branch
        Node cfRoot = randomAcg.getRoot();
        cfRoot.setHeight(cfRoot.getHeight() + 0.5);
        tree.updateCreatorHeight(cfRoot);
        tree.markUpToDate();
        tree.checkConsistency("after height update");

        MarginalNode rootNode = tree.getNodeForCFNode(cfRoot);
        if (rootNode != null)
            assertEquals(cfRoot.getHeight(), rootNode.getHeight(), 0.0);
        else
            assertNotNull(tree.getRoot());
    }
}
