package contactrees.test;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.spec.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.speciation.YuleModel;
import beast.base.spec.evolution.substitutionmodel.Frequencies;
import beast.base.spec.evolution.substitutionmodel.GTR;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.MarginalNode;
import contactrees.MarginalTree;
import contactrees.model.ConversionMovePrior;
import contactrees.model.FreezableClock;
import contactrees.model.ConversionPrior;
import contactrees.model.likelihood.CTreeLikelihood;
import contactrees.model.likelihood.SharedTransitionMatrixCache;
import contactrees.operators.ACGScaler;
import contactrees.operators.ACGTipDatesRandomWalker;
import contactrees.operators.AddRemoveConversion;
import contactrees.operators.AddRemoveConversionGibbs;
import contactrees.operators.CFConversionSwap;
import contactrees.operators.CFSubtreeExchange;
import contactrees.operators.CFSubtreeSlide;
import contactrees.operators.CFUniform;
import contactrees.operators.CFWilsonBalding;
import contactrees.operators.ConversionSplit;
import contactrees.operators.ConvertedEdgeFlip;
import contactrees.operators.ConvertedEdgeHop;
import contactrees.operators.ConvertedEdgeHopGibbs;
import contactrees.operators.ConvertedEdgeSlide;
import contactrees.operators.GibbsSampleMovesPerConversion;
import contactrees.operators.ResampleBorrowings;

/**
 * Runs a short MCMC with every ACG operator (including the Gibbs operators, which evaluate
 * likelihoods inside their proposals) and checks after every step that
 *  - every marginal tree matches a from-scratch sweep and has no pending proposal state, and
 *  - every cached tree likelihood equals a from-scratch recomputation.
 */
public class MarginalTreeCacheConsistencyTest extends MarginalTreeMappingTest {


    List<MarginalTree> trees = new ArrayList<>();
    List<TreeLikelihood> likelihoods = new ArrayList<>();
    long nChecks = 0;

    class ValidatingMCMC extends MCMC {
        @Override
        protected Operator propagateState(long sampleNr) {
            Operator operator = super.propagateState(sampleNr);
            checkCaches(sampleNr, operator);
            return operator;
        }
    }

    void checkCaches(long sampleNr, Operator operator) {
        for (int i = 0; i < trees.size(); i++) {
            MarginalTree tree = trees.get(i);
            TreeLikelihood likelihood = likelihoods.get(i);
            String context = "sample " + sampleNr + " after " + operator.getID();

            if (tree.isModifiedInProposal())
                throw new AssertionError(context + ": " + tree.getID() + " still has proposal state");
            tree.ensureUpToDate();
            if (tree.isModifiedInProposal())
                throw new AssertionError(context + ": " + tree.getID() + " was out of date");
            tree.checkConsistency(context);

            // Recompute the likelihood from scratch inside a store/restore window
            double cached = likelihood.getCurrentLogP();
            tree.makeAllDirty();
            likelihood.store();
            SharedTransitionMatrixCache.bypass = true;
            double fresh;
            try {
                fresh = likelihood.calculateLogP();
            } finally {
                SharedTransitionMatrixCache.bypass = false;
            }
            likelihood.restore();
            tree.setEverythingDirty(false);

            if (Math.abs(cached - fresh) > 1e-8 * Math.max(1.0, Math.abs(fresh)))
                throw new AssertionError(context + ": cached log likelihood of " + likelihood.getID()
                        + " is " + cached + " but recomputation gives " + fresh);
            nChecks++;
        }
    }

    @Test
    public void testCachesStayConsistentUnderAllOperators() throws Exception {
        runChain(20000, false, 7, false, false);
    }

    /**
     * Chain with an estimated GTR+G model, so that the likelihoods share transition matrices and
     * the substitution model changes (and is restored) during the run. Every cached likelihood is
     * compared against a recomputation that bypasses the shared matrices.
     */
    @Test
    public void testSharedTransitionMatricesStayConsistent() throws Exception {
        SharedTransitionMatrixCache.resetStatistics();
        runChain(10000, false, 13, true, false);
        for (TreeLikelihood likelihood : likelihoods)
            assertTrue(((CTreeLikelihood) likelihood).usesSharedTransitionMatrices());
        assertTrue("no matrix was ever shared", SharedTransitionMatrixCache.getHitCount() > 0);
    }

    /**
     * Chain with a relaxed clock (one estimated rate per CF branch) and a frozen taxon, so that the
     * marginal trees track the CF paths of their branches and recompute path-averaged rates. The
     * branch rates are scaled individually, so clock-only proposals are accepted and rejected.
     * checkConsistency compares paths and rates against a fresh sweep after every step.
     */
    @Test
    public void testRelaxedClockCachesStayConsistent() throws Exception {
        runChain(15000, false, 23, false, true);
        for (MarginalTree tree : trees)
            assertTrue(tree.getPathBranches(0) != null);
    }

    /**
     * The relaxed-clock chain in debug mode (checksums include the branch rates).
     */
    @Test
    public void testRelaxedClockRestoreChecksumsInDebugMode() throws Exception {
        String old = System.getProperty("beast.debug");
        System.setProperty("beast.debug", "true");
        try {
            runChain(3000, true, 19, false, true);
        } finally {
            if (old == null)
                System.clearProperty("beast.debug");
            else
                System.setProperty("beast.debug", old);
        }
    }

    /**
     * Shared and unshared transition matrices give identical likelihoods.
     */
    @Test
    public void testSharedMatricesMatchUnshared() {
        for (long seed = 1; seed <= 5; seed++) {
            simulateRandomACG(seed, 0.4, 0.5);
            Alignment alignment = getAlignment(10);
            SiteModel siteModel = gtrSiteModel(alignment, new RealScalarParam<>(3.0, PositiveReal.INSTANCE));
            for (Block block : randomBlockSet.getBlocks()) {
                MarginalTree tree = newMarginalTree(randomAcg, block);
                double[] logP = new double[2];
                for (int shared = 0; shared < 2; shared++) {
                    StrictClockModel clock = new StrictClockModel();
                    clock.initByName("clock.rate", new RealScalarParam<>(0.3, PositiveReal.INSTANCE));
                    CTreeLikelihood likelihood = new CTreeLikelihood();
                    likelihood.initByName("data", alignment, "tree", tree, "siteModel", siteModel,
                            "branchRateModel", clock, "shareTransitionMatrices", shared == 1);
                    assertTrue(likelihood.usesSharedTransitionMatrices() == (shared == 1));
                    logP[shared] = likelihood.calculateLogP();
                }
                assertTrue(logP[0] + " vs " + logP[1], logP[0] == logP[1]);
            }
        }
    }

    static SiteModel gtrSiteModel(Alignment alignment, RealScalarParam<PositiveReal> rateAG) {
        Frequencies freqs = new Frequencies();
        freqs.initByName("data", alignment);
        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs, "rateAG", rateAG);
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("gammaCategoryCount", 4,
                "shape", new RealScalarParam<>(0.5, PositiveReal.INSTANCE), "substModel", gtr);
        return siteModel;
    }

    /**
     * Same chain in BEAST's debug mode, which additionally checks the checksums of all calculation
     * nodes (including the marginal trees) after every rejected proposal and recomputes the
     * posterior from scratch.
     */
    @Test
    public void testRestoreChecksumsInDebugMode() throws Exception {
        String old = System.getProperty("beast.debug");
        System.setProperty("beast.debug", "true");
        try {
            runChain(3000, true, 11, false, false);
        } finally {
            if (old == null)
                System.clearProperty("beast.debug");
            else
                System.setProperty("beast.debug", old);
        }
    }

    void runChain(int chainLength, boolean debug, long seed, boolean gtr, boolean relaxed) throws Exception {
        trees.clear();
        likelihoods.clear();
        nChecks = 0;
        simulateRandomACG(seed, 0.4, 0.4);
        Randomizer.setSeed(seed * 1000 + 1);

        randomAcg.setID("acg");
        Alignment alignment = getAlignment(10);
        RealScalarParam<PositiveReal> rateAG = new RealScalarParam<>(3.0, PositiveReal.INSTANCE);
        rateAG.setID("rateAG");
        SiteModel siteModel;
        if (gtr) {
            siteModel = gtrSiteModel(alignment, rateAG);
        } else {
            JukesCantor jc = new JukesCantor();
            jc.initByName();
            siteModel = new SiteModel();
            siteModel.initByName("substModel", jc);
        }

        RealScalarParam<PositiveReal> clockRate = new RealScalarParam<>(0.2, PositiveReal.INSTANCE);
        clockRate.setID("clockRate");
        RealVectorParam<PositiveReal> branchRates = null;
        BranchRateModel clock;
        if (relaxed) {
            double[] rates = new double[randomAcg.getNodeCount() - 1];
            java.util.Random random = new java.util.Random(seed);
            for (int i = 0; i < rates.length; i++)
                rates[i] = 0.5 + random.nextDouble();
            branchRates = new RealVectorParam<>(rates, PositiveReal.INSTANCE);
            branchRates.setID("branchRates");
            UCRelaxedClockModel ucld = new UCRelaxedClockModel();
            ucld.initByName("rates", branchRates, "distr", new Uniform(), "tree", randomAcg, "clock.rate", clockRate);
            FreezableClock freezable = new FreezableClock();
            freezable.initByName("clock", ucld, "frozenTaxa", List.of(new Taxon("4")));
            clock = freezable;
        } else {
            StrictClockModel strict = new StrictClockModel();
            strict.initByName("clock.rate", clockRate);
            clock = strict;
        }

        // Likelihoods
        List<Distribution> likelihoodList = new ArrayList<>();
        for (Block block : randomBlockSet.getBlocks()) {
            // validate: every partial resweep is checked against a full sweep
            MarginalTree tree = new MarginalTree();
            tree.initByName("network", randomAcg, "block", block, "nodetype", MarginalNode.class.getName(),
                    "validate", true);
            tree.setID("mtree." + block.getID());
            CTreeLikelihood likelihood = new CTreeLikelihood();
            likelihood.initByName("data", alignment, "tree", tree, "siteModel", siteModel, "branchRateModel", clock);
            likelihood.setID("treeLikelihood." + block.getID());
            trees.add(tree);
            likelihoods.add(likelihood);
            likelihoodList.add(likelihood);
        }
        CompoundDistribution likelihood = new CompoundDistribution();
        likelihood.initByName("distribution", likelihoodList);
        likelihood.setID("likelihood");

        // Priors
        RealVectorParam<NonNegativeReal> conversionRates = new RealVectorParam<>(new double[] {0.4}, NonNegativeReal.INSTANCE);
        conversionRates.setID("conversionRates");
        RealScalarParam<UnitInterval> pMove = new RealScalarParam<>(0.4, UnitInterval.INSTANCE);
        pMove.setID("pMove");
        RealScalarParam<PositiveReal> birthRate = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
        birthRate.setID("birthRate");

        YuleModel yule = new YuleModel();
        yule.initByName("tree", randomAcg, "birthDiffRate", birthRate);
        ConversionPrior convPrior = new ConversionPrior();
        convPrior.initByName("network", randomAcg, "conversionRates", conversionRates);
        ConversionMovePrior movePrior = new ConversionMovePrior();
        movePrior.initByName("network", randomAcg, "blockSet", randomBlockSet, "pMove", pMove);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", List.of(yule, convPrior, movePrior));

        CompoundDistribution posterior = new CompoundDistribution();
        posterior.initByName("distribution", List.of(prior, likelihood));

        // State
        List<StateNode> stateNodes = new ArrayList<>();
        stateNodes.add(randomAcg);
        stateNodes.addAll(randomBlockSet.getBlocks());
        stateNodes.add(clockRate);
        if (gtr)
            stateNodes.add(rateAG);
        if (relaxed)
            stateNodes.add(branchRates);
        for (StateNode fixed : List.of(conversionRates, pMove, birthRate)) {
            fixed.isEstimatedInput.setValue(false, fixed);
            stateNodes.add(fixed);
        }
        State state = new State();
        state.initByName("stateNode", stateNodes);

        // Operators
        List<Operator> operators = new ArrayList<>();
        Object[] creationArgs = {"acg", randomAcg, "pMove", pMove, "blockSet", randomBlockSet,
                "conversionPrior", convPrior};
        operators.add(op(new ACGScaler(), "acg", randomAcg, "scaleFactor", 0.9));
        operators.add(op(new ACGScaler(), "acg", randomAcg, "scaleFactor", 0.95, "rootOnly", true));
        operators.add(op(new CFUniform(), creationArgs));
        operators.add(op(new CFSubtreeSlide(), creationArgs));
        operators.add(op(new CFSubtreeSlide(), with(creationArgs, "gibbsSampleBlockMoves", true, "treeLikelihood", likelihoods)));
        operators.add(op(new CFWilsonBalding(), with(creationArgs, "alpha", 0.1)));
        operators.add(op(new CFSubtreeExchange(), creationArgs));
        operators.add(op(new CFSubtreeExchange(), with(creationArgs, "isNarrow", false)));
        operators.add(op(new CFConversionSwap(), creationArgs));
        operators.add(op(new CFConversionSwap(), with(creationArgs, "resampleBorrowings", true, "treeLikelihood", likelihoods)));
        operators.add(op(new AddRemoveConversion(), creationArgs));
        operators.add(op(new AddRemoveConversionGibbs(), with(creationArgs, "treeLikelihood", likelihoods)));
        operators.add(op(new ConversionSplit(), creationArgs));
        operators.add(op(new ConvertedEdgeSlide(), "acg", randomAcg, "apertureSize", 0.2));
        operators.add(op(new ConvertedEdgeHop(), "acg", randomAcg));
        operators.add(op(new ConvertedEdgeFlip(), "acg", randomAcg));
        operators.add(op(new ConvertedEdgeHopGibbs(), with(creationArgs, "treeLikelihood", likelihoods)));
        operators.add(op(new ConvertedEdgeHopGibbs(), with(creationArgs, "treeLikelihood", likelihoods, "nClosestRelatives", 3)));
        operators.add(op(new GibbsSampleMovesPerConversion(), "acg", randomAcg, "pMove", pMove,
                "blockSet", randomBlockSet, "treeLikelihood", likelihoods));
        operators.add(op(new GibbsSampleMovesPerConversion(), "acg", randomAcg, "pMove", pMove,
                "blockSet", randomBlockSet, "treeLikelihood", likelihoods, "mcmcmc", true));
        operators.add(op(new ResampleBorrowings(), "acg", randomAcg, "pMove", pMove, "blockSet", randomBlockSet));
        operators.add(op(new ACGTipDatesRandomWalker(), "tree", randomAcg, "windowSize", 0.1));
        ScaleOperator clockScaler = new ScaleOperator();
        clockScaler.initByName("parameter", clockRate, "scaleFactor", 0.97, "weight", 1.0, "optimise", false);
        clockScaler.setID("ClockRateScaler");
        operators.add(clockScaler);
        if (relaxed) {
            ScaleOperator branchRateScaler = new ScaleOperator();
            branchRateScaler.initByName("parameter", branchRates, "scaleFactor", 0.7, "weight", 3.0, "optimise", false);
            branchRateScaler.setID("BranchRateScaler");
            operators.add(branchRateScaler);
        }
        if (gtr) {
            ScaleOperator rateScaler = new ScaleOperator();
            rateScaler.initByName("parameter", rateAG, "scaleFactor", 0.8, "weight", 3.0, "optimise", false);
            rateScaler.setID("RateAGScaler");
            operators.add(rateScaler);
        }

        Logger screenLog = new Logger();
        screenLog.initByName("logEvery", chainLength, "log", posterior);

        MCMC mcmc = new ValidatingMCMC();
        mcmc.initByName(
                "chainLength", (long) chainLength,
                "state", state,
                "distribution", posterior,
                "operator", operators,
                "logger", screenLog);
        mcmc.run();

        assertTrue(nChecks >= (long) chainLength * trees.size());
        long partial = 0;
        for (MarginalTree tree : trees)
            partial += tree.nPartialResweeps;
        assertTrue("no partial resweeps", partial > 0);
        // The short debug-mode run is too short for every operator to be accepted
        for (Operator operator : operators)
            assertTrue(operator.getID() + " was never accepted", debug ||
                    operator.get_m_nNrAccepted() > 0 || operator instanceof CFWilsonBalding
                    || operator instanceof ScaleOperator);
        for (Operator operator : operators)
            if (operator.getID().equals("RateAGScaler") || operator.getID().equals("BranchRateScaler"))
                assertTrue(operator.getID() + " must be both accepted and rejected",
                        operator.get_m_nNrAccepted() > 0 && operator.get_m_nNrRejected() > 0);
    }

    static Object[] with(Object[] base, Object... extra) {
        Object[] all = new Object[base.length + extra.length];
        System.arraycopy(base, 0, all, 0, base.length);
        System.arraycopy(extra, 0, all, base.length, extra.length);
        return all;
    }

    static Operator op(Operator operator, Object... args) {
        Object[] all = with(args, "weight", 1.0);
        operator.initByName(all);
        String id = operator.getClass().getSimpleName();
        for (int i = 0; i < args.length; i += 2) {
            if (args[i].equals("gibbsSampleBlockMoves") || args[i].equals("resampleBorrowings")
                    || args[i].equals("rootOnly") || args[i].equals("isNarrow") || args[i].equals("mcmcmc")
                    || args[i].equals("nClosestRelatives"))
                id += "." + args[i] + "=" + args[i + 1];
        }
        operator.setID(id);
        return operator;
    }
}
