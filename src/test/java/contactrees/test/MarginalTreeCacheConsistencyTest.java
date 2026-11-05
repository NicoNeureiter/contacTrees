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
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.speciation.YuleModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.MarginalTree;
import contactrees.model.ConversionMovePrior;
import contactrees.model.ConversionPrior;
import contactrees.model.likelihood.CTreeLikelihood;
import contactrees.operators.ACGScaler;
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
            double fresh = likelihood.calculateLogP();
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
        runChain(20000, false, 7);
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
            runChain(3000, true, 11);
        } finally {
            if (old == null)
                System.clearProperty("beast.debug");
            else
                System.setProperty("beast.debug", old);
        }
    }

    void runChain(int chainLength, boolean debug, long seed) throws Exception {
        trees.clear();
        likelihoods.clear();
        nChecks = 0;
        simulateRandomACG(seed, 0.4, 0.4);
        Randomizer.setSeed(seed * 1000 + 1);

        randomAcg.setID("acg");
        Alignment alignment = getAlignment(10);
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", jc);

        RealScalarParam<PositiveReal> clockRate = new RealScalarParam<>(0.2, PositiveReal.INSTANCE);
        clockRate.setID("clockRate");
        StrictClockModel clock = new StrictClockModel();
        clock.initByName("clock.rate", clockRate);

        // Likelihoods
        List<Distribution> likelihoodList = new ArrayList<>();
        for (Block block : randomBlockSet.getBlocks()) {
            MarginalTree tree = newMarginalTree(randomAcg, block);
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
        operators.add(op(new ACGScaler(), "acg", randomAcg, "scaleFactor", 0.8, "rootOnly", true));
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
        ScaleOperator clockScaler = new ScaleOperator();
        clockScaler.initByName("parameter", clockRate, "scaleFactor", 0.97, "weight", 1.0, "optimise", false);
        clockScaler.setID("ClockRateScaler");
        operators.add(clockScaler);

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
        // The short debug-mode run is too short for every operator to be accepted
        for (Operator operator : operators)
            assertTrue(operator.getID() + " was never accepted", debug ||
                    operator.get_m_nNrAccepted() > 0 || operator instanceof CFWilsonBalding
                    || operator instanceof ScaleOperator);
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
