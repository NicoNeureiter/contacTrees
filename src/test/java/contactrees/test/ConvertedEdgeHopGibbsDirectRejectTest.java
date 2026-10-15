package contactrees.test;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.inference.parameter.RealScalarParam;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalNode;
import contactrees.MarginalTree;
import contactrees.model.ConversionPrior;
import contactrees.operators.ConvertedEdgeHopGibbs;

/**
 * A ConvertedEdgeHopGibbs proposal that is rejected directly (returns -Infinity) must leave every
 * per-block likelihood consistent with the state MCMC restores afterwards.
 *
 * <p>MCMC's direct-reject path calls only state.restore(): no calculation node is stored or
 * restored. So an operator that has already recalculated marginal trees and evaluated
 * likelihoods (the Gibbs borrowing resampling) before returning -Infinity leaves those trees
 * built for the rejected proposal, and the next step that re-evaluates a likelihood uses them.
 *
 * <p>Setup: the caterpillar ((((1,2),3),4),5) with one conversion 4 -> 5 at height 0.5, and the
 * operator configured with sourceOnly and nClosestRelatives=1, so it always moves the arrival
 * end (leaf 5). Leaf 5's closest relatives at 0.5 are {1,2,3,4}; without the conversion's own
 * ends that leaves {1,2,3}, and none of those has leaf 5 among its own two closest relatives.
 * Every proposal therefore has no back-move and is rejected directly.
 */
public class ConvertedEdgeHopGibbsDirectRejectTest extends ContactreesTest {

    private static final String CATERPILLAR = "((((1:1.0,2:1.0)6:1.0,3:2.0)7:1.0,4:3.0)8:1.0,5:4.0)9:0.0;";

    private RealScalarParam<PositiveReal> mutationRate;
    private List<Block> blocks;
    private List<TreeLikelihood> treeLHs;
    private CompoundDistribution compound;
    private ConvertedEdgeHopGibbs operator;
    private State state;

    private Node leaf(String id) {
        for (Node node : acg.getExternalNodes())
            if (node.getID().equals(id))
                return node;
        throw new IllegalArgumentException("no leaf " + id);
    }

    private TreeLikelihood newTreeLikelihood(Block block, SiteModel siteModel) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg, "block", block, "nodetype", MarginalNode.class.getName());
        TreeLikelihood treeLH = new TreeLikelihood();
        treeLH.initByName("data", getAlignment(acg.getLeafNodeCount()), "tree", marginalTree, "siteModel", siteModel);
        return treeLH;
    }

    private SiteModel newSiteModel(RealScalarParam<PositiveReal> rate) {
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", jc, "mutationRate", rate);
        return siteModel;
    }

    private void setUpModel() {
        setupFromNewick(CATERPILLAR);
        Conversion conv = new Conversion(leaf("4"), leaf("5"), 0.5, acg, 1);
        acg.addConversion(conv);

        blocks = blockSet.getBlocks();
        // Half the blocks borrow over the conversion, so their marginal trees differ from the CF.
        for (int i = 0; i < blocks.size(); i += 2)
            blocks.get(i).addMove(conv);

        mutationRate = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
        SiteModel siteModel = newSiteModel(mutationRate);
        treeLHs = new ArrayList<>();
        for (Block block : blocks)
            treeLHs.add(newTreeLikelihood(block, siteModel));

        compound = new CompoundDistribution();
        compound.initByName("distribution", new ArrayList<>(treeLHs));

        RealScalarParam<NonNegativeReal> conversionRate = new RealScalarParam<>(0.1, NonNegativeReal.INSTANCE);
        ConversionPrior conversionPrior = new ConversionPrior();
        conversionPrior.initByName("network", acg, "conversionRate", conversionRate);

        operator = new ConvertedEdgeHopGibbs();
        operator.initByName("acg", acg, "blockSet", blockSet,
                "pMove", new RealScalarParam<>(0.5, UnitInterval.INSTANCE),
                "treeLikelihood", new ArrayList<>(treeLHs),
                "conversionRate", conversionRate, "conversionPrior", conversionPrior,
                "sourceOnly", true, "nClosestRelatives", 1, "weight", 1.0);

        List<StateNode> stateNodes = new ArrayList<>();
        stateNodes.add(acg);
        stateNodes.add(mutationRate);
        stateNodes.addAll(blocks);
        state = new State();
        state.initByName("stateNode", stateNodes);
        state.initialise();
        state.setPosterior(compound);

        // As MCMC.run() does before the first step.
        state.robustlyCalcPosterior(compound);
    }

    /** One MCMC step as MCMC.propagateState() runs it for a directly rejected proposal. */
    private void directlyRejectedStep(long sampleNr) {
        state.store(sampleNr);
        double logHR = operator.proposal();
        assertEquals("the setup must force a direct reject", Double.NEGATIVE_INFINITY, logHR, 0.0);
        // operator.requiresStateInitialisation() is true: no restoreCalculationNodes(), no
        // setEverythingDirty(false) -- only the state nodes are restored.
        state.restore();
    }

    /** One accepted MCMC step around an arbitrary state change, as MCMC.propagateState() runs it. */
    private void acceptedStep(long sampleNr, Runnable change) {
        state.store(sampleNr);
        change.run();
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
        compound.calculateLogP();
        state.acceptCalculationNodes();
        state.setEverythingDirty(false);
    }

    private void assertAllLikelihoodsMatchFreshRebuild(String phase) {
        SiteModel freshSiteModel = newSiteModel(
                new RealScalarParam<>(mutationRate.get(), PositiveReal.INSTANCE));
        for (int i = 0; i < blocks.size(); i++) {
            double fresh = newTreeLikelihood(blocks.get(i), freshSiteModel).calculateLogP();
            assertEquals(phase + ": block " + i + " likelihood must match a fresh rebuild",
                    fresh, treeLHs.get(i).getCurrentLogP(), 1e-8);
        }
    }

    @Test
    public void testSiteModelChangeAfterDirectRejectUsesRestoredTrees() {
        setUpModel();
        assertAllLikelihoodsMatchFreshRebuild("initial state");

        directlyRejectedStep(1);
        assertAllLikelihoodsMatchFreshRebuild("after the direct reject");

        // Touches only the likelihoods' site model, so no marginal tree is on this step's
        // calculation path: the likelihoods are re-evaluated on whatever trees they hold.
        acceptedStep(2, () -> mutationRate.set(1.7));
        assertAllLikelihoodsMatchFreshRebuild("after a site-model-only step");
    }

    @Test
    public void testClonalFrameChangeAfterDirectRejectUsesRestoredTrees() {
        setUpModel();
        directlyRejectedStep(1);

        // Touches the ACG, so every marginal tree is rebuilt on this step's calculation path
        // and each likelihood is updated incrementally from its cached partials.
        Node cherry = leaf("1").getParent();
        acceptedStep(2, () -> cherry.setHeight(cherry.getHeight() - 0.1));
        assertAllLikelihoodsMatchFreshRebuild("after a clonal-frame step");
    }
}
