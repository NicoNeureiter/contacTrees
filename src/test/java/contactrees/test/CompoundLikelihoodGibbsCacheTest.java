package contactrees.test;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.State;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalNode;
import contactrees.MarginalTree;
import contactrees.operators.BorrowingOperator;

/**
 * Reproduces (or refutes) a suspected cache-staleness bug one level ABOVE
 * {@link MarginalTreeGibbsCacheTest}: a plain BEAST {@code CompoundDistribution} wrapping several
 * per-block {@code TreeLikelihood}s, driven through a REAL {@link State} (store/
 * checkCalculationNodesDirtiness/calculateLogP/accept/restore), with the Gibbs choreography run
 * via the real {@link BorrowingOperator} methods ({@code flipBorrowing}/
 * {@code calcLogBorrowingPosterior}), exactly mirroring {@code sampleBlockMove}'s "keep" branch
 * (flip, recalc, do NOT revert).
 *
 * <p>Found empirically: BEAST's own periodic self-check (every 10000 samples, independent of
 * {@code -Dbeast.debug}) reported "Likelihood incorrectly calculated" for the real
 * {@code bantu-contactrees.xml} run (seed 12345), attributing it to
 * {@code AddRemoveConversionGibbs}. The breakdown showed every individual per-block
 * featureLikelihood matched its own stored value -- only the WRAPPING
 * {@code CompoundDistribution}'s own current vs. stored logP disagreed, by the full diff. This
 * test isolates whether a single Gibbs "keep" cycle (flip without revert) against a real
 * {@code State} leaves that wrapping compound's own store/restore bookkeeping out of sync with a
 * fresh rebuild.
 *
 * @author Nico Neureiter
 */
public class CompoundLikelihoodGibbsCacheTest extends ContactreesTest {

    private static class TestBorrowingOperator extends BorrowingOperator {
        @Override
        protected double doProposal() {
            throw new UnsupportedOperationException("not used directly in this test");
        }
    }

    private MarginalTree buildMarginalTreeWithRelaxedClock(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        double[] rates = {2.0, 0.5, 1.5, 0.8, 1.2, 2.5, 0.6, 1.8};
        UCRelaxedClockModel clock = new UCRelaxedClockModel();
        clock.initByName(
                "rates", new RealVectorParam<>(rates, PositiveReal.INSTANCE),
                "distr", new Uniform(),
                "tree", acg2);
        marginalTree.initByName(
                "network", acg2,
                "block", block,
                "nodetype", MarginalNode.class.getName(),
                "branchRateModel", clock);
        return marginalTree;
    }

    private TreeLikelihood buildTreeLikelihood(MarginalTree marginalTree) {
        Alignment alignment = getAlignment(acg2.getLeafNodeCount());
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", jc);

        TreeLikelihood treeLH = new TreeLikelihood();
        treeLH.initByName("data", alignment, "tree", marginalTree, "siteModel", siteModel);
        return treeLH;
    }

    private double freshCompoundLogP(Block block0, Block block1) {
        MarginalTree t0 = buildMarginalTreeWithRelaxedClock(block0);
        MarginalTree t1 = buildMarginalTreeWithRelaxedClock(block1);
        TreeLikelihood lh0 = buildTreeLikelihood(t0);
        TreeLikelihood lh1 = buildTreeLikelihood(t1);
        CompoundDistribution fresh = new CompoundDistribution();
        fresh.initByName("distribution", List.of(lh0, lh1));
        return fresh.calculateLogP();
    }

    @Test
    public void testCompoundStaysInSyncAfterGibbsKeepCycle() throws Exception {
        Block block0 = blockSet2.getBlocks().get(0);
        Block block1 = blockSet2.getBlocks().get(1);
        Conversion conv = conv2_1;

        MarginalTree marginalTree0 = buildMarginalTreeWithRelaxedClock(block0);
        MarginalTree marginalTree1 = buildMarginalTreeWithRelaxedClock(block1);
        TreeLikelihood treeLH0 = buildTreeLikelihood(marginalTree0);
        TreeLikelihood treeLH1 = buildTreeLikelihood(marginalTree1);

        CompoundDistribution compound = new CompoundDistribution();
        compound.initByName("distribution", List.of(treeLH0, treeLH1));

        TestBorrowingOperator op = new TestBorrowingOperator();
        RealScalar<UnitInterval> pMove = new RealScalarParam<>(0.3, UnitInterval.INSTANCE);
        op.initByName("acg", acg2, "blockSet", blockSet2, "pMove", pMove,
                "treeLikelihood", List.of(treeLH0, treeLH1), "weight", 1.0);

        State state = new State();
        state.initByName("stateNode", List.of(acg2, block0, block1));
        state.initialise();
        state.setPosterior(compound);

        // Establish the baseline, as MCMC would have it at the end of some earlier, settled
        // proposal: a real rebuild, then a real store/accept cycle.
        marginalTree0.recalculate();
        marginalTree1.recalculate();
        compound.calculateLogP();
        state.storeCalculationNodes();
        state.acceptCalculationNodes();

        // One outer MCMC-level proposal cycle, with a Gibbs "keep" (flip, recalc, no revert)
        // on block0 for the existing conversion -- mirrors sampleBlockMove()'s kept-flip branch,
        // using the real methods (flipBorrowing/calcLogBorrowingPosterior), not a reimplementation.
        state.store(1);
        treeLH0.store();
        op.calcLogBorrowingPosterior(conv, treeLH0, false);
        op.flipBorrowing(block0, conv);
        op.calcLogBorrowingPosterior(conv, treeLH0, true);
        // KEEP: no flip-back, no treeLH0.restore() -- matches the real "kept" code path.

        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
        double newLogP = compound.calculateLogP();

        double freshLogP = freshCompoundLogP(block0, block1);
        assertEquals("compound logP right after the kept flip must match a fresh rebuild",
                freshLogP, newLogP, 1e-6);

        // Outer-level ACCEPT (this flip is kept).
        state.acceptCalculationNodes();

        assertEquals("compound's OWN cached current logP must still match a fresh rebuild "
                + "immediately after accept", freshLogP, compound.getCurrentLogP(), 1e-6);

        // The critical check: after a REJECT of some LATER, unrelated proposal, restore() sets
        // compound.logP = compound.storedLogP. If storedLogP was never refreshed to the
        // just-accepted value above, this reveals it by going stale right here.
        state.store(2);
        // Touch a StateNode so this later proposal is non-trivial (mirrors some other operator
        // changing the ACG without touching these two blocks' marginal trees at all).
        block1.startEditing(null);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
        compound.calculateLogP();
        state.restore();
        state.restoreCalculationNodes();

        assertEquals("compound's cached logP after a later, unrelated reject must still match "
                + "the state accepted in the previous cycle, not some earlier stale value",
                freshLogP, compound.getCurrentLogP(), 1e-6);
    }
}
