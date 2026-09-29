package contactrees.test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.spec.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealVectorParam;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Fast (unit-test-scale) reproduction of Stage 2's Bug 2: a suspected cache-consistency bug in
 * the interaction between MarginalTree's per-node dirty marking and the Gibbs-sampling
 * choreography in {@code BorrowingOperator.sampleBlockMove()}/{@code getBorrowingsProbGibbs()}
 * (store -> recalculate(original) -> flip -> recalculate(flipped) -> maybe-revert-without-a-
 * third-recalculate -> treeLikelihood.restore()). See the plan file
 * ({@code in-coalre-we-recently-enchanted-feather.md}, "Stage 2 retry") and memory
 * {@code stage2-reverted-not-safe.md} for the full history: this exact scenario previously only
 * reproduced in a full empirical MCMC run (bantu-contactrees.xml), never at unit-test scale.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeGibbsCacheTest extends ContactreesTest {

    /**
     * Build a MarginalTree + TreeLikelihood pair over the given block, using the shared acg2
     * fixture and a fresh alignment/site model each time (so the two trees being compared never
     * share mutable state).
     */
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

    private MarginalTree buildMarginalTree(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName());
        return marginalTree;
    }

    /**
     * A real (non-rate-1) relaxed clock over acg2 (5 leaves, 8 non-root branches), so the
     * timeLength-based dirty criterion is genuinely exercised (timeLength != height for any
     * branch), rather than only incidentally through hasBranchRates always being true.
     */
    private UCRelaxedClockModel relaxedClock() {
        double[] rates = {2.0, 0.5, 1.5, 0.8, 1.2, 2.5, 0.6, 1.8};
        UCRelaxedClockModel clock = new UCRelaxedClockModel();
        clock.initByName(
                "rates", new RealVectorParam<>(rates, PositiveReal.INSTANCE),
                "distr", new Uniform(),
                "tree", acg2);
        return clock;
    }

    private MarginalTree buildMarginalTreeWithRelaxedClock(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName(
                "network", acg2,
                "block", block,
                "nodetype", MarginalNode.class.getName(),
                "branchRateModel", relaxedClock());
        return marginalTree;
    }

    /**
     * Independently, freshly rebuild a MarginalTree + TreeLikelihood over the block's CURRENT
     * state and return its logP -- the ground truth the Gibbs choreography's result is checked
     * against.
     */
    private double freshLogP(Block block) {
        MarginalTree freshTree = buildMarginalTree(block);
        TreeLikelihood freshLH = buildTreeLikelihood(freshTree);
        return freshLH.calculateLogP();
    }

    private double freshLogPRelaxed(Block block) {
        MarginalTree freshTree = buildMarginalTreeWithRelaxedClock(block);
        TreeLikelihood freshLH = buildTreeLikelihood(freshTree);
        return freshLH.calculateLogP();
    }

    /**
     * Mirrors BorrowingOperator.sampleBlockMove()'s exact choreography for a single Gibbs step
     * that REVERTS (the suspected failure mode): store -> recalc(original) -> flip ->
     * recalc(flipped) -> flip back (NO third recalculate) -> restore. Afterwards simulates
     * BEAST's own next dirty-check (as would happen when MCMC continues past this Gibbs step)
     * and compares against an independent fresh rebuild.
     */
    @Test
    public void testGibbsRevertMatchesFreshRebuild() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion conv = conv2_1;

        MarginalTree marginalTree = buildMarginalTree(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        treeLH.store();
        marginalTree.recalculate();
        treeLH.calculateLogP();

        block.addMove(conv);
        marginalTree.recalculate();
        treeLH.calculateLogP();

        // Revert: flip back, but do NOT call marginalTree.recalculate() again here -- matches
        // sampleBlockMove()'s real code path exactly.
        block.removeMove(conv);
        treeLH.restore();

        // What BEAST's own machinery does next, once this Gibbs step is done and the chain
        // continues: ask the tree if it's dirty, and if so, recompute.
        marginalTree.requiresRecalculation();
        double logPAfter = treeLH.calculateLogP();

        double logPExpected = freshLogP(block);

        assertEquals("Gibbs step that reverts must leave the TreeLikelihood matching a fresh "
                + "rebuild of the (unchanged) block state",
                logPExpected, logPAfter, 1e-6);
    }

    /**
     * Same choreography, but KEEPS the flip (no revert) -- the other real code path. Included as
     * a sanity check on the test harness itself: if this one also failed, the harness (not Bug 2)
     * would be the suspect.
     */
    @Test
    public void testGibbsKeepFlipMatchesFreshRebuild() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion conv = conv2_1;

        MarginalTree marginalTree = buildMarginalTree(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        treeLH.store();
        marginalTree.recalculate();
        treeLH.calculateLogP();

        block.addMove(conv);
        marginalTree.recalculate();
        double logPNew = treeLH.calculateLogP();

        // Keep the flip: no revert, no treeLH.restore() (matches sampleBlockMove()'s
        // not-reverted branch, which never calls restore() either -- treeLH.store()'s snapshot
        // is simply discarded).

        double logPExpected = freshLogP(block);

        assertEquals("Gibbs step that keeps the flip must leave the TreeLikelihood matching a "
                + "fresh rebuild of the (changed) block state",
                logPExpected, logPNew, 1e-6);
    }

    /**
     * Same as {@link #testGibbsRevertMatchesFreshRebuild()}, but with a real relaxed clock (rate
     * != 1 on every branch), so timeLength genuinely differs from height everywhere -- the
     * dirty criterion Stage 2 actually keys on.
     */
    @Test
    public void testGibbsRevertMatchesFreshRebuildWithRelaxedClock() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion conv = conv2_1;

        MarginalTree marginalTree = buildMarginalTreeWithRelaxedClock(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        treeLH.store();
        marginalTree.recalculate();
        treeLH.calculateLogP();

        block.addMove(conv);
        marginalTree.recalculate();
        treeLH.calculateLogP();

        block.removeMove(conv);
        treeLH.restore();

        marginalTree.requiresRecalculation();
        double logPAfter = treeLH.calculateLogP();

        double logPExpected = freshLogPRelaxed(block);

        assertEquals("Gibbs step that reverts (relaxed clock) must leave the TreeLikelihood "
                + "matching a fresh rebuild of the (unchanged) block state",
                logPExpected, logPAfter, 1e-6);
    }

    /**
     * Mirrors AddRemoveConversionGibbs.drawNewConversion() specifically: a brand-new conversion
     * (not yet on ANY block) is added to the ACG, then drawBorrowingsGibbs -> sampleBlockMove is
     * called with mtreeChanged=false for its FIRST posterior evaluation -- meaning
     * calcLogBorrowingPosterior's "old" branch does NOT call marginalTree.recalculate() at all,
     * it trusts treeLH.getCurrentLogP() (correct, since a conversion with zero borrowings can't
     * have changed any block's marginal tree yet). This is a materially different choreography
     * from the flip-an-existing-conversion tests above (only ONE recalculate() happens before
     * the possible revert, not two), and is the actual code path AddRemoveConversionGibbs uses
     * where Bug 2 was originally (empirically) observed.
     */
    @Test
    public void testGibbsAddNewConversionMatchesFreshRebuild() throws Exception {
        Block block = blockSet2.getBlocks().get(0);

        MarginalTree marginalTree = buildMarginalTreeWithRelaxedClock(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        // Establish baseline cache (this block carries no conversions yet).
        marginalTree.recalculate();
        treeLH.calculateLogP();

        // A brand-new conversion, added to the ACG but not yet on any block -- mirrors
        // ConversionCreationOperator.addNewConversion() + attachEdge().
        Conversion newConv = new Conversion(node2_2, node2_3, 0.5, acg2, 0);
        acg2.addConversion(newConv);

        // calcLogBorrowingPosterior(newConv, treeLH, mtreeChanged=false): "old" branch, no
        // recalculate() call at all.
        treeLH.store();
        double logPOld = treeLH.getCurrentLogP();

        // flipBorrowing: block doesn't have newConv yet, so this adds it.
        block.addMove(newConv);
        marginalTree.recalculate();
        double logPNew = treeLH.calculateLogP();

        // Revert (as if pRevert was drawn): remove the move again, no third recalculate().
        block.removeMove(newConv);
        treeLH.restore();
        assertEquals("treeLH.restore() should bring back exactly the pre-flip cached logP",
                logPOld, treeLH.getCurrentLogP(), 1e-9);

        marginalTree.requiresRecalculation();
        double logPAfter = treeLH.calculateLogP();

        double logPExpected = freshLogPRelaxed(block);

        assertEquals("Gibbs add-new-conversion step that reverts must leave the TreeLikelihood "
                + "matching a fresh rebuild (block still carries no conversions); logPNew=" + logPNew,
                logPExpected, logPAfter, 1e-6);
    }

    /**
     * Wraps the inner Gibbs revert choreography inside the OUTER MCMC-level store/reject cycle
     * that every real proposal goes through: {@code acg2.store()} before the proposal begins,
     * then (simulating the proposal being rejected, the common case) {@code acg2.restore()} --
     * which swaps in brand-new {@code Conversion} objects via {@code copyTo()} -- followed by
     * {@code block.restore()} and {@code marginalTree.restore()}. None of the tests above
     * exercise this outer layer at all; the original empirical Bug 2 was only ever seen inside a
     * full MCMC run, where this layer is always present. Block starts and ends this test with no
     * conversions, so its convIDsStored (never explicitly stored here, since Block.store() is
     * protected and unreachable from this package) legitimately still matches what it should
     * revert to.
     */
    @Test
    public void testGibbsThenOuterRejectMatchesFreshRebuild() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion conv = conv2_1;

        MarginalTree marginalTree = buildMarginalTreeWithRelaxedClock(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        marginalTree.recalculate();
        treeLH.calculateLogP();

        // Simulate the end of the PREVIOUS proposal (accepted): in real MCMC this happens via
        // State.acceptCalculationNodes() before the next proposal's store() -- it's what tells
        // MarginalTree.restore()'s lazy snapshot machinery that a fresh proposal is starting.
        marginalTree.accept();

        // Outer MCMC-level "store" before the proposal begins.
        acg2.store();

        // Inner Gibbs choreography (revert case), exactly as in
        // testGibbsRevertMatchesFreshRebuildWithRelaxedClock().
        treeLH.store();
        marginalTree.recalculate();
        treeLH.calculateLogP();

        block.addMove(conv);
        marginalTree.recalculate();
        treeLH.calculateLogP();

        block.removeMove(conv);
        treeLH.restore();

        marginalTree.requiresRecalculation();
        treeLH.calculateLogP();

        // Outer MCMC-level "reject": ConversionGraph.restore() swaps in fresh Conversion
        // objects (via copyTo()); Block.restore() swaps its move list back; MarginalTree.restore()
        // (-> refresh() -> recalculate()) is what BEAST's own machinery calls to bring the
        // derived tree back in sync afterward.
        acg2.restore();
        block.restore();
        marginalTree.restore();

        double logPAfterOuterReject = treeLH.calculateLogP();
        double logPExpected = freshLogPRelaxed(block);

        assertEquals("After an outer MCMC-level reject following a Gibbs revert step, the "
                + "TreeLikelihood must match a fresh rebuild",
                logPExpected, logPAfterOuterReject, 1e-6);
    }

    /**
     * A short sequence of several revert/keep Gibbs steps in a row, across different
     * conversions, closer to GibbsSampleMovesPerConversion's real usage pattern (many
     * recalculate() calls against the same MarginalTree/TreeLikelihood instance in one
     * proposal). Checked against a fresh rebuild after every single step, not just at the end,
     * so a failure is localized to the step that broke it.
     */
    @Test
    public void testGibbsSequenceMatchesFreshRebuildAfterEveryStep() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion[] convs = {conv2_1, conv2_2, conv2_3};
        // Deterministic revert/keep pattern covering both branches multiple times.
        boolean[] revert = {true, false, true, true, false, false};

        MarginalTree marginalTree = buildMarginalTree(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        // Establish the initial cached state, as the first Gibbs step's treeLH.store() would
        // see it.
        marginalTree.recalculate();
        treeLH.calculateLogP();

        for (int step = 0; step < revert.length; step++) {
            Conversion conv = convs[step % convs.length];

            treeLH.store();
            marginalTree.recalculate();
            treeLH.calculateLogP();

            boolean wasAffectedBefore = block.isAffected(conv);
            if (wasAffectedBefore) {
                block.removeMove(conv);
            } else {
                block.addMove(conv);
            }
            marginalTree.recalculate();
            double logPFlipped = treeLH.calculateLogP();

            double logPStep;
            if (revert[step]) {
                if (wasAffectedBefore) {
                    block.addMove(conv);
                } else {
                    block.removeMove(conv);
                }
                treeLH.restore();
                marginalTree.requiresRecalculation();
                logPStep = treeLH.calculateLogP();
            } else {
                logPStep = logPFlipped;
            }

            double logPExpected = freshLogP(block);
            assertEquals("Step " + step + " (revert=" + revert[step] + ", conv=" + conv.getID()
                    + ") must match a fresh rebuild",
                    logPExpected, logPStep, 1e-6);
        }
    }

    /**
     * A sequence of several full outer-level accept/reject cycles in a row (each ending with
     * marginalTree.accept() or marginalTree.restore()), checking that the store()/restore()
     * snapshot machinery's bookkeeping (prevBranchDeterminant, prevRootHeight, and the stable-
     * numbering fields) stays correct across repeated snapshot/restore/snapshot cycles, not just
     * a single one.
     */
    @Test
    public void testOuterAcceptRejectSequenceMatchesFreshRebuild() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        Conversion[] convs = {conv2_1, conv2_2, conv2_3};
        // Deterministic accept/reject pattern covering both branches multiple times.
        boolean[] reject = {true, false, true, true, false, false};

        MarginalTree marginalTree = buildMarginalTreeWithRelaxedClock(block);
        TreeLikelihood treeLH = buildTreeLikelihood(marginalTree);

        marginalTree.recalculate();
        treeLH.calculateLogP();
        marginalTree.accept();

        for (int step = 0; step < reject.length; step++) {
            Conversion conv = convs[step % convs.length];

            acg2.store();

            boolean wasAffectedBefore = block.isAffected(conv);
            if (wasAffectedBefore) {
                block.removeMove(conv);
            } else {
                block.addMove(conv);
            }
            marginalTree.recalculate();
            treeLH.calculateLogP();

            double logPStep;
            if (reject[step]) {
                acg2.restore();
                block.restore();
                marginalTree.restore();
                logPStep = treeLH.calculateLogP();
            } else {
                marginalTree.accept();
                logPStep = treeLH.getCurrentLogP();
            }
            // Mirrors MCMC.propagateState()'s state.setEverythingDirty(false) call after every
            // proposal (accept or reject): resets hasStartedEditing on the edited StateNodes so
            // the NEXT proposal's first edit calls store() again instead of skipping it (Tree.
            // setEverythingDirty(false) -> setSomethingIsDirty(false) -> hasStartedEditing =
            // false). Without this, acg2/block's OWN store()/restore() (not MarginalTree's) would
            // silently stop snapshotting after the first edit -- unrelated to this task's new
            // MarginalTree store()/restore(), but required for this test to mirror real MCMC.
            acg2.setEverythingDirty(false);
            block.setEverythingDirty(false);

            double logPExpected = freshLogPRelaxed(block);
            assertEquals("Step " + step + " (reject=" + reject[step] + ", conv=" + conv.getID()
                    + ") must match a fresh rebuild",
                    logPExpected, logPStep, 1e-6);

            marginalTree.assertMappingConsistent("step " + step);
        }
    }

    /**
     * The untouched-block fast path: a block that is never touched during a proposal has no
     * snapshot taken (storedHeight == null), so MarginalTree.restore() must fall back to
     * refresh() rather than crashing on null snapshot arrays.
     */
    @Test
    public void testUntouchedBlockRestoreFallsBackToRefresh() throws Exception {
        Block untouchedBlock = blockSet2.getBlocks().get(1);
        Block touchedBlock = blockSet2.getBlocks().get(0);
        Conversion conv = conv2_1;

        MarginalTree untouchedTree = buildMarginalTree(untouchedBlock);
        TreeLikelihood untouchedLH = buildTreeLikelihood(untouchedTree);
        MarginalTree touchedTree = buildMarginalTree(touchedBlock);
        TreeLikelihood touchedLH = buildTreeLikelihood(touchedTree);

        untouchedTree.recalculate();
        untouchedLH.calculateLogP();
        untouchedTree.accept();

        touchedTree.recalculate();
        touchedLH.calculateLogP();
        touchedTree.accept();

        acg2.store();

        // Only touchedBlock's conversion set changes; untouchedTree never gets recalculate()
        // called on it this proposal, so it never snapshots.
        touchedBlock.addMove(conv);
        touchedTree.recalculate();
        touchedLH.calculateLogP();

        // Reject: bring both back via their own restore().
        acg2.restore();
        touchedBlock.restore();
        touchedTree.restore();
        untouchedTree.restore();

        double logPTouchedAfter = touchedLH.calculateLogP();
        double logPUntouchedAfter = untouchedLH.calculateLogP();

        assertEquals("Touched block must match a fresh rebuild after reject",
                freshLogP(touchedBlock), logPTouchedAfter, 1e-6);
        assertEquals("Untouched block (fast-path refresh) must match a fresh rebuild after reject",
                freshLogP(untouchedBlock), logPUntouchedAfter, 1e-6);
    }
}
