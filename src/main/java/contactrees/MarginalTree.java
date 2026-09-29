package contactrees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import contactrees.CFEventList.Event;
import contactrees.CFEventList.EventType;


/**
 * A calculation node extracting the marginal tree for a block from the ACG.
 * The resulting marginal trees are the basis for the likelihood computation
 * of the corresponding block.
 *
 * @author Nico Neureiter
 */

public class MarginalTree extends Tree {

    public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph (network) from which the marginal tree is sampled.",
            Input.Validate.REQUIRED);
    public Input<Block> blockInput = new Input<>(
            "block",
            "The block object cinitArraysontaining the moves this marginal tree follows along the conversion graph.",
            Input.Validate.REQUIRED);
    public Input<BranchRateModel> branchRateModelInput = new Input<>(
            "branchRateModel",
            "The branch rate model for rolling out the branch lengths of the marginal tree (if not given, a strict clock is assumed).",
            new StrictClockModel());
    public Input<ArrayList<String>> frozenTaxaInput = new Input<>(
            "frozenTaxa",
            "Taxa for which the last branch should have a fixed branch rate of 0.",
            new ArrayList<>());
    public Input<Integer> validationLevelInput = new Input<>("validationLevel",
            "0 = off; 1 = check this marginal tree against a fresh rebuild after every proposal; 2 = also check the CF/conversion <-> marginal-node mapping.", 0);

    public ConversionGraph acg;
    public Block block;
    protected BranchRateModel branchRateModel;
    protected boolean hasBranchRates;
    protected ArrayList<String> frozenTaxa;
    protected int validationLevel;

    protected boolean outdated;
    protected long seenCfVersion = -1;
    protected long seenConvVersion = -1;
    protected long seenBlockVersion = -1;
    protected boolean manuallyUpdated = false;
    protected String lastBlockState;
    protected boolean changed = true;
    protected int nextNonLeafNr;
    ArrayList<MarginalNode> activeCFlineages = new ArrayList<>();

    /*
     * Stable marginal-node numbering and the CF/conversion <-> marginal-node mapping.
     *
     * Internal marginal-node numbers used to be reassigned in height/event order on every
     * rebuild (a fresh nextNonLeafNr counter reset each recalculate()), so the same physical
     * "creator" -- a CF coalescence or a conversion -- got a different Nr almost every time.
     * That made any cross-rebuild comparison keyed on Nr (e.g. "did this node's children
     * change?") unreliable. Below, each creator keeps a persistent Nr across rebuilds via
     * mark-and-sweep: claimed once, kept as long as it keeps producing a marginal coalescence,
     * released back to a free list the moment it stops (overshadowed CF node, or a conversion
     * with no effect), and reused by the next creator that needs a slot.
     *
     * Leaves are unaffected: they are already identity-mapped to the CF leaf Nr.
     */

    /** Creator kinds for {@link #creatorKind}. */
    public static final int CREATOR_LEAF = 0;
    public static final int CREATOR_CF_COALESCENCE = 1;
    public static final int CREATOR_CONVERSION = 2;

    /** CF node Nr -> marginal Nr, for CF nodes that currently coalesce in this marginal tree. */
    protected int[] cfCoalescenceNr;
    /** Conversion id -> marginal Nr, for conversions that currently coalesce in this marginal tree. */
    protected Map<Integer, Integer> convMarginalNr;
    /** Internal marginal Nrs released by a creator that stopped producing a coalescence. */
    protected int[] freeNrStack;
    protected int freeNrTop;
    /** Inverse mapping: marginal Nr -> creator kind / creator key (CF node Nr or conversion id). */
    protected int[] creatorKind;
    protected int[] creatorKey;
    /** CF node Nr -> the marginal Nr currently governing that point in the CF, or -1. */
    protected int[] cfBranchOccupant;
    /** Bumped every recalculate(); lets external caches know the mapping may have changed. */
    protected long mappingVersion = 0;

    /*
     * Per-CF-node occupancy history: a single append-only, height-ordered log shared by all CF
     * nodes, threaded per node via occPrev so a query can still walk one node's entries. Flat and
     * contiguous instead of one tiny row array per CF node, for cache locality -- recalculate()
     * is the hottest method in the program. Reset (occCount := 0, occLastOfNode filled with -1)
     * at the start of every recalculate() -- see recalculate()'s reset loop. Nothing but tests
     * reads this yet; it exists purely as groundwork for a future incremental rebuild.
     */
    /** Event heights, indexed by log position. */
    protected double[] occHeight;
    /** Occupant marginal Nrs (-1 = vacated), indexed by log position, parallel to {@link #occHeight}. */
    protected int[] occNr;
    /** Log position of the previous entry for the same CF node, or -1, indexed by log position. */
    protected int[] occPrev;
    /** Number of entries used in the log for the current rebuild. */
    protected int occCount;
    /** CF node Nr -> log position of its most recent entry, or -1. */
    protected int[] occLastOfNode;

    /**
     * The root's own branch-length determinant has no meaning (it has no parent branch), so its
     * dirtiness is tracked separately from MarginalNode.prevBranchDeterminant: just its height.
     * NaN means "never assessed yet" (always dirty).
     */
    protected double prevRootHeight = Double.NaN;

    /*
     * Real store()/restore(): a snapshot of this tree's state taken lazily on the first
     * recalculate() of a proposal (see recalculate()'s startEditing(null) call), restored in
     * place (no array swap) by restore(). See the plan file's "Real MarginalTree.store()/
     * restore()" section for the full design rationale.
     */
    protected boolean snapshotTakenThisProposal = false;

    /*
     * Lowest height changed during the current proposal, for a future incremental rebuild that
     * only replays the walk from there upward. POSITIVE_INFINITY doubles as both "nothing
     * reported yet" and "not eligible this proposal" (a CF/conversion/branch-rate change, or
     * makeOutdated(), makes the whole tree ineligible for height-localized replay). Cleared at
     * the same points as snapshotTakenThisProposal -- see accept()/restore(). Plumbing only:
     * nothing reads this yet to change the rebuild.
     */
    protected double lowestChangedHeight = Double.POSITIVE_INFINITY;

    protected double[] storedHeight;
    protected double[] storedTimeLength;
    protected double[] storedLastEventHeight;
    protected int[] storedLeftNr;
    protected int[] storedRightNr;
    protected int storedRootNr;
    protected double storedPrevRootHeight;

    protected int[] storedCfCoalescenceNr;
    protected Map<Integer, Integer> storedConvMarginalNr;
    protected int[] storedFreeNrStack;
    protected int storedFreeNrTop;
    protected int[] storedCreatorKind;
    protected int[] storedCreatorKey;
    protected int[] storedCfBranchOccupant;
    protected int storedNextNonLeafNr;
    protected long storedSeenCfVersion;
    protected long storedSeenConvVersion;
    protected long storedSeenBlockVersion;
    protected long storedMappingVersion;

    public void setManuallyUpdated() {
        manuallyUpdated = true;
    }

    boolean isInBeauti() {
        return ProgramStatus.name.equals("BEAUti");
    }

    // Diagnostics for single-threaded benchmark runs; not synchronized.
    public static long nRecalculate = 0;
    public static long nFilthyNodesMarked = 0;
    public static long nNodesChecked = 0;
    public static long nRestore = 0;
    public static long nRestoreUntouched = 0;
    public static long nOutdatedChecks = 0;
    public static long nOutdatedTrue = 0;
    /** recalculate() calls where lowestChangedHeight was usable (finite) vs. not (infinite/ineligible). */
    public static long nRecalculateWithHeight = 0;
    public static long nRecalculateIneligible = 0;

    public static void resetCounters() {
        nRecalculate = 0;
        nFilthyNodesMarked = 0;
        nNodesChecked = 0;
        nRestore = 0;
        nRestoreUntouched = 0;
        nOutdatedChecks = 0;
        nOutdatedTrue = 0;
        nRecalculateWithHeight = 0;
        nRecalculateIneligible = 0;
    }

    public static String counterSummary() {
        long nHeightTotal = nRecalculateWithHeight + nRecalculateIneligible;
        return "MarginalTree counters: recalculate=" + nRecalculate
                + ", filthyNodesMarked=" + nFilthyNodesMarked
                + ", nodesChecked=" + nNodesChecked
                + ", filthyFraction=" + (nNodesChecked == 0 ? 0.0 : (double) nFilthyNodesMarked / nNodesChecked)
                + ", restore=" + nRestore
                + ", restoreUntouched=" + nRestoreUntouched
                + ", outdatedChecks=" + nOutdatedChecks
                + ", outdatedTrue=" + nOutdatedTrue
                + ", recalculateWithHeight=" + nRecalculateWithHeight
                + ", recalculateIneligible=" + nRecalculateIneligible
                + ", recalculateWithHeightFraction=" + (nHeightTotal == 0 ? 0.0 : (double) nRecalculateWithHeight / nHeightTotal);
    }

    public void setBranchRateModel(Base brm) {
        branchRateModelInput.setValue(brm, this);
        branchRateModel = brm;
        hasBranchRates = (branchRateModel != null);
    }

    @Override
    public void initAndValidate() {

        acg = networkInput.get();
        block = blockInput.get();
        // After block is set: the registry indexes trees by block.
        if (!isInBeauti())
            acg.registerMarginalTree(this);
        frozenTaxa = frozenTaxaInput.get();
        validationLevel = validationLevelInput.get();

        branchRateModel = branchRateModelInput.get();
//        hasBranchRates = !(branchRateModel instanceof StrictClockModel);
        hasBranchRates = (branchRateModel != null);

        // Initialize to clonal frame of acg
        String beastID = ID;
        assignFrom(acg);
        setID(beastID);

        super.initAndValidate();
        activeCFlineages = new ArrayList<>();
        for (int i=0; i<acg.getNodeCount(); i++)
            activeCFlineages.add(null);

        // One-time setup of the stable-numbering state, now that m_nodes/nodeCount are known.
        int nLeafs = acg.getLeafNodeCount();
        nextNonLeafNr = nLeafs;
        cfCoalescenceNr = newIntArrayFilled(acg.getNodeCount(), -1);
        cfBranchOccupant = newIntArrayFilled(acg.getNodeCount(), -1);
        convMarginalNr = new HashMap<>();
        freeNrStack = new int[Math.max(nLeafs - 1, 0)];
        freeNrTop = 0;
        creatorKind = new int[m_nodes.length];
        creatorKey = new int[m_nodes.length];

        int initialOccCapacity = 256;
        occHeight = new double[initialOccCapacity];
        occNr = new int[initialOccCapacity];
        occPrev = new int[initialOccCapacity];
        occCount = 0;
        occLastOfNode = newIntArrayFilled(acg.getNodeCount(), -1);

        recalculate();
        makeOutdated();
    }

    private static int[] newIntArrayFilled(int size, int fillValue) {
        int[] arr = new int[size];
        Arrays.fill(arr, fillValue);
        return arr;
    }

    @Override
    public boolean requiresRecalculation() {
        outdated = checkOutdated();
        // OR in, don't overwrite: a direct recalculate() call earlier this proposal (the
        // Gibbs-sampling code in BorrowingOperator calls it directly, bypassing this method)
        // may already have set changed=true to signal a real change to TreeLikelihood, which
        // reads only somethingIsDirty()/changed, never re-derives dirtiness itself. Overwriting
        // that with a fresh (possibly false) checkOutdated() result here would hide the change.
        changed = changed || outdated;

        if (outdated) {
            recalculate();

            if (manuallyUpdated) {
                manuallyUpdated = false;
                return false;
            } else
                return true;


        } else {
            return false;
        }
    }

    public boolean isValidating() {
        return validationLevel > 0;
    }

    /**
     * Verify that this marginal tree matches a rebuild from the current ACG and block state.
     * A mismatch means the tree was stale, i.e. the dirty-tracking logic failed to invalidate it.
     *
     * Calls {@link #recalculate()}, which mutates the tree into the correct state as a side
     * effect; that self-healing is intentional, not a bug.
     *
     * @param phase description of the calling context, used in the failure message
     */
    public void assertConsistent(String phase) {
        int n = getNodeCount();
        String oldNewick = getRoot().toSortedNewick(new int[1], false);
        double[] oldHeights = new double[n];
        double[] oldTimeLengths = new double[n];
        for (int i = 0; i < n; i++) {
            Node node = getNode(i);
            oldHeights[i] = node.getHeight();
            oldTimeLengths[i] = ((MarginalNode) node).timeLength;
        }

        recalculate();

        String newNewick = getRoot().toSortedNewick(new int[1], false);
        if (!oldNewick.equals(newNewick)) {
            throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): stale tree detected"
                    + " for marginal tree " + getID() + " (block " + block.getID() + "):"
                    + " sorted Newick mismatch."
                    + " before=" + truncate(oldNewick) + " after=" + truncate(newNewick));
        }

        for (int i = 0; i < n; i++) {
            Node node = getNode(i);
            double newHeight = node.getHeight();
            if (Math.abs(oldHeights[i] - newHeight) > 1e-12) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): stale tree detected"
                        + " for marginal tree " + getID() + " (block " + block.getID() + "):"
                        + " height mismatch at node " + i
                        + " before=" + oldHeights[i] + " after=" + newHeight);
            }

            double newTimeLength = ((MarginalNode) node).timeLength;
            if (Math.abs(oldTimeLengths[i] - newTimeLength) > 1e-12) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): stale tree detected"
                        + " for marginal tree " + getID() + " (block " + block.getID() + "):"
                        + " timeLength mismatch at node " + i
                        + " before=" + oldTimeLengths[i] + " after=" + newTimeLength);
            }
        }

        if (validationLevel >= 2) {
            assertMappingConsistent(phase);
        }
    }

    /**
     * Verify the stable-numbering invariants (validationLevel &gt;= 2): every internal Nr is
     * occupied exactly once by either a claimed creator or the free stack, and the forward/
     * inverse mapping round-trips in both directions.
     */
    public void assertMappingConsistent(String phase) {
        int nLeafs = acg.getLeafNodeCount();
        int nInternal = Math.max(nLeafs - 1, 0);
        boolean[] seen = new boolean[nInternal];

        for (int cfNr = 0; cfNr < cfCoalescenceNr.length; cfNr++) {
            int nr = cfCoalescenceNr[cfNr];
            if (nr < 0)
                continue;
            markOccupied(seen, nr, nLeafs, phase, "cfCoalescenceNr[" + cfNr + "]");
            if (creatorKind[nr] != CREATOR_CF_COALESCENCE || creatorKey[nr] != cfNr) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): "
                        + "creator round-trip failed for CF node " + cfNr + " -> marginal Nr " + nr);
            }
        }

        for (Map.Entry<Integer, Integer> e : convMarginalNr.entrySet()) {
            int nr = e.getValue();
            markOccupied(seen, nr, nLeafs, phase, "convMarginalNr[" + e.getKey() + "]");
            if (creatorKind[nr] != CREATOR_CONVERSION || creatorKey[nr] != e.getKey()) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): "
                        + "creator round-trip failed for conversion " + e.getKey() + " -> marginal Nr " + nr);
            }
        }

        for (int i = 0; i < freeNrTop; i++) {
            markOccupied(seen, freeNrStack[i], nLeafs, phase, "freeNrStack[" + i + "]");
        }

        for (int i = 0; i < nInternal; i++) {
            if (!seen[i]) {
                throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): "
                        + "internal Nr " + (i + nLeafs) + " is neither claimed nor free");
            }
        }
    }

    private void markOccupied(boolean[] seen, int nr, int nLeafs, String phase, String source) {
        int idx = nr - nLeafs;
        if (idx < 0 || idx >= seen.length || seen[idx]) {
            throw new IllegalStateException("MarginalTree.assertConsistent(" + phase + "): "
                    + "internal Nr " + nr + " (from " + source + ") is out of range or double-claimed");
        }
        seen[idx] = true;
    }

    private static String truncate(String s) {
        return s.length() > 500 ? s.substring(0, 500) + "..." : s;
    }

    public void makeOutdated() {
        outdated = true;
        setEverythingDirty(true);
        manuallyUpdated = false;
        makeIneligibleForIncrementalUpdate();
    }

    /**
     * Report that something changed at the given height during the current proposal (called by
     * operators that flip a block's borrowing of a single, unmoved conversion). Keeps the
     * running minimum across the proposal; has no effect once the tree has been made ineligible.
     */
    public void noteChangedAtHeight(double height) {
        if (height < lowestChangedHeight)
            lowestChangedHeight = height;
    }

    /**
     * @return the lowest height reported via {@link #noteChangedAtHeight(double)} so far this
     *         proposal, or {@link Double#POSITIVE_INFINITY} if nothing was reported, or if the
     *         tree is not eligible for a height-localized incremental rebuild this proposal.
     */
    public double getLowestChangedHeight() {
        return lowestChangedHeight;
    }

    /**
     * Force this tree out of eligibility for a height-localized incremental rebuild for the rest
     * of the current proposal: something changed that is not a plain block-move-set toggle on an
     * unmoved conversion (CF topology/height, conversion re-attachment/height, branch-rate-model
     * dirtiness, or a manual makeOutdated()).
     */
    public void makeIneligibleForIncrementalUpdate() {
        lowestChangedHeight = Double.POSITIVE_INFINITY;
    }

    @Override
    public boolean somethingIsDirty() {
        return this.changed;
    }

    public boolean checkOutdated() {
        nOutdatedChecks++;

        boolean result;
        if (outdated) {
            // Forced by makeOutdated(), which already made this tree ineligible.
            result = true;
        } else if (branchRateModel instanceof CalculationNode brm && brm.somethingIsDirty()) {
            result = true;
            makeIneligibleForIncrementalUpdate();
        } else if (acg.getCFVersion() != seenCfVersion) {
            result = true;
            makeIneligibleForIncrementalUpdate();
        } else if (block.getMoveSetVersion() != seenBlockVersion) {
            result = true;
        } else if (acg.getConvVersion() != seenConvVersion) {
            // Some conversion changed somewhere in the ACG; check precisely whether it was
            // one of THIS block's own conversions before concluding we're outdated.
            result = false;
            for (Conversion conv : getBlockConversions()) {
                if (conv.getVersion() > seenConvVersion) {
                    result = true;
                    makeIneligibleForIncrementalUpdate();
                    break;
                }
            }
            if (!result) {
                // Nothing relevant changed -- re-baseline so we don't re-scan on every
                // subsequent check until convVersion moves again.
                seenConvVersion = acg.getConvVersion();
            }
        } else {
            result = false;
        }

        if (result)
            nOutdatedTrue++;

        return result;
    }


    /**
     * Claim a stable marginal Nr for a creator (a CF node or a conversion) that does not
     * currently have one: reuse a freed slot if one is available, otherwise extend the range.
     */
    private int claimMarginalNr(int kind, int key) {
        int nr = (freeNrTop > 0) ? freeNrStack[--freeNrTop] : nextNonLeafNr++;
        if (nr >= m_nodes.length) {
            throw new IllegalStateException("MarginalTree exceeded node capacity while claiming a new"
                    + " internal node number; nr=" + nr + ", capacity=" + m_nodes.length);
        }
        creatorKind[nr] = kind;
        creatorKey[nr] = key;
        return nr;
    }

    /**
     * Record the active lineage at a CF node and mirror it into {@link #cfBranchOccupant}, so
     * the coverage map ("which CF nodes does this marginal Nr currently govern?") comes for
     * free from the existing bookkeeping.
     */
    private void setActiveLineage(int cfNodeNr, MarginalNode marginalNode, double height) {
        activeCFlineages.set(cfNodeNr, marginalNode);
        int nr = (marginalNode == null) ? -1 : marginalNode.getNr();
        cfBranchOccupant[cfNodeNr] = nr;
        recordOccupancy(cfNodeNr, height, nr);
    }

    /**
     * Append one occupancy-change entry to the shared log, growing it (doubling) if full, and
     * thread it onto the CF node's chain via occPrev/occLastOfNode. The log is reset (occCount
     * := 0) at the start of every recalculate(), not reallocated, so this only grows the log's
     * capacity, never its count of live entries at reset time.
     */
    private void recordOccupancy(int cfNodeNr, double height, int marginalNr) {
        if (occCount == occHeight.length) {
            int newCap = occHeight.length * 2;
            occHeight = Arrays.copyOf(occHeight, newCap);
            occNr = Arrays.copyOf(occNr, newCap);
            occPrev = Arrays.copyOf(occPrev, newCap);
        }
        occHeight[occCount] = height;
        occNr[occCount] = marginalNr;
        occPrev[occCount] = occLastOfNode[cfNodeNr];
        occLastOfNode[cfNodeNr] = occCount++;
    }

    /**
     * @param cfNodeNr a CF node Nr
     * @param height a height along that CF node's parent branch
     * @return the marginal Nr occupying the branch at that height during the most recent
     *         recalculate(), or -1 if none did (before the branch's first occupancy, or after a
     *         vacate with no later occupant, or the CF node Nr is unknown to this height).
     *         Entries are recorded in height order, so a node's chain is non-increasing in
     *         height; this walks it newest-to-oldest for the first entry at or below
     *         {@code height}, identical semantics to a backward scan over that node's history.
     */
    public int getOccupantAtHeight(int cfNodeNr, double height) {
        for (int i = occLastOfNode[cfNodeNr]; i >= 0; i = occPrev[i]) {
            if (occHeight[i] <= height)
                return occNr[i];
        }
        return -1;
    }

    /**
     * @param cfNodeNr a CF node Nr
     * @return the number of occupancy-history entries recorded for that CF node during the most
     *         recent recalculate(). Exposed for tests verifying the history is reset, not
     *         accumulated, across rebuilds.
     */
    public int getOccupancyHistoryLength(int cfNodeNr) {
        int len = 0;
        for (int i = occLastOfNode[cfNodeNr]; i >= 0; i = occPrev[i])
            len++;
        return len;
    }

    /**
     * Determine, without touching any MarginalNode, exactly which CF nodes and conversions will
     * produce a marginal coalescence in the upcoming rebuild.
     *
     * This mirrors the active-lineage transitions of the real walk below using plain booleans
     * (keep the two in sync). It exists because claims and releases cannot safely be interleaved
     * with the real walk: a conversion that suppresses a CF node's coalescence is processed
     * BEFORE that CF node's own event (conversions are inserted into the walk by height, and a
     * suppressing conversion is necessarily below the CF node it suppresses), so if release were
     * deferred until the CF node's event was reached, a transient third slot would be needed even
     * though the final tree never has more than n-1 internal nodes at once -- overflowing the
     * fixed 2n-1 node capacity. Precomputing the firing set lets recalculate() reconcile all
     * releases before any new claims, so capacity is never transiently exceeded.
     */
    private void computeFiringSets(List<Event> cfEvents, List<Conversion> convs,
            java.util.Set<Integer> firingCfNodes, java.util.Set<Integer> firingConvIds) {
        boolean[] active = new boolean[acg.getNodeCount()];
        int iConv = 0;
        for (int iEvent = 0; iEvent < cfEvents.size(); iEvent++) {
            Event event = cfEvents.get(iEvent);
            Node node = event.getNode();

            if (event.getType() == EventType.SAMPLE) {
                active[node.getNr()] = true;
            } else {
                Node left = node.getChild(0);
                Node right = node.getChild(1);
                if (active[left.getNr()] && active[right.getNr()]) {
                    firingCfNodes.add(node.getNr());
                    active[left.getNr()] = false;
                    active[right.getNr()] = false;
                    active[node.getNr()] = true;
                } else if (active[left.getNr()]) {
                    active[left.getNr()] = false;
                    active[node.getNr()] = true;
                } else if (active[right.getNr()]) {
                    active[right.getNr()] = false;
                    active[node.getNr()] = true;
                }
            }

            while (iConv < convs.size() &&
                    (event.node.isRoot() || isBelowNextCFEvent(convs.get(iConv), cfEvents.get(iEvent + 1)))) {
                Conversion conv = convs.get(iConv++);
                Node node1 = conv.getNode1();
                Node node2 = conv.getNode2();
                if (active[node1.getNr()] && active[node2.getNr()]) {
                    firingConvIds.add(conv.getID());
                    active[node1.getNr()] = false;
                    active[node2.getNr()] = false;
                    active[node2.getNr()] = true;
                } else if (active[node1.getNr()]) {
                    active[node1.getNr()] = false;
                    active[node2.getNr()] = true;
                }
            }
        }
    }

    /**
     * Snapshot this tree's live state into the stored* fields, called lazily by recalculate()
     * on the first call of a new proposal. m_nodes only ever changes as a result of
     * recalculate()'s own logic, so this state -- right before the first recalculate() of a
     * proposal -- is exactly the correct baseline for restore() to revert to.
     */
    private void snapshotThisProposal() {
        int n = m_nodes.length;
        if (storedHeight == null || storedHeight.length != n) {
            storedHeight = new double[n];
            storedTimeLength = new double[n];
            storedLastEventHeight = new double[n];
            storedLeftNr = new int[n];
            storedRightNr = new int[n];
        }

        for (int nr = 0; nr < n; nr++) {
            MarginalNode node = (MarginalNode) m_nodes[nr];
            storedHeight[nr] = node.getHeight();
            storedTimeLength[nr] = node.timeLength;
            storedLastEventHeight[nr] = node.lastEventHeight;
            storedLeftNr[nr] = node.isLeaf() ? -1 : node.getLeft().getNr();
            storedRightNr[nr] = node.isLeaf() ? -1 : node.getRight().getNr();
        }
        storedRootNr = root.getNr();
        storedPrevRootHeight = prevRootHeight;

        storedCfCoalescenceNr = cfCoalescenceNr.clone();
        storedConvMarginalNr = new HashMap<>(convMarginalNr);
        storedFreeNrStack = freeNrStack.clone();
        storedFreeNrTop = freeNrTop;
        storedCreatorKind = creatorKind.clone();
        storedCreatorKey = creatorKey.clone();
        storedCfBranchOccupant = cfBranchOccupant.clone();
        storedNextNonLeafNr = nextNonLeafNr;
        storedSeenCfVersion = seenCfVersion;
        storedSeenConvVersion = seenConvVersion;
        storedSeenBlockVersion = seenBlockVersion;
        storedMappingVersion = mappingVersion;
    }

    public void recalculate() {
        nRecalculate++;
        if (lowestChangedHeight < Double.POSITIVE_INFINITY)
            nRecalculateWithHeight++;
        else
            nRecalculateIneligible++;
        startEditing(null);

        if (!snapshotTakenThisProposal) {
            snapshotThisProposal();
            snapshotTakenThisProposal = true;
        }

        if (cfCoalescenceNr.length != acg.getNodeCount()) {
            throw new IllegalStateException("MarginalTree: CF node count changed between rebuilds ("
                    + cfCoalescenceNr.length + " -> " + acg.getNodeCount()
                    + "); the stable marginal-node numbering assumes it is fixed for the life of the run.");
        }

        List<Event> cfEvents = acg.getCFEvents();
        activeCFlineages.replaceAll((oldNode) -> {
            return null;
        });
        occCount = 0;
        Arrays.fill(occLastOfNode, -1);
        List<Conversion> convs = getBlockConversions();
        convs.sort((c1, c2) -> {
            if (c1.height < c2.height)
                return -1;

            if (c1.height > c2.height)
                return 1;

            return 0;
        });

        // Only a change confined to this block's move set is height-localizable. Validated here
        // rather than in checkOutdated() for two reasons: convs is already materialised, so the
        // scan is free; and the Gibbs operators call recalculate() DIRECTLY, bypassing
        // checkOutdated() entirely, so this is the only place guaranteed to see every rebuild.
        for (int i = 0; i < convs.size(); i++) {
            if (convs.get(i).getVersion() > seenConvVersion) {
                makeIneligibleForIncrementalUpdate();
                break;
            }
        }

        // Precompute which creators fire this round, then reconcile ALL releases before ANY
        // new claim -- see computeFiringSets() for why that order matters.
        java.util.Set<Integer> firingCfNodes = new java.util.HashSet<>();
        java.util.Set<Integer> firingConvIds = new java.util.HashSet<>();
        computeFiringSets(cfEvents, convs, firingCfNodes, firingConvIds);

        for (int cfNr = 0; cfNr < cfCoalescenceNr.length; cfNr++) {
            if (cfCoalescenceNr[cfNr] >= 0 && !firingCfNodes.contains(cfNr)) {
                freeNrStack[freeNrTop++] = cfCoalescenceNr[cfNr];
                cfCoalescenceNr[cfNr] = -1;
            }
        }
        convMarginalNr.entrySet().removeIf(entry -> {
            if (!firingConvIds.contains(entry.getKey())) {
                freeNrStack[freeNrTop++] = entry.getValue();
                return true;
            }
            return false;
        });
        for (int cfNr : firingCfNodes) {
            if (cfCoalescenceNr[cfNr] < 0) {
                cfCoalescenceNr[cfNr] = claimMarginalNr(CREATOR_CF_COALESCENCE, cfNr);
            }
        }
        for (int convId : firingConvIds) {
            convMarginalNr.computeIfAbsent(convId, id -> claimMarginalNr(CREATOR_CONVERSION, id));
        }

        int iConv = 0;

        for (int iEvent = 0; iEvent < cfEvents.size(); iEvent++) {
            Event event = cfEvents.get(iEvent);
            Node node = event.getNode();

            // Process the current CF-event
            switch (event.getType()) {
                case SAMPLE:
                    MarginalNode marginalLeaf = registerLeafNode(node);
                    setActiveLineage(node.getNr(), marginalLeaf, node.getHeight());
                    break;

                case COALESCENCE:
                    Node left = node.getChild(0);
                    Node right = node.getChild(1);

                    if (activeCFlineages.get(left.getNr()) != null && activeCFlineages.get(right.getNr()) != null) {
                        // Create a new marginal node at the coalescence event, using the stable
                        // marginal Nr already reconciled for this CF node above.
                        int nr = cfCoalescenceNr[node.getNr()];
                        MarginalNode marginalNode = registerNode(node, activeCFlineages, nr);

                        // Remove the old and add the new marginal node to the active lineages.
                        setActiveLineage(left.getNr(), null, node.getHeight());
                        setActiveLineage(right.getNr(), null, node.getHeight());
                        setActiveLineage(node.getNr(), marginalNode, node.getHeight());
                    } else {
                        // Only one side is active -> no coalescence in marginal tree (i.e. no marginal node)

                        if (activeCFlineages.get(left.getNr()) != null) {
                            MarginalNode marginalLeft = activeCFlineages.get(left.getNr());
                            updateTimeLength(node.getHeight(), left, marginalLeft);

                            setActiveLineage(left.getNr(), null, node.getHeight());
                            setActiveLineage(node.getNr(), marginalLeft, node.getHeight());
                            break;
                        }

                        if (activeCFlineages.get(right.getNr()) != null) {
                            MarginalNode marginalRight = activeCFlineages.get(right.getNr());
                            updateTimeLength(node.getHeight(), right, marginalRight);

                            setActiveLineage(right.getNr(), null, node.getHeight());
                            setActiveLineage(node.getNr(), marginalRight, node.getHeight());
                            break;
                        }
                    }
                    break;
            }

            // Process all conversion below the next CF-event
            while (iConv < convs.size() &&
                    (event.node.isRoot() || isBelowNextCFEvent(convs.get(iConv), cfEvents.get(iEvent + 1)))) {

                Conversion conv = convs.get(iConv++);
                Node node1 = conv.getNode1();
                Node node2 = conv.getNode2();

                if ((activeCFlineages.get(node1.getNr()) != null) && (activeCFlineages.get(node2.getNr()) != null)) {
                    // Both lineages at the conversion are active --> coalescence in the marginal tree

                    // Create a MarginalNode at the point of the conversion, using the stable
                    // marginal Nr already reconciled for this conversion above.
                    int nr = convMarginalNr.get(conv.getID());
                    MarginalNode convNode = registerNode(conv, activeCFlineages, nr);

                    // Remove child lineages and add new one
                    setActiveLineage(node1.getNr(), null, conv.getHeight());
                    setActiveLineage(node2.getNr(), null, conv.getHeight());
                    setActiveLineage(node2.getNr(), convNode, conv.getHeight());
                } else {
                    // node1 or node2 already moved away (overshadowed by another conversion)

                    if (activeCFlineages.get(node1.getNr()) != null) {
                        // node1 passes conversion, but node2 branched away --> CF lineage of node2 is continued by node1
                        MarginalNode margNode1 = activeCFlineages.get(node1.getNr());
                        updateTimeLength(conv.getHeight(), node1, margNode1);

                        setActiveLineage(node1.getNr(), null, conv.getHeight());
                        setActiveLineage(node2.getNr(), margNode1, conv.getHeight());
                    }
                    // else: node1 already branched away --> conversion has no effect
                }
            }
        }

        // A single active CF lineage (the root) should remain:;
        MarginalNode newRoot = activeCFlineages.get(acg.getRoot().getNr());
        // With sequential numbering the root always landed on the highest Nr, which was also
        // the slot assignFrom(acg)'s initial clone gave the CF root -- already parentless. Stable
        // numbering can now put the root on ANY internal slot, including one that still carries a
        // stale parent pointer from whatever node originally occupied it, so this must be cleared
        // explicitly: setRootOnly() (Tree.java) does not do it, and isRoot()/getLength() rely on it.
        newRoot.setParent(null);
        setRootOnly(newRoot);

        if (hasBranchRates) {
            rollOutTime(newRoot, newRoot.getHeight());
        }

        // Per-node dirty marking, keyed on what TreeLikelihood actually consumes: a node's
        // transition matrix depends on its own branch length, which in rolled-out coordinates
        // is timeLength (or, if there's no branch-rate model, the raw height directly -- either
        // way, "the branch determinant"). This runs once, after the whole walk (every node's
        // timeLength/height is only final at this point), and is keyed against the LAST
        // COMPUTED state, not the last accepted one: recalculate() legitimately runs several
        // times per proposal in the Gibbs operators, and each call must be judged against
        // whatever this same tree last actually looked like, not against some more distant
        // "committed" baseline.
        for (int nr = 0; nr < m_nodes.length; nr++) {
            MarginalNode node = (MarginalNode) m_nodes[nr];
            if (node == newRoot)
                continue;
            nNodesChecked++;
            double determinant = hasBranchRates ? node.timeLength : node.getHeight();
            if (Double.isNaN(node.prevBranchDeterminant) || determinant != node.prevBranchDeterminant) {
                nFilthyNodesMarked++;
                node.makeDirty(Tree.IS_FILTHY);
            }
            node.prevBranchDeterminant = determinant;
        }
        nNodesChecked++;
        if (Double.isNaN(prevRootHeight) || newRoot.getHeight() != prevRootHeight) {
            nFilthyNodesMarked++;
            newRoot.makeDirty(Tree.IS_FILTHY);
        }
        prevRootHeight = newRoot.getHeight();

        outdated = false;
        seenCfVersion = acg.getCFVersion();
        seenConvVersion = acg.getConvVersion();
        seenBlockVersion = block.getMoveSetVersion();
        mappingVersion++;

        // Some operators (the Gibbs-sampling borrowing code in BorrowingOperator) call
        // recalculate() directly, bypassing requiresRecalculation(). somethingIsDirty() is
        // the ONLY way TreeLikelihood learns this tree changed (its own requiresRecalculation()
        // reads treeInput.get().somethingIsDirty(), never re-derives it), so a direct call must
        // still flip it -- otherwise the change is invisible to TreeLikelihood and it silently
        // keeps a stale cached logP. Always safe: the next requiresRecalculation() call
        // overwrites this with a fresh checkOutdated() result regardless.
        changed = true;
//        block.updatedMarginalTree();
//        lastBlockState = block.toString();

    }

    /**
     * Decide whether a conversion should still be processed before the next CF event.
     *
     * A conversion exactly at the height of the next CF event only belongs below that event
     * if the event is the coalescence of one of the conversion's own lineages.
     */
    private boolean isBelowNextCFEvent(Conversion conv, Event nextEvent) {
        double convHeight = conv.getHeight();
        double nextHeight = nextEvent.getHeight();

        if (convHeight < nextHeight)
            return true;

        if (convHeight > nextHeight)
            return false;

        if (nextEvent.getType() != EventType.COALESCENCE)
            return false;

        Node nextNode = nextEvent.getNode();
        return conv.getNode1().getParent() == nextNode || conv.getNode2().getParent() == nextNode;
    }

    public MarginalNode registerLeafNode(Node node) {
        MarginalNode marginalNode = (MarginalNode) m_nodes[node.getNr()];
        // Dirtiness is decided by the post-walk pass in recalculate(), keyed on the branch
        // determinant (timeLength, accumulated below and elsewhere as this leaf's lineage
        // passes through pass-through CF events) -- not here, and not by height: once
        // hasBranchRates rolls heights out, a leaf's raw CF height is no longer what
        // TreeLikelihood's branch length actually depends on.
        marginalNode.setHeight(node.getHeight());
        marginalNode.lastEventHeight = node.getHeight();
        marginalNode.setID(node.getID());
        marginalNode.timeLength = 0;

        return marginalNode;
    }

    public MarginalNode registerNode(Node node, ArrayList<MarginalNode> activeCFlineages, int nodeNr) {
        Node left = node.getChild(0);
        Node right = node.getChild(1);
        double height = node.getHeight();
        MarginalNode newLeft = activeCFlineages.get(left.getNr());
        MarginalNode newRight = activeCFlineages.get(right.getNr());

        if (nodeNr >= m_nodes.length) {
            throw binaryMergeCapacityError(nodeNr, newLeft, newRight,
                "coalescence at node " + node.getNr() + " height=" + height);
        }

        // Take the old marginal node from the given index
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];
        Node oldLeft = marginalNode.getLeft();
        Node oldRight = marginalNode.getRight();

        // Mark the node as filthy if a child slot now holds a different creator than last time
        // (stable numbering, per-creator, makes this Nr comparison meaningful across rebuilds --
        // see the class-level comment on cfCoalescenceNr/convMarginalNr). The node's OWN branch
        // determinant (timeLength) is checked by the post-walk pass in recalculate(), not here.
        if (oldLeft.getNr() != newLeft.getNr() || oldRight.getNr() != newRight.getNr()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }

        // Update the marginal node and meta information
        marginalNode.update(height, newLeft, newRight);
        marginalNode.timeLength = 0;
        marginalNode.lastEventHeight = height;

        // Update time-length of both children
        if (hasBranchRates) {
            updateTimeLength(height, left, newLeft);
            updateTimeLength(height, right, newRight);
        }

        return marginalNode;
    }

    public MarginalNode registerNode(Conversion conv, ArrayList<MarginalNode> activeCFlineages, int nodeNr) {
        Node left = conv.getNode1();
        Node right = conv.getNode2();
        double height = conv.getHeight();
        MarginalNode marginalLeft = activeCFlineages.get(left.getNr());
        MarginalNode marginalRight = activeCFlineages.get(right.getNr());

        if (nodeNr >= m_nodes.length) {
            throw binaryMergeCapacityError(nodeNr, marginalLeft, marginalRight,
                "conversion " + conv.getID() + " at height=" + height
                    + " between nodes " + left.getNr() + " and " + right.getNr());
        }

        // Take the old marginal node from the given index
        MarginalNode marginalNode = (MarginalNode) m_nodes[nodeNr];

        // Mark the node as filthy if a child slot now holds a different creator than last time
        // (see the analogous comment in the CF-coalescence overload above). The node's OWN
        // branch determinant (timeLength) is checked by the post-walk pass in recalculate().
        if (marginalNode.getLeft().getNr() != marginalLeft.getNr()
                || marginalNode.getRight().getNr() != marginalRight.getNr()) {
            nFilthyNodesMarked++;
            marginalNode.makeDirty(Tree.IS_FILTHY);
        }

        marginalNode.update(height, marginalLeft, marginalRight);
        marginalNode.timeLength = 0;
        marginalNode.lastEventHeight = height;

        // Update time-length of both children
        if (hasBranchRates) {
            updateTimeLength(height, left, marginalLeft);
            updateTimeLength(height, right, marginalRight);
        }

        return marginalNode;
    }

    /**
     * Build a descriptive error for a merge that would write past the end of the node array,
     * instead of letting it surface as an ArrayIndexOutOfBoundsException. Callers check the
     * capacity themselves so this (and the context string) is only ever built on the failure
     * path, not on every registerNode() call.
     */
    private IllegalStateException binaryMergeCapacityError(int nodeNr, MarginalNode left, MarginalNode right, String context) {
        String leftId = left == null ? "null" : Integer.toString(left.getNr());
        String rightId = right == null ? "null" : Integer.toString(right.getNr());
        return new IllegalStateException("MarginalTree exceeded node capacity while registering " + context
                + "; nodeNr=" + nodeNr + ", capacity=" + m_nodes.length
                + ", leftMarginal=" + leftId + ", rightMarginal=" + rightId
                + ", sameLineage=" + (left != null && left == right));
    }

    void updateTimeLength(double parentHeight, Node child, MarginalNode marginalChild) {
        marginalChild.timeLength += (parentHeight - marginalChild.lastEventHeight) * branchRateModel.getRateForBranch(child);
        marginalChild.lastEventHeight = parentHeight;
    }

    /**
     * Assign new heights to the nodes in the tree by
     * rolling out the branch time-lengths from the
     * root down (starting at rootHeight).
     *
     * @param root
     * @param rootHeight
     */
    void rollOutTime(MarginalNode root, double rootHeight) {
        root.setHeight(rootHeight);

        List<Node> children = root.getChildren();
        if (children.size() > 0) {
            MarginalNode leftChild = (MarginalNode) children.get(0);
            rollOutTime(leftChild, rootHeight - leftChild.timeLength);
        }
        if (children.size() > 1) {
            MarginalNode rightChild = (MarginalNode) children.get(1);
            rollOutTime(rightChild, rootHeight - rightChild.timeLength);
        }

////        for (Node childAsNode : root.getChildren()) {
//        for (int i=0; i< root.getChildren().size(); i++) {
//            MarginalNode child = (MarginalNode) root.getChildren().get(i);
//            rollOutTime(child, rootHeight - child.timeLength);
//        }
    }

    @Override
    public void accept() {
        super.accept();
        setEverythingDirty(false);
        // Reset the sticky "changed" signal (see requiresRecalculation()) now that this
        // whole proposal cycle -- and any TreeLikelihood read of it -- is done.
        changed = false;
        // Next proposal's first recalculate() must snapshot again.
        snapshotTakenThisProposal = false;
        // Next proposal starts with nothing reported.
        lowestChangedHeight = Double.POSITIVE_INFINITY;
    }

    @Override
    public void store() {}

    /**
     * Restore this tree to the state it had right before the first recalculate() of the
     * rejected proposal, in place (no m_nodes/m_storedNodes swap -- see the plan file's "Real
     * MarginalTree.store()/restore()" section for why).
     *
     * A tree that saw no recalculate() at all this proposal needs no restoring and no rebuild:
     * see the untouched fast path below.
     */
    @Override
    public void restore() {
        nRestore++;

        if (!snapshotTakenThisProposal) {
            // Untouched this proposal: m_nodes changes ONLY via recalculate(), and none ran, so
            // this tree still matches the state its inputs have just been restored to (MCMC calls
            // state.restore() before state.restoreCalculationNodes(), so the ACG and Blocks are
            // already back to the accepted state). Nothing to restore, and nothing to rebuild.
            //
            // The stored arrays must NOT be applied here: they hold the snapshot of some EARLIER
            // proposal, and if that one was accepted, applying it would revert an accepted change.
            //
            // Re-baseline the seen-version stamps instead. The ACG bumps its version counters
            // conservatively in its own restore(), so leaving the old stamps in place would make
            // the next checkOutdated() demand a full rebuild of a tree already known to be
            // correct -- which is exactly the waste this fast path exists to remove. A later
            // genuine change still bumps the counters past these values, so it is still detected.
            nRestoreUntouched++;
            seenCfVersion = acg.getCFVersion();
            seenConvVersion = acg.getConvVersion();
            seenBlockVersion = block.getMoveSetVersion();
            outdated = false;
            manuallyUpdated = false;
            return;
        }

        int n = m_nodes.length;
        for (int nr = 0; nr < n; nr++) {
            MarginalNode node = (MarginalNode) m_nodes[nr];
            boolean isLeaf = node.isLeaf();
            boolean isRoot = (nr == storedRootNr);

            int preLeftNr = isLeaf ? -1 : node.getLeft().getNr();
            int preRightNr = isLeaf ? -1 : node.getRight().getNr();
            boolean childrenChanged = !isLeaf
                    && (preLeftNr != storedLeftNr[nr] || preRightNr != storedRightNr[nr]);

            boolean dirty;
            if (isRoot) {
                // The root's own branch determinant is meaningless (no parent branch); its
                // dirtiness is judged by height alone, mirroring recalculate()'s own
                // prevRootHeight handling.
                dirty = childrenChanged || node.getHeight() != storedPrevRootHeight;
            } else {
                double preDeterminant = hasBranchRates ? node.timeLength : node.getHeight();
                double restoredDeterminant = hasBranchRates ? storedTimeLength[nr] : storedHeight[nr];
                dirty = childrenChanged || preDeterminant != restoredDeterminant;
            }
            nNodesChecked++;
            if (dirty) {
                nFilthyNodesMarked++;
                node.makeDirty(Tree.IS_FILTHY);
            }

            if (isLeaf) {
                node.setHeight(storedHeight[nr]);
            } else {
                MarginalNode left = (MarginalNode) m_nodes[storedLeftNr[nr]];
                MarginalNode right = (MarginalNode) m_nodes[storedRightNr[nr]];
                node.update(storedHeight[nr], left, right);
            }
            node.timeLength = storedTimeLength[nr];
            node.lastEventHeight = storedLastEventHeight[nr];
            node.prevBranchDeterminant = hasBranchRates ? storedTimeLength[nr] : storedHeight[nr];
        }

        MarginalNode restoredRoot = (MarginalNode) m_nodes[storedRootNr];
        restoredRoot.setParent(null);
        setRootOnly(restoredRoot);
        prevRootHeight = storedPrevRootHeight;

        cfCoalescenceNr = storedCfCoalescenceNr;
        convMarginalNr = storedConvMarginalNr;
        freeNrStack = storedFreeNrStack;
        freeNrTop = storedFreeNrTop;
        creatorKind = storedCreatorKind;
        creatorKey = storedCreatorKey;
        cfBranchOccupant = storedCfBranchOccupant;
        nextNonLeafNr = storedNextNonLeafNr;
        seenCfVersion = storedSeenCfVersion;
        seenConvVersion = storedSeenConvVersion;
        seenBlockVersion = storedSeenBlockVersion;
        mappingVersion = storedMappingVersion;

        postCache = null;
        outdated = false;
        changed = true;
        manuallyUpdated = false;
        snapshotTakenThisProposal = false;
        // Next proposal starts with nothing reported.
        lowestChangedHeight = Double.POSITIVE_INFINITY;
    }

    /**
     * Rebuild this marginal tree against the current ACG and block state.
     *
     * This is NOT an undo: use it from operators that mutate the ACG mid-proposal and need the
     * derived tree brought up to date.
     */
    public void refresh() {
        postCache = null;

        manuallyUpdated = false;
        recalculate();
//        setEverythingDirty(false);
    }

    protected void initArraysSlim() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new MarginalNode[nodeCount];
        listMarginalNodes((MarginalNode) root, m_nodes);
        postCache = null;
    }

    /**
     * convert tree to array representation *
     */
    void listMarginalNodes(final MarginalNode node, final Node[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);

        for (final Node child : node.getChildren()) {
            listMarginalNodes((MarginalNode) child, nodes);
        }
    }

    ArrayList<Conversion> _blockConvs = new ArrayList<>();

    /**
     * Obtain the list of conversions which affect this block, i.e. the ones defining this marginal tree.
     *
     * @return List of relevant conversions.
     */
    public List<Conversion> getBlockConversions() {
        _blockConvs.clear();
        ConversionList convList = acg.getConversions();
        List<Integer> cIDs = block.getConversionIDs();

        for (int i=0; i < cIDs.size(); i++) {
            Conversion c = convList.get(cIDs.get(i));
            validateConversionRef(c);
            _blockConvs.add(c);
        }

        return _blockConvs;
    }

    /**
     * Reject conversions that a block refers to but that are missing or degenerate.
     */
    private void validateConversionRef(Conversion conversion) {
        if (conversion == null) {
            throw new IllegalStateException(
                    "MarginalTree encountered a missing conversion referenced by block " + block.getID());
        }

        if (!conversion.isValid()) {
            throw new IllegalStateException("MarginalTree encountered invalid conversion " + conversion.getID()
                    + " referenced by block " + block.getID());
        }
    }

    /*
     * Public CF/conversion <-> marginal-node mapping API.
     */

    /**
     * @param cfNodeNr a node Nr on the clonal frame
     * @return the marginal Nr that CF node currently maps to, or -1 if it has no marginal
     *         coalescence right now (a leaf always maps to the same Nr).
     */
    public int getMarginalNr(int cfNodeNr) {
        if (cfNodeNr < acg.getLeafNodeCount())
            return cfNodeNr;
        return cfCoalescenceNr[cfNodeNr];
    }

    /**
     * @param conv a conversion
     * @return the marginal Nr that conversion currently maps to, or -1 if it has no marginal
     *         coalescence right now (e.g. it is overshadowed, or not carried by this block).
     */
    public int getMarginalNr(Conversion conv) {
        Integer nr = convMarginalNr.get(conv.getID());
        return (nr == null) ? -1 : nr;
    }

    public MarginalNode getMarginalNode(int cfNodeNr) {
        int nr = getMarginalNr(cfNodeNr);
        return (nr < 0) ? null : (MarginalNode) m_nodes[nr];
    }

    public MarginalNode getMarginalNode(Conversion conv) {
        int nr = getMarginalNr(conv);
        return (nr < 0) ? null : (MarginalNode) m_nodes[nr];
    }

    /**
     * @param marginalNr a Nr in this marginal tree
     * @return {@link #CREATOR_LEAF}, {@link #CREATOR_CF_COALESCENCE} or {@link #CREATOR_CONVERSION}
     */
    public int getCreatorKind(int marginalNr) {
        return (marginalNr < acg.getLeafNodeCount()) ? CREATOR_LEAF : creatorKind[marginalNr];
    }

    /**
     * @param marginalNr a Nr whose creator kind is {@link #CREATOR_LEAF} or {@link #CREATOR_CF_COALESCENCE}
     * @return the CF node Nr that created it
     */
    public int getCreatorCFNodeNr(int marginalNr) {
        return (marginalNr < acg.getLeafNodeCount()) ? marginalNr : creatorKey[marginalNr];
    }

    /**
     * @param marginalNr a Nr whose creator kind is {@link #CREATOR_CONVERSION}
     * @return the conversion that created it, or null if that kind doesn't match
     */
    public Conversion getCreatorConversion(int marginalNr) {
        if (getCreatorKind(marginalNr) != CREATOR_CONVERSION)
            return null;
        int convId = creatorKey[marginalNr];
        for (Conversion conv : getBlockConversions()) {
            if (conv.getID() == convId)
                return conv;
        }
        return null;
    }

    /**
     * @return all CF node Nrs currently governed by the given marginal Nr (i.e. the marginal
     *         branch at marginalNr is the ancestor of these CF lineages), including any
     *         pass-through CF coalescences it covers.
     */
    public int[] collectCoveredCFNodes(int marginalNr) {
        int count = 0;
        for (int occupant : cfBranchOccupant) {
            if (occupant == marginalNr)
                count++;
        }
        int[] result = new int[count];
        int i = 0;
        for (int cfNr = 0; cfNr < cfBranchOccupant.length; cfNr++) {
            if (cfBranchOccupant[cfNr] == marginalNr)
                result[i++] = cfNr;
        }
        return result;
    }

    /** Bumped every recalculate(); lets external caches know the mapping may have changed. */
    public long getMappingVersion() {
        return mappingVersion;
    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Tree tree = (Tree) other;

        final MarginalNode[] nodes = new MarginalNode[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = newNode();
        }

        setID(tree.getID());

        root = (MarginalNode) nodes[tree.getRoot().getNr()];
        root.assignFrom(nodes, tree.getRoot());
        root.setParent(null);

        nodeCount = tree.getNodeCount();
        internalNodeCount = tree.getInternalNodeCount();
        leafNodeCount = tree.getLeafNodeCount();

        initArraysSlim();
    }

    @Override
    protected MarginalNode newNode() {
        return new MarginalNode();
    }

    /**
     * The marginal tree in Newick format.
     */
    @Override
    public String toString() {
        return root.toString();
    }


    /**
     * Checksum over everything {@link #restore()} is responsible for putting back, so that BEAST's
     * own reject validation ({@code MCMC.validateReject()} under {@code -Dbeast.debug=true})
     * independently catches an incomplete restore on every rejected proposal.
     *
     * This used to return 0 -- opting out of that check -- on the grounds that a marginal tree is
     * a lean calculation node whose intermediate state is arbitrary. That no longer holds: now
     * that restore() is a real undo rather than a full rebuild, the state it restores to is
     * exactly the state at the start of the proposal, which is precisely when
     * {@code CalculationNode.storeChecksum()} is called.
     *
     * Deliberately excluded: the seen-version stamps and {@link #mappingVersion}. The untouched
     * fast path in restore() re-baselines the stamps to the just-restored inputs' current values
     * rather than reverting them, by design -- including them here would report that as a failure.
     * Also excluded are the dirty/outdated flags, which describe freshness rather than tree state.
     */
    @Override
    public int getChecksum() {
        int checksum = 1;
        for (int nr = 0; nr < m_nodes.length; nr++) {
            MarginalNode node = (MarginalNode) m_nodes[nr];
            checksum = 31 * checksum + Double.hashCode(node.getHeight());
            checksum = 31 * checksum + Double.hashCode(node.timeLength);
            checksum = 31 * checksum + Double.hashCode(node.lastEventHeight);
            checksum = 31 * checksum + (node.isLeaf() ? -1 : node.getLeft().getNr());
            checksum = 31 * checksum + (node.isLeaf() ? -1 : node.getRight().getNr());
        }
        checksum = 31 * checksum + root.getNr();

        // The stable-numbering bookkeeping: a mapping restored incorrectly corrupts the NEXT
        // rebuild rather than this state, so no likelihood comparison would catch it.
        checksum = 31 * checksum + Arrays.hashCode(cfCoalescenceNr);
        checksum = 31 * checksum + Arrays.hashCode(cfBranchOccupant);
        checksum = 31 * checksum + Arrays.hashCode(creatorKind);
        checksum = 31 * checksum + Arrays.hashCode(creatorKey);
        checksum = 31 * checksum + nextNonLeafNr;
        // Only the live prefix of the free stack: entries above freeNrTop are dead leftovers.
        checksum = 31 * checksum + freeNrTop;
        for (int i = 0; i < freeNrTop; i++)
            checksum = 31 * checksum + freeNrStack[i];
        // Map.hashCode() is the sum of entry hashes, so iteration order is irrelevant.
        checksum = 31 * checksum + convMarginalNr.hashCode();

        return checksum;
    }

    /**
     * Skip the checksum comparison in the one case where BEAST never gave this tree a chance to
     * restore itself.
     *
     * On the DIRECT-reject path (the operator returned -infinity), {@code MCMC.propagateState()}
     * calls {@code state.restore()} but then guards {@code state.restoreCalculationNodes()} behind
     * {@code !operator.requiresStateInitialisation()} (MCMC.java:574). That method defaults to
     * {@code true} in {@code Operator} and no contacTrees operator overrides it, so calculation
     * nodes are never restored on that path: this tree is deliberately left holding the rejected
     * proposal's state. It self-heals before next use, because {@code ConversionGraph.restore()}
     * conservatively bumps the version counters past the stamps recorded during the proposal,
     * which makes {@link #checkOutdated()} force a rebuild.
     *
     * {@code snapshotTakenThisProposal} is exactly the signal for that case: a rebuild happened
     * this proposal and neither {@link #restore()} nor {@link #accept()} has cleared it since.
     * Whenever restore() genuinely ran, the flag is clear and the comparison goes ahead.
     */
    @Override
    public boolean matchesOldChecksum() {
        if (snapshotTakenThisProposal)
            return true;
        return super.matchesOldChecksum();
    }

}

