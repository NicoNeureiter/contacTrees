package contactrees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import contactrees.CFEventList.Event;
import contactrees.CFEventList.EventType;


/**
 * The marginal tree of one block: the tree that the block's sites evolve along, given the
 * clonal frame (CF) of the conversion graph and the conversions the block moves over.
 *
 * <h2>Creators and stable node numbers</h2>
 * Every internal marginal node is produced by exactly one <i>creator</i>: either a CF
 * coalescence whose two child branches both carry the block, or a conversion the block
 * moves over whose two endpoint branches are both occupied at its height. Leaves are
 * identified with the CF leaves (same node numbers). A creator keeps its marginal node
 * number for as long as it keeps firing, so cached partials of untouched subtrees stay valid
 * across updates. The forward mapping (creator to node number) and the inverse mapping are kept
 * here, per marginal tree.
 *
 * <h2>Heights and branch rates</h2>
 * Marginal node heights are the creator heights (no rolling out of branch rates). Branch
 * rates are supplied to the likelihood by {@link contactrees.model.MarginalClockModel}.
 * For a clock with branch-dependent rates, the tree also records the <i>path</i> of every
 * marginal branch: the CF branches it runs along (bottom to top) and the events where it
 * changes CF branch (a CF node it passes through, or a conversion over which it jumps to an
 * empty branch). Events are stored by identity (CF node nr or conversion ID), not by time, so
 * height-only updates leave the paths valid. The rate of a marginal branch is the time-weighted
 * mean of the CF rates along its path; rates are recomputed from the paths whenever the
 * conversion graph, the block or the clock changed, and a node whose rate changed is flagged dirty.
 *
 * <h2>Updates and caching</h2>
 * The structure is (re)derived by {@link #resweep()}: a sweep over the CF events and the block's
 * conversions computes the firing creators, which are then reconciled with the current tree.
 * Only node numbers whose children or height changed are written and flagged dirty. The first
 * modification within an MCMC proposal takes a snapshot, which {@link #restore()} reinstates
 * exactly; {@link #store()} is a no-op because BEAST calls it after the proposal.
 * Trees notice changes of the conversion graph or block through version numbers
 * ({@link ConversionGraph#getVersion()}, {@link Block#getVersion()}); a tree that does not
 * reflect the current versions resweeps itself on demand ({@link #ensureUpToDate()}).
 *
 * @author Nico Neureiter
 */
public class MarginalTree extends Tree {

    public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph (network) from which the marginal tree is derived.",
            Input.Validate.REQUIRED);
    public Input<Block> blockInput = new Input<>(
            "block",
            "The block object containing the moves this marginal tree follows along the conversion graph.",
            Input.Validate.REQUIRED);
    public Input<BranchRateModel> branchRateModelInput = new Input<>(
            "branchRateModel",
            "No longer supported: marginal trees keep their true heights. Put the clock model "
            + "on the CTreeLikelihood instead.");
    public Input<Boolean> validateInput = new Input<>(
            "validate",
            "Check the marginal tree against a from-scratch sweep whenever it is used (slow, for debugging).",
            false);

    /** Creator kinds stored in {@link #creatorKind}. */
    public static final byte CREATOR_NONE = 0;
    public static final byte CREATOR_LEAF = 1;
    public static final byte CREATOR_CF = 2;
    public static final byte CREATOR_CONVERSION = 3;

    public ConversionGraph acg;
    public Block block;
    protected boolean validate;

    protected int nLeaves;
    protected int nCFNodes;

    /* Forward and inverse creator mapping. */
    /** CF node nr -> marginal node nr (leaves map to themselves, -1 if the CF node does not coalesce here). */
    protected int[] cfToMarg;
    /** Conversion ID -> marginal node nr, for conversions that create a node in this tree. */
    protected HashMap<Integer, Integer> convToMarg = new HashMap<>();
    /** Marginal node nr -> creator kind / creator key (CF node nr or conversion ID). */
    protected byte[] creatorKind;
    protected int[] creatorKey;
    /** Internal node numbers released during a reconcile and not yet reclaimed. */
    protected int[] freeStack;
    protected int freeTop;

    /** Versions of the conversion graph and block that this tree currently reflects. */
    protected long reflectedAcgVersion = -1, reflectedBlockVersion = -1;

    /* Proposal bookkeeping. */
    /** Structure or mapping changed in the current proposal (a snapshot exists). */
    protected boolean touched = false;
    /** The reflected versions changed in the current proposal (their old values are saved). */
    protected boolean versionsTouched = false;
    /** All nodes are flagged dirty (initialisation or setEverythingDirty(true)). */
    protected boolean allDirty = false;
    /** Suppresses snapshots and proposal registration while initialising. */
    protected boolean initialising = false;
    /**
     * True once BEAST stored this tree as a calculation node, i.e. the tree is on the path to the
     * posterior and BEAST accepts/restores it (and recomputes its likelihood) after each proposal.
     * Unmanaged trees (e.g. used only by Gibbs operators) never clear their dirty flags, so their
     * likelihood recomputes the right nodes when it is next evaluated.
     */
    protected boolean managed = false;
    /** Node flags were set and not yet cleared (beyond the current proposal). */
    protected boolean flagsPending = false;
    /** Registered in the conversion graph's list of trees modified in the current proposal. */
    boolean inModifiedList = false;

    /* Snapshot of the pre-proposal state. */
    protected int[] sParent, sLeft, sRight, sCreatorKey, sCfToMarg;
    protected byte[] sCreatorKind;
    protected double[] sHeight;
    protected int sRoot;
    protected long sAcgVersion, sBlockVersion;

    /* Scratch space for the sweep. One record per firing creator, in firing order. */
    protected int[] activeLineage;
    protected byte[] recKind;
    protected int[] recKey, recChildA, recChildB, recNr;
    protected double[] recHeight;
    protected int nRecords;
    protected boolean[] keepNr;
    protected final ArrayList<Conversion> sortedConvs = new ArrayList<>();

    /* Branch paths and rates, only used with a branch-dependent clock (see setBranchRateClock). */
    /** The clock on the clonal frame, or null if the clock is branch independent. */
    protected BranchRateModel rateClock;
    /** The clock changed since the rates were last computed. */
    protected boolean clockStale = false;
    /** Per node nr: the CF branches (CF node nrs) occupied by the marginal branch above it, bottom to top. */
    protected int[][] pathBranches;
    /** Per node nr: the events where that marginal branch changes CF branch (creatorId-encoded). */
    protected long[][] pathEvents;
    /** Per node nr: the rate of the marginal branch above it. */
    protected double[] rates;
    protected int[][] sPathBranches;
    protected long[][] sPathEvents;
    protected double[] sRates;
    /* Sweep scratch for the paths: per lineage a linked list of events in an event pool. */
    protected int[] linStart, linHead, linTail;
    protected long[] evId = new long[0];
    protected int[] evBranch = new int[0], evNext = new int[0];
    protected int nEvents;

    // Partial sweeps: lineages swept (paths to check), departures per CF branch below the frontier
    protected int[] sweptLineages;
    protected int nSwept, nSweptLineagesStarted, rootLineage;
    protected int[] depHead, depTail, depNext = new int[0];

    /** Counters for diagnostics and tests. */
    public long nResweeps = 0, nPartialResweeps = 0, nNodesWritten = 0, nRateUpdates = 0;

    @Override
    public void initAndValidate() {
        acg = networkInput.get();
        block = blockInput.get();
        validate = validateInput.get();

        if (branchRateModelInput.get() != null) {
            throw new IllegalArgumentException("MarginalTree " + getID() + ": the branchRateModel input "
                    + "is no longer supported. Marginal trees keep their true heights; put the clock "
                    + "model on the CTreeLikelihood (it is wrapped in a MarginalClockModel).");
        }

        // Initialise the node objects as a copy of the clonal frame
        String beastID = ID;
        assignFrom(acg);
        setID(beastID);

        super.initAndValidate();

        nLeaves = acg.getLeafNodeCount();
        nCFNodes = acg.getNodeCount();
        if (getNodeCount() != nCFNodes)
            throw new IllegalStateException("MarginalTree " + getID() + ": node count does not match the conversion graph.");

        initMapping();

        acg.registerMarginalTree(this);
        block.registerMarginalTree(this);

        initialising = true;
        try {
            resweep();
        } finally {
            initialising = false;
        }
        touched = false;
        versionsTouched = false;
        makeAllDirty();
    }

    protected void initMapping() {
        int nNodes = getNodeCount();
        cfToMarg = new int[nCFNodes];
        Arrays.fill(cfToMarg, -1);
        creatorKind = new byte[nNodes];
        creatorKey = new int[nNodes];
        for (int i = 0; i < nLeaves; i++) {
            cfToMarg[i] = i;
            creatorKind[i] = CREATOR_LEAF;
            creatorKey[i] = i;
        }
        convToMarg.clear();

        // All internal node numbers start out free
        freeStack = new int[Math.max(nNodes - nLeaves, 0)];
        freeTop = 0;
        for (int nr = nNodes - 1; nr >= nLeaves; nr--)
            freeStack[freeTop++] = nr;

        sParent = new int[nNodes];
        sLeft = new int[nNodes];
        sRight = new int[nNodes];
        sHeight = new double[nNodes];
        sCreatorKind = new byte[nNodes];
        sCreatorKey = new int[nNodes];
        sCfToMarg = new int[nCFNodes];

        int nInternal = Math.max(nLeaves - 1, 0);
        activeLineage = new int[nCFNodes];
        recKind = new byte[nInternal];
        recKey = new int[nInternal];
        recChildA = new int[nInternal];
        recChildB = new int[nInternal];
        recNr = new int[nInternal];
        recHeight = new double[nInternal];
        keepNr = new boolean[nNodes];
        sweptLineages = new int[nCFNodes];
        depHead = new int[nCFNodes];
        depTail = new int[nCFNodes];
    }

    /*
     * Version tracking and lazy fallback
     */

    /**
     * @return true iff this tree reflects the current versions of the conversion graph and the block.
     */
    public boolean isUpToDate() {
        return reflectedAcgVersion == acg.getVersion() && reflectedBlockVersion == block.getVersion();
    }

    /**
     * Bring the tree up to date with the conversion graph and the block, resweeping if necessary.
     */
    public void ensureUpToDate() {
        if (!isUpToDate())
            resweep(lowestChange());
        else if (clockStale)
            updateRates(Double.NEGATIVE_INFINITY);
    }

    /**
     * @return a height below which the conversion graph and the block are known to be unchanged
     * since the state this tree reflects (-infinity if unknown).
     */
    protected double lowestChange() {
        if (acg.getVersion() != reflectedAcgVersion)
            return Double.NEGATIVE_INFINITY;
        return block.getLowestChangeSince(reflectedBlockVersion);
    }

    /**
     * Declare that this tree reflects the current conversion graph and block, without changing it.
     * Operators call this for trees whose structure they know to be unaffected by their changes.
     */
    public void markUpToDate() {
        markUpToDate(Double.NEGATIVE_INFINITY);
    }

    /**
     * Like {@link #markUpToDate()}, for changes that left everything below the given height
     * unchanged (only the rates of branches reaching that height are recomputed).
     */
    public void markUpToDate(double lowestChange) {
        boolean changed = setReflectedVersions();
        // With a branch-dependent clock, rates depend on heights along the paths (e.g. of CF
        // nodes the block only passes through), so they are recomputed after any change.
        if (rateClock != null && (changed || clockStale))
            updateRates(clockStale ? Double.NEGATIVE_INFINITY : lowestChange);
    }

    /**
     * Declare that the changes since the state this tree reflects did not touch anything the tree
     * depends on (e.g. only conversions its block does not follow): structure and rates stay.
     */
    public void markUnaffected() {
        markUpToDate(Double.POSITIVE_INFINITY);
    }

    /**
     * @return true iff the reflected versions changed.
     */
    private boolean setReflectedVersions() {
        long acgVersion = acg.getVersion();
        long blockVersion = block.getVersion();
        if (reflectedAcgVersion == acgVersion && reflectedBlockVersion == blockVersion)
            return false;
        beginVersionChange();
        reflectedAcgVersion = acgVersion;
        reflectedBlockVersion = blockVersion;
        return true;
    }

    /*
     * The sweep
     */

    /**
     * Recompute the structure from the conversion graph and the block, and reconcile it with the
     * current tree: creators keep their node numbers, only changed nodes are written and flagged.
     */
    public void resweep() {
        resweep(Double.NEGATIVE_INFINITY);
    }

    /**
     * Resweep, given that the conversion graph and the block are unchanged below the height
     * lowestChange (compared to the state this tree reflects). Only the part of the tree above
     * that height is swept again.
     */
    protected void resweep(double lowestChange) {
        nResweeps++;
        double from = lowestChange;
        if (from == Double.NEGATIVE_INFINITY || !computeStructure(from)) {
            from = Double.NEGATIVE_INFINITY;
            computeStructure(from);
        } else {
            nPartialResweeps++;
        }
        reconcile(from);
        setReflectedVersions();
        if (rateClock != null)
            updateRates(clockStale ? Double.NEGATIVE_INFINITY : from);
        if (validate && from > Double.NEGATIVE_INFINITY)
            checkConsistency("partial resweep from " + from);
    }

    /**
     * Full sweep (see {@link #computeStructure(double)}).
     */
    protected void computeStructure() {
        computeStructure(Double.NEGATIVE_INFINITY);
    }

    /**
     * Sweep through the CF events and the block's conversions (from the leaves to the root) and
     * record every firing creator in the scratch arrays. Lineages are identified by the node at
     * their bottom: leaf nr for leaves, nLeaves + record index for new internal nodes, and
     * nCFNodes + node nr for existing nodes on the frontier of a partial sweep.
     *
     * A partial sweep (from > -infinity) keeps the current tree below the height {@code from}:
     * it starts from the marginal branches crossing that height and records only the creators
     * at or above it. This requires the conversion graph and the block to be unchanged below
     * that height.
     *
     * @return false if a partial sweep is not possible (exact ties around the frontier)
     */
    protected boolean computeStructure(double from) {
        List<Event> cfEvents = acg.getCFEvents();
        List<Conversion> convs = getSortedBlockConversions();

        Arrays.fill(activeLineage, -1);
        nRecords = 0;
        nEvents = 0;
        nSwept = 0;
        nSweptLineagesStarted = 0;
        rootLineage = -1;
        boolean paths = (rateClock != null);
        int iConv = 0;
        int iEvent = 0;

        if (from > Double.NEGATIVE_INFINITY) {
            while (iConv < convs.size() && convs.get(iConv).getHeight() < from)
                iConv++;
            iEvent = firstEventAtOrAbove(cfEvents, from);
            if (!startFromFrontier(convs, iConv, from, paths))
                return false;
            // Conversions at or above the frontier that come before the first CF event swept
            if (iEvent < cfEvents.size())
                iConv = sweepConversions(convs, iConv, cfEvents.get(iEvent), paths);
        }

        for (; iEvent < cfEvents.size(); iEvent++) {
            Event event = cfEvents.get(iEvent);
            Node node = event.getNode();
            int nodeNr = node.getNr();

            if (event.getType() == EventType.SAMPLE) {
                activeLineage[nodeNr] = nodeNr;
                lineageStarted(nodeNr);
                if (paths)
                    startPath(nodeNr, nodeNr);
            } else {
                int left = node.getChild(0).getNr();
                int right = node.getChild(1).getNr();
                int lineageLeft = activeLineage[left];
                int lineageRight = activeLineage[right];

                if (lineageLeft >= 0 && lineageRight >= 0) {
                    // Both children carry the block: the CF node is a marginal coalescence
                    int lineage = addRecord(CREATOR_CF, nodeNr, node.getHeight(), lineageLeft, lineageRight);
                    activeLineage[left] = -1;
                    activeLineage[right] = -1;
                    activeLineage[nodeNr] = lineage;
                    lineageStarted(lineage);
                    if (paths)
                        startPath(lineage, nodeNr);
                } else if (lineageLeft >= 0) {
                    activeLineage[left] = -1;
                    activeLineage[nodeNr] = lineageLeft;
                    if (paths)
                        addPathEvent(lineageLeft, creatorId(CREATOR_CF, nodeNr), nodeNr);
                } else if (lineageRight >= 0) {
                    activeLineage[right] = -1;
                    activeLineage[nodeNr] = lineageRight;
                    if (paths)
                        addPathEvent(lineageRight, creatorId(CREATOR_CF, nodeNr), nodeNr);
                }
            }

            // Process all conversions below the next CF event
            iConv = sweepConversions(convs, iConv, node.isRoot() ? null : cfEvents.get(iEvent + 1), paths);
        }

        if (from == Double.NEGATIVE_INFINITY) {
            if (nRecords != Math.max(nLeaves - 1, 0)) {
                throw new IllegalStateException("MarginalTree " + getID() + ": sweep produced " + nRecords
                        + " internal nodes, expected " + (nLeaves - 1) + ".");
            }
        } else if (nRecords != nSweptLineagesStarted - nRecords - 1) {
            // Lineages started: frontier + leaves above it + one per record; all merge into one
            return false;
        }
        rootLineage = (nRecords > 0) ? nLeaves + nRecords - 1 : rootLineage;
        return true;
    }

    /**
     * Process the conversions from index iConv on that come before the given CF event (all
     * remaining ones if it is null).
     *
     * @return the index of the first conversion not processed
     */
    private int sweepConversions(List<Conversion> convs, int iConv, Event nextEvent, boolean paths) {
        while (iConv < convs.size() && (nextEvent == null || isBelowNextCFEvent(convs.get(iConv), nextEvent))) {
            Conversion conv = convs.get(iConv++);
            int nr1 = conv.getNode1().getNr();
            int nr2 = conv.getNode2().getNr();
            int lineage1 = activeLineage[nr1];
            int lineage2 = activeLineage[nr2];

            if (lineage1 >= 0 && lineage2 >= 0) {
                // Both lineages are active: coalescence at the conversion, continuing on node2
                int lineage = addRecord(CREATOR_CONVERSION, conv.getID(), conv.getHeight(), lineage1, lineage2);
                activeLineage[nr1] = -1;
                activeLineage[nr2] = lineage;
                lineageStarted(lineage);
                if (paths)
                    startPath(lineage, nr2);
            } else if (lineage1 >= 0) {
                // The lineage of node1 jumps to the (empty) branch of node2
                activeLineage[nr1] = -1;
                activeLineage[nr2] = lineage1;
                if (paths)
                    addPathEvent(lineage1, creatorId(CREATOR_CONVERSION, conv.getID()), nr2);
            }
            // else: node1 already left this branch, the conversion has no effect
        }
        return iConv;
    }

    private void lineageStarted(int lineage) {
        sweptLineages[nSwept++] = lineage;
        nSweptLineagesStarted++;
        rootLineage = lineage;
    }

    /**
     * @return the index of the first CF event at or above the given height
     */
    private static int firstEventAtOrAbove(List<Event> cfEvents, double height) {
        int lo = 0, hi = cfEvents.size();
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (cfEvents.get(mid).getHeight() < height)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo;
    }

    /**
     * Start a partial sweep: activate the lineages of the marginal branches that cross the height
     * {@code from} on the CF branches they occupy there. These are found by following each branch
     * up from its bottom node: below its parent, the lineage passes through every CF node it meets
     * and jumps at every conversion of the block that departs from its current CF branch (if
     * either was a coalescence, it would be the parent). Exact ties make the order of events
     * ambiguous; then the partial sweep is abandoned.
     *
     * @param iFrom index of the first conversion at or above {@code from} in the sorted list
     * @return false if a partial sweep is not possible
     */
    private boolean startFromFrontier(List<Conversion> convs, int iFrom, double from, boolean paths) {
        // The block's conversions below the frontier, grouped by the CF branch they depart from
        if (depNext.length < convs.size())
            depNext = new int[Math.max(convs.size(), 2 * depNext.length)];
        Arrays.fill(depHead, -1);
        for (int k = 0; k < iFrom; k++) {
            int b = convs.get(k).getNode1().getNr();
            depNext[k] = -1;
            if (depHead[b] < 0)
                depHead[b] = k;
            else
                depNext[depTail[b]] = k;
            depTail[b] = k;
        }

        ConversionList allConvs = acg.getConversions();
        int nNodes = getNodeCount();
        for (int nr = 0; nr < nNodes; nr++) {
            Node node = m_nodes[nr];
            double t = node.getHeight();
            Node parent = node.getParent();
            if (!(t < from) || (parent != null && parent.getHeight() < from))
                continue;
            if (creatorKind[nr] == CREATOR_NONE)
                return false;

            int lineage = (nr < nLeaves) ? nr : nCFNodes + nr;
            int b;
            if (creatorKind[nr] == CREATOR_CONVERSION) {
                Conversion creator = allConvs.get(creatorKey[nr]);
                if (creator == null)
                    return false;
                b = creator.getNode2().getNr();
            } else {
                b = creatorKey[nr];
            }
            if (paths)
                startPath(lineage, b);

            int dep = depHead[b];
            while (true) {
                while (dep >= 0 && convs.get(dep).getHeight() < t)
                    dep = depNext[dep];
                double tDep = (dep >= 0) ? convs.get(dep).getHeight() : Double.POSITIVE_INFINITY;
                Node cfNode = acg.getNode(b);
                double tTop = cfNode.isRoot() ? Double.POSITIVE_INFINITY : cfNode.getParent().getHeight();
                double next = Math.min(tDep, tTop);
                if (next >= from)
                    break;
                if (tDep == t || tDep == tTop)
                    return false;

                long eventId;
                if (tDep < tTop) {
                    Conversion conv = convs.get(dep);
                    eventId = creatorId(CREATOR_CONVERSION, conv.getID());
                    b = conv.getNode2().getNr();
                } else {
                    b = cfNode.getParent().getNr();
                    eventId = creatorId(CREATOR_CF, b);
                }
                t = next;
                if (paths)
                    addPathEvent(lineage, eventId, b);
                dep = depHead[b];
            }

            if (activeLineage[b] >= 0)
                return false;
            activeLineage[b] = lineage;
            lineageStarted(lineage);
        }
        return true;
    }

    private void startPath(int lineage, int branch) {
        linStart[lineage] = branch;
        linHead[lineage] = -1;
        linTail[lineage] = -1;
    }

    private void addPathEvent(int lineage, long eventId, int newBranch) {
        if (nEvents >= evId.length) {
            int capacity = Math.max(2 * evId.length, nCFNodes + 16);
            evId = Arrays.copyOf(evId, capacity);
            evBranch = Arrays.copyOf(evBranch, capacity);
            evNext = Arrays.copyOf(evNext, capacity);
        }
        int e = nEvents++;
        evId[e] = eventId;
        evBranch[e] = newBranch;
        evNext[e] = -1;
        if (linTail[lineage] < 0)
            linHead[lineage] = e;
        else
            evNext[linTail[lineage]] = e;
        linTail[lineage] = e;
    }

    private int pathEventCount(int lineage) {
        int n = 0;
        for (int e = linHead[lineage]; e >= 0; e = evNext[e])
            n++;
        return n;
    }

    /**
     * @return true iff the stored path matches the swept path of the given lineage.
     */
    private boolean pathMatches(int[] branches, long[] events, int lineage) {
        if (branches == null || branches[0] != linStart[lineage] || events.length != pathEventCount(lineage))
            return false;
        int i = 0;
        for (int e = linHead[lineage]; e >= 0; e = evNext[e], i++) {
            if (events[i] != evId[e] || branches[i + 1] != evBranch[e])
                return false;
        }
        return true;
    }

    private int addRecord(byte kind, int key, double height, int childA, int childB) {
        if (nRecords >= recKind.length) {
            throw new IllegalStateException("MarginalTree " + getID() + " exceeded its node capacity while "
                    + "registering " + (kind == CREATOR_CF ? "CF node " : "conversion ") + key
                    + " at height " + height + ".");
        }
        int i = nRecords++;
        recKind[i] = kind;
        recKey[i] = key;
        recHeight[i] = height;
        recChildA[i] = childA;
        recChildB[i] = childB;
        return nLeaves + i;
    }

    /**
     * Decide whether a conversion should still be processed before the next CF event.
     *
     * A conversion exactly at the height of the next CF event only belongs below that event
     * if the event is the coalescence of one of the conversion's own lineages.
     */
    private static boolean isBelowNextCFEvent(Conversion conv, Event nextEvent) {
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

    /**
     * The block's conversions sorted by height (ties broken by conversion ID, so that the result
     * does not depend on the order of moves in the block).
     */
    protected List<Conversion> getSortedBlockConversions() {
        sortedConvs.clear();
        ConversionList convList = acg.getConversions();
        List<Integer> cIDs = block.getConversionIDs();
        for (int i = 0; i < cIDs.size(); i++) {
            Conversion c = convList.get(cIDs.get(i));
            validateConversionRef(c);
            sortedConvs.add(c);
        }
        sortedConvs.sort((c1, c2) -> {
            int cmp = Double.compare(c1.getHeight(), c2.getHeight());
            return (cmp != 0) ? cmp : Integer.compare(c1.getID(), c2.getID());
        });
        return sortedConvs;
    }

    /**
     * Obtain the list of conversions which affect this block, i.e. the ones defining this marginal tree.
     *
     * @return List of relevant conversions (sorted by height).
     */
    public List<Conversion> getBlockConversions() {
        return new ArrayList<>(getSortedBlockConversions());
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
     * Reconciling the swept structure with the current tree
     */

    private int lookupNr(byte kind, int key) {
        if (kind == CREATOR_CF)
            return cfToMarg[key];
        Integer nr = convToMarg.get(key);
        return (nr == null) ? -1 : nr;
    }

    private int lineageToNr(int lineage) {
        if (lineage < nLeaves)
            return lineage;
        if (lineage >= nCFNodes)
            return lineage - nCFNodes;
        return recNr[lineage - nLeaves];
    }

    /**
     * Apply the structure in the scratch records to the tree (see {@link #reconcile(double)}).
     */
    protected void reconcile() {
        reconcile(Double.NEGATIVE_INFINITY);
    }

    /**
     * Apply the structure in the scratch records to the tree: release creators that no longer fire,
     * claim node numbers for new creators, then write (and flag) only what changed. After a
     * partial sweep from the given height, the nodes below it are kept as they are.
     */
    protected void reconcile(double from) {
        int nNodes = getNodeCount();

        // 1. Look up existing node numbers of the firing creators
        Arrays.fill(keepNr, false);
        for (int i = 0; i < nRecords; i++) {
            int nr = lookupNr(recKind[i], recKey[i]);
            recNr[i] = nr;
            if (nr >= 0)
                keepNr[nr] = true;
        }

        // 2. Release all internal nodes whose creator stopped firing
        for (int nr = nLeaves; nr < nNodes; nr++) {
            if (creatorKind[nr] != CREATOR_NONE && !keepNr[nr] && !(m_nodes[nr].getHeight() < from)) {
                beginModification();
                if (creatorKind[nr] == CREATOR_CF)
                    cfToMarg[creatorKey[nr]] = -1;
                else
                    convToMarg.remove(creatorKey[nr]);
                creatorKind[nr] = CREATOR_NONE;
                freeStack[freeTop++] = nr;
            }
        }

        // 3. Claim node numbers for new creators (CF creators prefer their own node number)
        for (int i = 0; i < nRecords; i++) {
            if (recNr[i] >= 0)
                continue;
            beginModification();
            int nr = claimFreeNr(recKind[i] == CREATOR_CF ? recKey[i] : -1);
            recNr[i] = nr;
            creatorKind[nr] = recKind[i];
            creatorKey[nr] = recKey[i];
            if (recKind[i] == CREATOR_CF)
                cfToMarg[recKey[i]] = nr;
            else
                convToMarg.put(recKey[i], nr);
            keepNr[nr] = false;   // marks "newly claimed" for step 4
        }
        assert freeTop == 0 : "Free node numbers left after reconcile";

        // 4. Leaves: heights follow the CF leaves (tip dates)
        for (int i = 0; i < nLeaves; i++) {
            Node leaf = m_nodes[i];
            double h = acg.getNode(i).getHeight();
            if (leaf.getHeight() != h) {
                beginModification();
                ((MarginalNode) leaf).setHeight(h);
                leaf.makeDirty(Tree.IS_DIRTY);
                nNodesWritten++;
            }
        }

        // 5. Internal nodes: write children and heights where they differ
        for (int i = 0; i < nRecords; i++) {
            MarginalNode node = (MarginalNode) m_nodes[recNr[i]];
            MarginalNode childA = (MarginalNode) m_nodes[lineageToNr(recChildA[i])];
            MarginalNode childB = (MarginalNode) m_nodes[lineageToNr(recChildB[i])];
            boolean isNew = !keepNr[recNr[i]];

            Node oldLeft = node.getLeft();
            Node oldRight = node.getRight();
            boolean sameChildren = (oldLeft == childA && oldRight == childB)
                    || (oldLeft == childB && oldRight == childA);
            boolean sameHeight = node.getHeight() == recHeight[i];

            if (isNew || !sameChildren) {
                beginModification();
                if (!(oldLeft == childB && oldRight == childA))
                    node.update(recHeight[i], childA, childB);
                node.setHeight(recHeight[i]);
                childA.setParent(node);
                childB.setParent(node);
                node.makeDirty(Tree.IS_FILTHY);
                childA.makeDirty(Tree.IS_DIRTY);
                childB.makeDirty(Tree.IS_DIRTY);
                nNodesWritten++;
            } else if (!sameHeight) {
                beginModification();
                node.setHeight(recHeight[i]);
                node.makeDirty(Tree.IS_DIRTY);
                childA.makeDirty(Tree.IS_DIRTY);
                childB.makeDirty(Tree.IS_DIRTY);
                nNodesWritten++;
            }
        }

        // 6. Paths of the marginal branches (only with a branch-dependent clock). A changed
        //    path does not flag the node itself: updateRates() does if the rate changed.
        if (rateClock != null) {
            for (int s = 0; s < nSwept; s++) {
                int lineage = sweptLineages[s];
                int nr = lineageToNr(lineage);
                if (lineage == rootLineage) {
                    if (pathBranches[nr] != null) {
                        beginModification();
                        pathBranches[nr] = null;
                        pathEvents[nr] = null;
                    }
                } else if (!pathMatches(pathBranches[nr], pathEvents[nr], lineage)) {
                    beginModification();
                    int n = pathEventCount(lineage);
                    int[] branches = new int[n + 1];
                    long[] events = new long[n];
                    branches[0] = linStart[lineage];
                    int i = 0;
                    for (int e = linHead[lineage]; e >= 0; e = evNext[e], i++) {
                        events[i] = evId[e];
                        branches[i + 1] = evBranch[e];
                    }
                    // Paths are replaced, never changed in place (the snapshot shares them)
                    pathBranches[nr] = branches;
                    pathEvents[nr] = events;
                }
            }
        }

        // 7. Root: the last firing creator (or the only lineage)
        Node newRoot = m_nodes[lineageToNr(rootLineage)];
        if (newRoot != root || newRoot.getParent() != null) {
            beginModification();
            newRoot.setParent(null);
            setRootOnly(newRoot);
            newRoot.makeDirty(Tree.IS_FILTHY);
        }
        postCache = null;
    }

    private int claimFreeNr(int preferred) {
        if (freeTop == 0) {
            throw new IllegalStateException("MarginalTree " + getID() + " ran out of internal node numbers.");
        }
        if (preferred >= 0) {
            for (int j = 0; j < freeTop; j++) {
                if (freeStack[j] == preferred) {
                    freeStack[j] = freeStack[--freeTop];
                    return preferred;
                }
            }
        }
        return freeStack[--freeTop];
    }

    /*
     * Direct updates (used by operators that know the topology is unchanged)
     */

    /**
     * Set the height of the marginal node created by the given CF node (if it creates one here)
     * to the CF node's current height. The caller guarantees that the change does not alter the
     * order of events along the branches of the CF node (i.e. the topology is unchanged).
     */
    public void updateCreatorHeight(Node cfNode) {
        int nr = cfToMarg[cfNode.getNr()];
        if (nr >= 0)
            setNodeHeight(nr, cfNode.getHeight());
    }

    /**
     * Set the height of the marginal node created by the given conversion (if it creates one
     * here) to the conversion's current height. Same precondition as {@link #updateCreatorHeight(Node)}.
     */
    public void updateCreatorHeight(Conversion conv) {
        Integer nr = convToMarg.get(conv.getID());
        if (nr != null)
            setNodeHeight(nr, conv.getHeight());
    }

    /**
     * Set the heights of all marginal nodes to the current heights of their creators. The caller
     * guarantees that the order of events along every CF branch is unchanged (e.g. scaling all
     * heights), so the topology stays the same.
     */
    public void updateAllHeights() {
        ConversionList convs = acg.getConversions();
        for (int nr = 0; nr < getNodeCount(); nr++) {
            double height;
            switch (creatorKind[nr]) {
                case CREATOR_LEAF:
                case CREATOR_CF:
                    height = acg.getNode(creatorKey[nr]).getHeight();
                    break;
                case CREATOR_CONVERSION:
                    height = convs.get(creatorKey[nr]).getHeight();
                    break;
                default:
                    throw new IllegalStateException("MarginalTree " + getID() + ": node " + nr + " has no creator.");
            }
            setNodeHeight(nr, height);
        }
    }

    private void setNodeHeight(int nr, double height) {
        Node node = m_nodes[nr];
        if (node.getHeight() == height)
            return;
        beginModification();
        ((MarginalNode) node).setHeight(height);
        node.makeDirty(Tree.IS_DIRTY);
        for (Node child : node.getChildren())
            child.makeDirty(Tree.IS_DIRTY);
        postCache = null;
        nNodesWritten++;
    }

    /*
     * Branch-dependent clocks: paths and rates
     */

    /**
     * Use a clock with branch-dependent rates on the clonal frame: from now on the tree records
     * the paths of its marginal branches and provides their rates ({@link #getBranchRate(int)}).
     */
    public void setBranchRateClock(BranchRateModel clock) {
        if (rateClock == clock)
            return;
        if (rateClock != null)
            throw new IllegalStateException("MarginalTree " + getID() + " already uses clock " + rateClock
                    + "; one marginal tree cannot be used with two different clocks.");
        rateClock = clock;

        int nNodes = getNodeCount();
        pathBranches = new int[nNodes][];
        pathEvents = new long[nNodes][];
        rates = new double[nNodes];
        sPathBranches = new int[nNodes][];
        sPathEvents = new long[nNodes][];
        sRates = new double[nNodes];
        // Lineage ids: leaves and new records below nCFNodes, frontier nodes above (see computeStructure)
        linStart = new int[2 * nCFNodes];
        linHead = new int[2 * nCFNodes];
        linTail = new int[2 * nCFNodes];

        initialising = true;
        try {
            computeStructure();
            reconcile();
            setReflectedVersions();
            updateRates(Double.NEGATIVE_INFINITY);
        } finally {
            initialising = false;
        }
        touched = false;
        versionsTouched = false;
        makeAllDirty();
    }

    /**
     * Notify the tree that the clock changed (its rates are recomputed on the next use).
     */
    public void invalidateRates() {
        if (rateClock != null)
            clockStale = true;
    }

    /**
     * @return the rate of the marginal branch above the node with the given number (only
     * available with a branch-dependent clock; call {@link #ensureUpToDate()} first).
     */
    public double getBranchRate(int nr) {
        return rates[nr];
    }

    /**
     * Recompute the branch rates from the paths and the current heights and clock. Nodes whose
     * rate changed are flagged dirty.
     *
     * @param lowestChange nothing below this height changed since the rates were last computed
     *        (and the clock did not change): branches ending below it keep their rates.
     */
    protected void updateRates(double lowestChange) {
        assert !clockStale || lowestChange == Double.NEGATIVE_INFINITY;
        clockStale = false;
        nRateUpdates++;
        ConversionList convs = acg.getConversions();
        int rootNr = root.getNr();
        for (int nr = 0; nr < getNodeCount(); nr++) {
            Node parent = m_nodes[nr].getParent();
            if (parent != null && parent.getHeight() < lowestChange)
                continue;
            double rate = (nr == rootNr) ? 1.0 : computeRate(nr, convs);
            if (Double.doubleToRawLongBits(rate) != Double.doubleToRawLongBits(rates[nr])) {
                beginModification();
                rates[nr] = rate;
                m_nodes[nr].makeDirty(Tree.IS_DIRTY);
            }
        }
    }

    /**
     * Time-weighted mean of the CF rates along the path of the marginal branch above node nr.
     * The terms are summed in path order, so equal paths give bit-identical rates in all trees.
     */
    protected double computeRate(int nr, ConversionList convs) {
        int[] branches = pathBranches[nr];
        long[] events = pathEvents[nr];
        double bottom = m_nodes[nr].getHeight();
        double top = m_nodes[nr].getParent().getHeight();

        double t = bottom;
        double sum = 0.0;
        for (int i = 0; i < events.length; i++) {
            double tEvent = eventHeight(events[i], convs);
            sum += cfRate(branches[i]) * (tEvent - t);
            t = tEvent;
        }
        sum += cfRate(branches[events.length]) * (top - t);

        double length = top - bottom;
        return (length > 0.0) ? sum / length : cfRate(branches[0]);
    }

    private double cfRate(int cfNr) {
        return rateClock.getRateForBranch(acg.getNode(cfNr));
    }

    private double eventHeight(long eventId, ConversionList convs) {
        byte kind = (byte) (eventId >>> 32);
        int key = (int) eventId;
        return (kind == CREATOR_CF) ? acg.getNode(key).getHeight() : convs.get(key).getHeight();
    }

    /** @return the CF branches of the path above node nr (bottom to top), or null for the root. */
    public int[] getPathBranches(int nr) {
        return pathBranches == null ? null : pathBranches[nr];
    }

    /*
     * Mapping queries
     */

    /** @return the marginal node created by the given CF node, or null. */
    public MarginalNode getNodeForCFNode(Node cfNode) {
        int nr = cfToMarg[cfNode.getNr()];
        return (nr < 0) ? null : (MarginalNode) m_nodes[nr];
    }

    /** @return the marginal node created by the given conversion, or null. */
    public MarginalNode getNodeForConversion(Conversion conv) {
        Integer nr = convToMarg.get(conv.getID());
        return (nr == null) ? null : (MarginalNode) m_nodes[nr];
    }

    /** @return the creator kind of the marginal node with the given number. */
    public byte getCreatorKind(int nr) {
        return creatorKind[nr];
    }

    /** @return the creator key (CF node nr or conversion ID) of the marginal node with the given number. */
    public int getCreatorKey(int nr) {
        return creatorKey[nr];
    }

    /*
     * Snapshot / store / restore / accept
     */

    /**
     * Called before the first change of the structure or mapping in the current proposal.
     */
    protected void beginModification() {
        if (touched || initialising)
            return;
        takeSnapshot();
        touched = true;
        acg.markTreeModified(this);
    }

    /**
     * Called before the reflected versions change in the current proposal.
     */
    protected void beginVersionChange() {
        if (versionsTouched || initialising)
            return;
        sAcgVersion = reflectedAcgVersion;
        sBlockVersion = reflectedBlockVersion;
        versionsTouched = true;
        acg.markTreeModified(this);
    }

    protected void takeSnapshot() {
        int nNodes = getNodeCount();
        for (int nr = 0; nr < nNodes; nr++) {
            Node node = m_nodes[nr];
            sHeight[nr] = node.getHeight();
            sParent[nr] = (node.getParent() == null) ? -1 : node.getParent().getNr();
            if (nr >= nLeaves) {
                sLeft[nr] = node.getLeft().getNr();
                sRight[nr] = node.getRight().getNr();
            }
        }
        System.arraycopy(creatorKind, 0, sCreatorKind, 0, nNodes);
        System.arraycopy(creatorKey, 0, sCreatorKey, 0, nNodes);
        System.arraycopy(cfToMarg, 0, sCfToMarg, 0, nCFNodes);
        // convToMarg is not copied: it is the inverse of the creator arrays and rebuilt on restore
        sRoot = root.getNr();
        if (rateClock != null) {
            System.arraycopy(pathBranches, 0, sPathBranches, 0, nNodes);
            System.arraycopy(pathEvents, 0, sPathEvents, 0, nNodes);
            System.arraycopy(rates, 0, sRates, 0, nNodes);
        }
    }

    protected void restoreSnapshot() {
        int nNodes = getNodeCount();
        for (int nr = 0; nr < nNodes; nr++) {
            MarginalNode node = (MarginalNode) m_nodes[nr];
            node.setHeight(sHeight[nr]);
            node.setParent(sParent[nr] < 0 ? null : m_nodes[sParent[nr]]);
            if (nr >= nLeaves)
                node.update(sHeight[nr], (MarginalNode) m_nodes[sLeft[nr]], (MarginalNode) m_nodes[sRight[nr]]);
        }
        // update() sets the children's parent pointers; reinstate the stored ones
        for (int nr = 0; nr < nNodes; nr++)
            ((MarginalNode) m_nodes[nr]).setParent(sParent[nr] < 0 ? null : m_nodes[sParent[nr]]);

        System.arraycopy(sCreatorKind, 0, creatorKind, 0, nNodes);
        System.arraycopy(sCreatorKey, 0, creatorKey, 0, nNodes);
        System.arraycopy(sCfToMarg, 0, cfToMarg, 0, nCFNodes);
        convToMarg.clear();
        for (int nr = nLeaves; nr < nNodes; nr++)
            if (creatorKind[nr] == CREATOR_CONVERSION)
                convToMarg.put(creatorKey[nr], nr);
        setRootOnly(m_nodes[sRoot]);
        freeTop = 0;
        if (rateClock != null) {
            System.arraycopy(sPathBranches, 0, pathBranches, 0, nNodes);
            System.arraycopy(sPathEvents, 0, pathEvents, 0, nNodes);
            System.arraycopy(sRates, 0, rates, 0, nNodes);
        }
        postCache = null;
    }

    /**
     * @return true iff the structure of this tree changed in the current proposal.
     */
    public boolean hasStructuralChanges() {
        return touched;
    }

    /**
     * @return true iff anything (structure or reflected versions) changed in the current proposal.
     */
    public boolean isModifiedInProposal() {
        return touched || versionsTouched;
    }

    /**
     * Flag every node as filthy, e.g. after initialisation or when the whole state is marked dirty.
     */
    public void makeAllDirty() {
        allDirty = true;
        super.setEverythingDirty(true);
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        if (isDirty) {
            makeAllDirty();
        } else {
            allDirty = false;
            flagsPending = false;
            super.setEverythingDirty(false);
        }
    }

    @Override
    public boolean requiresRecalculation() {
        ensureUpToDate();
        if (validate)
            checkConsistency("requiresRecalculation");
        return touched || allDirty || flagsPending;
    }

    @Override
    public boolean somethingIsDirty() {
        ensureUpToDate();
        return touched || allDirty || flagsPending;
    }

    @Override
    public void store() {
        // Snapshots are taken lazily at the first modification of a proposal
        // (BEAST calls store() of calculation nodes only after the proposal).
        managed = true;
    }

    @Override
    public void restore() {
        boolean hadFlags = touched;
        if (touched)
            restoreSnapshot();
        if (versionsTouched) {
            reflectedAcgVersion = sAcgVersion;
            reflectedBlockVersion = sBlockVersion;
        }
        touched = false;
        versionsTouched = false;
        if (managed) {
            clearFlags();
        } else if (hadFlags) {
            flagsPending = true;
        }
        if (validate && isUpToDate())
            checkConsistency("restore");
    }

    @Override
    protected void accept() {
        boolean hadFlags = touched;
        touched = false;
        versionsTouched = false;
        if (managed) {
            clearFlags();
        } else if (hadFlags) {
            flagsPending = true;
        }
    }

    /**
     * Accept the current proposal. Used by the clock wrapper: after a clock-only proposal, BEAST
     * accepts or restores the clock wrapper but not the tree, whose rates changed.
     */
    public void acceptProposal() {
        accept();
    }

    /**
     * Keep the changes of the current proposal without clearing the dirty flags. Used for trees
     * that BEAST does not accept or restore (see {@link #managed}); if the proposal was in fact
     * rejected, the version check resweeps the tree when it is next used.
     */
    public void commitProposal() {
        if (touched)
            flagsPending = true;
        touched = false;
        versionsTouched = false;
    }

    private void clearFlags() {
        if (touched || allDirty || flagsPending)
            super.setEverythingDirty(false);
        allDirty = false;
        flagsPending = false;
    }

    /**
     * @return true iff the cached value of a (BEAST-managed) likelihood on this tree is valid for
     * the current tree: the tree is up to date and no node is flagged dirty.
     */
    public boolean isLikelihoodCacheCurrent() {
        return managed && isUpToDate() && !touched && !allDirty && !flagsPending;
    }

    /**
     * @return true iff the tree is up to date and neither changed in the current proposal nor
     * flagged dirty, so that {@link #revertToCheckpoint()} can return to the current state.
     */
    public boolean isCheckpoint() {
        return isUpToDate() && !clockStale && !touched && !allDirty && !flagsPending;
    }

    /**
     * Return to the state of a checkpoint (see {@link #isCheckpoint()}) after edits that were
     * undone, e.g. a Gibbs flip that was flipped back. The caller guarantees that the conversion
     * graph and the block have the same content as at the checkpoint (only their versions
     * differ). The tree was unchanged in the proposal up to the checkpoint, so the snapshot holds
     * the checkpoint state; restoring it and clearing the flags lets the likelihood keep its
     * cached value.
     */
    public void revertToCheckpoint() {
        if (touched) {
            restoreSnapshot();
            touched = false;
            super.setEverythingDirty(false);
        }
        setReflectedVersions();
        if (validate)
            checkConsistency("revertToCheckpoint");
    }

    @Override
    public void startEditing(Operator operator) {
        // Marginal trees are not part of the state: no StateNode bookkeeping.
        postCache = null;
    }

    /*
     * Validation
     */

    /**
     * Compare the tree with a from-scratch sweep of the current conversion graph and block.
     * Throws an IllegalStateException describing the first mismatch.
     */
    public void checkConsistency(String context) {
        computeStructure();

        Map<Long, Integer> recByCreator = new HashMap<>();
        for (int i = 0; i < nRecords; i++)
            recByCreator.put(creatorId(recKind[i], recKey[i]), i);

        for (int i = 0; i < nLeaves; i++) {
            if (m_nodes[i].getHeight() != acg.getNode(i).getHeight())
                throw inconsistency(context, "leaf " + i + " has height " + m_nodes[i].getHeight()
                        + " but CF leaf has " + acg.getNode(i).getHeight());
        }

        int nInternal = 0;
        for (int nr = nLeaves; nr < getNodeCount(); nr++) {
            if (creatorKind[nr] == CREATOR_NONE)
                throw inconsistency(context, "internal node " + nr + " has no creator");
            nInternal++;
            Integer rec = recByCreator.get(creatorId(creatorKind[nr], creatorKey[nr]));
            if (rec == null)
                throw inconsistency(context, "node " + nr + " (" + creatorString(nr) + ") does not fire");
            if (lookupNr(creatorKind[nr], creatorKey[nr]) != nr)
                throw inconsistency(context, "forward mapping of " + creatorString(nr) + " does not point to node " + nr);
            Node node = m_nodes[nr];
            if (node.getHeight() != recHeight[rec])
                throw inconsistency(context, "node " + nr + " (" + creatorString(nr) + ") has height "
                        + node.getHeight() + ", expected " + recHeight[rec]);
            long a = lineageCreatorId(recChildA[rec]);
            long b = lineageCreatorId(recChildB[rec]);
            long l = creatorId(creatorKind[node.getLeft().getNr()], creatorKey[node.getLeft().getNr()]);
            long r = creatorId(creatorKind[node.getRight().getNr()], creatorKey[node.getRight().getNr()]);
            if (!((a == l && b == r) || (a == r && b == l)))
                throw inconsistency(context, "node " + nr + " (" + creatorString(nr) + ") has wrong children");
            if (node.getLeft().getParent() != node || node.getRight().getParent() != node)
                throw inconsistency(context, "parent pointers below node " + nr + " are inconsistent");
        }
        if (nInternal != nRecords)
            throw inconsistency(context, "tree has " + nInternal + " internal nodes, sweep has " + nRecords);

        Node expectedRoot = (nRecords > 0) ? null : m_nodes[0];
        if (nRecords > 0) {
            int rootRec = nRecords - 1;
            if (creatorId(creatorKind[root.getNr()], creatorKey[root.getNr()]) != creatorId(recKind[rootRec], recKey[rootRec]))
                throw inconsistency(context, "root is " + creatorString(root.getNr()));
        } else if (root != expectedRoot) {
            throw inconsistency(context, "root of single-leaf tree is not the leaf");
        }
        if (root.getParent() != null)
            throw inconsistency(context, "root has a parent");

        if (rateClock != null) {
            int rootLineage = (nRecords > 0) ? nLeaves + nRecords - 1 : 0;
            for (int lineage = 0; lineage < nLeaves + nRecords; lineage++) {
                int nr = (lineage < nLeaves) ? lineage
                        : lookupNr(recKind[lineage - nLeaves], recKey[lineage - nLeaves]);
                if (lineage == rootLineage)
                    continue;
                if (!pathMatches(pathBranches[nr], pathEvents[nr], lineage))
                    throw inconsistency(context, "path above node " + nr + " (" + creatorString(nr) + ") is "
                            + Arrays.toString(pathBranches[nr]) + " " + Arrays.toString(pathEvents[nr])
                            + ", sweep gives start " + linStart[lineage] + " with " + pathEventCount(lineage) + " events");
            }
            if (!clockStale) {
                ConversionList convs = acg.getConversions();
                for (int nr = 0; nr < getNodeCount(); nr++) {
                    double expected = (nr == root.getNr()) ? 1.0 : computeRate(nr, convs);
                    if (rates[nr] != expected)
                        throw inconsistency(context, "rate above node " + nr + " (" + creatorString(nr) + ") is "
                                + rates[nr] + ", expected " + expected);
                }
            }
        }
    }

    private long lineageCreatorId(int lineage) {
        if (lineage < nLeaves)
            return creatorId(CREATOR_LEAF, lineage);
        int rec = lineage - nLeaves;
        return creatorId(recKind[rec], recKey[rec]);
    }

    private static long creatorId(byte kind, int key) {
        return ((long) kind << 32) | (key & 0xffffffffL);
    }

    private String creatorString(int nr) {
        switch (creatorKind[nr]) {
            case CREATOR_LEAF: return "leaf " + creatorKey[nr];
            case CREATOR_CF: return "CF node " + creatorKey[nr];
            case CREATOR_CONVERSION: return "conversion " + creatorKey[nr];
            default: return "no creator";
        }
    }

    private IllegalStateException inconsistency(String context, String message) {
        return new IllegalStateException("MarginalTree " + getID() + " (block " + block.getID()
                + ") inconsistent at " + context + ": " + message);
    }

    /**
     * Checksum over the structure and the creator mapping, so that BEAST's debug mode
     * (validateReject) detects restore errors.
     */
    @Override
    public int getChecksum() {
        int h = root.getNr();
        for (int nr = 0; nr < getNodeCount(); nr++) {
            Node node = m_nodes[nr];
            h = 31 * h + (node.getParent() == null ? -1 : node.getParent().getNr());
            h = 31 * h + Double.hashCode(node.getHeight());
            h = 31 * h + creatorKind[nr];
            h = 31 * h + creatorKey[nr];
            if (rateClock != null)
                h = 31 * h + Double.hashCode(rates[nr]);
        }
        return h;
    }

    /*
     * Tree plumbing
     */

    protected void initArraysSlim() {
        m_nodes = new MarginalNode[nodeCount];
        listMarginalNodes((MarginalNode) root, m_nodes);
        postCache = null;
    }

    void listMarginalNodes(final MarginalNode node, final Node[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);

        for (final Node child : node.getChildren()) {
            listMarginalNodes((MarginalNode) child, nodes);
        }
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

        root = nodes[tree.getRoot().getNr()];
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

}
