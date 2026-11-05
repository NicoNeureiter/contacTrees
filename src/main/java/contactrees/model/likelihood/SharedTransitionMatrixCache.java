package contactrees.model.likelihood;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.spec.evolution.substitutionmodel.BasicGeneralSubstitutionModel;

/**
 * Transition probability matrices shared by all tree likelihoods that use the same substitution
 * model.
 *
 * The marginal trees have real heights, so a branch that appears in several marginal trees (same
 * creators at both ends, same clock) has exactly the same length in all of them. Each distinct
 * distance (branch length times rate) therefore needs its transition matrix only once.
 *
 * The cache is content-addressed: an entry is fully determined by the distance and the current
 * eigendecomposition of the substitution model. It is never stored or restored. All entries are
 * dropped whenever the eigendecomposition object of the model changes, which for
 * {@link BasicGeneralSubstitutionModel} happens exactly when the rate matrix is recomputed
 * (parameter change) or swapped back (restore).
 *
 * Each thread has its own cache instance, so multi-threaded likelihood evaluation needs no locks.
 * The cache is direct-mapped: a slot holds one distance and a colliding distance overwrites it.
 */
public class SharedTransitionMatrixCache {

    static final int LOG2_CAPACITY = 14;
    static final int CAPACITY = 1 << LOG2_CAPACITY;

    // One ThreadLocal of caches per substitution model
    /** Bypass the cache and compute every matrix (for validation in tests). */
    public static boolean bypass = false;

    static final Map<SubstitutionModel, ThreadLocal<SharedTransitionMatrixCache>> REGISTRY = new IdentityHashMap<>();
    static final List<SharedTransitionMatrixCache> ALL_CACHES = new ArrayList<>();

    static {
        if (Boolean.getBoolean("contactrees.matrixCacheStats"))
            Runtime.getRuntime().addShutdownHook(new Thread(() -> System.err.println(statistics())));
    }

    /**
     * @return the cache of the given model for the calling thread
     */
    static ThreadLocal<SharedTransitionMatrixCache> forModel(BasicGeneralSubstitutionModel model) {
        synchronized (REGISTRY) {
            return REGISTRY.computeIfAbsent(model, m -> ThreadLocal.withInitial(() -> {
                SharedTransitionMatrixCache cache = new SharedTransitionMatrixCache(model);
                synchronized (REGISTRY) {
                    ALL_CACHES.add(cache);
                }
                return cache;
            }));
        }
    }

    /**
     * Sharing is only valid if the matrices depend on nothing but the distance and the
     * eigendecomposition, i.e. if the model uses the transition probability code of
     * BasicGeneralSubstitutionModel unchanged.
     */
    static boolean supports(SubstitutionModel model) {
        if (!(model instanceof BasicGeneralSubstitutionModel))
            return false;
        try {
            Class<?> c = model.getClass();
            Method full = c.getMethod("getTransitionProbabilities", Node.class, double.class, double.class,
                    double.class, double[].class, boolean.class);
            Method normalized = c.getMethod("getTransitionProbabilities", Node.class, double.class, double.class,
                    double.class, double[].class);
            Method eigen = c.getMethod("getEigenDecomposition", Node.class);
            return full.getDeclaringClass() == BasicGeneralSubstitutionModel.class
                    && normalized.getDeclaringClass() == BasicGeneralSubstitutionModel.class
                    && eigen.getDeclaringClass() == BasicGeneralSubstitutionModel.class
                    && !model.canReturnComplexDiagonalization();
        } catch (NoSuchMethodException e) {
            return false;
        }
    }

    public static String statistics() {
        long hits = 0, misses = 0, resets = 0;
        synchronized (REGISTRY) {
            for (SharedTransitionMatrixCache cache : ALL_CACHES) {
                hits += cache.hits;
                misses += cache.misses;
                resets += cache.resets;
            }
        }
        long total = hits + misses;
        return String.format("Transition matrix cache: %d requests, %d hits (%.1f%%), %d computed, %d resets",
                total, hits, total == 0 ? 0.0 : 100.0 * hits / total, misses, resets);
    }

    public static long getHitCount() {
        long hits = 0;
        synchronized (REGISTRY) {
            for (SharedTransitionMatrixCache cache : ALL_CACHES)
                hits += cache.hits;
        }
        return hits;
    }

    public static void resetStatistics() {
        synchronized (REGISTRY) {
            for (SharedTransitionMatrixCache cache : ALL_CACHES)
                cache.hits = cache.misses = cache.resets = 0;
        }
    }

    final BasicGeneralSubstitutionModel model;
    final int matrixSize;
    final long[] keys = new long[CAPACITY];
    final int[] generations = new int[CAPACITY];
    final double[] matrices;
    int generation = 1;
    EigenDecomposition eigenDecomposition;

    long hits, misses, resets;

    SharedTransitionMatrixCache(BasicGeneralSubstitutionModel model) {
        this.model = model;
        int stateCount = model.getStateCount();
        this.matrixSize = stateCount * stateCount;
        this.matrices = new double[CAPACITY * matrixSize];
    }

    /**
     * Drop all entries if the substitution model changed since the last call. Called before each
     * likelihood evaluation; the model cannot change during one.
     */
    void sync() {
        EigenDecomposition current = model.getEigenDecomposition(null);
        if (current != eigenDecomposition) {
            eigenDecomposition = current;
            generation++;
            resets++;
        }
    }

    void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // Same expression as BasicGeneralSubstitutionModel, which uses nothing else
        double distance = (startTime - endTime) * rate;
        long bits = Double.doubleToRawLongBits(distance);
        int slot = (int) ((bits * 0x9E3779B97F4A7C15L) >>> (64 - LOG2_CAPACITY));
        int offset = slot * matrixSize;

        if (generations[slot] == generation && keys[slot] == bits) {
            System.arraycopy(matrices, offset, matrix, 0, matrixSize);
            hits++;
            return;
        }

        model.getTransitionProbabilities(node, startTime, endTime, rate, matrix);
        System.arraycopy(matrix, 0, matrices, offset, matrixSize);
        keys[slot] = bits;
        generations[slot] = generation;
        misses++;
    }

    /**
     * The view of the substitution model given to one tree likelihood. It routes transition
     * probabilities through the cache of the evaluating thread.
     */
    static class CachingSubstitutionModel implements SubstitutionModel {
        final BasicGeneralSubstitutionModel model;
        final ThreadLocal<SharedTransitionMatrixCache> caches;
        SharedTransitionMatrixCache active;

        CachingSubstitutionModel(BasicGeneralSubstitutionModel model) {
            this.model = model;
            this.caches = forModel(model);
        }

        /** Select and sync the cache of the calling thread; call before each evaluation. */
        void beginEvaluation() {
            active = caches.get();
            active.sync();
        }

        @Override
        public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
            if (bypass) {
                model.getTransitionProbabilities(node, startTime, endTime, rate, matrix);
                return;
            }
            if (active == null)
                beginEvaluation();
            active.getTransitionProbabilities(node, startTime, endTime, rate, matrix);
        }

        @Override
        public double[] getRateMatrix(Node node) {
            return model.getRateMatrix(node);
        }

        @Override
        public double[] getFrequencies() {
            return model.getFrequencies();
        }

        @Override
        public int getStateCount() {
            return model.getStateCount();
        }

        @Override
        public EigenDecomposition getEigenDecomposition(Node node) {
            return model.getEigenDecomposition(node);
        }

        @Override
        public boolean canReturnComplexDiagonalization() {
            return model.canReturnComplexDiagonalization();
        }

        @Override
        public boolean canHandleDataType(DataType dataType) {
            return model.canHandleDataType(dataType);
        }
    }
}
