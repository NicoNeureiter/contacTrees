package contactrees.test;

import contactrees.test.util.Counter;
import contactrees.test.util.DefaultHashMap;
import org.junit.Assert;
import beast.base.inference.Distribution;
import java.util.*;
import java.util.function.Function;


/**
 * Generic detailed-balance test harness.
 *
 * <p>Ported from {@code targetedbeast.DetailedBalanceTest} (io.github.nicfel:targetedbeast),
 * by the same author, adapted for contacTrees (contacTrees cannot depend on targetedbeast
 * itself: it is pinned to a different BEAST version and has no commons-math3 dependency).
 *
 * <p>Operators that move between "classes" of states (where the class is defined
 * by some quantisation/grouping function) should satisfy detailed balance
 * between those classes under their target distribution:
 * {@code p_i * q_ij = p_j * q_ji}. This class draws samples from a prior, applies
 * a proposal, and accumulates the Metropolis-Hastings acceptance mass flowing
 * between classes in both directions, then asserts the forward and backward flows
 * agree within a statistically estimated tolerance.
 *
 * <p>The harness is deliberately agnostic to:
 * <ul>
 *   <li>the operator being tested,</li>
 *   <li>the prior used to draw and score states,</li>
 *   <li>the quantisation method that maps a state to its class.</li>
 * </ul>
 * A subclass supplies these via {@link #getStateMappers()} and a {@link Trial}
 * per operator (passed to {@link #testDetailedBalance(Trial)}), and defines the
 * actual JUnit {@code @Test} methods.
 *
 * @param <S> the type of state being sampled and grouped (e.g. a tree or network)
 */
public abstract class DetailedBalanceTest<S> {

    /**
     * A named quantisation: maps a state to a string label identifying the class it
     * belongs to. Subclasses provide these to specify what classes to test balance over.
     * The classifier is a plain function, so mappers can be written as lambdas.
     */
    public static final class StateMapper<S> {
        private final String name;
        private final Function<S, String> classifier;

        public StateMapper(String name, Function<S, String> classifier) {
            this.name = name;
            this.classifier = classifier;
        }

        public String name() {
            return name;
        }

        public String classify(S state) {
            return classifier.apply(state);
        }
    }

    /**
     * Drives a single operator over the course of one detailed-balance test.
     *
     * <p>Created once per operator and reused across all samples, so that
     * expensive operator/prior set-up happens only once. Each iteration the
     * harness calls {@link #nextSample()} to (re)draw a state, evaluates
     * {@link #prior()} before and after a {@link #proposal()}.
     */
    public interface Trial<S> {
        /** Draw or refresh the state to test for the next iteration. */
        S nextSample() throws Exception;

        /** The prior used to score the current state (may be constant across iterations). */
        Distribution prior();

        /** Perform the proposal on the current state and return the log Hastings ratio. */
        double proposal();
    }

    /** Number of samples to draw per test. */
    protected abstract int getNumSamples();

    /** Number of standard deviations of slack allowed when comparing forward/backward flow. */
    protected double getFlowSigmaMultiplier() {
        return 4.0;
    }

    /** The named quantisation functions to verify detailed balance over. */
    protected abstract List<StateMapper<S>> getStateMappers();

    /**
     * Hook invoked after every proposal with its acceptance probability.
     * Default is a no-op; subclasses may override to gather extra diagnostics.
     */
    protected void onProposal(Trial<S> trial, double pAccept) {
    }

    /**
     * Hook invoked once after all samples have been drawn, before the balance
     * assertions. Default is a no-op; subclasses may override to add extra
     * sanity checks (e.g. that a secondary counter behaved as expected).
     */
    protected void afterSampling() {
    }

    protected void testDetailedBalance(Trial<S> trial) throws Exception {
        List<StateMapper<S>> mappers = getStateMappers();

        // One accumulator per quantisation, gathering its own statistics.
        Map<StateMapper<S>, BalanceAccumulator> accumulators = new LinkedHashMap<>();
        for (StateMapper<S> mapper : mappers) {
            accumulators.put(mapper, new BalanceAccumulator());
        }

        // Draw samples, propose, and feed the observed class transitions to each accumulator.
        int numSamples = getNumSamples();
        double totalAcceptanceMass = 0.0;
        for (int i = 0; i < numSamples; i++) {
            S state = trial.nextSample();
            Distribution prior = trial.prior();

            // The class each mapper assigns to the state before the proposal.
            Map<StateMapper<S>, String> beforeKeys = new LinkedHashMap<>();
            for (StateMapper<S> mapper : mappers) {
                String beforeKey = mapper.classify(state);
                beforeKeys.put(mapper, beforeKey);
                accumulators.get(mapper).recordState(beforeKey);
            }

            double logPBefore = prior.calculateLogP();
            double logHR = trial.proposal();
            double logPAfter = prior.calculateLogP();
            double pAccept = Math.min(1, Math.exp(logPAfter - logPBefore + logHR));

            totalAcceptanceMass += pAccept;

            onProposal(trial, pAccept);

            // The class each mapper assigns after the proposal, plus its acceptance mass.
            for (StateMapper<S> mapper : mappers) {
                String afterKey = mapper.classify(state);
                accumulators.get(mapper).recordTransition(beforeKeys.get(mapper), afterKey, pAccept);
            }
        }

        Assert.assertTrue(
            "Operator produced zero acceptance mass across " + numSamples + " proposals",
            totalAcceptanceMass > 0.0);

        afterSampling();

        // Each accumulator checks detailed balance over the transitions it observed.
        for (StateMapper<S> mapper : mappers) {
            accumulators.get(mapper).assertBalanced(mapper.name());
        }
    }

    /**
     * Gathers, for a single quantisation, the per-class state counts and the
     * Metropolis-Hastings acceptance mass flowing between classes, and asserts
     * detailed balance over what it has seen.
     */
    protected class BalanceAccumulator {
        private record Pair(String from, String to) {}

        private final Counter<String> groupCounter = new Counter<>();
        private final Counter<Pair> proposalCounter = new Counter<>();
        private final DefaultHashMap<Pair, Double> flow = new DefaultHashMap<>(0.0);
        private final DefaultHashMap<Pair, Double> flowSquares = new DefaultHashMap<>(0.0);

        /** Record the class a sampled state belongs to (before any proposal). */
        void recordState(String groupKey) {
            groupCounter.increment(groupKey);
        }

        /** Record a proposed transition between classes and its acceptance mass. */
        void recordTransition(String beforeKey, String afterKey, double pAccept) {
            Pair transitionKey = new Pair(beforeKey, afterKey);
            proposalCounter.increment(transitionKey);
            flow.put(transitionKey, flow.get(transitionKey) + pAccept);
            // Sum of squared per-sample contributions, for the empirical flow variance.
            flowSquares.put(transitionKey, flowSquares.get(transitionKey) + pAccept * pAccept);
        }

        /**
         * Assert detailed balance between classes: {@code p_i * q_ij = p_j * q_ji}.
         * Forward and backward acceptance flow must agree within a tolerance derived
         * from their estimated sampling variance.
         */
        void assertBalanced(String mapperName) {
            System.out.println("\n=== Detailed balance test for mapper: " + mapperName + " ===");
            System.out.println(groupCounter);
            System.err.println(proposalCounter);

            double flowSigmaMultiplier = getFlowSigmaMultiplier();
            String firstFailure = null;
            for (String fromGroup : groupCounter.keySet()) {
                for (String toGroup : groupCounter.keySet()) {
                    // Check each unordered pair once, skipping trivially symmetric self-loops.
                    if (fromGroup.compareTo(toGroup) <= 0)
                        continue;

                    Pair keyForward = new Pair(fromGroup, toGroup);
                    Pair keyBackward = new Pair(toGroup, fromGroup);

                    double transitionsForward = flow.get(keyForward);
                    double transitionsBackward = flow.get(keyBackward);

                    // Skip low count/high variance groups
                    if (Math.min(groupCounter.getCount(fromGroup), groupCounter.getCount(toGroup)) < 20)
                        continue;

                    // Skip low count/high variance transitions
                    if (Math.min(proposalCounter.getCount(keyForward), proposalCounter.getCount(keyBackward)) < 20)
                        continue;

                    double forwVariance = estimateFlowVariance(flowSquares.get(keyForward), transitionsForward);
                    double backVariance = estimateFlowVariance(flowSquares.get(keyBackward), transitionsBackward);
                    double sigma = Math.sqrt(forwVariance + backVariance);
                    double tolerance = flowSigmaMultiplier * sigma;
                    double diff = Math.abs(transitionsForward - transitionsBackward);

                    if (Math.max(transitionsForward, transitionsBackward) > 0)
                        System.out.println(String.format(
                            "%-16s:  %8.4f <-> %8.4f    (tol=%8.4f     diff=%8.4f   nsigma=%6.2f)%s",
                            fromGroup + "->" + toGroup, transitionsForward, transitionsBackward,
                            tolerance, diff, diff / sigma, diff > tolerance ? "  <<< FAIL" : ""));
                    if (diff > tolerance && firstFailure == null) {
                        firstFailure = String.format(
                            "Detailed balance [%s] %s->%s expected:<%s> but was:<%s>",
                            mapperName, fromGroup, toGroup, transitionsForward, transitionsBackward);
                    }
                }
            }
            Assert.assertNull(firstFailure, firstFailure);
        }
    }

    /**
     * Empirical variance of a transition flow.
     *
     * <p>The flow {@code F = sum_s X_s} is a sum over the {@code N} iid samples of the
     * per-sample contribution {@code X_s = 1{transition observed} * pAccept_s}, so
     * {@code Var(F) = N * Var(X)}, estimated by the (Bessel-corrected) sample variance:
     * {@code N/(N-1) * (sum_s X_s^2 - F^2/N)}. No assumptions are made about how
     * {@code pAccept} is distributed within a transition class; its spread is measured
     * directly through the sum of squares.
     *
     * @param sumSquaredFlow {@code sum_s X_s^2}, the sum of squared acceptance masses
     * @param totalFlow      {@code F}, the total acceptance mass for the transition
     */
    protected double estimateFlowVariance(final double sumSquaredFlow, final double totalFlow) {
        final int n = getNumSamples();
        return (sumSquaredFlow - totalFlow * totalFlow / n) * n / (n - 1.0);
    }

}
