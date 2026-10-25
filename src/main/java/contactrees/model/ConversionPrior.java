/**
 *
 */
package contactrees.model;

import java.util.Arrays;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;
import java.util.function.Supplier;

import org.apache.commons.statistics.distribution.PoissonDistribution;

import beast.base.inference.CalculationNode;
import beast.base.inference.Distribution;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.CFEventList;
import contactrees.CFEventList.Event;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.util.Util;

/**
 *
 *
 * @author Nico Neureiter
 */
public class ConversionPrior extends Distribution {

    final public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph containing the conversion edges.",
            Input.Validate.REQUIRED);
    
    final public Input<RealScalar<PositiveReal>> originInput = new Input<>(
    		"origin",
    		"The height of the origin of the underlying transmission process."
    		+ "This is used for measuring the boundaries of epochs in forward time if multiple conversion rates are specified."
    		+ "If specified, it must not be less than the height of the tree."
    		+ "If not specified, the height of the tree's root (i.e., the MRCA of all taxa) is used instead.");

    final public Input<RealVector<NonNegativeReal>> conversionRatesInput = new Input<>(
            "conversionRates",
            "Epoch-dependent rates at which a pair of lineages will get in contact and form a conversion."
            + "If a single rate is specified, then a single epoch containing the entire tree will be assumed."
            + "The conversion rates should be specified in forward time, starting from the epoch farthest in the past and moving towards the epoch ending in the present."
            + "If more than one rate is specified, then the times at which the epochs change can be specified in the conversionRateChangeTimes input.");
    
    final public Input<RealVector<NonNegativeReal>> expectedConversionsInput = new Input<>(
            "expectedConversions",
            "The expected number of conversions in each epoch."
            + "If a single number is specified, then a single epoch containing the entire tree will be assumed.",
            Validate.XOR, conversionRatesInput);
    
    final public Input<RealVector<PositiveReal>> conversionRateChangeTimesInput = new Input<>(
            "conversionRateChangeTimes",
            "The times specifying when conversion rate changes occur."
            + "The times can be specified in forward time (from the origin/root to the present) or in backward time depending on the reverseChangeTimes input."
            + "They can be specified in absolute heights or in proportions of the tree's total height (or the origin height, if specified) depending on the relativeChangeTimes input."
            + "They do not have to be specified in order, and they can be estimated during the analysis."
            + "If the conversionRates input consists of a single value, this input is ignored."
            + "If specified, this input's dimension must be one less than the dimension of the conversionRates input."
            + "If not specified, the change times are distributed evenly over the height of the tree (or of the origin, if specified).");
    
    final public Input<Boolean> relativeChangeTimesInput =
            new Input<Boolean>("relativeChangeTimes", "True if the conversion rate change times are specified relative to tree height (or the origin height, if specified).", false);
    
    final public Input<Boolean> reverseChangeTimesInput =
            new Input<Boolean>("reverseChangeTimes", "True if the time arrays are given in backwards time (measured in time before the present)." +
                    "Note that the rate array must still be given in forward time (from the epoch farthest in the past to the epoch ending in the present).", false);
    
    final public Input<Integer> lowerCCBoundInput = new Input<>("lowerConvCountBound",
            "Lower bound on conversion count.", 0);

    final public Input<Integer> upperCCBoundInput = new Input<>("upperConvCountBound",
            "Upper bound on conversion count.", Integer.MAX_VALUE);

    final public Input<Boolean> linearContactGrowthInput = new Input<>(
            "linearContactGrowth",
            "Contact process is applied per lineage, i.e. the expected number of contact edges grows linearly with lineages.",
            false);

    final public Input<Distribution> clonalFramePriorInput = new Input<>(
            "clonalFramePrior",
            "Tree prior of the clonal frame, only used for direct simulation (sample()). "
            + "Only needed if it can not be found automatically, i.e. if it does not list the network "
            + "as an argument (e.g. a coalescent prior on TreeIntervals).");

    ConversionGraph acg;
    // Data structure containing the conversion objects
    PriorityQueue<Conversion> convQueue;
    // An internal array of heights at which the conversion rate changes
    protected double[] times;
    // An internal array of total clonal frame lengths for each epoch
    protected double[] clonalFrameLengths;
    // An internal array of conversion rate values for each epoch
    protected double[] rates;
    // An internal array of expected conversions for each epoch
    protected double[] expectedConversions;
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        acg = networkInput.get();
        convQueue = new PriorityQueue<>();
        
        Distribution cfPrior = clonalFramePriorInput.get();
        if (cfPrior instanceof ConversionPrior)
            throw new IllegalArgumentException("clonalFramePrior of " + getID() + " must be a tree prior, "
                    + "not a conversion prior (" + cfPrior.getID() + ").");

        // If multiple conversion rates and their change times are specified, 
        // then ensure that there is one fewer change time than there are rates:
		if (conversionRateChangeTimesInput.get() != null && conversionRatesInput.get() != null && conversionRateChangeTimesInput.get().size() != conversionRatesInput.get().size() - 1) {
			throw new RuntimeException("The specified conversionRateChangeTimes input should have dimension one less than that of conversionRates, but it does not.");
		}
		// If multiple expected numbers of conversions and conversion rate change times are specified, 
		// then ensure that there is one fewer change time than there are expected numbers of conversions:
		if (conversionRateChangeTimesInput.get() != null && expectedConversionsInput.get() != null && conversionRateChangeTimesInput.get().size() != expectedConversionsInput.get().size() - 1) {
			throw new RuntimeException("The specified conversionRateChangeTimes input should have dimension one less than that of expectedConversions, but it does not.");
		}
        // Then attempt to initialize the internal arrays, and report any invalid assignments:
		times = conversionRatesInput.get() != null ? new double[conversionRatesInput.get().size() - 1] : new double[expectedConversionsInput.get().size() - 1];
        Double timesLogP = updateTimes();
        if (timesLogP.isInfinite()) {
        	throw new RuntimeException("The times inputs are not valid. Please ensure that the origin height, if specified, is not less than the tree height and that all change times, if specified in forward time, do not exceed the origin/root height (or 1, if they are specified relative to this height).");
        }
        clonalFrameLengths = new double[times.length + 1];
        Double clonalFrameLengthsLogP = updateClonalFrameLengths();
        if (clonalFrameLengthsLogP.isInfinite()) {
        	throw new RuntimeException("The network input is producing clonal frame lengths that are not valid in one or more epochs. Please check the network input.");
        }
		rates = new double[times.length + 1];
		expectedConversions = new double[times.length + 1];
		Double ratesAndExpectedConversionsLogP = updateRatesAndExpectedConversions();
        if (ratesAndExpectedConversionsLogP.isInfinite()) {
        	throw new RuntimeException("The rates input or expected conversions input is not valid. Please ensure that all specified values are non-negative and that the expected number of conversions is not positive in epochs where contact edges cannot be added to the network (e.g., if at most one lineage in the network occurs in those epochs).");
        }
    }
    
    /**
     * Updates the internal change time arrays of this prior.
     * If any of the current values for the change time inputs are invalid 
     * (typically due to a bad proposal by an operator),
     * a log-likelihood of negative infinity is returned.
     * @return 0 if the current time inputs are valid, and negative infinity otherwise.
     */
    protected Double updateTimes() {
    	// Get the tree's height, as it will be used in subsequent checks and updates:
    	double maxTime = originInput.get() != null ? originInput.get().get() : networkInput.get().getRoot().getHeight();
    	// If the origin height has been specified and it is currently less than the tree's height,
    	// the inputs are invalid:
    	if (maxTime < networkInput.get().getRoot().getHeight()) {
    		return Double.NEGATIVE_INFINITY;
    	}
    	// Otherwise, populate the times array based on which inputs have been specified:
		if (conversionRateChangeTimesInput.get() == null) {
			// If no change times are specified, then use equidistant change times:
			double intervalWidth = maxTime / (times.length + 1);
            for (int i = 1; i <= times.length; i++) {
            	times[i - 1] = intervalWidth * i;
            }
            return 0.0;
		}
		// Otherwise, ensure that the specified change times are valid:
		Boolean relative = relativeChangeTimesInput.get();
		Boolean reverse = reverseChangeTimesInput.get();
		if (!reverse) {
			// If change times are specified in forward time, ensure that they are all strictly less than the tree height:
			if (relative) {
				for (int i = 0; i < conversionRateChangeTimesInput.get().size(); i++) {
		    		if (conversionRateChangeTimesInput.get().get(i) >= 1.0) {
		    			return Double.NEGATIVE_INFINITY;
		    		}
		    	}
			} else {
				for (int i = 0; i < conversionRateChangeTimesInput.get().size(); i++) {
		    		if (conversionRateChangeTimesInput.get().get(i) >= maxTime) {
		    			return Double.NEGATIVE_INFINITY;
		    		}
		    	}
			}
		}
		// If they are, then convert them to heights, add them to the times array, 
		// and sort them in ascending order:
		ArrayList<Double> conversionRateChangeTimes = new ArrayList<Double>();
		for (int i = 0; i < conversionRateChangeTimesInput.get().size(); i++) {
    		double t = conversionRateChangeTimesInput.get().get(i) * (relative ? maxTime : 1);
    		if (!reverse) {
    			t = maxTime - t;
    		}
    		conversionRateChangeTimes.add(t);
    	}
		conversionRateChangeTimes.sort(null);
		// Then populate the times array based on the options:
		for (int i = 0; i < conversionRateChangeTimes.size(); i++) {
			times[i] = conversionRateChangeTimes.get(i);
		}
		return 0.0;
    }
    
    /**
     * Returns the index of the epoch that contains 
     * (or whose lowest height equals) the given time.
     * @param t The time in question.
     * @return The index of the epoch in the times array
     * that contains (or whose lowest height equals) the given time.
     */
    public int index(double t) {
    	if (times.length == 0 || t >= times[times.length - 1])
            return times.length;
        int epoch = Arrays.binarySearch(times, t);
        // The binarySearch algorithm returns a positive index 
        // if the value is present in the array
        // and a negative index (-(insertion point) - 1) if it is absent:
        if (epoch < 0) {
            epoch = -epoch - 1;
        }
        return epoch;
    }
    
    /**
     * Updates the internal array of total clonal frame lengths by epoch for this prior
     * using this prior's time arrays and the input network's events list.
     * The total lengths will depend on whether this prior assumes linear or quadratic contact growth.
     * A log-likelihood of 0 is returned if everything is valid and negative infinity if anything is invalid.
     * @return 0 if the current time inputs are valid, and negative infinity otherwise.
     */
    protected Double updateClonalFrameLengths() {
    	// First, reinitialize the internal array:
    	for (int i = 0; i < clonalFrameLengths.length; i++) {
    		clonalFrameLengths[i] = 0.0;
    	}
    	// Then proceed through all sampling and coalescence events 
    	// and add their clonal frame lengths to the internal array:
        List<Event> events = acg.getCFEvents();
        int eStartEpoch = 0;
        int eEndEpoch = 0;
        for (int i = 0; i < events.size() - 1; i++) {
            Event eStart = events.get(i);
            Event eEnd = events.get(i+1);
            // The current event starts at the previous event's end time,
            // so its epoch is the previous event's end epoch, 
            // and we just have to update the end epoch:
            eStartEpoch = eEndEpoch;
            while (eEndEpoch < times.length && eEnd.getHeight() > times[eEndEpoch]) {
            	eEndEpoch++;
            }
            // If this event and the next take place at the same time,
            // then nothing can happen between them:
            if (eStart.getHeight() == eEnd.getHeight()) continue;
            // By definition, between any two events, 
            // there will be a constant number of lineages in the tree:
            int k = eStart.getLineageCount();
            // If there is only one, then no contact events are possible, so we move on:
            if (k == 1) continue;
            // Otherwise, for every epoch that contains part of this interval,
            // we calculate the expected number of conversions 
            // in the part of the interval contained in that epoch:
            for (int epoch = eStartEpoch; epoch <= eEndEpoch; epoch++) {
            	double start = (epoch == 0 || eStart.getHeight() > times[epoch - 1]) ? eStart.getHeight() : times[epoch - 1];
            	double end = (epoch >= times.length || eEnd.getHeight() < times[epoch]) ? eEnd.getHeight() : times[epoch];
            	double dt = end - start;
            	double length = linearContactGrowthInput.get() ? dt * k : dt * k * (k-1);
            	clonalFrameLengths[epoch] += length;
            }
        }
        return 0.0;
    }
    
    /**
     * Updates the internal rate and expected conversion arrays of this prior.
     * This method assumes that the internal arrays of epoch start times and clonal frame lengths by epoch 
     * have been updated already.
     * If any of the current values are invalid (typically due to a bad proposal by an operator),
     * a log-likelihood of negative infinity is returned.
     * @return 0 if the current rate inputs are valid, and negative infinity otherwise.
     */
    protected Double updateRatesAndExpectedConversions() {
    	// Proceed based on which input has been specified:
    	if (conversionRatesInput.get() != null) {
    		// If the conversion rates by epoch have been specified, then validate them,
    		// copy them into the rates array directly,
    		// and calculate the expected conversions using them:
    		for (int i = 0; i < rates.length; i++) {
        		if (conversionRatesInput.get().get(i) < 0.0) {
            		return Double.NEGATIVE_INFINITY;
            	}
        	}
        	for (int i = 0; i < rates.length; i++) {
        		// The rates are always specified from the highest epoch to the lowest, 
        		// so we populate the internal array in reverse order:
        		rates[i] = conversionRatesInput.get().get(rates.length - i - 1);
        		expectedConversions[i] = rates[i] * clonalFrameLengths[i];
        	}
        	return 0.0;
    	}
        // Otherwise, the expected numbers of conversions by epoch have been specified;
    	// validate them, copy them into the expected conversions array directly,
    	// and calculate the expected conversion rates using them:
    	for (int i = 0; i < expectedConversions.length; i++) {
    		if (expectedConversionsInput.get().get(i) < 0.0) {
        		return Double.NEGATIVE_INFINITY;
        	}
    		if (expectedConversionsInput.get().get(i) != 0.0 && clonalFrameLengths[i] == 0.0) {
        		return Double.NEGATIVE_INFINITY;
        	}
    	}
    	for (int i = 0; i < expectedConversions.length; i++) {
    		expectedConversions[i] = expectedConversionsInput.get().get(i);
    		// The above check ensures that the expected number of conversions never exceeds zero 
    		// when the total clonal frame length is zero,
    		// so a division by zero should only occur when the expected number of conversions is also zero;
    		// in this case, we can safely assume that the conversion rate is zero:
    		rates[i] = expectedConversions[i] == 0.0 ? 0.0 : expectedConversions[i] / clonalFrameLengths[i];
    	}
    	return 0.0;
    }
    
    /**
     * Returns an array containing this prior's epoch start times.
     * @return A list containing this prior's epoch start times.
     */
    public double[] getTimes() {
    	return times;
    }
    
    /**
     * Returns an array containing this prior's clonal frame lengths by epoch.
     * @return An array containing this prior's clonal frame lengths by epoch.
     */
    public double[] getClonalFrameLengths() {
    	return clonalFrameLengths;
    }
    
    /**
     * Returns an array containing this prior's conversion rates by epoch.
     * @return An array containing this prior's conversion rates by epoch.
     */
    public double[] getRates() {
    	return rates;
    }
    
    /**
     * Returns an array containing this prior's expected numbers of conversions by epoch.
     * @return An array containing this prior's expected numbers of conversions by epoch.
     */
    public double[] getExpectedConversions() {
    	return expectedConversions;
    }

    @Override
    public double calculateLogP() {
    	// First, update the internal arrays and exit early 
    	// if any of their corresponding inputs are currently invalid:
    	Double timesLogP = updateTimes();
        if (timesLogP.isInfinite()) {
        	logP = Double.NEGATIVE_INFINITY;
        	return logP;
        }
        Double clonalFrameLengthsLogP = updateClonalFrameLengths();
        if (clonalFrameLengthsLogP.isInfinite()) {
        	logP = Double.NEGATIVE_INFINITY;
        	return logP;
        }
    	Double ratesAndExpectedConversionsLogP = updateRatesAndExpectedConversions();
        if (ratesAndExpectedConversionsLogP.isInfinite()) {
        	logP = Double.NEGATIVE_INFINITY;
        	return logP;
        }
        // If all arrays are valid, then proceed as usual:
        logP = 0.0;
        // Check whether conversion count exceeds bounds.
        if (acg.getConvCount() < lowerCCBoundInput.get() || acg.getConvCount() > upperCCBoundInput.get()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }
        // Reject conversions outside the valid range of their edges (e.g. after a tip date move)
        for (Conversion conv : acg.getConversions()) {
            if (!conv.isValid()) {
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }
        // We assume a Poisson prior on the number of conversions per event interval
        // (between sampling/coalescent events)
        for (int i = 0; i < clonalFrameLengths.length; i++) {
        	// For a Poisson-distributed random variable X,
        	// the probability of no event happening in time t
        	// is e^-E(X) = e^-(conversion rate * t).
        	// In our case, t is the total clonal frame length in the current epoch,
        	// so E(X) is the expected number of conversions in this epoch:
        	logP -= expectedConversions[i];
        }
        // To account for the conversions that did occur, 
        // we incorporate the rates in the epochs at which they did occur, 
        // proceeding in increasing order of the conversions' heights:
        convQueue.clear();
        convQueue.addAll(acg.getConversions().getConversions());
        int epoch = 0;
        while (convQueue.peek() != null) {
        	Conversion conv = convQueue.poll();
        	double convHeight = conv.getHeight();
        	// Update the current epoch to the one containing this conversion:
        	while ((epoch < times.length) && (convHeight >= times[epoch])) {
        		epoch++;
        	}
        	double conversionRate = rates[epoch];
        	// If the conversion rate for this epoch is 0, 
        	// then the probability of this conversion occurring is 0, and we can exit early:
        	if (conversionRate == 0.0) {
        		logP = Double.NEGATIVE_INFINITY;
                return logP;
        	}
            logP += Math.log(conversionRate);
        }
        assert convQueue.isEmpty();  // We iterated through all conversions end removed them from `convs`

        // Correct for probability mass outside the specified bounds on number of conversions.
        if (lowerCCBoundInput.get() > 0 || upperCCBoundInput.get() < Integer.MAX_VALUE) {
        	// The expected value of an inhomogeneous Poisson process that is piecewise constant
        	// is just the sum of the expected values of the pieces.
        	// In our case, it is the sum of the expected number of conversions in each epoch:
            double poissonMean = 0.0;
            for (int i = 0; i < expectedConversions.length; i++) {
            	poissonMean += expectedConversions[i];
            }
            PoissonDistribution poi = PoissonDistribution.of(poissonMean);
            int lower = lowerCCBoundInput.get();
            int upper = upperCCBoundInput.get();
            double pUpper = poi.cumulativeProbability(upper);
            double pBelowLower = lower > 0 ? poi.cumulativeProbability(lower - 1) : 0.0;
            double pInBounds = pUpper - pBelowLower;
            if (pInBounds <= 0.0) {
                throw new RuntimeException("Error computing modification to ARG prior density required by conversion number constraint.");
            }
            logP -= Math.log(pInBounds);
        }

        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        // For now we use the safe version (always recalculate)
        return true;
        // TODO: Use the version below when sure that dirty logic is fine in ACG.
        // JJM: If uncommented, the code below should also check 
        // if any of the other inputs (origin, epoch times, conversion rates or expected conversions) have changed.
        // return networkInput.get().somethingIsDirty();
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        arguments.add(acg.getID());
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        if (conversionRatesInput.get() != null)
            maybeAddInputToConditions(conditions, conversionRatesInput);
        else
            maybeAddInputToConditions(conditions, expectedConversionsInput);
        return conditions;
    }

    /**
     * Add an input to the `conditions` list if it is a BEASTInterface.
     * Otherwise, we assume it is constant and doesn't need to be sampled.
     */
    void maybeAddInputToConditions(List<String> conditions, Input<?> input) {
        if (input.get() instanceof BEASTInterface beastInterface) {
            conditions.add(beastInterface.getID());
        }
    }

    /**
     * Sample the network: first the clonal frame (by triggering its tree prior, see
     * {@link #getClonalFramePrior}), then the conversions given the clonal frame.
     * If the network is not estimated and has no tree prior, the clonal frame stays fixed.
     *
     * The tree prior has to be part of the simulated distribution, since DirectSimulator
     * only resets the sampledFlag of distributions reachable from there.
     */
    @Override
    public void sample(State state, Random random) {
        if (sampledFlag)
            return;
        sampledFlag = true;

        checkNetworkConditions();

        Distribution cfPrior = getClonalFramePrior();
        if (cfPrior != null)
            cfPrior.sample(state, random);

        // Conditions are sampled after the clonal frame, since the conversion rate may depend on it
        sampleConditions(state, random);
        generateConversions();
    }

    /**
     * Find the tree prior of the clonal frame. This is the explicitly specified clonalFramePrior, if
     * any. Otherwise, it is the other distribution which has the network as a direct argument.
     * Tree priors connected only indirectly (e.g. a coalescent on TreeIntervals) have to be
     * specified explicitly.
     *
     * @return the clonal frame prior, or null if there is none and the network is not estimated
     *         (fixed clonal frame).
     * @throws IllegalStateException if none is found and the network is estimated, or if it is ambiguous.
     */
    private Distribution getClonalFramePrior() {
        if (clonalFramePriorInput.get() != null)
            return clonalFramePriorInput.get();

        List<Distribution> candidates = new ArrayList<>();
        for (BEASTInterface output : acg.getOutputs()) {
            if (output == this || !(output instanceof Distribution))
                continue;
            Distribution distr = (Distribution) output;
            if (!safeIDs(distr::getArguments).contains(acg.getID()))
                continue;
            if (distr instanceof ConversionPrior)
                throw new IllegalStateException("Network " + acg.getID() + " is governed by multiple "
                        + "conversion priors (" + getID() + ", " + distr.getID() + ").");
            candidates.add(distr);
        }

        if (candidates.size() > 1)
            throw new IllegalStateException("Ambiguous clonal frame prior for network " + acg.getID()
                    + ": multiple distributions have it as an argument. Specify the tree prior via clonalFramePrior.");

        if (candidates.isEmpty()) {
            if (!acg.isEstimatedInput.get())
                return null;
            throw new IllegalStateException("No tree prior found for the clonal frame of " + acg.getID()
                    + ". Specify it via clonalFramePrior, or set estimate=\"false\" on the network "
                    + "to keep the clonal frame fixed during simulation.");
        }
        return candidates.get(0);
    }

    /**
     * Ensure that no other distribution is conditioned on the network, directly or through
     * a CalculationNode depending on it. Distribution.sampleConditions() would resolve such a
     * condition by sampling an arbitrary distribution governing the network (possibly the
     * tree prior instead of this conversion prior), making the simulation depend on the order
     * of distributions. Distributions conditioned on the network should call
     * {@link #sampleNetwork} instead.
     */
    @SuppressWarnings("deprecation")  // legacy GenericTreeLikelihood
    private void checkNetworkConditions() {
        Set<BEASTInterface> visited = new HashSet<>();
        Deque<BEASTInterface> queue = new ArrayDeque<>();
        queue.add(acg);
        while (!queue.isEmpty()) {
            BEASTInterface node = queue.poll();
            if (!visited.add(node))
                continue;

            for (BEASTInterface output : node.getOutputs()) {
                // Tree likelihoods never list their tree as a condition. Skipping them also avoids
                // core TreeLikelihood.getConditions(), which throws a NullPointerException when the
                // likelihood is delegated to BEAGLE (its m_siteModel is never set in that case).
                if (output instanceof beast.base.spec.evolution.likelihood.GenericTreeLikelihood
                        || output instanceof beast.base.evolution.likelihood.GenericTreeLikelihood)
                    continue;

                if (output instanceof Distribution && output != this && node.getID() != null
                        && safeIDs(((Distribution) output)::getConditions).contains(node.getID()))
                    throw new IllegalStateException("Distribution " + output.getID() + " lists " + node.getID()
                            + " as a condition, which depends on the network " + acg.getID() + ". This makes "
                            + "direct simulation order dependent. Use ConversionPrior.sampleNetwork() instead.");

                if (output instanceof CalculationNode)
                    queue.add(output);
            }
        }
    }

    /**
     * Sample the network by triggering the conversion prior
     * governing it. To be used by distributions conditioned on the network, instead of listing
     * the network in getConditions() (see {@link #checkNetworkConditions}).
     */
    public static void sampleNetwork(ConversionGraph acg, State state, Random random) {
        List<Distribution> networkPriors = new ArrayList<>();
        List<Distribution> owners = new ArrayList<>();
        for (BEASTInterface output : acg.getOutputs()) {
            if (!(output instanceof Distribution))
                continue;
            Distribution distr = (Distribution) output;
            if (distr instanceof ConversionPrior)
                networkPriors.add(distr);
            else if (safeIDs(distr::getArguments).contains(acg.getID()))
                owners.add(distr);
        }

        // Without a conversion prior, fall back to the tree prior (if any)
        if (networkPriors.isEmpty())
            networkPriors = owners;

        if (networkPriors.size() > 1)
            throw new IllegalStateException("Ambiguous prior for network " + acg.getID()
                    + ": multiple distributions have it as an argument.");

        if (networkPriors.size() == 1)
            networkPriors.get(0).sample(state, random);
    }

    /** Get list of IDs from getArguments() or getConditions(), which may be null or unsupported. */
    private static List<String> safeIDs(Supplier<List<String>> getter) {
        try {
            List<String> ids = getter.get();
            return ids == null ? Collections.emptyList() : ids;
        } catch (UnsupportedOperationException e) {
            return Collections.emptyList();
        }
    }

    /**
     * Sample conversion edges according to the current conversion rates
     * in the current clonal frame.
     */
    private void generateConversions() {
        acg.removeAllConversions();

        // The total expected number of conversions 
        // is the sum of each epoch's expected number of conversions:
        int nConv = 0;
        for (int epoch = 0; epoch < expectedConversions.length; epoch++) {
        	// Draw number of conversions from the Poisson distribution for this epoch:
        	double nConvMean = expectedConversions[epoch];
        	nConv += (int) Randomizer.nextPoisson(nConvMean);
        }
        // Then add the new edges to the conversion graph,
        // (they will be distributed according to the epochs' expected numbers of conversions 
        // in the attachEdge method):
        for (int i = 0; i < nConv; i++) {
            Conversion conv = acg.addNewConversion();
            attachEdge(conv);
        }
    }

    /**
     * Attach the conversion edge at a random point in the clonal
     * frame, selecting points of departure and arrival.
     * The location of the edge is distributed according to the conversion rate 
     * and number of active lineages at any given time.
     *
     * @param conv recombination to associate
     * @return the log probability density of the attachment points
     */
    public double attachEdge(Conversion conv) {
        double logQ = 0.0;
        // As an edge case, if all conversion rates are zero, then no edge can be attached:
        double rateSum = 0.0;
        for (int i = 0; i < rates.length; i++) {
        	rateSum += rates[i];
        }
        if (rateSum == 0.0) {
        	return Double.NEGATIVE_INFINITY;
        }
        // Otherwise, populate lists of starting heights, lineage counts, 
        // and weighted volumes for candidate intervals, 
        // where the intervals are divided at event boundaries and epoch boundaries:
        ArrayList<Double> intervalStarts = new ArrayList<Double>();
        ArrayList<Integer> intervalLineageCounts = new ArrayList<Integer>();
        ArrayList<Double> intervalVolumes = new ArrayList<Double>();
        CFEventList cfEventList = acg.getCFEventList();
        List<Event> cfEvents = cfEventList.getCFEvents();
        int eStartEpoch = 0;
        int eEndEpoch = 0;
        for (int i = 0; i < cfEvents.size() - 1; i++) {
            Event eStart = cfEvents.get(i);
            Event eEnd = cfEvents.get(i+1);
            // The current event starts at the previous end's end time,
            // so its epoch is the previous event's end epoch:
            eStartEpoch = eEndEpoch;
            while (eEndEpoch < times.length && eEnd.getHeight() > times[eEndEpoch]) {
            	eEndEpoch++;
            }
            // By definition, between any two events, 
            // there will be a constant number of lineages in the tree:
            int k = eStart.getLineageCount();
            // For every epoch that contains part of this interval,
            // we calculate the necessary values for the part of the interval contained in that epoch:
            for (int epoch = eStartEpoch; epoch <= eEndEpoch; epoch++) {
            	double start = (epoch == 0 || eStart.getHeight() > times[epoch - 1]) ? eStart.getHeight() : times[epoch - 1];
            	double end = (epoch >= times.length || eEnd.getHeight() < times[epoch]) ? eEnd.getHeight() : times[epoch];
            	double dt = end - start;
            	double length = linearContactGrowthInput.get() ? dt * k : dt * k * (k-1);
            	intervalStarts.add(start);
            	intervalLineageCounts.add(k);
            	intervalVolumes.add(length*rates[epoch]);
            }
        }
        // Finally, add a final interval covering everything above the root,
        // assigning it a volume of 0:
        intervalStarts.add(cfEvents.get(cfEvents.size() - 1).getHeight());
        intervalLineageCounts.add(1);
        intervalVolumes.add(0.0);
        // Sample an interval according to its weight:
        int iInterval = Util.sampleCategorical(intervalVolumes);
        double pInterval = intervalVolumes.get(iInterval) / Util.sum(intervalVolumes);
        logQ += Math.log(pInterval);
        // Then uniformly sample a height within the interval
        // (note that intervalStarts.get(iInterval + 1) should not throw an IndexOutOfBounds exception, 
        // because the single lineage above the root has a volume of zero and should never get sampled):
        double height = Randomizer.uniform(intervalStarts.get(iInterval), intervalStarts.get(iInterval + 1));
        logQ -= Math.log(intervalStarts.get(iInterval + 1) - intervalStarts.get(iInterval));
        // Then uniformly sample the source lineage at this height:
        Set<Node> activeLineages = acg.getLineagesAtHeight(height);
        if (activeLineages.size() <= 1)
            return Double.POSITIVE_INFINITY;
        conv.setHeight(height);
        Node node1 = Util.sampleFrom(activeLineages);
        logQ -= Math.log(activeLineages.size());
        conv.setNode1(node1);
        assert node1.getHeight() < height;

        // Choose destination lineage (given the height and node1)
        activeLineages.remove(node1);
        Node node2 = Util.sampleFrom(activeLineages);
        logQ -= Math.log(activeLineages.size());
        conv.setNode2(node2);

        assert conv.isValid();

        return logQ;
    }

    /**
     * Get the probability density of the attachment points of the given edge.
     *
     * @param conv recombination to associate
     * @return the log probability density of the attachment points
     */
    public double getEdgeAttachmentProb(Conversion conv) {
        double logQ = 0.0;
        double height = conv.getHeight();
        // As an edge case, if all conversion rates are zero, then no edge can be attached:
        double rateSum = 0.0;
        for (int i = 0; i < rates.length; i++) {
        	rateSum += rates[i];
        }
        if (rateSum == 0.0) {
        	return Double.NEGATIVE_INFINITY;
        }
        // Otherwise, find the height and volume of the interval containing this conversion,
        // checking all intervals partitioned by event and epoch boundaries:
        double totalVolume = 0.0;
        double intervalVolume = 0.0;
        double intervalLength = 0.0;
        CFEventList cfEventList = acg.getCFEventList();
        List<Event> cfEvents = cfEventList.getCFEvents();
        int eStartEpoch = 0;
        int eEndEpoch = 0;
        int hits = 0;
        for (int i = 0; i < cfEvents.size() - 1; i++) {
            Event eStart = cfEvents.get(i);
            Event eEnd = cfEvents.get(i+1);
            // The current event starts at the previous event's end time,
            // so its epoch is the previous event's end epoch, 
            // and we just have to update the end epoch:
            eStartEpoch = eEndEpoch;
            while (eEndEpoch < times.length && eEnd.getHeight() > times[eEndEpoch]) {
            	eEndEpoch++;
            }
            // By definition, between any two events, 
            // there will be a constant number of lineages in the tree:
            int k = eStart.getLineageCount();
            // If there is only one, then we can skip it
            // (because the edge in question could not be attached here):
            if (k == 1) continue;
            // For every epoch that contains part of this interval,
            // we calculate the necessary values for the part of the interval contained in that epoch:
            for (int epoch = eStartEpoch; epoch <= eEndEpoch; epoch++) {
            	double start = (epoch == 0 || eStart.getHeight() > times[epoch - 1]) ? eStart.getHeight() : times[epoch - 1];
            	double end = (epoch >= times.length || eEnd.getHeight() < times[epoch]) ? eEnd.getHeight() : times[epoch];
            	double dt = end - start;
            	double length = linearContactGrowthInput.get() ? dt * k : dt * k * (k-1);
            	if (start <= height && height < end) {
            		// Probability of choosing this interval (normalization over intervals is done further down)
                    intervalVolume = rates[epoch] * length;
                    // Probability density of height within interval
                    intervalLength = dt;
                    // Just to make sure we find one and only one active interval
                    hits++;
            	}
            	totalVolume += rates[epoch] * length;
            }
        }
        assert hits == 1;

        double pInterval = intervalVolume / totalVolume;
        logQ += Math.log(pInterval);

        // Prob. of height within interval
        logQ -= Math.log(intervalLength);

        Set<Node> activeLineages = acg.getLineagesAtHeight(height);
        if (activeLineages.size() <= 1)
            return Double.POSITIVE_INFINITY;

        // Prob. of source lineage (given the height)
        int k = activeLineages.size();
        logQ -= Math.log(k);
        // Prob. of destination lineage (given the height and node1)
        logQ -= Math.log(k - 1);
        return logQ;
    }
}
