/**
 *
 */
package contactrees.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;

import org.apache.commons.statistics.distribution.PoissonDistribution;

import beast.base.inference.Distribution;
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
            "conversionRate",
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

    // The network containing conversion edges
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
        // If the conversion rates input is a calculation node 
        // (typically based on the expected number of conversions per epoch),
        // then issue a warning if its options do not match those of this prior:
        if (conversionRatesInput.get() != null && conversionRatesInput.get() instanceof ConversionRate) {
        	ConversionRate rate = (ConversionRate) conversionRatesInput.get();
        	if (rate.linearContactGrowthInput.get() != linearContactGrowthInput.get()) {
        		rate.linearContactGrowthInput.setValue(linearContactGrowthInput.get(), rate);
        		Log.warning("==========================================================================");
        		Log.warning("==========================================================================");        		
        		Log.warning("LinearContactGrowth differs between " + getID()  + " and " + rate.getID());
        		Log.warning("Setting linearContactGrowth to " + linearContactGrowthInput.get() + " on conversion rate parameter " + rate.getID());
        		Log.warning("==========================================================================");
        		Log.warning("==========================================================================");
        	}
        }
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
		times = conversionRatesInput.get() != null ? new double[conversionRatesInput.get().size()] : new double[expectedConversionsInput.get().size()];
        Double timesLogP = updateTimes();
        if (timesLogP.isInfinite()) {
        	throw new RuntimeException("The times inputs are not valid. Please ensure that the origin height, if specified, is not less than the tree height and that all change times, if specified in forward time, do not exceed the origin/root height (or 1, if they are specified relative to this height).");
        }
        clonalFrameLengths = new double[times.length];
        Double clonalFrameLengthsLogP = updateClonalFrameLengths();
        if (clonalFrameLengthsLogP.isInfinite()) {
        	throw new RuntimeException("The network input is producing clonal frame lengths that are not valid in one or more epochs. Please check the network input.");
        }
		rates = new double[times.length];
		expectedConversions = new double[times.length];
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
			double intervalWidth = maxTime / times.length;
            for (int i = 0; i < times.length; i++) {
            	times[i] = intervalWidth * i;
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
		// If they are, then convert them to heights, add them to the times array, starting with a height of 0, 
		// and sort them in ascending order:
		ArrayList<Double> conversionRateChangeTimes = new ArrayList<Double>();
		conversionRateChangeTimes.add(0.0);
		for (int i = 0; i < conversionRateChangeTimesInput.get().size(); i++) {
    		double t = conversionRateChangeTimesInput.get().get(i) * (relative ? maxTime : 1);
    		if (reverse) {
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
     * Returns the index of the epoch that contains (or ends at) the given time.
     * @param t The time in question.
     * @return The index of this time or the next smallest time in the times array.
     */
    public int index(double t) {
        if (t >= times[times.length - 1])
            return times.length - 1;
        int epoch = Arrays.binarySearch(times, t);
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
    	for (int i = 0; i < times.length; i++) {
    		clonalFrameLengths[i] = 0.0;
    	}
    	// Then proceed through all sampling and coalescence events 
    	// and add their clonal frame lengths to the internal array:
        List<CFEventList.Event> events = acg.getCFEvents();
        for (int i = 0; i < events.size() - 1; i++) {
            CFEventList.Event eStart = events.get(i);
            int eStartEpoch = index(eStart.getHeight());
            CFEventList.Event eEnd = events.get(i+1);
            int eEndEpoch = index(eEnd.getHeight());
            // By definition, between any two events, 
            // there will be a constant number of lineages in the tree:
            int k = eStart.getLineageCount();
            // If there is only one, then no contact events are possible, so we move on:
            if (k == 1) continue;
            // Otherwise, for every epoch that contains part of this interval,
            // we calculate the expected number of conversions 
            // in the part of the interval contained in that epoch:
            for (int epoch = eStartEpoch; epoch <= eEndEpoch; epoch++) {
            	double start = times[epoch] > eStart.getHeight() ? times[epoch] : eStart.getHeight();
            	double end = (epoch < times.length - 1 && times[epoch + 1] < eEnd.getHeight()) ? times[epoch + 1] : eEnd.getHeight();
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
    		for (int i = 0; i < times.length; i++) {
        		if (conversionRatesInput.get().get(i) < 0.0) {
            		return Double.NEGATIVE_INFINITY;
            	}
        	}
        	for (int i = 0; i < times.length; i++) {
        		rates[i] = conversionRatesInput.get().get(i);
        		expectedConversions[i] = rates[i] * clonalFrameLengths[i];
        	}
        	return 0.0;
    	}
        // Otherwise, the expected numbers of conversions by epoch have been specified;
    	// validate them, copy them into the expected conversions array directly,
    	// and calculate the expected conversion rates using them:
    	for (int i = 0; i < times.length; i++) {
    		if (expectedConversionsInput.get().get(i) < 0.0) {
        		return Double.NEGATIVE_INFINITY;
        	}
    		if (expectedConversionsInput.get().get(i) != 0.0 && clonalFrameLengths[i] == 0.0) {
        		return Double.NEGATIVE_INFINITY;
        	}
    	}
    	for (int i = 0; i < times.length; i++) {
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
     * Returns a list containing this prior's epoch start times.
     * @return A list containing this prior's epoch start times.
     */
    public List<Double> getTimes() {
    	List<Double> timesList = new ArrayList<Double>();
    	for (int i = 0; i < times.length; i++) {
    		timesList.add(times[i]);
    	}
    	return timesList;
    }
    
    /**
     * Returns a list containing this prior's clonal frame lengths by epoch.
     * @return A list containing this prior's clonal frame lengths by epoch.
     */
    public List<Double> getClonalFrameLengths() {
    	List<Double> clonalFrameLengthsList = new ArrayList<Double>();
    	for (int i = 0; i < clonalFrameLengths.length; i++) {
    		clonalFrameLengthsList.add(clonalFrameLengths[i]);
    	}
    	return clonalFrameLengthsList;
    }
    
    /**
     * Returns a list containing this prior's conversion rates by epoch.
     * @return A list containing this prior's conversion rates by epoch.
     */
    public List<Double> getRates() {
    	List<Double> ratesList = new ArrayList<Double>();
    	for (int i = 0; i < rates.length; i++) {
    		ratesList.add(rates[i]);
    	}
    	return ratesList;
    }
    
    /**
     * Returns a list containing this prior's expected numbers of conversions by epoch.
     * @return A list containing this prior's expected numbers of conversions by epoch.
     */
    public List<Double> getExpectedConversions() {
    	List<Double> expectedConversionsList = new ArrayList<Double>();
    	for (int i = 0; i < expectedConversions.length; i++) {
    		expectedConversionsList.add(expectedConversions[i]);
    	}
    	return expectedConversionsList;
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
        // Poisson prior on the number of conversions per event interval (between sampling/coalescent events)
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
        	while ((epoch < times.length - 1) && (convHeight >= times[epoch + 1])) {
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
            logP -= pInBounds;
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
        throw new UnsupportedOperationException("For direct simulation, use ACGDistribution.");
    }

    @Override
    public List<String> getConditions() {
        throw new UnsupportedOperationException("For direct simulation, use ACGDistribution.");
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("For direct simulation, use ACGDistribution.");
    }

    // /**
    //  * Sample conversion edges according to the current convRate
    //  * in the current clonal frame.
    //  */
    // private void generateConversions(double expectedConversions) {
    //     acg.removeAllConversions();

    //     // Draw number of conversions:
    //     double nConvMean = expectedConversions;
    //     int nConv = (int) Randomizer.nextPoisson(nConvMean);

    //     // Generate conversions:
    //     for (int i=0; i<nConv; i++) {
    //         Conversion conv = acg.addNewConversion();
    //         attachEdge(conv);
    //     }
    // }

    /**
     * Attach the conversion edge at a random point in the clonal
     * frame, selecting points of departure and arrival.
     *
     * @param conv recombination to associate
     * @return the log probability density of the attachment points
     */
    public double attachEdge(Conversion conv) {
        double logQ = 0.0;
        CFEventList cfEventList = acg.getCFEventList();
        List<Event> cfEvents = cfEventList.getCFEvents();

        // Choose event interval
        double[] intervalVolumes = cfEventList.getIntervalVolumes(!linearContactGrowthInput.get());
        int iEvent = Util.sampleCategorical(intervalVolumes);
        Event event = cfEvents.get(iEvent);
        double pInterval = intervalVolumes[iEvent] / Util.sum(intervalVolumes);
        logQ += Math.log(pInterval);

        // Choose height within interval
        double height = Randomizer.uniform(event.getHeight(), cfEvents.get(iEvent+1).getHeight());
        logQ -= Math.log(cfEvents.get(iEvent+1).getHeight() - event.getHeight());
        Set<Node> activeLineages = acg.getLineagesAtHeight(height);
        if (activeLineages.size() <= 1)
            return Double.POSITIVE_INFINITY;
        conv.setHeight(height);

        // Choose source lineage (given the height)
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
        if (linearContactGrowthInput.get()) {
            double logQ = 0.0;
            double height = conv.getHeight();
            List<CFEventList.Event> events = acg.getCFEvents();

            // Prob. of event interval
            double totalVolume = 0.0;
            double intervalVolume = 0.0;
            double intervalLength = 0.0;
            int hits = 0;
            for (int i = 0; i<events.size()-1; i++) {
                CFEventList.Event eStart = events.get(i);
                CFEventList.Event eEnd = events.get(i+1);
                double dt = eEnd.getHeight() - eStart.getHeight();
                int k = eStart.getLineageCount();
                if (k == 1) continue;
                if (eStart.getHeight() <= height && height < eEnd.getHeight()) {
                    // Probability of choosing this interval (normalization over intervals is done further down)
                    intervalVolume = k * dt;
                    // Probability density of height within interval
                    intervalLength = dt;
                    // Just to make sure we find one and only one active interval
                    hits++;
                }
                totalVolume += (k * dt);
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
        } else {
            return -Math.log(acg.getClonalFramePairedLength());
        }
    }
}
