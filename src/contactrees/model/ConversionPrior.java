/**
 *
 */
package contactrees.model;

import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
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

    final public Input<RealParameter> conversionRateInput = new Input<>(
            "conversionRate",
            "The rate at which a pair of lineages will get in contact and form a conversion.");

    final public Input<RealParameter> expectedConversionsInput = new Input<>(
            "expectedConversions",
            "The expected number of conversions in the whole tree.",
            Input.Validate.XOR, conversionRateInput);

    final public Input<Integer> lowerCCBoundInput = new Input<>("lowerConvCountBound",
            "Lower bound on conversion count.", 0);

    final public Input<Integer> upperCCBoundInput = new Input<>("upperConvCountBound",
            "Upper bound on conversion count.", Integer.MAX_VALUE);

    final public Input<BooleanParameter> linearContactGrowthInput = new Input<>(
            "linearContactGrowth",
            "Contact process is applied per lineage, i.e. the expected number of contact edges grows linearly with lineages.",
            new BooleanParameter(new Boolean[] {false}));

    ConversionGraph acg;
    PriorityQueue<Conversion> convQueue;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        acg = networkInput.get();
        convQueue = new PriorityQueue<>();
    }

    protected double getExpectedConversions() {
        if (expectedConversionsInput.get() != null){
            return expectedConversionsInput.get().getValue();
        } else {
            double convRate = conversionRateInput.get().getValue();
            if (linearContactGrowthInput.get().getValue()) {
                return convRate * acg.getClonalFrameLength();
            } else {
                return convRate * acg.getClonalFramePairedLength();
            }
        }
    }

    protected double getConversionRate() {
        if (conversionRateInput.get() != null){
            return conversionRateInput.get().getValue();
        } else {
            double eConv = expectedConversionsInput.get().getValue();
            if (linearContactGrowthInput.get().getValue()) {
                return eConv / acg.getClonalFrameLength();
            } else {
                return eConv / acg.getClonalFramePairedLength();
            }

        }
    }

    @Override
    public double calculateLogP() {
        logP = 0;
        double convRate = getConversionRate();

        // Handle some corner cases
        assert convRate >= 0.0;
        if (convRate == 0.0) {
            if (acg.getConvCount() == 0)
                logP += 0.0;
            else
                logP += Double.NEGATIVE_INFINITY;
            return logP;
        }

        // Check whether conversion count exceeds bounds.
        if (acg.getConvCount() < lowerCCBoundInput.get()
                || acg.getConvCount() > upperCCBoundInput.get()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        convQueue.clear();
        convQueue.addAll(acg.getConversions().getConversions());

        // Poisson prior on the number of conversions per event interval (between sampling/coalescent events)
        double localConvRate = convRate;
        List<CFEventList.Event> events = acg.getCFEvents();
        for (int i = 0; i<events.size()-1; i++) {
            CFEventList.Event eStart= events.get(i);
            CFEventList.Event eEnd = events.get(i+1);
            double dt = eEnd.getHeight() - eStart.getHeight();
            int k = eStart.getLineageCount();
            if (k == 1) continue;

            double waitingtime = dt * k * (k-1);

            if (linearContactGrowthInput.get().getValue())
                localConvRate = convRate / (k-1);

            // Add probability for waiting time in the current interval
            logP -= localConvRate * waitingtime;

            // Add event probabilities for the contact edges in this interval
            while ((convQueue.peek() != null) && (convQueue.peek().getHeight() < eEnd.getHeight())) {
                logP += Math.log(localConvRate);
                convQueue.poll();
            }
        }
        assert convQueue.isEmpty();  // We iterated through all conversions end removed them from `convs`

        // Correct for probability mass outside the specified bounds on number of conversions.
        if (lowerCCBoundInput.get()>0 || upperCCBoundInput.get()<Integer.MAX_VALUE) {
            double poissonMean = getExpectedConversions();
            try {
                logP -= new PoissonDistributionImpl(poissonMean)
                        .cumulativeProbability(
                                lowerCCBoundInput.get(),
                                upperCCBoundInput.get());
            } catch (MathException e) {
                throw new RuntimeException("Error computing modification to ARG " +
                        "prior density required by conversion number constraint.");
            }
        }

        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        // For now we use the safe version (always recalculate)
        return true;
        // TODO: Use the version below when sure that dirty logic is fine in ACG.
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

    /**
     * Sample conversion edges according to the current convRate
     * in the current clonal frame.
     */
    private void generateConversions(double expectedConversions) {
        acg.removeAllConversions();

        // Draw number of conversions:
        double nConvMean = expectedConversions;
        int nConv = (int) Randomizer.nextPoisson(nConvMean);

        // Generate conversions:
        for (int i=0; i<nConv; i++) {
            Conversion conv = acg.addNewConversion();
            attachEdge(conv);
        }
    }

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
        double[] intervalVolumes = cfEventList.getIntervalVolumes(!linearContactGrowthInput.get().getValue());
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
        if (linearContactGrowthInput.get().getValue()) {
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
