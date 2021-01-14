/**
 *
 */
package contactrees.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.util.Randomizer;
import contactrees.CFEventList;
import contactrees.CFEventList.Event;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.util.Util;

/**
 * @author Nico Neureiter
 */
public class ACGDistribution extends Distribution {

	final public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The conversion graph containing the conversion edges.",
			Input.Validate.REQUIRED);

	final public Input<TreeDistribution> cfModelInput = new Input<>(
	        "cfModel",
	        "The tree prior of the clonal frame.",
	        Input.Validate.REQUIRED);

	final public Input<RealParameter> conversionRateInput = new Input<>(
			"conversionRate",
			"The rate at which a pair of lineages will get in contact and form a conversion."
			);

	final public Input<RealParameter> expectedConversionsInput = new Input<>(
	        "expectedConversions",
	        "The expected number of conversions in the whole phylogeny (alternative parameterization to conversionRate).",
	        Input.Validate.XOR, conversionRateInput);

	final public Input<Integer> lowerCCBoundInput = new Input<>("lowerConvCountBound",
            "Lower bound on conversion count.", 0);

	final public Input<Integer> upperCCBoundInput = new Input<>("upperConvCountBound",
            "Upper bound on conversion count.", Integer.MAX_VALUE);

	ConversionGraph acg;
	TreeDistribution cfModel;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		acg = networkInput.get();
		cfModel = cfModelInput.get();
	}

	protected double getExpectedConversions() {
	    if (conversionRateInput.get() != null) {
	        double convRate = conversionRateInput.get().getValue();
	        return convRate * acg.getClonalFramePairedLength();
        } else {
            return expectedConversionsInput.get().getValue();
        }
	}

	@Override
	public double calculateLogP() {
		logP = cfModel.calculateLogP();

        // Check whether conversion count exceeds bounds.
        if (acg.getConvCount()<lowerCCBoundInput.get()
                || acg.getConvCount()>upperCCBoundInput.get()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // Poisson prior on the number of conversions
        double poissonMean = getExpectedConversions();
        logP += -poissonMean + acg.getConvCount() * Math.log(poissonMean);

        assert poissonMean >= 0.0;
        if (poissonMean == 0) {
        	if (acg.getConvCount() == 0) {
        	    logP = 0.0;
            } else {
                logP = Double.NEGATIVE_INFINITY;
    		}
        	return logP;
        }

        // Probability density of each conversion placement
		for (Conversion conv : acg.getConversions())
			logP += calculateConversionLogP(conv);

        // Correct for probability mass outside the specified bounds (on number of conversions)
        if (lowerCCBoundInput.get()>0 || upperCCBoundInput.get()<Integer.MAX_VALUE) {
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

	public double calculateConversionLogP(Conversion conv) {
		ConversionGraph acg = networkInput.get();

		// For now, a uniform distribution over time and pairs of lineages
		return -Math.log(acg.getClonalFramePairedLength());
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
        List<String> arguments = new ArrayList<>();
        arguments.add(networkInput.get().getID());

        return arguments;
	}

	@Override
	public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(conversionRateInput.get().getID());

        return conditions;
	}

	@Override
	public void sample(State state, Random random) {
        if (sampledFlag)
            return;
        sampledFlag = true;

        // Cause conditional parameters to be sampled

        sampleConditions(state, random);
        cfModel.sample(state, random);

        acg.assignFromFragile((ConversionGraph) cfModel.treeInput.get());

        generateConversions(getExpectedConversions());
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
            associateConversionWithCF(conv);
        }
    }

    /**
     * Associates recombination with the clonal frame, selecting
     * points of departure and arrival.
     *
     * @param conv recombination to associate
     */
    private void associateConversionWithCF(Conversion conv) {
        CFEventList cfEventList = acg.getCFEventList();
        List<Event> cfEvents = cfEventList.getCFEvents();

        // Choose event interval
        double[] intervalVolumes = cfEventList.getIntervalVolumes();
        int iEvent = Util.sampleCategorical(intervalVolumes);
        Event event = cfEvents.get(iEvent);

        // Choose height within interval
        double height = Randomizer.uniform(event.getHeight(), cfEvents.get(iEvent+1).getHeight());
        conv.setHeight(height);

        // Choose source lineage (given the height)
        Set<Node> activeLineages = acg.getLineagesAtHeight(height);
        Node node1 = Util.sampleFrom(activeLineages);
        conv.setNode1(node1);
        assert node1.getHeight() < height;

        // Choose destination lineage (given the height and node1)
        activeLineages.remove(node1);
        Node node2 = Util.sampleFrom(activeLineages);
        conv.setNode2(node2);

        assert conv.isValid();
    }

    @Override
    public void printDebug() {
        System.out.println(getID() + "  \t" + isDirtyCalculation() + "\t\t" + logP);
    }

}
