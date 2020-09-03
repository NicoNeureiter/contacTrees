/**
 * 
 */
package contactrees.model;

import java.util.List; 
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import contactrees.CFEventList;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.util.Util;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;

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
			"The rate at which a pair of lineages will get in contact and form a conversion.",
			Input.Validate.REQUIRED);

	final public Input<Integer> lowerCCBoundInput = new Input<>("lowerConvCountBound",
            "Lower bound on conversion count.", 0);

	final public Input<Integer> upperCCBoundInput = new Input<>("upperConvCountBound",
            "Upper bound on conversion count.", Integer.MAX_VALUE);
	
	ConversionGraph acg;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		acg = networkInput.get();
	}
	
	@Override
	public double calculateLogP() {
		logP = 0.0;
		double convRate = conversionRateInput.get().getValue();
		
        // Check whether conversion count exceeds bounds.
        if (acg.getConvCount()<lowerCCBoundInput.get()
                || acg.getConvCount()>upperCCBoundInput.get()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // Poisson prior on the number of conversions
        double poissonMean = convRate * acg.getClonalFramePairedLength();
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
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public List<String> getConditions() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public void sample(State state, Random random) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

}
