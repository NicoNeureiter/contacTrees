package contactrees.model;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;


public class ConversionPriorLogger extends BEASTObject implements Loggable {
	
	protected ConversionPrior conversionPrior;

    final public Input<ConversionPrior> conversionPriorInput = new Input<>(
            "conversionPrior",
            "The conversion prior whose internal values will be logged.",
            Input.Validate.REQUIRED);
    
    final public Input<Boolean> logTimesInput = new Input<>(
            "logTimes",
            "Log all epoch change times (measured backwards from the present).",
            false);
    
    final public Input<Boolean> logRatesInput = new Input<>(
            "logRates",
            "Log conversion rates for all epochs (from the lowest epoch to the highest).",
            false);
    
    final public Input<Boolean> logExpectedConversionsInput = new Input<>(
            "logExpectedConversions",
            "Log expected number of conversions for all epochs (from the lowest epoch to the highest).",
            false);
    
    final public Input<Boolean> logProbabilityInput = new Input<>(
            "logProbability",
            "Log the log-probability of the conversion prior.",
            true);

    @Override
	public void initAndValidate() {
		conversionPrior = conversionPriorInput.get();
	}

    @Override
    public void init(PrintStream out) {
    	if (logTimesInput.get()) {
    		double[] times = conversionPrior.getTimes();
	    	for (int i = 0; i < times.length; i++) {
	    		out.print("ConversionChangeTimes_" + i + "\t");
	    	}
    	}
    	if (logRatesInput.get()) {
    		double[] rates = conversionPrior.getRates();
	    	for (int i = 0; i < rates.length; i++) {
	    		out.print("ConversionRates_" + i + "\t");
	    	}
    	}
    	if (logExpectedConversionsInput.get()) {
    		double[] expectedConversions = conversionPrior.getExpectedConversions();
	    	for (int i = 0; i < expectedConversions.length; i++) {
	    		out.print("ExpectedConversions_" + i + "\t");
	    	}
    	}
    	if (logProbabilityInput.get()) {
	    	out.print("ConversionPrior" + "\t");
    	}
    }

    @Override
    public void log(long sample, PrintStream out) {
    	if (logTimesInput.get()) {
    		double[] times = conversionPrior.getTimes();
    		for (int i = 0; i < times.length; i++) {
        		out.print(times[i] + "\t");
        	}
    	}
    	if (logRatesInput.get()) {
    		double[] rates = conversionPrior.getRates();
        	for (int i = 0; i < rates.length; i++) {
        		out.print(rates[i] + "\t");
        	}
    	}
    	if (logExpectedConversionsInput.get()) {
    		double[] expectedConversions = conversionPrior.getExpectedConversions();
	    	for (int i = 0; i < expectedConversions.length; i++) {
	    		out.print(expectedConversions[i] + "\t");
	    	}
    	}
    	if (logProbabilityInput.get()) {
    		double logP = conversionPrior.calculateLogP();
	    	out.print(logP + "\t");
    	}
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }
}
