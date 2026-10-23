package contactrees.model;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import contactrees.ConversionGraph;


public class ConversionRateLogger extends BEASTObject implements Loggable {
	
	protected ConversionPrior conversionPrior;

    final public Input<ConversionPrior> conversionPriorInput = new Input<>(
            "conversionPrior",
            "The conversion prior whose rates will be logged.",
            Input.Validate.REQUIRED);

    @Override
	public void initAndValidate() {
		conversionPrior = conversionPriorInput.get();
	}

    @Override
    public void init(PrintStream out) {
    	double[] rates = conversionPrior.getRates();
    	for (int i = 0; i < rates.length; i++) {
    		out.print("ConversionPriorRates_" + i + "\t");
    	}
    }

    @Override
    public void log(long sample, PrintStream out) {
    	double[] rates = conversionPrior.getRates();
    	for (int i = 0; i < rates.length; i++) {
    		out.print(rates.length + "\t");
    	}
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }
}
