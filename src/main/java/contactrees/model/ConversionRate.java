package contactrees.model;

import java.io.PrintStream;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.type.RealScalar;
import contactrees.ConversionGraph;


public class ConversionRate extends CalculationNode implements RealScalar<NonNegativeReal>, Loggable {

    final public Input<RealScalar<NonNegativeReal>> expectedConversionsInput = new Input<>(
            "expectedConversions",
            "The expected number of conversions in the whole phylogeny.",
            Input.Validate.REQUIRED);

    final public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph containing, composed of a tree (clonal frame) and conversion edges.",
            Input.Validate.REQUIRED);

    final public Input<Boolean> linearContactGrowthInput = new Input<>(
            "linearContactGrowth",
            "Contact process is applied per lineage, i.e. the expected number of contact edges grows linearly with lineages.",
            false);
    @Override
    public void initAndValidate() {
    }

    @Override
    public double get() {
        double expectedConversions = expectedConversionsInput.get().get();
        return expectationToRate(expectedConversions);
    }

    @Override
    public Double getLower() {
        return expectationToRate(getExpectation().getLower());
    }

    @Override
    public Double getUpper() {
        return expectationToRate(getExpectation().getUpper());
    }

    protected double expectationToRate(Double expectation) {
        if (linearContactGrowthInput.get()) {
            return expectation / getTreeLength();
        } else {
            return expectation / getPairedTreeLength();
        }
    }

    protected double rateToExpectation(Double rate) {
        if (linearContactGrowthInput.get()) {
            return rate * getTreeLength();
        } else {
            return rate * getPairedTreeLength();
        }
    }

    protected double getPairedTreeLength() {
        return networkInput.get().getClonalFramePairedLength();
    }

    protected double getTreeLength() {
        return networkInput.get().getClonalFrameLength();
    }

    protected RealScalar<NonNegativeReal> getExpectation() {
        return expectedConversionsInput.get();
    }


    @Override
    public void store() {}

    @Override
    public void restore() {}

    @Override
    public String toString() {
        return Double.toString(get());
    }

    @Override
    public NonNegativeReal getDomain() {
        return NonNegativeReal.INSTANCE;
    }

    @Override
    public void init(PrintStream out) {
        out.print(getID() + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(get() + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

}
