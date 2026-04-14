package contactrees.model;

import org.w3c.dom.Node;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import contactrees.ConversionGraph;


public class ConversionRate extends RealScalarParam<NonNegativeReal> {

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

    public ConversionRate() {
        valuesInput.setRule(Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
//        setInputValue("value", expectationToRate(getExpectation().getValue()));
    }

    @Override
    public double get() {
        double expectedConversions = expectedConversionsInput.get().get();
        return expectationToRate(expectedConversions);
    }

    @Override
    public void set(double value) {
        RealScalar<?> expectedConversions = getExpectation();
        if (expectedConversions instanceof RealScalarParam expParam)
            expParam.set(rateToExpectation(value));
    }

    @Override
    public Double getLower() {
        return expectationToRate(getExpectation().getLower());
    }

    @Override
    public Double getUpper() {
        return expectationToRate(getExpectation().getUpper());
    }

    /*
     * Convenience methods.
     */

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


    /*
     * Override Parameter/StateNode methods
     */


    @Override
    public void setEverythingDirty(final boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
    public void store() {}

    @Override
    public void restore() {}

    @Override
    public void fromXML(final Node node) {}

    @Override
    public void assignFromFragile(final StateNode other) {}

    @Override
    public ConversionRate copy() {
        try {
            return (ConversionRate) this.clone();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public String toString() {
        return Double.toString(get());
    }

}
