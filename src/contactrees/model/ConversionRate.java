package contactrees.model;

import org.w3c.dom.Node;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import contactrees.ConversionGraph;


public class ConversionRate extends RealParameter {

    final public Input<RealParameter> expectedConversionsInput = new Input<>(
            "expectedConversions",
            "The expected number of conversions in the whole phylogeny.",
            Input.Validate.REQUIRED);

    final public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph containing, composed of a tree (clonal frame) and conversion edges.",
            Input.Validate.REQUIRED);

    final public Input<BooleanParameter> linearContactGrowthInput = new Input<>(
            "linearContactGrowth",
            "Contact process is applied per lineage, i.e. the expected number of contact edges grows linearly with lineages.",
            new BooleanParameter(new Boolean[] {false}));

    public ConversionRate() {
        valuesInput.setRule(Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
//        setInputValue("value", expectationToRate(getExpectation().getValue()));
    }

    @Override
    public Double getValue() {
        double expectedConversions = expectedConversionsInput.get().getValue();
        return expectationToRate(expectedConversions);
    }

    @Override
    public int getDimension() {
        return getExpectation().getDimension();
    }

    @Override
    public double getArrayValue(int dim) {
        return expectationToRate(getExpectation().getArrayValue(dim));
    }

    @Override
    public Double getValue(int i) {
        return expectationToRate(getExpectation().getValue(i));
    }

    @Override
    public void setValue(int i, Double value) {
        getExpectation().setValue(i, rateToExpectation(value));
    }

    @Override
    public void setValue(Double value) {
        getExpectation().setValue(rateToExpectation(value));
    }

    @Override
    public Double getLower() {
        return expectationToRate(getExpectation().getLower());
    }

    @Override
    public void setLower(Double lower) {
        getExpectation().setLower(rateToExpectation(lower));
    }

    @Override
    public Double getUpper() {
        return expectationToRate(getExpectation().getUpper());
    }

    @Override
    public void setUpper(Double upper) {
        getExpectation().setUpper(rateToExpectation(upper));
    }

    @Override
    public Double[] getValues() {
        return getExpectation().getValues();
    }

    @Override
    public int getMinorDimension1() {
        return getExpectation().getMinorDimension1();
    }

    @Override
    public int getMinorDimension2() {
        return getExpectation().getMinorDimension2();
    }

    @Override
    public Double getMatrixValue(int i, int j) {
        return getExpectation().getMatrixValue(i, j);
    }

    @Override
    public void swap(int i, int j) {
        getExpectation().swap(i, j);
    }

    /*
     * Convenience methods.
     */

    protected double expectationToRate(Double expectation) {
        if (linearContactGrowthInput.get().getValue()) {
            return expectation / getTreeLength();
        } else {
            return expectation / getPairedTreeLength();
        }
    }

    protected double rateToExpectation(Double rate) {
        if (linearContactGrowthInput.get().getValue()) {
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

    protected RealParameter getExpectation() {
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
        return Double.toString(getValue());
    }

}
