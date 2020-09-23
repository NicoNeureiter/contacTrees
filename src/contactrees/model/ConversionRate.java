package contactrees.model;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
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
        return expectation / getPairedTreeLength();
    }

    protected double rateToExpectation(Double rate) {
        return rate * getPairedTreeLength();
    }
    
    protected double getPairedTreeLength() {
        return networkInput.get().getClonalFramePairedLength();
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
