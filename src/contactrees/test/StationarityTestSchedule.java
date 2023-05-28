/**
 * 
 */
package contactrees.test;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Formatter;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.util.Randomizer;

/**
 * OperatorSchedule which is designed to test the stationarity conditions of a single operator (without requiring ergodicity). 
 * The schedule works as follows:
 *  - In the beginning and after each logged sample: simulate the state directly using the simulatorOperator.
 *  - For the remaining steps (to the next logged sample): apply the operator(s) given in operatorsInput.  
 *  
 * @author Nico Neureiter
 */
public class StationarityTestSchedule extends OperatorSchedule {

    final public Input<Integer> logIntervalInput = new Input<>(
            "logInterval",
            "The interval at which the logger logs and at which new direct simulator samples need to be generated.",
            Input.Validate.REQUIRED);

    final public Input<Operator> simulatorOperatorInput = new Input<>(
            "simulatorOperator",
            "The operator generating a sample from the direct simulator.",
            Input.Validate.REQUIRED);
    
	ArrayList<Operator> opMap;
	Iterator<Operator> opIterator;
	int iSample, logInterval;
	Operator simulatorOperator;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		buildOpMap();
		iSample = 0;
		logInterval = logIntervalInput.get();
	}
	
	public void buildOpMap() {
		operators.clear();
		operators.addAll(operatorsInput.get());
        for (Operator o : operators) {
        	o.setOperatorSchedule(this);
        }
		
        simulatorOperator = simulatorOperatorInput.get();
        simulatorOperator.setOperatorSchedule(this);
        
		System.out.println("Number of operators: " + operators.size());		
	}
	
    /**
     * Select operators cycling in the order they were added with a frequency proportional to their weight.
     * @return
     */
    public Operator selectOperator() {
        Operator selectedOperator;
        
        if (iSample % logInterval == 1)
            selectedOperator = simulatorOperator;
        else
            selectedOperator = super.selectOperator();

//        System.out.println("Step " + iSample + ": " + selectedOperator.getName());
        
        iSample += 1;
    	
        return selectedOperator;
    }
    
    @Override
    public void showOperatorRates(final PrintStream out) {}


}
