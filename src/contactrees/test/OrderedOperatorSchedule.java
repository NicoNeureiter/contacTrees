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

import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.util.Randomizer;

/**
 * 
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class OrderedOperatorSchedule extends OperatorSchedule {

	List<Operator> opMap;
	Iterator<Operator> opIterator;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		buildOpMap();
	}
	
	public void buildOpMap() {
		operators.clear();
		operators.addAll(operatorsInput.get());
        for (Operator o : operators) {
        	o.setOperatorSchedule(this);
        }
		
		opMap = new LinkedList<Operator>();
		for (Operator p : operators) {
			System.out.println(p.getName() + " : " + p.getWeight());
			for (int pRepeat=0; pRepeat<p.getWeight(); pRepeat++) {
				opMap.add(p);
			}
		}
		System.out.println("Number of operators: " + opMap.size());
		opIterator = opMap.iterator();
		
	}
	
    /**
     * Select operators cycling in the order they were added with a frequency proportional to their weight.
     * @return
     */
    public Operator selectOperator() {
    	if (!opIterator.hasNext()) {
    		if (opMap.isEmpty()) {
    			buildOpMap();
    		}
    		opIterator = opMap.iterator();
    	}
    	
    	assert opIterator.hasNext();
        return opIterator.next();
    }
    
    @Override
    public void showOperatorRates(final PrintStream out) {}


}
