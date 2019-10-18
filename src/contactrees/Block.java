/**
 * 
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.lang.IllegalArgumentException;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("A block which is affected by a marginal tree and defined through a set of conversions wihtin a network.")
public class Block extends StateNode { 

	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The network on which the conversion moves are defined.",
			Input.Validate.REQUIRED);
	
	ConversionGraph acg;
	List<Conversion> convs;

	@Override
	public void initAndValidate() {
		acg = networkInput.get();
		convs = new ArrayList<Conversion>();
	}

    /**
     * StateNode implementation *
     */
	@Override
	protected void store() {
		// TODO Auto-generated method stub

	}

	@Override
	public void restore() {
		// TODO Auto-generated method stub

	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
		// TODO Auto-generated method stub

	}

	@Override
	public StateNode copy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void assignTo(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFrom(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void fromXML(final org.w3c.dom.Node node) {
		// TODO Auto-generated method stub

	}

	@Override
	public int scale(double scale) {
		throw new IllegalArgumentException("Marginal tree moves can not be scaled");
	}
	
	/**
	 * Loggable implementation *
	 */
	@Override
	public void init(PrintStream out) {
		// TODO Auto-generated method stub
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
	}
	
	/**
	 * Function implementation *
	 */
	
	@Override
	public int getDimension() {
		// TODO ?  
		return convs.size();
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO ?
		return convs.get(dim).height;
	}
}
