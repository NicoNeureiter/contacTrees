/**
 * 
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;

/**
 * 
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class ACGStatsLogger extends BEASTObject implements Loggable {


	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The conversion graph to be logged.",
			Input.Validate.REQUIRED);
	public Input<BlockSet> blockSetInput = new Input<>(
			"blockSet",
			"The moves each local tree takes along the conversion graph.",
			Input.Validate.REQUIRED);
	
	protected ConversionGraph acg;
	protected BlockSet blockSet;
	
	@Override
	public void initAndValidate() {
		blockSet = blockSetInput.get();
		acg = networkInput.get();
	}

    @Override
    public void init(final PrintStream out) {
        out.print("convCount\trootHeight\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	int nConv = acg.getConvCount();
    	double height = acg.getRoot().getHeight();
    	
        out.print(nConv + "\t" + height + "\t");
    }


    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }

    /*
     * TESTING INTERFACE
     */
    static public ACGStatsLogger getACGStatsLogger(ConversionGraph acg, BlockSet blockSet) {
    	ACGStatsLogger acgStatsLogger = new ACGStatsLogger();
    	acgStatsLogger.initAndValidate();
    	acgStatsLogger.acg = acg;
    	acgStatsLogger.blockSet = blockSet;
    	return acgStatsLogger;
    }
    
    public ConversionGraph getACG() {
    	return acg;
    }
}
