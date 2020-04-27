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
        out.print("rootHeight\ttreeLength\tpairedTreeLength\tconvCount\tmeanConvHeight\tmoveCount\tmovesPerConv\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	double height = acg.getRoot().getHeight();
    	double treeLength = acg.getClonalFrameLength();
    	double pairedTreeLength = acg.getClonalFramePairedLength();
    	int nConv = acg.getConvCount();
    	int moveCount = blockSet.countMoves();
    	double meanConvHeight = 0;
    	double movesPerConv = 0;
        if (nConv > 0) {
            meanConvHeight = meanConvHeight(acg);
            movesPerConv = (double) moveCount / nConv;
        }
        
        out.print(height + "\t" + treeLength + "\t" + pairedTreeLength + "\t" + nConv + "\t" + meanConvHeight + "\t" + moveCount + "\t" + movesPerConv + "\t");
    }
	
	protected double meanConvHeight(ConversionGraph acg) {
		double sum = 0.0;
		for (Conversion conv : acg.getConversions()) {
			sum += conv.getHeight();
		}
		return sum / acg.getConvCount();
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
