/**
 *
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;

/**
 *
 *
 * @author Nico Neureiter
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
	public Input<Boolean> logGeneFlowInput = new Input<>(
	        "logGeneFlow",
	        "Include logs about the gene flow (percentage of blocks borrowed) between different clades of the tree.",
	        false);

	protected ConversionGraph acg;
	protected BlockSet blockSet;
	final protected String[] columnNames = {"rootHeight", "treeLength", "pairedTreeLength", "convCount", "meanConvHeight", "moveCount", "movesPerConv"};

	@Override
	public void initAndValidate() {
		blockSet = blockSetInput.get();
		acg = networkInput.get();
	}

    @Override
    public void init(final PrintStream out) {
        for (String column : columnNames)
            out.print(column + "\t");

        if (logGeneFlowInput.get())
            for (CladePair cladePair : getCladePairs(acg))
                out.print(cladePair.toString() + "\t");

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

        if (logGeneFlowInput.get()) {
            HashMultiset<String> geneFlows = HashMultiset.create();
            HashMultimap<Node, Integer> branchedOut = HashMultimap.create();
            for (Conversion conv : acg.getConversions().asSortedArray()) {
                CladePair cladePair = new CladePair(conv);
                for (int iBlock : blockSet.getAffectedBlockIDs(conv)) {
                    if (!branchedOut.containsEntry(cladePair.target, iBlock)) {
                        branchedOut.put(cladePair.target, iBlock);
                        geneFlows.add(cladePair.toString());
                    }
                }
            }
        }
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
     * CLADE PAIR CLASS FOR CONVENIENCE
     */

    class CladePair {

        Node source, target;

        public CladePair(Node target, Node source) {
            this.target = target;
            this.source = source;
        }

        public CladePair(Conversion conv) {
            this(conv.node1, conv.node2);
        }

        @Override
        public String toString() {
            return "flow_" + source.getNr() + "_to_" + target.getNr();
        }
    }

    public ArrayList<CladePair> getCladePairs(ConversionGraph acg) {
        ArrayList<CladePair> cladePairs = new ArrayList<>();

        for (Node source : acg.getNodesAsArray()) {
            for (Node target : acg.getNodesAsArray()) {
                if (source == target)
                    continue;

                if (source.isLeaf() && target.isLeaf())
                    cladePairs.add(new CladePair(target, source));
            }
        }

        return cladePairs;
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
