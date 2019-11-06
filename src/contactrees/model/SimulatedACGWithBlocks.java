/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package contactrees.model;

import contactrees.ACGWithMetaDataLogger;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.CFEventList;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.CFEventList.Event;
import contactrees.util.Util;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.Binomial;
import beast.util.Randomizer;
import feast.nexus.NexusBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import feast.nexus.TreesBlock;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.distribution.BinomialDistribution;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulates an ARG under the full ClonalOrigin model - can be used"
    + " for chain initialization or for sampler validation.")
public class SimulatedACGWithBlocks extends BEASTObject {

	public Input<ConversionGraph> networkInput = new Input<>(
			"network",
			"The conversion graph to be logged.",
			Input.Validate.REQUIRED);
    
    public Input<BlockSet> blockSetInput = new Input<>(
    		"blockSet",
    		"The set of blocks evolving along the ACG.",
    		Input.Validate.REQUIRED);
	
    public Input<Double> conversionRateInput = new Input<>(
            "conversionRate",
            "Conversion rate parameter.",
            Input.Validate.REQUIRED);
    
    public Input<Double> moveProbInput = new Input<>(
            "moveProb",
            "Probability of a block moving over a conversion edge.",
            Input.Validate.REQUIRED);

    public Input<PopulationFunction> popFuncInput = new Input<>(
            "populationModel",
            "Demographic model to use.",
            Input.Validate.REQUIRED);
    
    public Input<Tree> clonalFrameInput = new Input<>(
            "clonalFrame",
            "Optional tree specifying fixed clonal frame.");
    
    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated ARG is additionally written to this file.");
    

    private ConversionGraph acg;
    private BlockSet blockSet;
    private Double conversionRate, moveProb;
    private PopulationFunction popFunc;

    
    @Override
    public void initAndValidate() {
    	blockSet = blockSetInput.get();
        acg = networkInput.get();
        conversionRate = conversionRateInput.get();
        moveProb = moveProbInput.get();
        popFunc = popFuncInput.get();
        
        acg.initAndValidate();
        blockSet.initAndValidate();
        

        if (clonalFrameInput.get() == null) {
            simulateClonalFrame();
        } else {
            acg.assignFromWithoutID(clonalFrameInput.get());
        }
        
        // Need to do this here as this sets the tree object that the nodes
        // point to, so without it they point to the dummy tree created by
        // super.initAndValidate().
        acg.initTreeArrays();
        
        // Generate recombinations
        generateConversions();
        
        // Write output file
        if (outputFileNameInput.get() != null) {

            NexusBuilder nexusBuilder = new NexusBuilder();
            
            nexusBuilder.append(new TaxaBlock(acg.m_taxonset.get()));
            
            nexusBuilder.append((new TreesBlock() {
                @Override
                public String getTreeString(Tree tree) {
                	ACGWithMetaDataLogger acgLogger = ACGWithMetaDataLogger.getACGWMDLogger((ConversionGraph) tree, blockSet);
                    return acgLogger.getExtendedNewick();
                }
            }).addTree(acg, "simulatedARG"));
            
            nexusBuilder.append(new NexusBlock() {

                @Override
                public String getBlockName() {
                    return "contactrees";
                }

                @Override
                public List<String> getBlockLines() {
                    List<String> lines = new ArrayList<>();

                    String blockSetLine = "blockSet";
                    for (Block block: blockSet.getBlocks())
                    	blockSetLine += " " + block.getID();
                    lines.add(blockSetLine);

                    lines.add("clonalframe_labeled " + acg.getRoot().toNewick());
                    lines.add("clonalframe_numbered " + acg.getRoot().toShortNewick(true));
                    for (Conversion conv : acg.getConversions()) {
                        lines.add("conversion"
                        		+ " node2=" + conv.getNode1().getNr()
                                + " node3=" + conv.getNode2().getNr()
                                + " affectedBlocks=" + blockSet.getAffectedBlocks(conv));
                    }

                    return lines;
                }
            });

            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                nexusBuilder.write(pstream);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Use coalescent model to simulate clonal frame.
     */
    private void simulateClonalFrame() {
        // Initialize leaf nodes
        List<Node> leafNodes = new ArrayList<>();
        for (int i=0; i<acg.m_taxonset.get().getTaxonCount(); i++) {
            Node leaf = new Node();
            leaf.setNr(i);
            leaf.setID(acg.m_taxonset.get().getTaxonId(i));
                        
            if (acg.hasDateTrait())
                leaf.setHeight(acg.getDateTrait().getValue(leaf.getID()));
            else
                leaf.setHeight(0.0);
            
            leafNodes.add(leaf);
        }
        
        // Create and sort list of inactive nodes
        List<Node> inactiveNodes = new ArrayList<>(leafNodes);
        Collections.sort(inactiveNodes, (Node n1, Node n2) -> {
            if (n1.getHeight()<n2.getHeight())
                return -1;
            
            if (n1.getHeight()>n2.getHeight())
                return 1;
            
            return 0;
        });
        
        List<Node> activeNodes = new ArrayList<>();
        
        double tau = 0.0;
        int nextNr = leafNodes.size();
        while (true) {
            
            // Calculate coalescence propensity
            int k = activeNodes.size();
            double chi = 0.5*k*(k-1);
            
            // Draw scaled coalescent time
            if (chi>0.0)
                tau += Randomizer.nextExponential(chi);
            else
                tau = Double.POSITIVE_INFINITY;
            
            // Convert to real time
            double t = popFunc.getInverseIntensity(tau);
            
            // If new time takes us past next sample time, insert that sample
            if (!inactiveNodes.isEmpty() && t>inactiveNodes.get(0).getHeight()) {
                Node nextActive = inactiveNodes.remove(0);
                activeNodes.add(nextActive);
                tau = popFunc.getIntensity(nextActive.getHeight());
                continue;
            }
            
            // Coalesce random pair of active nodes.
            Node node1 = activeNodes.remove(Randomizer.nextInt(k));
            Node node2 = activeNodes.remove(Randomizer.nextInt(k-1));
            
            Node parent = new Node();
            parent.addChild(node1);
            parent.addChild(node2);
            parent.setHeight(t);
            parent.setNr(nextNr++);
            
            activeNodes.add(parent);
            
            if (inactiveNodes.isEmpty() && activeNodes.size()<2)
                break;
        }
        
        // Remaining active node is root
        acg.setRoot(activeNodes.get(0));
    }
    
    private void generateConversions() {

        // Draw number of conversions:
    	double nConvMean = conversionRate * acg.getClonalFramePairedLength();
        int nConv = (int) Randomizer.nextPoisson(nConvMean);

        // Generate conversions:
        for (int i=0; i<nConv; i++) {
        	
            Conversion conv = acg.addNewConversion();
            associateConversionWithCF(conv);
            
            // Choose affected blocks:
            int nAffected = sampleBinomial(blockSet.getBlockCount(), moveProb);

            assert blockSet.getAffectedBlocks(conv).isEmpty();
            
            int[] shuffledBlockIdxs = Randomizer.shuffled(nAffected);
            for (int j=0; j<nAffected; j++) {
            	int blockIdx = shuffledBlockIdxs[j];
            	blockSet.getBlocks().get(blockIdx).addMove(conv);
            }
        }
    }
    
    private int sampleBinomial(int n, double p) {
    	int sum = 0;
    	for (int i =0; i<n; i++) 
    		if (Randomizer.nextDouble() < p) 
    			sum += 1;
    	return sum;
    }
    
    
    /**
     * Associates recombination with the clonal frame, selecting points of
     * departure and arrival.
     * 
     * @param conv recombination to associate
     */
    private void associateConversionWithCF(Conversion conv) {
    	CFEventList cfEventList = acg.getCFEventList();
    	List<Event> cfEvents = cfEventList.getCFEvents();
    	
        // Choose event interval
    	double[] intervalVolumes = cfEventList.getIntervalVolumes();
    	int iEvent = Util.sampleCategorical(intervalVolumes);
        Event event = cfEvents.get(iEvent);
        
    	// Choose height within interval
        double height = Randomizer.uniform(event.getHeight(), cfEvents.get(iEvent+1).getHeight());
    	conv.setHeight(height);
    	
    	// Choose source lineage (given the height)
    	Set<Node> activeLineages = acg.getLineagesAtHeight(height);
    	Node node1 = Util.sampleFrom(activeLineages);
    	conv.setNode1(node1);
    	assert node1.getHeight() < height;
        
        // Choose destination lineage (given the height and node1)
        activeLineages.remove(node1);
        Node node2 = Util.sampleFrom(activeLineages);
        conv.setNode2(node2);
        
        // Some validity checks...
        assert conv.getNode1() != null;
        assert conv.getNode2() != null;
        assert conv.getNode1() != conv.getNode2();
        assert !conv.getNode1().isRoot();
        assert !conv.getNode2().isRoot();
        assert conv.isValid();
    }

    public void rescaleCF(ConversionGraph acg, double targetHeight) {
    	double scale = targetHeight / acg.getRoot().getHeight();
    	for (Node node : acg.getNodesAsArray()) {
    		node.setHeight(scale*node.getHeight());
    	}
    }

    public ConversionGraph getACG() {
    	return acg;
    }
    
    public BlockSet getBlockSet() {
    	return blockSet;
    }
    
}
