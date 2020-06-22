/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

package contactrees.acgannotator;

import contactrees.ACGWithBlocks;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.acgannotator.ACGCladeSystem.BitSetPair;
import beast.app.treeannotator.CladeSystem;
import beast.evolution.tree.Node;

import java.util.*;
import java.util.function.BiFunction;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.Range;

/**
 * Adds conversion summary tools to CladeSystem.
 * 
 * Nico Neureiter
 */
public class ACGCladeSystem extends CladeSystem {

    protected ArrayListMultimap<BitSetPair, ConversionWithBlocks> conversionLists = ArrayListMultimap.create();
    protected List<Multiset<BitSetPair>> geneFlow = new ArrayList<>();
    protected BitSet[] bitSets;
    protected int nBlocks = -1;

    protected int acgIndex = 1;

    public ACGCladeSystem() { }

    public ACGCladeSystem(ConversionGraph acg) {
        add(acg, true);
    }

    /**
     * Assemble list of bitSets for this ACG.
     */
    public BitSet[] getBitSets(ConversionGraph acg) {

        if (bitSets == null)
            bitSets = new BitSet[acg.getNodeCount()];

        applyToClades(acg.getRoot(), (cladeNode, bits) -> {
            bitSets[cladeNode.getNr()] = bits;
            return null;
        });

        return bitSets;
    }

    /**
     * Add conversions described on provided acg to the internal list
     * for later summary.
     *
     * @param acg conversion graph from which to extract conversions
     */
    public void collectConversions(ACGWithBlocks acg) {
        getBitSets(acg);
        BlockSet blockSet = acg.blockSet;
        nBlocks = blockSet.getBlockCount();
        Multiset<BitSetPair> geneFlowSample = HashMultiset.create();
                
        for (Conversion conv : acg.getConversions())  {
            conv.acgIndex = acgIndex;
            BitSetPair bsPair = new BitSetPair(conv);

            conversionLists.get(bsPair).add(
                    new ConversionWithBlocks(conv, blockSet.getAffectedBlockIDs(conv))
            );

            // Record gene flow
            geneFlowSample.add(bsPair, blockSet.countAffectedBlocks(conv));
//            if (!geneFlowSample.containsKey(bsPair.from))
//                geneFlowSample.put(bsPair.from, new HashMap<>());
//            long oldFlow = geneFlowSample.get(bsPair.from).getOrDefault(bsPair.to, 0L);
//            geneFlowSample.get(bsPair.from).put(bsPair.to, oldFlow + blockSet.countAffectedBlocks(conv));
        }

        geneFlow.add(geneFlowSample);

        acgIndex += 1;
    }

    /**
     * Determine contiguous regions on specified locus where the fraction of
     * ACGs having a conversion active is greater than the given threshold.
     *
     * @param cladePair representing the source and destination clades
     * @param Total number of samples
     * @param threshold minimum fraction of sampled conversions included
     * @return List of regions
     */
    public List<ConversionSummary> getConversionSummaries(BitSetPair cladePair,
                                                          int nACGs,
                                                          double threshold) {

        List<ConversionSummary> convSummaryList = new ArrayList<>(nBlocks);
        for (int i=0; i<nBlocks; i++) {
            convSummaryList.add(new ConversionSummary());
        }
        
        // Return empty list if on conversions meet the criteria.
        if (!conversionLists.containsKey(cladePair))
            return convSummaryList;

//        int thresholdCount = (int)Math.ceil(nACGs*threshold);
//        ConversionSummary conversionSummary = null;

        for (ConversionWithBlocks convWB : conversionLists.get(cladePair)) {
            for (int blockID : convWB.blocks) {
                convSummaryList.get(blockID).addConv(convWB.conversion);
            }
        }

        return convSummaryList;
    }

    /**
     * @return list of multisets specifying gene flow between clades in each sample.
     */
    public List<Multiset<BitSetPair>> getGeneFlowMap() {
        return geneFlow;
    }
    
    /**
     * @return a multiset specifying the gene flow between clades.
     */
    public Multiset<BitSetPair> getGeneFlow() {
        Multiset<BitSetPair> totalGeneFlow = HashMultiset.create();
        for (Multiset<BitSetPair> geneFlowSample : geneFlow) {
            for (BitSetPair cladePair : geneFlowSample) {
                totalGeneFlow.add(cladePair, geneFlowSample.count(cladePair));
            }
        }
        return totalGeneFlow;
    }
    
    public BitSetPair createCladePair(BitSet from, BitSet to) {
        return new BitSetPair(from, to);
    }

    /**
     * Apply a function to each sub-clade.
     *
     * @param node MRCA of clade
     * @param function function to apply. Given sub-clade parent node
     *                 and bitset as arguments.
     * @return BitSet representing clade.
     */
    public BitSet applyToClades(Node node, BiFunction<Node, BitSet, Void> function) {
        BitSet bits = new BitSet();

        if (node.isLeaf()) {
            bits.set(2 * getTaxonIndex(node));
        } else {
            for (Node child : node.getChildren())
                bits.or(applyToClades(child, function));
        }

        function.apply(node, bits);

        return bits;
    }

    public List<BitSetPair> listCladePairs(ConversionGraph acg){
        getBitSets(acg);
        List <BitSetPair> cladePairs = new ArrayList<>();
        
        for (BitSet from : bitSets) {
            for (BitSet to : bitSets) {
                if (from == to) {
                    continue;
                }
                cladePairs.add(createCladePair(from, to));
            }
        }
        
        return cladePairs;
    }


    /**
     * Class representing an ordered pair of BitSets.
     */
    protected class BitSetPair {
        public BitSet from, to;

        public BitSetPair(BitSet from, BitSet to) {
            this.from = from;
            this.to = to;
        }

        public BitSetPair(Conversion conv) {
            this.from = bitSets[conv.getNode1().getNr()];
            this.to = bitSets[conv.getNode2().getNr()];
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            BitSetPair that = (BitSetPair) o;

            return from.equals(that.from) && to.equals(that.to);

        }

        @Override
        public int hashCode() {
            int result = from.hashCode();
            result = 31 * result + to.hashCode();
            return result;
        }

        @Override
        public String toString() {
            return from.toString() + " -> " + to.toString();
        }
    }

    /**
     * Class representing a summary of similar conversions between two
     * points in the summarized clonal frame.
     */
    public class ConversionSummary {
        
        List<Double> heights = new ArrayList<>();

        public int nIncludedACGs = 0;
        
        ConversionGraph acg;

        /**
         * Add metrics associated with given conversion to summary.
         *
         * @param conv conversion
         */
        public void addConv(Conversion conv) {
            heights.add(conv.getHeight());
        }

        /**
         * Add metrics associated with each of the conversions in the
         * given list to the summary.
         *
         * @param convs list of conversions
         */
        public void addConvs(List<Conversion> convs) {
            for (Conversion conv : convs)
                addConv(conv);
        }

        /**
         * @return number of conversions included in summary.
         */
        public int summarizedConvCount() {
            return heights.size();
        }

        public List<Double> getHeights() {
            return heights;
        }
    }
    
    public class ConversionWithBlocks {
        Conversion conversion;
        List<Integer> blocks;
        
        public ConversionWithBlocks(Conversion conversion, List<Integer> blocks) {
            this.conversion = conversion;
            this.blocks = blocks;
        }
    }
}
