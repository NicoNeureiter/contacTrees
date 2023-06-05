package contactrees.operators;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.util.Pair;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import contactrees.Block;
import contactrees.Conversion;

/**
 * Split: Split one conversion into two conversions between the same branches.
 * Merge: Combine two conversions between the same branches into one.
 *
 * @author Nico Neureiter
 */
@Description("Operator which splits or merges conversions")
public class ConversionSplit extends ConversionCreationOperator {

    final public Input<Boolean> flipInput = new Input<>(
            "flip",
            "Flip the the duplicated conversion edge to point in the oposite direction (compared to the original one).",
            false
            );

    @Override
    public double proposal() {

        if (Randomizer.nextBoolean())
            return splitProposal();
        else
            return mergeProposal();
    }

    protected double splitProposal() {
        double logHGF = 0;
        int nConv = acg.getConvCount();

        // The split move requires at least one existing conversion
        if (nConv == 0)
            return Double.NEGATIVE_INFINITY;

        // Chose the conversion to split
        Conversion conv = acg.getConversions().getRandomConversion();
        // Create the duplicate
        Conversion convNew = acg.addDuplicateConversion(conv);
        // Update Hastings ratio:
        //        logHGF -= Math.log(1.0/nConv);
        //        logHGF += Math.log(1.0/((nConv+1)*nConv));
        logHGF -= Math.log(nConv+1);



        if (flipInput.get())
            flipDirection(convNew);

        // Sample new height uniformly on the current lineages
        double hMin = Math.max(conv.getNode1().getHeight(),
                               conv.getNode2().getHeight());
        double hMax = Math.min(conv.getNode1().getParent().getHeight(),
                               conv.getNode2().getParent().getHeight());
        double heightNew = Randomizer.uniform(hMin, hMax);
        convNew.setHeight(heightNew);

        // Update Hastings ratio for sampled height
        logHGF -= Math.log(1.0/(hMax - hMin));

        assert hMin <= conv.getHeight() && conv.getHeight() <= hMax :
            "Conversion height must be in validHeightRange";


        // Choose which blocks-moves stay on conv and which are moved to convNew
        double pDoubleMove = pMoveInput.get().getValue();
        List<Block> blocksToRemove = new LinkedList<Block>();
        for (Block block : blockSet.getAffectedBlocks(conv)) {
            if (Randomizer.nextDouble() < pDoubleMove) {
                // Put ´block´ on both edges with probability ´pDoubleMove´
                blockSet.addBlockMove(convNew, block);

                logHGF -= Math.log(pDoubleMove);

            } else {
                // Otherwise decide whether ´block´ is moved to convNew or stays on conv
                if (Randomizer.nextBoolean()) {
                    blockSet.addBlockMove(convNew, block);
                    blocksToRemove.add(block);
                }

                logHGF -= Math.log(1 - pDoubleMove);
                logHGF -= Math.log(0.5);
            }
        }

        for (Block block : blocksToRemove)
            blockSet.removeBlockMove(conv, block);

        assert acg.getConvCount() == nConv + 1;
        assert !acg.isInvalid();

        return logHGF;
    }

    protected double mergeProposal() {
        double logHGF = 0;
        int nConv = acg.getConvCount();

        // The merge move requires at least two existing conversions
        if (nConv < 2)
            return Double.NEGATIVE_INFINITY;

        // Get all candidate conversion pairs (with same start/end nodes)
        ArrayList<Pair<Conversion, Conversion>> convPairs = getConversionPairs();
        if (convPairs.isEmpty())
            return Double.NEGATIVE_INFINITY;

        int iPair = Randomizer.nextInt(convPairs.size());
        Conversion conv1 = convPairs.get(iPair).getFirst();
        Conversion conv2 = convPairs.get(iPair).getSecond();

        // Update
        logHGF -= Math.log(1.0 / convPairs.size());

//        // Choose two random conversions
//        Conversion conv1 = acg.getConversions().getRandomConversion();
//        Conversion conv2;
//        do {
//            conv2 = acg.getConversions().getRandomConversion();
//        } while (conv2 == conv1);

        if (flipInput.get())
            flipDirection(conv2);

        // Operator fails if conversions don't share the same nodes
        if (conv1.getNode1() != conv2.getNode1() || conv1.getNode2() != conv2.getNode2()) {
            assert false; // NN: This should never happen, since we are only consider conversion pairs with same nodes
            return Double.NEGATIVE_INFINITY;
        }

//        // Update Hastings ratio
//        //        logHGF -= Math.log(1.0/(nConv*(nConv-1)));
//        //        logHGF += Math.log(1.0/(nConv-1));
//        logHGF += Math.log(nConv);


        // Update Hastings ratio with back-probability for height
        double hMin = Math.max(conv1.getNode1().getHeight(),
                               conv1.getNode2().getHeight());
        double hMax = Math.min(conv1.getNode1().getParent().getHeight(),
                               conv1.getNode2().getParent().getHeight());
        logHGF += Math.log(1.0/(hMax - hMin));




        // Update Hastings ratio with back-probability for sampled block moves
        double pDoubleMove = pMoveInput.get().getValue();

        for (Block block : blockSet.getAffectedBlocks(conv1)) {
            if (block.isAffected(conv2)) {
                // Blocks affected by conv1 and conv2  ->  reverse move is a duplication move
                logHGF += Math.log(pDoubleMove);
            } else {
                // Blocks affected by conv1 but not conv2  ->  reverse move is no duplication
                logHGF += Math.log(1 - pDoubleMove);
                logHGF += Math.log(0.5);
            }
        }

        for (Block block : blockSet.getAffectedBlocks(conv2)) {
            if (!block.isAffected(conv1)) {
                // Blocks affected by conv1 but not conv2  ->  reverse move is no duplication
                logHGF += Math.log(1 - pDoubleMove);
                logHGF += Math.log(0.5);

                blockSet.addBlockMove(conv1, block);
            }
//            blockSet.removeBlockMove(conv2, block);
        }

        // Remove conv2
        blockSet.removeConversion(conv2);
        acg.removeConversion(conv2);

        assert acg.getConvCount() == nConv - 1;
        assert !acg.isInvalid();

        return logHGF;
    }

    private ArrayList<Pair<Conversion, Conversion>> getConversionPairs() {
        int n = acg.getNodeCount();
        ArrayList<ArrayList<ArrayList<Conversion>>> convsByNodes = new ArrayList<>(n);

        // Initialize convsByNodes
        for (int i=0; i<n; i++) {
            convsByNodes.add(i, new ArrayList<>(n));
            for (int j=0; j<n; j++) {
                convsByNodes.get(i)
                            .add(j, new ArrayList<>());
            }
        }

        ArrayList<Pair<Conversion, Conversion>> conversionPairs = new ArrayList<>();

        // Collect all pairs of conversion with same nodes
        for (Conversion conv : acg.getConversions()) {
            ArrayList<Conversion> convs = convsByNodes.get(conv.getNode1().getNr())
                                                      .get(conv.getNode2().getNr());

            ArrayList<Conversion> others;
            if (flipInput.get())
                others = convsByNodes.get(conv.getNode2().getNr())
                                     .get(conv.getNode1().getNr());
            else
                others = convs;

            for (Conversion other : others) {
                conversionPairs.add(new Pair<Conversion, Conversion>(conv, other));
                conversionPairs.add(new Pair<Conversion, Conversion>(other, conv));
            }

            convs.add(conv);
        }

        return conversionPairs;
    }

    /**
     * Flip direction (i.e. swap node1 and node2) of the given conversion
     * @param conv
     */
    private void flipDirection(Conversion conv) {
        Node node1 = conv.getNode1();
        Node node2 = conv.getNode2();
        conv.setNode1(node2);
        conv.setNode2(node1);
    }

}
