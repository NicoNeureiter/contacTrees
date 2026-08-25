package contactrees.test;

import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.spec.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.evolution.tree.TreeParser;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import contactrees.MarginalNode;
import contactrees.MarginalTree;

/**
 * Unit test for marginal tree traversal.
 *
 * @author Nico Neureiter
 */
public class MarginalTreeTest extends ContactreesTest {

    public MarginalTreeTest() {
    }

    @Test
    public void testNonOverlapping() throws Exception {

        // Test all marginals against truth
        // (I have eyeballed each of these trees and claim that they are correct.)
        String[] correctNewickStrings = {
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
            "(1:2.5,(2:0.5,3:0.5)4:2.0)5:0.5",  // Conv 1
            "((1:1.0,2:1.0)4:0.5,3:1.5)5:1.5",  // Conv 2
            "(1:1.5,(2:0.5,3:0.5)4:1.0)5:1.5",  // Conv 1 & 2
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
            "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5",  // No conv
        };

        List<Block> blocks = blockSet.getBlocks();
        blocks.get(1).addMove(conv1);
        blocks.get(2).addMove(conv2);
        blocks.get(3).addMove(conv1);
        blocks.get(3).addMove(conv2);

        for (int b=0; b<N_BLOCKS; b++) {
            MarginalTree marginalTree = new MarginalTree();
            marginalTree.initByName("network", acg, "block", blocks.get(b), "nodetype", MarginalNode.class.getName());

            String newickStr = correctNewickStrings[b];
            Tree correctTree = new TreeParser(newickStr, false, true, false, 1);
            assertTrue(treesEquivalent(marginalTree, correctTree, 1e-15));

            equalLikelihood(correctTree, marginalTree);
        }
    }

    @Test
    public void testNonOverlapping_2() throws Exception {


        // Test all marginals against truth
        // (I have eyeballed each of these trees and claim that they are correct.)
        String[] correctNewickStrings = {
            "(((1:1.0,2:1.0)6:1.5,3:2.5)8:1.0,(4:1.5,5:1.5)7:2.0)9:0.5;",  // No conv
            "((1:2.5,(2:0.5,3:0.5)6:2.0)8:1.0,(4:1.5,5:1.5)7:2.0)9:0.5",  // Conv 1
            "(((1:1.0,2:1.0)6:0.5,3:1.5)8:2.0,(4:1.5,5:1.5)7:2.0)9:0.5",  // Conv 2
            "((1:1.0,2:1.0)6:2.5,((3:1.0,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 3
            "((1:1.5,(2:0.5,3:0.5)6:1.0)8:2.0,(4:1.5,5:1.5)7:2.0)9:0.5",  // Conv 1 & 2
            "(1:3.5,(((2:0.5,3:0.5)6:0.5,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 1 & 3
            "((1:1.0,2:1.0)6:2.5,((3:1.0,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 2 & 3
            "(1:3.5,(((2:0.5,3:0.5)6:0.5,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 1 & 2 & 3
        };

        List<Block> blocks = blockSet2.getBlocks();

        blocks.get(1).addMove(conv2_1);

        blocks.get(2).addMove(conv2_2);

        blocks.get(3).addMove(conv2_3);

        blocks.get(4).addMove(conv2_1);
        blocks.get(4).addMove(conv2_2);

        blocks.get(5).addMove(conv2_1);
        blocks.get(5).addMove(conv2_3);

        blocks.get(6).addMove(conv2_2);
        blocks.get(6).addMove(conv2_3);

        blocks.get(7).addMove(conv2_1);
        blocks.get(7).addMove(conv2_2);
        blocks.get(7).addMove(conv2_3);

        for (int b=0; b<N_BLOCKS; b++) {
            System.out.println(b);
            MarginalTree marginalTree = new MarginalTree();
            marginalTree.initByName("network", acg2, "block", blocks.get(b), "nodetype", MarginalNode.class.getName());

            String newickStr = correctNewickStrings[b];
            Tree correctTree = new TreeParser(newickStr, false, true, false, 1);

            assertTrue(treesEquivalent(marginalTree, correctTree, 1e-15));
            equalLikelihood(correctTree, marginalTree);
        }
    }

    @Test
    public void testNonOverlapping_3() throws Exception {


        // Test all marginals against truth
        // (I have eyeballed each of these trees and claim that they are correct.)
        String[] correctNewickStrings = {
            "(((1:1.0,2:1.0)6:1.5,3:2.5)8:1.0,(4:1.5,5:1.5)7:2.0)9:0.5;",  // No conv
            "((1:2.5,(2:0.5,3:0.5)6:2.0)8:1.0,(4:1.5,5:1.5)7:2.0)9:0.5",  // Conv 1
            "(((1:1.0,2:1.0)6:0.5,3:1.5)8:2.0,(4:1.5,5:1.5)7:2.0)9:0.5",  // Conv 2
            "((1:1.0,2:1.0)6:2.5,((3:1.0,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 3
            "((1:1.5,(2:0.5,3:0.5)6:1.0)8:2.0,(4:1.5,5:1.5)7:2.0)9:0.5",  // Conv 1 & 2
            "(1:3.5,(((2:0.5,3:0.5)6:0.5,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 1 & 3
            "((1:1.0,2:1.0)6:2.5,((3:1.0,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 2 & 3
            "(1:3.5,(((2:0.5,3:0.5)6:0.5,4:1.0)8:0.5,5:1.5)7:2.0)9:0.5",  // Conv 1 & 2 & 3
        };

        List<Block> blocks = blockSet2.getBlocks();
        List<List<Conversion>> convsByBlock = new ArrayList<>();

        convsByBlock.add(new ArrayList<>());

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(1).add(conv2_1);

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(2).add(conv2_2);

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(3).add(conv2_3);

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(4).add(conv2_1);
        convsByBlock.get(4).add(conv2_2);

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(5).add(conv2_1);
        convsByBlock.get(5).add(conv2_3);

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(6).add(conv2_2);
        convsByBlock.get(6).add(conv2_3);

        convsByBlock.add(new ArrayList<>());
        convsByBlock.get(7).add(conv2_1);
        convsByBlock.get(7).add(conv2_2);
        convsByBlock.get(7).add(conv2_3);

        for (int b=0; b<N_BLOCKS; b++) {
            System.out.println(b);
            Block block = blocks.get(b);

            MarginalTree marginalTree = new MarginalTree();
            marginalTree.initByName("network", acg2, "block", blocks.get(b), "nodetype", MarginalNode.class.getName());

            // Add conversions to block
            for (Conversion conv : convsByBlock.get(b))
                block.addMove(conv);

            // Update marginal tree
            marginalTree.requiresRecalculation();

            //System.out.println(marginalTree + ";");
            String newickStr = correctNewickStrings[b];
            System.out.println(marginalTree.toString());
            System.out.println(newickStr);
            Tree correctTree = new TreeParser(newickStr, false, true, false, 1);

            assertTrue(treesEquivalent(marginalTree, correctTree, 1e-15));
            equalLikelihood(correctTree, marginalTree);
        }
    }

    @Test
    public void testBranchRates() throws Exception {
        double[] rates = {2., 2., 1., 1.};
        UCRelaxedClockModel clock = new UCRelaxedClockModel();
        clock.initByName(
                "rates", new RealVectorParam<>(rates, PositiveReal.INSTANCE),
                "distr", new Uniform(),
                "tree", acg
                );

        // Test all marginals against truth
        // (I have eyeballed each of these trees and claim that they are correct.)
        String[] correctNewickStrings = {
            "((1:2.0,2:2.0)4:1.5,3:2.5)5:0.0",  // No conv
            "(1:3.5,(2:1.0,3:0.5)4:2.0)5:0.0",  // Conv 1
            "((1:2.0,2:2.0)4:0.5,3:1.5)5:0.0",  // Conv 2
            "(1:2.5,(2:1.0,3:0.5)4:1.0)5:0.0",  // Conv 1 & 2
            "((1:2.0,2:2.0)4:1.5,3:2.5)5:0.0",  // No conv
            "((1:2.0,2:2.0)4:1.5,3:2.5)5:0.0",  // No conv
            "((1:2.0,2:2.0)4:1.5,3:2.5)5:0.0",  // No conv
            "((1:2.0,2:2.0)4:1.5,3:2.5)5:0.0",  // No conv
            "((1:2.0,2:2.0)4:1.5,3:2.5)5:0.0",  // No conv
        };

        List<Block> blocks = blockSet.getBlocks();
        blocks.get(1).addMove(conv1);
        blocks.get(2).addMove(conv2);
        blocks.get(3).addMove(conv1);
        blocks.get(3).addMove(conv2);

        for (int b=0; b<N_BLOCKS; b++) {
            MarginalTree marginalTree = new MarginalTree();
            marginalTree.initByName(
                    "network", acg,
                    "block", blocks.get(b),
                    "nodetype", MarginalNode.class.getName(),
                    "branchRateModel", clock);

            String newickStr = correctNewickStrings[b];
            Tree correctTree = new TreeParser(newickStr, false, true, false, 1);

            assertTrue(treesEquivalentShifted(marginalTree, correctTree, 1e-15));
            equalLikelihood(correctTree, marginalTree);
        }
    }

    /**
     * A conversion sitting exactly at the height of a CF coalescence must be processed
     * below that coalescence, i.e. give the same marginal tree as one placed just below it.
     */
    @Test
    public void testConversionAtCoalescenceHeightMatchesLimitFromBelow() throws Exception {
        String newick = "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5;";

        MarginalTree exactMarginalTree = marginalTreeWithSingleConversion(newick, 0.0);
        MarginalTree belowMarginalTree = marginalTreeWithSingleConversion(newick, -1e-12);

        assertTrue(treesEquivalent(exactMarginalTree, belowMarginalTree, 1e-8));
    }

    /**
     * Build a marginal tree over a single block whose only conversion attaches node3 to node2
     * at the height of node4, offset by the given amount.
     */
    private MarginalTree marginalTreeWithSingleConversion(String newick, double heightOffset) {
        ConversionGraph localAcg = getACGFromNewick(newick);
        BlockSet localBlockSet = getBlockSet(1, localAcg);

        Node localRoot = localAcg.getRoot();
        Node localNode4 = localRoot.getLeft();
        Node localNode3 = localRoot.getRight();
        Node localNode2 = localNode4.getRight();

        Conversion conv = new Conversion(localNode3, localNode2, localNode4.getHeight() + heightOffset, localAcg, 1);
        localAcg.addConversion(conv);
        localBlockSet.getBlocks().get(0).addMove(conv);

        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName(
                "network", localAcg,
                "block", localBlockSet.getBlocks().get(0),
                "nodetype", MarginalNode.class.getName());
        return marginalTree;
    }

    /**
     * A block referring to a degenerate (self-loop) conversion must fail with a clear error
     * rather than silently producing a broken marginal tree.
     */
    @Test
    public void testDegenerateLoopConversionThrowsClearError() throws Exception {
        String newick = "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5;";

        ConversionGraph loopAcg = getACGFromNewick(newick);
        BlockSet loopBlockSet = getBlockSet(1, loopAcg);
        Node loopNode4 = loopAcg.getRoot().getLeft();
        Conversion loopConversion = new Conversion(loopNode4, loopNode4, 2.0, loopAcg, 1);
        loopAcg.addConversion(loopConversion);
        loopBlockSet.getBlocks().get(0).addMove(loopConversion);

        RuntimeException exception = assertThrows(RuntimeException.class, () -> {
            MarginalTree loopMarginalTree = new MarginalTree();
            loopMarginalTree.initByName(
                    "network", loopAcg,
                    "block", loopBlockSet.getBlocks().get(0),
                    "nodetype", MarginalNode.class.getName());
        });

        assertTrue(exception.getMessage().contains("invalid conversion"));
    }

    public boolean equalLikelihood(Tree correctTree, MarginalTree derivedTree) {
        Alignment alignment = getAlignment(correctTree.getLeafNodeCount());

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "substModel", jc);

        // Likelihood
        TreeLikelihood correctLikelihood = new TreeLikelihood();
        correctLikelihood.initByName(
                "data", alignment,
                "tree", correctTree,
                "siteModel", siteModel);
        TreeLikelihood derivedLikelihood = new TreeLikelihood();
        derivedLikelihood.initByName(
                "data", alignment,
                "tree", derivedTree,
                "siteModel", siteModel);

        correctTree.setEverythingDirty(true);
        derivedTree.setEverythingDirty(true);

        double logPtrue = correctLikelihood.calculateLogP();
        double logP = derivedLikelihood.calculateLogP();

        double diff = Math.abs(logPtrue - logP);

        assertTrue(diff < 0.001);

        return true;
    }

}
