package contactrees.test;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Tree;
import beast.math.distributions.Uniform;
import beast.util.TreeParser;
import contactrees.Block;
import contactrees.Conversion;
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
        Double[] rates = {2., 2., 1., 1.};
        UCRelaxedClockModel clock = new UCRelaxedClockModel();
        clock.initByName(
                "rates", new RealParameter(rates),
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
