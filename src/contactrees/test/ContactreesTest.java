/**
 * 
 */
package contactrees.test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Before;

import beast.core.Logger;
import beast.core.MCMC;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contactrees.ACGWithMetaDataLogger;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;

/**
 * 
 * 
 * @author Nico Neureiter <n.neureiter.research@gmail.com>
 */
public class ContactreesTest {

	protected static final double EPS = 1E-8;
	
//	List<ConversionGraph> acgs;
//	List<TestBlockSet> blockSets;
//	List<TestACGWithMetaDataLogger> acgLoggers;

	protected ConversionGraph acg, acg2;
	protected BlockSet blockSet, blockSet2;
	protected ACGWithMetaDataLogger acgLogger, acgLogger2;
    protected Node node1, node2, node3, node4, root;
    protected Node root2, node2_1, node2_2, node2_3, node2_4, node2_5, node2_6, node2_7, node2_8;
	protected Conversion conv1, conv2, conv2_1, conv2_2, conv2_3;
	public int N_BLOCKS = 8;
	
	@Before
	public void setup() {
		buildACG_1();
		buildACG_2();
		assertAllValid();
	}
	
	public void assertAllValid() {
	    assertFalse(acg.isInvalid());
	    assertFalse(acg2.isInvalid());
	}	
    
    public void buildACG_1() {
        //
        //      Height             Topology
        //
        //       3.0                   |
        //       2.5              _____5_____
        //       2.0             |           |
        //       1.5             |  <------  |
        //       1.0          ___4___        |
        //       0.5         |       |  -->  |
        //       0.0         1       2       3
        //
        String newick = "((1:1.0,2:1.0)4:1.5,3:2.5)5:0.5;";
        setupFromNewick(newick);

        // Assign nodes
        root = acg.getRoot();
        node4 = root.getLeft();
        node3 = root.getRight();
        node2 = node4.getRight();
        node1 = node4.getLeft();
        
        // Add conversion
        conv1 = new Conversion(node2, node3, 0.5, acg, 1);
        conv2 = new Conversion(node3, node4, 1.5, acg, 2);
        acg.addConversion(conv1);
        acg.addConversion(conv2); 
    }
    
    public void buildACG_2() {
        //
        //      Height             Topology
        //       
        //       4.0                       |
        //       3.5                 ______9______
        //       3.0                |             |
        //       2.5             ___8____         |
        //       2.0            |        |        |
        //       1.5            | <----- |      __7__
        //       1.0          __6__      | --> |     |
        //       0.5         |     | --> |     |     |
        //       0.0         1     2     3     4     5
        //
        String newick = "(((1:1.0,2:1.0)6:1.5,3:2.5)8:1.0,(4:1.5,5:1.5)7:2.0)9:0.5;";
        

        acg2 = getACGFromNewick(newick);
        blockSet2 = getBlockSet(N_BLOCKS, acg2);
        acgLogger2 = new ACGWithMetaDataLogger(acg2, blockSet2);       

        // Assign nodes
        root2 = acg2.getRoot();
        
        node2_8 = root2.getLeft();
        node2_7 = root2.getRight();
        
        node2_6 = node2_8.getLeft();
        node2_3 = node2_8.getRight();
        
        node2_4 = node2_7.getLeft();
        node2_5 = node2_7.getRight();
        
        node2_1 = node2_6.getLeft();
        node2_2 = node2_6.getRight();
        
        // Add conversion
        conv2_1 = new Conversion(node2_2, node2_3, 0.5, acg2, 1);
        conv2_2 = new Conversion(node2_3, node2_6, 1.5, acg2, 2);
        conv2_3 = new Conversion(node2_3, node2_4, 1.0, acg2, 3);
        
        acg2.addConversion(conv2_1);
        acg2.addConversion(conv2_2);
        acg2.addConversion(conv2_3); 
    }
	
	public void setupFromNewick(String newick) {
	    acg = getACGFromNewick(newick);
		blockSet = getBlockSet(N_BLOCKS, acg);
		acgLogger = new ACGWithMetaDataLogger(acg, blockSet);		
	}
	
	public ConversionGraph getACGFromNewick(String newick) {
        ConversionGraph acg = new ConversionGraph();
        acg.initAndValidate();
        acg.fromExtendedNewick(newick);
        return acg;
	}
	
	/**
	 * Create a BlockSet with nBlocks new blocks.
	 * 
	 * @param nBlocks
	 * @param acg
	 * @return The new BlockSet
	 */
	public BlockSet getBlockSet(int nBlocks, ConversionGraph acg) {

        List<Block> blocks = new ArrayList<Block>();
        for (int i=0; i<nBlocks; i++) {
            Block b = new Block();
            b.initAndValidate();
            blocks.add(b);
        }
        
        blockSet = new BlockSet();
        blockSet.initByName("network", acg, "block", blocks);
        
        return blockSet;
	}
	

    /**
     * @return An Alignment object containing primate data for unit tests.
     * @throws Exception 
     */
    public Alignment getAlignment(int nLeafs) {
        List<Sequence> allSequences = new ArrayList<>();
        allSequences.add(new Sequence("1","AAGTTTCATTGGAGCCACCARTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTATTTTGCCTAGCAAATACAAACTACGAACGAGTCCACAGTCGAACAATAGCACTAGCCCGTGGCCTTCAAACCCTATTACCTCTTGCAGCAACATGATGACTCCTCGCCAGCTTAACCAACCTGGCCCTTCCCCCAACAATTAATTTAATCGGTGAACTGTCCGTAATAATAGCAGCATTTTCATGGTCACACCTAACTATTATCTTAGTAGGCCTTAACACCCTTATCACCGCCCTATATTCCCTATATATACTAATCATAACTCAACGAGGAAAATACACATATCATATCAACAATATCATGCCCCCTTTCACCCGAGAAAATACATTAATAATCATACACCTATTTCCCTTAATCCTACTATCTACCAACCCCAAAGTAATTATAGGAACCATGTACTGTAAATATAGTTTAAACAAAACATTAGATTGTGAGTCTAATAATAGAAGCCCAAAGATTTCTTATTTACCAAGAAAGTA-TGCAAGAACTGCTAACTCATGCCTCCATATATAACAATGTGGCTTTCTT-ACTTTTAAAGGATAGAAGTAATCCATCGGTCTTAGGAACCGAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATTTATTTTCATCCTCCATTTTACTATCACTTACACTCTTAATTACCCCATTTATTATTACAACAACTAAAAAATATGAAACACATGCATACCCTTACTACGTAAAAAACTCTATCGCCTGCGCATTTATAACAAGCCTAGTCCCAATGCTCATATTTCTATACACAAATCAAGAAATAATCATTTCCAACTGACATTGAATAACGATTCATACTATCAAATTATGCCTAAGCTT"));
        allSequences.add(new Sequence("2","AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTATTCTGTCTAGCCAACTCTAACTACGAACGAATCCATAGCCGTACAATACTACTAGCACGAGGGATCCAAACCATTCTCCCTCTTATAGCCACCTGATGACTACTCGCCAGCCTAACTAACCTAGCCCTACCCACCTCTATCAATTTAATTGGCGAACTATTCGTCACTATAGCATCCTTCTCATGATCAAACATTACAATTATCTTAATAGGCTTAAATATGCTCATCACCGCTCTCTATTCCCTCTATATATTAACTACTACACAACGAGGAAAACTCACATATCATTCGCACAACCTAAACCCATCCTTTACACGAGAAAACACCCTTATATCCATACACATACTCCCCCTTCTCCTATTTACCTTAAACCCCAAAATTATTCTAGGACCCACGTACTGTAAATATAGTTTAAA-AAAACACTAGATTGTGAATCCAGAAATAGAAGCTCAAAC-CTTCTTATTTACCGAGAAAGTAATGTATGAACTGCTAACTCTGCACTCCGTATATAAAAATACGGCTATCTCAACTTTTAAAGGATAGAAGTAATCCATTGGCCTTAGGAGCCAAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATCTATTATCCTCTTTCACCCTTGTCACACTGATTATCCTAACTTTACCTATCATTATAAACGTTACAAACATATACAAAAACTACCCCTATGCACCATACGTAAAATCTTCTATTGCATGTGCCTTCATCACTAGCCTCATCCCAACTATATTATTTATCTCCTCAGGACAAGAAACAATCATTTCCAACTGACATTGAATAACAATCCAAACCCTAAAACTATCTATTAGCTT"));
        allSequences.add(new Sequence("3","AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTA-CGACCCCTTATTTACCGAGAAAGCT-CACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT"));
        allSequences.add(new Sequence("4","AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATTCTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCTAACCTCGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAAGCCCTCATTCACACGAGAAAATACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCTATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCA-CGACCCCTTATTTACCGAGAAAGCT-TATAAGAACTGCTAATTCATATCCCCATGCCTGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACCACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCGTCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCTATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT"));
        allSequences.add(new Sequence("5","AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGGACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCTCGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTCCTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTCCCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAAACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCTATCCCTCAACCCCGATATTATCACCGGGTTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGATAACAGAGGCTCA-CAACCCCTTATTTACCGAGAAAGCT-CGTAAGAGCTGCTAACTCATACCCCCGTGCTTGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACCACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCATCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCTATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT"));
        allSequences.add(new Sequence("6","AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCATGGACTCACATCCTCCCTACTGTTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGGCCTTCAAACTCTACTCCCCCTAATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCTTGCCCTACCACCCACCATCAACCTTCTAGGAGAACTCTCCGTACTAATAGCCATATTCTCTTGATCTAACATCACCATCCTACTAACAGGACTCAACATACTAATCACAACCCTATACTCTCTCTATATATTCACCACAACACAACGAGGTACACCCACACACCACATCAACAACATAAAACCTTCTTTCACACGCGAAAATACCCTCATGCTCATACACCTATCCCCCATCCTCCTCTTATCCCTCAACCCCAGCATCATCGCTGGGTTCGCCTACTGTAAATATAGTTTAACCAAAACATTAGATTGTGAATCTAATAATAGGGCCCCA-CAACCCCTTATTTACCGAGAAAGCT-CACAAGAACTGCTAACTCTCACT-CCATGTGTGACAACATGGCTTTCTCAGCTTTTAAAGGATAACAGCTATCCCTTGGTCTTAGGATCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACAGCCATGTTTACCACCATAACTGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCTACCCTCATTAACCCCAACAAAAAAAACCCATACCCCCACTATGTAAAAACGGCCATCGCATCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGACAAGAAACCATCGTCACAAACTGATGCTGAACAACCACCCAGACACTACAACTCTCACTAAGCTT"));
        allSequences.add(new Sequence("7","AAGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATTCTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGGGCTCCAAGCCTTACTCCCACTGATAGCYTTCTGATGACTCGCAGCAAGCCTCGCTAACCTCGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTCCTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTCCCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAAACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCTAACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAACATTAGATTGTGAATCTAACAATAGAGGCTCG-AAACCTCTTGCTTACCGAGAAAGCC-CACAAGAACTGCTAACTCACTATCCCATGTATGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAGCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCCACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCCTCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACCATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT"));
        allSequences.add(new Sequence("8","AAGCTTTTCCGGCGCAACCATCCTTATGATCGCTCACGGACTCACCTCTTCCATATATTTCTGCCTAGCCAATTCAAACTATGAACGCACTCACAACCGTACCATACTACTGTCCCGAGGACTTCAAATCCTACTTCCACTAACAGCCTTTTGATGATTAACAGCAAGCCTTACTAACCTTGCCCTACCCCCCACTATCAATCTACTAGGTGAACTCTTTGTAATCGCAACCTCATTCTCCTGATCCCATATCACCATTATGCTAACAGGACTTAACATATTAATTACGGCCCTCTACTCTCTCCACATATTCACTACAACACAACGAGGAACACTCACACATCACATAATCAACATAAAGCCCCCCTTCACACGAGAAAACACATTAATATTCATACACCTCGCTCCAATTATCCTTCTATCCCTCAACCCCAACATCATCCTGGGGTTTACCTCCTGTAGATATAGTTTAACTAAAACACTAGATTGTGAATCTAACCATAGAGACTCA-CCACCTCTTATTTACCGAGAAAACT-CGCAAGGACTGCTAACCCATGTACCCGTACCTAAAATTACGGTTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGACCTTAGGAGTCAAAAACATTGGTGCAACTCCAAATAAAAGTAATAATCATGCACACCCCCATCATTATAACAACCCTTATCTCCCTAACTCTCCCAATTTTTGCCACCCTCATCAACCCTTACAAAAAACGTCCATACCCAGATTACGTAAAAACAACCGTAATATATGCTTTCATCATCAGCCTCCCCTCAACAACTTTATTCATCTTCTCAAACCAAGAAACAACCATTTGGAGCTGACATTGAATAATGACCCAAACACTAGACCTAACGCTAAGCTT"));
        allSequences.add(new Sequence("9","AAGCTTTTCTGGCGCAACCATCCTCATGATTGCTCACGGACTCACCTCTTCCATATATTTCTGCCTAGCCAATTCAAACTATGAACGCACTCACAACCGTACCATACTACTGTCCCGGGGACTTCAAATCCTACTTCCACTAACAGCTTTCTGATGATTAACAGCAAGCCTTACTAACCTTGCCCTACCCCCCACTATCAACCTACTAGGTGAACTCTTTGTAATCGCGACCTCATTCTCCTGGTCCCATATCACCATTATATTAACAGGATTTAACATACTAATTACGGCCCTCTACTCCCTCCACATATTCACCACAACACAACGAGGAGCACTCACACATCACATAATCAACATAAAACCCCCCTTCACACGAGAAAACATATTAATATTCATACACCTCGCTCCAATCATCCTCCTATCTCTCAACCCCAACATCATCCTGGGGTTTACTTCCTGTAGATATAGTTTAACTAAAACATTAGATTGTGAATCTAACCATAGAGACTTA-CCACCTCTTATTTACCGAGAAAACT-CGCGAGGACTGCTAACCCATGTATCCGTACCTAAAATTACGGTTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGACCTTAGGAGTCAAAAATATTGGTGCAACTCCAAATAAAAGTAATAATCATGCACACCCCTATCATAATAACAACCCTTATCTCCCTAACTCTCCCAATTTTTGCCACCCTCATCAACCCTTACAAAAAACGTCCATACCCAGATTACGTAAAAACAACCGTAATATATGCTTTCATCATCAGCCTCCCCTCAACAACTTTATTCATCTTCTCAAACCAAGAAACAACCATTTGAAGCTGACATTGAATAATAACCCAAACACTAGACCTAACACTAAGCTT"));
        allSequences.add(new Sequence("10","AAGCTTCTCCGGCGCAACCACCCTTATAATCGCCCACGGGCTCACCTCTTCCATGTATTTCTGCTTGGCCAATTCAAACTATGAGCGCACTCATAMCCGTACCATACTACTATCCCGAGGACTTCAAATTCTACTTCCATTGACAGCCTTCTGATGACTCACAGCAAGCCTTACTAACCTTGCCCTACCCCCCACTATTAATCTACTAGGCGAACTCTTTGTAATCACAACTTCATTTTCCTGATCCCATATCACCATTGTGTTAACGGGCCTTAATATACTAATCACAGCCCTCTACTCTCTCCACATGTTCATTACAGTACAACGAGGAACACTCACACACCACATAATCAATATAAAACCCCCCTTCACACGAGAAAACATATTAATATTCATACACCTCGCTCCAATTATCCTTCTATCTCTCAACCCCAACATCATCCTGGGGTTTACCTCCTGTAAATATAGTTTAACTAAAACATTAGATTGTGAATCTAACTATAGAGGCCTA-CCACTTCTTATTTACCGAGAAAACT-CGCAAGGACTGCTAATCCATGCCTCCGTACTTAAAACTACGGTTTCCTCAACTTTTAAAGGATAACAGCTATCCATTGACCTTAGGAGTCAAAAACATTGGTGCAACTCCAAATAAAAGTAATAATCATGCACACCCCCATCATAATAACAACCCTCATCTCCCTGACCCTTCCAATTTTTGCCACCCTCACCAACCCCTATAAAAAACGTTCATACCCAGACTACGTAAAAACAACCGTAATATATGCTTTTATTACCAGTCTCCCCTCAACAACCCTATTCATCCTCTCAAACCAAGAAACAACCATTTGGAGTTGACATTGAATAACAACCCAAACATTAGACCTAACACTAAGCTT"));
        allSequences.add(new Sequence("11","AAGCTTCTCCGGTGCAACTATCCTTATAGTTGCCCATGGACTCACCTCTTCCATATACTTCTGCTTGGCCAACTCAAACTACGAACGCACCCACAGCCGVATCATACTACTATCCCGAGGACTCCAAATCCTACTCCCACTAACAGCCTTCTGATGATTCACAGCAAGCCTTACTAATCTTGCTCTACCCTCCACTATTAATCTACTGGGCGAACTCTTCGTAATCGCAACCTCATTTTCCTGATCCCACATCACCATCATACTAACAGGACTGAACATACTAATTACAGCCCTCTACTCTCTTCACATATTCACCACAACACAACGAGGAGCGCTCACACACCACATAATTAACATAAAACCACCTTTCACACGAGAAAACATATTAATACTCATACACCTCGCTCCAATTATTCTTCTATCTCTTAACCCCAACATCATTCTAGGATTTACTTCCTGTAAATATAGTTTAATTAAAACATTAGACTGTGAATCTAACTATAGAAGCTTA-CCACTTCTTATTTACCGAGAAAACT-TGCAAGGACCGCTAATCCACACCTCCGTACTTAAAACTACGGTTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGCCTTAGGAGTCAAAAATATTGGTGCAACTCCAAATAAAAGTAATAATCATGTATACCCCCATCATAATAACAACTCTCATCTCCCTAACTCTTCCAATTTTCGCTACCCTTATCAACCCCAACAAAAAACACCTATATCCAAACTACGTAAAAACAGCCGTAATATATGCTTTCATTACCAGCCTCTCTTCAACAACTTTATATATATTCTTAAACCAAGAAACAATCATCTGAAGCTGGCACTGAATAATAACCCAAACACTAAGCCTAACATTAAGCTT"));
        allSequences.add(new Sequence("12","AAGCTTCACCGGCGCAATGATCCTAATAATCGCTCACGGGTTTACTTCGTHTATGCTATTCTGCCTAGCAAACTCAAATTACGAACGAATTCACAGCCGAACAA-BACH-TTACTCGAGGGCTCCAAACACTATTCCCGCTTATAGGCCTCTGATGACTCCTAGCAAATCTCGCTAACCTCGCCCTACCCACAGCTATTAATCTAGTAGGAGAATTACTCACAATCGTATCTTCCTTCTCTTGATCCAACTTTACTATTATATTCACAGGACTTAATATACTAATTACAGCACTCTACTCACTTCATATGTATGCCTCTACACAGCGAGGTCCACTTACATACAGCACCAGCAATATAAAACCAATATTTACACGAGAAAATACGCTAATATTTATACATATAACACCAATCCTCCTCCTTACCTTGAGCCCCAAGGTAATTATAGGACCCTCACCTTGTAATTATAGTTTAGCTAAAACATTAGATTGTGAATCTAATAATAGAAGAATA-TAACTTCTTAATTACCGAGAAAGTG-CGCAAGAACTGCTAATTCATGCTCCCAAGACTAACAACTTGGCTTCCTCAACTTTTAAAGGATAGTAGTTATCCATTGGTCTTAGGAGCCAAAAACATTGGTGCAACTCCAAATAAAAGTAATA---ATACACTTCTCCATCACTCTAATAACACTAATTAGCCTACTAGCGCCAATCCTAGCTACCCTCATTAACCCTAACAAAAGCACACTATACCCGTACTACGTAAAACTAGCCATCATCTACGCCCTCATTACCAGTACCTTATCTATAATATTCTTTATCCTTACAGGCCAAGAATCAATAATTTCAAACTGACACTGAATAACTATCCAAACCATCAAACTATCCCTAAGCTT"));
        
        List<Sequence> sequences = new ArrayList<>();
        for (int i=0; i<nLeafs; i++) {
            sequences.add(allSequences.get(i));
        }
        
        return new Alignment(sequences, "nucleotide");
    }
    
    /**
     * @return An Alignment object containing primate data for unit tests.
     * @throws Exception 
     */
    public Alignment getAlignment() {
        return getAlignment(3);
    }
    
    public void shift(Node node, double offset) {
        node.setHeight(node.getHeight() + offset);
        for (Node child : node.getChildren()) {
            shift(child, offset);
        }
    }
    
    public boolean treesEquivalentShifted(Tree treeA, Tree treeB, double tolerance) {
        shift(treeA.getRoot(), -treeA.getRoot().getHeight());
        shift(treeB.getRoot(), -treeB.getRoot().getHeight());
        return treesEquivalent(treeA, treeB, tolerance);
    }
    
    /**
     * Tests whether treeA and treeB are equivalent.  That is, whether they
     * have the same node heights (to within tolerance) and clade sets.
     * 
     * @param treeA
     * @param treeB
     * @param tolerance
     * @return 
     */
    public boolean treesEquivalent(Tree treeA, Tree treeB, double tolerance) {
        return treesEquivalent(treeA.getRoot(), treeB.getRoot(), tolerance);
    }

    /**
     * Tests whether trees below rootA and rootB are equivalent. That is,
     * whether they have the same node heights (to within tolerance) and clade
     * sets.
     *
     * @param rootA
     * @param rootB
     * @param tolerance
     * @return
     */
    public boolean treesEquivalent(Node rootA, Node rootB, double tolerance) {

        // Early exit if trees are different sizes
        if (rootA.getLeafNodeCount() != rootB.getLeafNodeCount())
            return false;
        
        Map<Clade, Double> cladeHeightsA = getCladeHeights(rootA);
        Map<Clade, Double> cladeHeightsB = getCladeHeights(rootB);
        
        if (!cladeHeightsA.keySet().containsAll(cladeHeightsB.keySet()))
            return false;
        
        for (Clade clade : cladeHeightsA.keySet()) {
            if (Math.abs(cladeHeightsA.get(clade)-cladeHeightsB.get(clade))>tolerance)
                return false;
        }
        
        return true;
    }

    /**
     * Method to test whether the topologies of two trees are equivalent.
     * 
     * @param treeA
     * @param treeB
     * @return true iff topologies are equivalent.
     */
    public boolean topologiesEquivalent(Tree treeA, Tree treeB) {
        return topologiesEquivalent(treeA.getRoot(), treeB.getRoot());
    }

    /**
     * Method to test whether the topologies of two trees with the
     * given root nodes are equivalent.
     * 
     * @param rootA
     * @param rootB
     * @return true iff topologies are equivalent.
     */
    public boolean topologiesEquivalent(Node rootA, Node rootB) {
        
        // Early exit if trees are different sizes
        if (rootA.getLeafNodeCount() != rootB.getLeafNodeCount())
            return false;
        
        Map<Clade, Double> cladeHeightsA = getCladeHeights(rootA);
        Map<Clade, Double> cladeHeightsB = getCladeHeights(rootB);
        
        return cladeHeightsA.keySet().containsAll(cladeHeightsB.keySet());
    }
    
    /**
     * Retrieve clades and heights of clade MRCAs from tree.
     * 
     * @param root
     * @return Map from clades to corresponding MRCA heights.
     */
    private Map<Clade, Double> getCladeHeights(Node root) {
        Map<Clade, Double> cladeHeights = new HashMap<>();
        
        cladeHeights.put(new Clade(root), root.getHeight());
        for (Node node : root.getAllChildNodesAndSelf())
            if (!node.isLeaf())
                cladeHeights.put(new Clade(node), node.getHeight());
        
        return cladeHeights;
    }
    
    /**
     * Convenience clade class.
     */
    private class Clade extends BitSet {

        /**
         * Construct clade from leaves below node.
         * 
         * @param node 
         */
        public Clade(Node node) {
            for (Node leaf : node.getAllLeafNodes())
                set(leaf.getNr());
        }
    }

    /**
     * Remove screen log (if it exists) from given runnable, if that
     * runnable is of dynamic type MCMC.  This prevents needlessly
     * verbose test log output.
     *
     * @param runnable from which to remove log.
     */
    public void disableScreenLog(beast.core.Runnable runnable) {
        if (runnable instanceof MCMC) {
            MCMC mcmc = (MCMC)runnable;

            Logger screenLog = null;
            for (Logger logger : mcmc.loggersInput.get()) {
                if (logger.fileNameInput.get() == null) {
                    screenLog = logger;
                    break;
                }
            }

            if (screenLog != null)
                mcmc.loggersInput.get().remove(screenLog);
        }
    }
	
}
