/**
 * 
 */
package contactrees.test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;

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
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class ContactreesTest {

	protected static final double EPS = 1E-8;
	
//	List<ConversionGraph> acgs;
//	List<TestBlockSet> blockSets;
//	List<TestACGWithMetaDataLogger> acgLoggers;

	protected ConversionGraph acg;
	protected BlockSet blockSet;
	protected ACGWithMetaDataLogger acgLogger;
	protected Node node1, node2, node3, node4, root;
	
	@Before
	public void setup() {
		buildACG_1();
		assertAllValid();
	}
	
	public void assertAllValid() {
		assertFalse(acg.isInvalid());
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
		acg.addConversion(new Conversion(node2, node3, 0.5, acg, 0));
		acg.addConversion(new Conversion(node3, node4, 1.5, acg, 0)); 
	}
	
	public void setupFromNewick(String newick) {
		// Build conversion graph from tree 
		acg = new ConversionGraph();
		acg.initAndValidate();
		acg.fromExtendedNewick(newick);

		// Initialize BlockSet and ACG logger
//		blockSet = BlockSet.getBlockSet(acg);
		List<Block> blocks = new ArrayList<Block>();
		for (int i=0; i<3; i++) {
		    Block b = new Block();
		    b.initAndValidate();
		    blocks.add(b);
		}
		
		blockSet = new BlockSet();
		blockSet.initByName("network", acg, "block", blocks);
		
		acgLogger = ACGWithMetaDataLogger.getACGWMDLogger(acg, blockSet);		
	}

	public void printConversions() {
		System.out.println(String.format(
				"ACG contains %d conversions:", 
				acg.getConvCount()));
		for (Conversion c : acg.getConversions())
			System.out.println("\t" + c);	
	}
	
}
