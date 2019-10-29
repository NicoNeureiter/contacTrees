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
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;

/**
 * 
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class ContactreesTest {

//	List<ConversionGraph> acgs;
//	List<TestBlockSet> blockSets;
//	List<TestACGWithMetaDataLogger> acgLoggers;

	ConversionGraph acg;
	BlockSet blockSet;
	ACGWithMetaDataLogger acgLogger;
	Node node0, node1, node2, node3, root;
	
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
		node3 = root.getLeft();
		node2 = root.getRight();
		node1 = node3.getRight();
		node0 = node3.getLeft();
		
		// Add conversion
		acg.addConversion(new Conversion(node1, node2, 0.5, acg, 0));
		acg.addConversion(new Conversion(node2, node3, 1.5, acg, 0)); 
	}
	
	public void setupFromNewick(String newick) {
		// Build conversion graph from tree 
		acg = new ConversionGraph();
		acg.initAndValidate();
		acg.fromExtendedNewick(newick);

		// Initialize BlockSet and ACG logger
		blockSet = BlockSet.getBlockSet(acg);
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
