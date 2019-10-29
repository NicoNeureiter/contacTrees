/**
 * 
 */
package contactrees.test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import javax.naming.directory.InvalidAttributesException;

import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import contactrees.ACGWithMetaDataLogger;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.CFEventList;
import contactrees.Conversion;
import contactrees.ConversionGraph;



/**
 * Test basic ConversionGraph functionality. 
 * 
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
public class  ConversionGraphTest extends ContactreesTest {

	static final double EPS = 1E-8;
	
	@Test
	public void testExtendedNewickParsing( ) {
		String NEWICK = "((1:1.0,(2:0.5, #1337[&conv=1415,blocks={1,2}]:0.00000001):0.5)4:1.5,(3:0.5)#1337:2.0)5:0.5;";
		setupFromNewick(NEWICK);
		 
		assertEquals(acg.getConvCount(), 1);
		
		Conversion conv = acg.getConversions().iterator().next();
		assertEquals(conv.getHeight(), 0.5, 1E-6);
		assertTrue(conv.getNode1().isLeaf());
		assertTrue(conv.getNode2().isLeaf());
		assertTrue(conv.getNode1().getParent().isRoot());
		assertEquals(conv.getNode2().getParent().getHeight(), 1.0, 1E-6);
	}
	
	@Test
	public void testRestore() {
		Node[] nodesBefore;
		String newickBefore, newickAfter;
		
		for (ModificationType mod : ModificationType.values()) {
			// Prep
			acg.store();
			newickBefore = acgLogger.getExtendedNewick();
			
			// Add/remove conversion
			try {
				acgModification(mod);
			} catch (InvalidAttributesException e) {
				assert false : "This should not happen (acg should not be empty here).";
			}					
			
			// "something changed"
			newickAfter = acgLogger.getExtendedNewick();
			assertFalse(newickBefore.equals(newickAfter));
			
			// Restore
			acg.restore();
			for (Block b : blockSet) {
				b.restore();
			}

			// Back to "nothing changed"
			newickAfter = acgLogger.getExtendedNewick();
			assertTrue(newickBefore.equals(newickAfter));
		}
	}
	
	public enum ModificationType {ADD_CONV, REMOVE_CONV}
	
	public void acgModification (ModificationType modification) throws InvalidAttributesException {
		switch (modification) {
		case ADD_CONV:
			Conversion conv = new Conversion(node3, node2, 2.0, acg, 0);
			acg.addConversion(conv);
			break;
		case REMOVE_CONV:		
			conv = acg.getConversions().getRandomConversion();
			acg.removeConversion(conv);
			blockSet.removeConversion(conv);
			break;
		}
	}
	
	@Test
	public void TestCFLength() {
		double height = acg.getRoot().getHeight();;
		assertEquals(2.5, height, EPS);
		assertEquals(6.0, acg.getClonalFrameLength(), EPS);
		assertEquals(9.0, acg.getClonalFramePairedLength(), EPS);
		
		setupFromNewick("(1:2.0,2:2.0):0.0;");
		height = acg.getRoot().getHeight();
		assertEquals(2.0, height, EPS);
		assertEquals(4.0, acg.getClonalFrameLength(), EPS);
		assertEquals(4.0, acg.getClonalFramePairedLength(), EPS);
		
		setupFromNewick("((1:2.0,2:2.0):1.0,3:3.0):0.0;");
		height = acg.getRoot().getHeight();
		assertEquals(3.0, height, EPS);
		assertEquals(8.0, acg.getClonalFrameLength(), EPS);
		assertEquals(14.0, acg.getClonalFramePairedLength(), EPS);
	}
	
	
	@Test
	public void TestLineagesAtHeight() {
		HashSet<Node> lineages;
		
		assertEquals(3, acg.countLineagesAtHeight(0.5));
		lineages = acg.getLineagesAtHeight(0.5);
		assertEquals(3, lineages.size());
		assertTrue(lineages.contains(node0));
		assertTrue(lineages.contains(node1));
		assertTrue(lineages.contains(node2));
		assertFalse(lineages.contains(node3));
		
		assertEquals(2, acg.countLineagesAtHeight(1.5));
		lineages = acg.getLineagesAtHeight(1.5);
		assertEquals(2, lineages.size());
		assertTrue(lineages.contains(node3));
		assertTrue(lineages.contains(node2));

		assertEquals(1, acg.countLineagesAtHeight(2.75));
		lineages = acg.getLineagesAtHeight(2.75);
		assertEquals(1, lineages.size());
		assertTrue(lineages.contains(root));
		
//		System.out.println(acg.countLineagesAtHeight(1.0));
//		System.out.println(acg.getLineagesAtHeight(1.0).size());
//		for (CFEventList.Event e : acg.getCFEvents()) {
//			System.out.println(e + ": " + e.getHeight() + ", " + e.getNode() + ", " + e.getLineageCount());
//		}
	}
	
	
//	@Test
//	public void test() {
//		String newickBefore, newickAfter;
//		assertEquals(3, acg.getLeafNodeCount());
//		
//		Conversion conv;
//		
//		scaleACG(2.0);
//		assertAllValid();
//
//		scaleACG(0.5);
//		assertAllValid();
//		
//		conv = new Conversion(node2, node3, 1.5, acg, 0);
//		acg.addConversion(conv);
//		assertAllValid();
//
//		scaleACG(2.0);
//		System.out.println(acgLogger.getExtendedNewick());
//		assertAllValid();
//		
//		printConversions();
//		
//		newickBefore = acgLogger.getExtendedNewick();
//		acg.store();
//		newickAfter = acgLogger.getExtendedNewick();
//		assertEquals(newickAfter, newickBefore);
//		
//		acg.removeConversion(conv);
//		scaleACG(2.0);
//		assertAllValid();
//		
//		acg.restore();
//		newickAfter = acgLogger.getExtendedNewick();
//		printConversions();
//		assertEquals(newickAfter, newickBefore);
//		assertAllValid();
//		
//	}	
//	
//	public double scaleACG(double f) {
////        double scaleParam = 1.2;
//        boolean rootOnly = false;
//		
//        // Keep track of number of positively scaled elements minus
//        // negatively scaled elements.
//        int count = 0;
//
//        // Choose scaling factor:
////        double f = scaleParam + Randomizer.nextDouble()*(1.0/scaleParam - scaleParam);
////        double f = Randomizer.uniform(scaleParam, 1.0/scaleParam);
//        
//        // Scale clonal frame:
//        if (rootOnly) {
//            acg.getRoot().setHeight(acg.getRoot().getHeight()*f);
//            count += 1;
//        } else {
//            for (Node node : acg.getInternalNodes()) {
//            	System.out.print("Rescale node ");
//            	System.out.println(node.getNr());
//            	System.out.print("Old height ");
//            	System.out.println(node.getHeight());
//                node.setHeight(node.getHeight()*f);
//            	System.out.print("New height ");
//            	System.out.println(node.getHeight());
//                count += 1;
//            }
//        }       
//        
//        // Scale conversion edges:
//        for (Conversion conv : acg.getConversions()) {
//        	Node node1 = conv.getNode1();
//        	Node node2 = conv.getNode2();
//        	
//        	boolean isRootChild = node1.getParent().isRoot(); 
//        	
//        	if (rootOnly) {
//        		if (node2.getParent().isRoot() != isRootChild) {
//        			// One node gets rescaled, the other does not -> Illegal move
//        			return Double.NEGATIVE_INFINITY;
//        		}
//        	}
//        	
//            if (!rootOnly || isRootChild) {
//                conv.setHeight(conv.getHeight() * f);
//                count += 1;
//            }
//
//            if (conv.getHeight() < node1.getHeight()
//                    || conv.getHeight() < node2.getHeight())
//            	// Conversion below its node -> Illegal move
//                return Double.NEGATIVE_INFINITY;
//        }
//        
//        // Check for illegal node heights:
//        if (rootOnly) {
//            for (Node node : acg.getRoot().getChildren()) {
//                if (node.getHeight()>node.getParent().getHeight())
//                    return Double.NEGATIVE_INFINITY;
//            }
//        } else {
//            for (Node node : acg.getExternalNodes()) {
//                if (node.getHeight()>node.getParent().getHeight())
//                    return Double.NEGATIVE_INFINITY;
//            }
//        }
//        
//        System.out.println(f);
//
//        assert !acg.isInvalid() : "ACGScaler produced invalid state.";
//        
//        // Return log of Hastings ratio:
//        return (count-2)*Math.log(f);
//    }
}
