/**
 *
 */
package contactrees.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashSet;

import javax.naming.directory.InvalidAttributesException;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import contactrees.Block;
import contactrees.Conversion;
import contactrees.ConversionGraph;



/**
 * Test basic ConversionGraph functionality.
 *
 * @author Nico Neureiter
 */
public class  ConversionGraphTest extends ContactreesTest {

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

	/**
	 * Restoring an ACG from a state file must keep every node number, since node-indexed state
	 * (e.g. relaxed clock rate categories) relies on them. Internal numbers here are deliberately
	 * not in post-order, as after a few MCMC moves.
	 */
	@Test
	public void testStateFileRoundTripKeepsNodeNumbers() throws Exception {
		// Post-order numbering would be 6 for (1,2), 7 for its parent and 8 for (4,5)
		ConversionGraph original = new ConversionGraph();
		original.initAndValidate();
		original.fromExtendedNewick("(((1:0.4,2:0.4)8:0.6,3:1.0)6:0.7,(4:0.8,5:0.8)7:0.9)9:0.0;", true, 1);
		assertEquals(original.getNode(7).getHeight(), 0.4, 1e-9);
		assertEquals(original.getNode(5).getHeight(), 1.0, 1e-9);
		original.addConversion(new Conversion(original.getNode(2), original.getNode(3), 0.5, original, 1));
		original.addConversion(new Conversion(original.getNode(7), original.getNode(3), 0.6, original, 2));

		ConversionGraph restored = restoreFromXML(original);
		assertEquals(original.getNodeCount(), restored.getNodeCount());
		for (int nr = 0; nr < original.getNodeCount(); nr++) {
			Node n0 = original.getNode(nr), n1 = restored.getNode(nr);
			assertEquals(nr, n1.getNr());
			assertEquals(n0.getHeight(), n1.getHeight(), 1e-9);
			assertEquals(n0.isRoot() ? -1 : n0.getParent().getNr(), n1.isRoot() ? -1 : n1.getParent().getNr());
		}
		assertEquals(original.getConvCount(), restored.getConvCount());
		for (Conversion c0 : original.getConversions()) {
			Conversion c1 = restored.getConversions().get(c0.getID());
			assertEquals(c0.getNode1().getNr(), c1.getNode1().getNr());
			assertEquals(c0.getNode2().getNr(), c1.getNode2().getNr());
			assertEquals(c0.getHeight(), c1.getHeight(), 1e-9);
		}
	}

	/** A root that is not the last node is still moved to the last slot, keeping numbers distinct. */
	@Test
	public void testNumberedParsingMovesRootLast() {
		ConversionGraph acg = new ConversionGraph();
		acg.initAndValidate();
		acg.fromExtendedNewick("((1:1.0,2:1.0)7:1.0,(3:0.5,4:0.5)6:1.5)5:0.0;", true, 1);
		assertEquals(6, acg.getRoot().getNr());
		HashSet<Integer> nrs = new HashSet<>();
		for (int nr = 0; nr < acg.getNodeCount(); nr++) {
			assertEquals(nr, acg.getNode(nr).getNr());
			nrs.add(acg.getNode(nr).getNr());
		}
		assertEquals(7, nrs.size());
	}

	/** Without valid internal labels the parser still numbers internal nodes itself. */
	@Test
	public void testNumberedParsingFallsBackWithoutLabels() {
		ConversionGraph acg = new ConversionGraph();
		acg.initAndValidate();
		acg.fromExtendedNewick("((1:1.0,2:1.0):1.0,(3:0.5,4:0.5):1.5):0.0;", true, 1);
		assertEquals(6, acg.getRoot().getNr());
		HashSet<Integer> nrs = new HashSet<>();
		for (int nr = 0; nr < acg.getNodeCount(); nr++)
			nrs.add(acg.getNode(nr).getNr());
		assertEquals(7, nrs.size());
	}

	static ConversionGraph restoreFromXML(ConversionGraph acg) throws Exception {
		String xml = acg.toXML();
		org.w3c.dom.Document doc = javax.xml.parsers.DocumentBuilderFactory.newInstance().newDocumentBuilder()
				.parse(new java.io.ByteArrayInputStream(xml.getBytes()));
		ConversionGraph restored = new ConversionGraph();
		restored.initAndValidate();
		restored.fromXML(doc.getDocumentElement());
		return restored;
	}

	@Test
	public void testRestore() {
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
			Conversion conv = new Conversion(node4, node3, 2.0, acg, 0);
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
		assertTrue(lineages.contains(node1));
		assertTrue(lineages.contains(node2));
		assertTrue(lineages.contains(node3));
		assertFalse(lineages.contains(node4));

		assertEquals(2, acg.countLineagesAtHeight(1.5));
		lineages = acg.getLineagesAtHeight(1.5);
		assertEquals(2, lineages.size());
		assertTrue(lineages.contains(node4));
		assertTrue(lineages.contains(node3));

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
//		conv = new Conversion(node3, node4, 1.5, acg, 0);
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
//        	Node node2 = conv.getNode1();
//        	Node node3 = conv.getNode2();
//
//        	boolean isRootChild = node2.getParent().isRoot();
//
//        	if (rootOnly) {
//        		if (node3.getParent().isRoot() != isRootChild) {
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
//            if (conv.getHeight() < node2.getHeight()
//                    || conv.getHeight() < node3.getHeight())
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
