package contactrees.test;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import contactrees.Conversion;
import contactrees.ConversionGraph;
import junit.framework.TestCase;

public class TreeTest extends TestCase {
	final static double EPSILON = 1e-10;

	@Test
	public void testTreeScaling() {
        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";

        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
        ConversionGraph acg = new ConversionGraph();
        acg.initAndValidate();
        acg.assignFrom(treeParser);
        
        Node [] node = acg.getNodesAsArray();
        assertEquals(0.0, node[0].getHeight(), EPSILON);
        assertEquals(0.0, node[1].getHeight(), EPSILON);
        // leaf node, not scaled
        assertEquals(0.5, node[2].getHeight(), EPSILON);
        assertEquals(0.5, node[3].getHeight(), EPSILON);
        // internal nodes, all scaled
        assertEquals(1.0, node[4].getHeight(), EPSILON);
        assertEquals(1.5, node[5].getHeight(), EPSILON);
        assertEquals(2.0, node[6].getHeight(), EPSILON);
        
                assertEquals(3.0 * Math.log(2.0), acg.scale(2.0), EPSILON);
        
        // leaf node
        node = acg.getNodesAsArray();
        assertEquals(0.0, node[0].getHeight(), EPSILON);
        assertEquals(0.0, node[1].getHeight(), EPSILON);
        // leaf node, not scaled
        assertEquals(0.5, node[2].getHeight(), EPSILON);
        assertEquals(0.5, node[3].getHeight(), EPSILON);
                // internal margins scaled
        assertEquals(2.0, node[4].getHeight(), EPSILON);
                assertEquals(2.5, node[5].getHeight(), EPSILON);
                assertEquals(3.5, node[6].getHeight(), EPSILON);
		
	}

        @Test
        public void testTreeScalingAlsoScalesConversions() {
                String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";

                TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
                ConversionGraph acg = new ConversionGraph();
                acg.initAndValidate();
                acg.assignFrom(treeParser);

                Node[] node = acg.getNodesAsArray();
                Conversion conv1 = new Conversion(node[4], node[5], 1.75, acg, 1);
                acg.addConversion(conv1);

                double initialScalableValue = acg.getScalableValue();
                assertEquals(2.75, initialScalableValue, EPSILON);

                assertEquals(4.0 * Math.log(2.0), acg.scale(2.0), EPSILON);
                assertEquals(5.5, acg.getScalableValue(), EPSILON);
                assertEquals(3.0, conv1.getHeight(), EPSILON);
                assertFalse(acg.isInvalid());
        }

        @Test
        public void testSetScalableValueMatchesScaleForConversions() {
                String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";

                TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
                ConversionGraph acg = new ConversionGraph();
                acg.initAndValidate();
                acg.assignFrom(treeParser);

                Node[] node = acg.getNodesAsArray();
                acg.addConversion(new Conversion(node[4], node[5], 1.75, acg, 1));

                ConversionGraph scaled = acg.copy();
                ConversionGraph set = acg.copy();
                double targetValue = acg.getScalableValue() * 1.7;

                scaled.scale(1.7);
                set.setScalableValue(targetValue);

                assertEquals(targetValue, scaled.getScalableValue(), EPSILON);
                assertEquals(targetValue, set.getScalableValue(), EPSILON);

                Node[] scaledNodes = scaled.getNodesAsArray();
                Node[] setNodes = set.getNodesAsArray();
                for (int i = 0; i < scaledNodes.length; i++) {
                        assertEquals(scaledNodes[i].getHeight(), setNodes[i].getHeight(), EPSILON);
                }

                Conversion scaledConv = scaled.getConversions().iterator().next();
                Conversion setConv = set.getConversions().iterator().next();
                assertEquals(scaledConv.getHeight(), setConv.getHeight(), EPSILON);
                assertEquals(scaledConv.getNode1().getNr(), setConv.getNode1().getNr());
                assertEquals(scaledConv.getNode2().getNr(), setConv.getNode2().getNr());
        }
}
