/**
 *
 */
package contactrees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contactrees.util.parsers.ExtendedNewickBaseVisitor;
import contactrees.util.parsers.ExtendedNewickLexer;
import contactrees.util.parsers.ExtendedNewickParser;

/**
 * @author Nico Neureiter
 */
@Description("Conversion graph based around the clonal frame.")
public class ConversionGraph extends Tree {

    public Input<String> newickInput = new Input<>(
            "newick",
            "Initialise ARG from extended Newick representation.");
    public Input<Boolean> dropNewickConvsInput = new Input<>(
            "dropNewickConvs",
            "Drop the conversions specified in the extended Newick tree.",
            false);


    /**
     * List of conversion edges on graph (and a copy for restore).
     */
    protected ConversionList convs, storedConvs;

    /**
     * Clonal frame event list.
     */
    protected CFEventList cfEventList;

    @Override
    public void initAndValidate() {
//        System.setProperty("java.only", "true");

        super.initAndValidate();

        // Initialise conversions lists
        convs = new ConversionList(this);
        storedConvs = new ConversionList(this);


        if (newickInput.get() != null) {
            fromExtendedNewick(newickInput.get());
        }

        cfEventList = new CFEventList(this);
    }

    /**
     * Add conversion event to graph.
     *
     * @param conv conversion events to add
     */
    public void addConversion(Conversion conv) {
        startEditing(null);

        conv.setConversionGraph(this);
        convs.add(conv);
    }

    /**
     * Create a new conversion, add it to the graph and return it.
     *
     * @return The newly created conversion
     */
    public Conversion addNewConversion() {
        startEditing(null);

        Conversion conv = convs.addNewConversion();
        conv.setConversionGraph(this);
        return conv;
    }

    /**
     * Create a new conversion, add it to the graph and return it.
     *
     * @return The newly created conversion
     */
    public Conversion addDuplicateConversion(Conversion original) {
        startEditing(null);

        Conversion copy = convs.duplicateConversion(original);
        copy.setConversionGraph(this);
        return copy;
    }

    /**
     * Remove a conversion edge from the graph.
     *
     * @param conv conversion to remove.
     */
    public void removeConversion(Conversion conv) {
        startEditing(null);
        convs.remove(conv);
    }

    /**
     * Retrieve list of conversions associated with given locus.
     *
     * @param locus locus with which conversions are associated
     * @return List of conversions.
     */
    public ConversionList getConversions() {
        return convs;
    }

    /**
     * Remove all conversion edges from the graph.
     */
    public void removeAllConversions() {
        startEditing(null);
        convs.clear();
    }

    /**
     * Obtain total number of conversion events.
     *
     * @return Number of conversions.
     */
    public int getConvCount() {
        return convs.size();
    }

    /**
     * Obtain ordered list of events that make up the clonal frame.  Used
     * for ACG probability density calculations and for various state proposal
     * operators.
     *
     * @return List of events.
     */
    public List<CFEventList.Event> getCFEvents() {
        return cfEventList.getCFEvents();
    }

    /**
     * Obtain the CFEventList object, containting the ordered list of events
     * that make up the clonal frame and some convenience methods.
     *
     * @return cfEventList.
     */
    public CFEventList getCFEventList() {
    	cfEventList.updateEvents();
        return cfEventList;
    }

    /**
     * @return Total length of all edges in clonal frame.
     */
    public double getClonalFrameLength() {
        double length = 0.0;
        for (Node node : m_nodes) {
            if (node.isRoot())
                continue;
            length += node.getLength();
        }

        return length;
    }

    /**
     * @return Total length of all pairs of edges in clonal frame.
     */
    public double getClonalFramePairedLength() {
    	double length = 0.0;

        List<CFEventList.Event> events = getCFEvents();
        for (int i = 0; i<events.size()-1; i++) {
        	CFEventList.Event eStart= events.get(i);
        	CFEventList.Event eEnd = events.get(i+1);
        	double dt = eEnd.t - eStart.t;
        	int k = eStart.lineages;
        	length += dt * k * (k-1);
        }

        return length;
    }

    /**
     * Obtain the set of lineages active at the specified height.
     * @param height
     * @return Set of active lineages
     */
    public HashSet<Node> getLineagesAtHeight(double height) {
    	HashSet<Node> lineages = new HashSet<>();

        List<CFEventList.Event> events = getCFEvents();

        for (int i=0; i<events.size(); i++) {
        	CFEventList.Event event = events.get(i);
        	if (i < events.size()-1)
        		assert event.t <= events.get(i+1).t;
        	if (i < events.size()-1)
        		assert event.getNode().getHeight() <= events.get(i+1).getNode().getHeight();

        	Node node = event.getNode();
        	assert event.t == node.getHeight();

        	if (event.getHeight() > height) {
        		// We found the interval of "height"
        		break;
        	}

        	// Add new node (for COALESCENCE and SAMPLE events)
        	lineages.add(node);
        	// Remove children on COALESCENCE events
        	if (event.getType() == CFEventList.EventType.COALESCENCE) {
        		for (Node child : node.getChildren()) {
            		lineages.remove(child);
        		}
        	} else {
        		assert node.isLeaf();
        	}
        }

        return lineages;
    }

    public int countLineagesAtHeight(double height) {
    	getCFEvents();
    	return cfEventList.getEventAtHeight(height).lineages;
    }

    /**
     * Check validity of conversions.  Useful for probability densities
     * over the ACG to decide whether to return 0 based on an unphysical
     * state.
     *
     * @return true if all conversions are valid w.r.t. clonal frame.
     */
    public boolean isInvalid() {
    	// A node is the root iff it is at the last index of m_nodes array
        for (int i=0; i<m_nodes.length; i++) {
        	Node node = m_nodes[i];
    		Node parent = node.getParent();
        	if (node.isRoot())
        		assert i == m_nodes.length - 1;
        	if (!node.isRoot())
        		assert i != m_nodes.length - 1;
        	if ((node.isRoot()) != (i == m_nodes.length-1))
        		return true;
        	if (!node.isRoot())
        		assert m_nodes[parent.getNr()] == parent;
        	if (!node.isRoot())
        		assert node.getLength() > 0;
        }

        for (Conversion conv : convs) {
            if (!conv.isValid()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Create a copy of this conversion graph.
     *
     * @return The new copy of this conversion graph
     */
    @Override
    public ConversionGraph copy() {
        ConversionGraph acg = new ConversionGraph();

        // From Tree.copy():
        acg.setID(getID());
        acg.index = index;
        acg.root = root.copy();
        acg.nodeCount = nodeCount;
        acg.internalNodeCount = internalNodeCount;
        acg.leafNodeCount = leafNodeCount;

        acg.initArrays();
        acg.m_taxonset.setValue(m_taxonset.get(), acg);

        acg.convs = new ConversionList(this);
        acg.storedConvs = new ConversionList(this);

        for (Conversion conv : convs) {
            Conversion convCopy = conv.getCopy();
            convCopy.setConversionGraph(acg);
            convCopy.setNode1(acg.m_nodes[conv.getNode1().getNr()]);
            convCopy.setNode2(acg.m_nodes[conv.getNode2().getNr()]);
            acg.convs.convs.put(convCopy.getID(), convCopy);
        }
        for (Conversion conv : storedConvs) {
            Conversion convCopy = conv.getCopy();
            convCopy.setConversionGraph(acg);
            convCopy.setNode1(acg.m_nodes[conv.getNode1().getNr()]);
            convCopy.setNode2(acg.m_nodes[conv.getNode2().getNr()]);
            acg.storedConvs.convs.put(convCopy.getID(), convCopy);
        }

        return acg;
    }

    private void generalAssignFrom(StateNode other, boolean fragileAssignment) {
        if (fragileAssignment)
            super.assignFromFragile(other);
        else
            super.assignFrom(other);

        if (other instanceof ConversionGraph) {
            ConversionGraph acg = (ConversionGraph)other;

            convs.clear();
            storedConvs.clear();

            for (Conversion conv : acg.convs) {
                Conversion convCopy = conv.getCopy();
                convCopy.setConversionGraph(this);
                convCopy.setNode1(m_nodes[conv.getNode1().getNr()]);
                convCopy.setNode2(m_nodes[conv.getNode2().getNr()]);
                convs.convs.put(convCopy.getID(), convCopy);
            }

            if (cfEventList == null)
                cfEventList = new CFEventList(this);
            else
                cfEventList.makeDirty();
        }

//        nodeCount = m_nodes.length;
//        initArrays();
        assert !isInvalid();
    }

    /**
     * Use another StateNode to configure this ACG.  If the other StateNode
     * is merely a tree, only the clonal frame is configured.
     *
     * @param other StateNode used to configure ACG
     */
    @Override
    public void assignFrom(StateNode other) {
        generalAssignFrom(other, false);
    }

    @Override
    public void assignFromFragile(StateNode other) {
        generalAssignFrom(other, true);
    }

    /*
    * StateNode implementation
    */

    @Override
    public void store() {
        super.store();

        // Copy the conversion list
        storedConvs = convs.copy();

        // Change copied node references from m_nodes to m_storedNodes
        for (Conversion conv : storedConvs) {
        	Conversion original = convs.get(conv.id);
            conv.setNode1(m_storedNodes[conv.getNode1().getNr()]);
            conv.setNode2(m_storedNodes[conv.getNode2().getNr()]);
            conv.newickMetaDataBottom = original.newickMetaDataBottom;
            conv.newickMetaDataMiddle = original.newickMetaDataMiddle;
            conv.newickMetaDataTop = original.newickMetaDataTop;
        }

    }

    @Override
    public void restore() {
    	super.restore();

        // Swap conversions with storedConversions
        ConversionList tmp = storedConvs;
        storedConvs = convs;
        convs = tmp;

        cfEventList.makeDirty();

        assert !isInvalid();
    }

    @Override
    public void startEditing(Operator operator) {
        if (state != null)
            super.startEditing(operator);

        if (cfEventList != null)
            cfEventList.makeDirty();
    }

    /**
     * @return true iff clonal frame is dirty
     */
    public boolean clonalFrameIsDirty() {
        for (Node node : getNodesAsArray())
            if (node.isDirty() > Tree.IS_CLEAN)
                return true;

        return false;
    }


    /*
    * Parsing functions
    */

    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        fromExtendedNewick(node.getTextContent().replaceAll("&amp", "&"));
    }

    /**
     * Read in an ACG from a string in extended newick format.  Assumes
     * that the network is stored with exactly the same metadata as written
     * by the getExtendedNewick() method.
     *
     * @param string extended newick representation of ACG
     */
    public void fromExtendedNewick(String string) {
        fromExtendedNewick(string, false, taxaTranslationOffset);
    }
    /**
     * Read in an ACG from a string in extended newick format.  Assumes
     * that the network is stored with exactly the same metadata as written
     * by the getExtendedNewick() method.
     *
     * @param string extended newick representation of ACG
     * @param numbered true indicates that the ACG is numbered.
     */
    public void fromExtendedNewick(String string, boolean numbered, int nodeNumberoffset) {
        // Spin up ANTLR
        CharStream input = CharStreams.fromString(string);
        ExtendedNewickLexer lexer = new ExtendedNewickLexer(input);
        CommonTokenStream tokens = new CommonTokenStream(lexer);
        ExtendedNewickParser parser = new ExtendedNewickParser(tokens);
        ParseTree parseTree = parser.tree();

        Map<String, Conversion> convIDMap = new HashMap<>();
        Node root = new ExtendedNewickBaseVisitor<Node>() {

            /**
             * Convert branch lengths to node heights for all nodes in clade.
             *
             * @param node clade parent
             * @return minimum height assigned in clade.
             */
            private double branchLengthsToHeights(Node node) {
                if (node.isRoot())
                    node.setHeight(0.0);
                else
                    node.setHeight(node.getParent().getHeight() - node.getHeight());

                double minHeight = node.getHeight();

                for (Node child : node.getChildren()) {
                    minHeight = Math.min(minHeight, branchLengthsToHeights(child));
                }

                return minHeight;
            }

            /**
             * Remove height offset from all nodes in clade
             * @param node parent of clade
             * @param offset offset to remove
             */
            private void removeOffset(Node node, double offset) {
                node.setHeight(node.getHeight() - offset);

                for (Node child : node.getChildren())
                    removeOffset(child, offset);
            }

            /**
             * Get the next clonal frame node (at or below the given node).
             * @param node
             * @return Clonal frame node.
             */
            private Node getTrueNode(Node node) {
                if (node.isLeaf()) {
                    assert !convIDMap.containsKey(node.getID());
                    return node;
                }

                if (convIDMap.containsKey(node.getID()))
                    return getTrueNode(node.getChild(0));

                int hybridIdx = -1;
                int nonHybridIdx = -1;
                for (int i=0; i<node.getChildCount(); i++) {
                    if (node.getChild(i).isLeaf() && convIDMap.containsKey(node.getChild(i).getID()))
                        hybridIdx = i;
                    else
                        nonHybridIdx = i;
                }

                if (hybridIdx>0)
                    return getTrueNode(node.getChild(nonHybridIdx));

                return node;
            }

            /**
             * Traverse the newly constructed tree looking for
             * hybrid nodes and using these to set the heights of
             * Conversion objects.
             *
             * @param node parent of clade
             */
            private void findConversionAttachments(Node node) {
                if (convIDMap.containsKey(node.getID())) {
                    Conversion conv = convIDMap.get(node.getID());
                    if (node.isLeaf()) {
                        conv.setHeight(node.getParent().getHeight());
                        conv.setNode2(getTrueNode(node.getParent()));
                    } else
                        conv.setNode1(getTrueNode(node));
                }

                for (Node child : node.getChildren())
                    findConversionAttachments(child);
            }

            /**
             * Remove all conversion-associated nodes, leaving only
             * the clonal frame.
             *
             * @param node parent of clade
             * @return new parent of same clade
             */
            private Node stripHybridNodes(Node node) {
                Node trueNode = getTrueNode(node);
                List<Node> trueChildren = new ArrayList<>();

                for (Node child : trueNode.getChildren()) {
                    trueChildren.add(stripHybridNodes(child));
                }

                trueNode.removeAllChildren(false);
                for (Node trueChild : trueChildren)
                    trueNode.addChild(trueChild);

                return trueNode;
            }

            private int numberInternalNodes(Node node, int nextNr) {
                if (node.isLeaf())
                    return nextNr;

                for (Node child : node.getChildren())
                    nextNr = numberInternalNodes(child, nextNr);

                node.setNr(nextNr);

                return nextNr + 1;
            }


            @Override
            public Node visitTree(ExtendedNewickParser.TreeContext ctx) {
                Node root =  visitNode(ctx.node());

                double minHeight = branchLengthsToHeights(root);
                removeOffset(root, minHeight);

                findConversionAttachments(root);

                root = stripHybridNodes(root);
                root.setParent(null);

                if (!numbered)
                    numberInternalNodes(root, root.getAllLeafNodes().size());

                return root;
            }

            @Override
            public Node visitNode(ExtendedNewickParser.NodeContext ctx) {
                Node node = new Node();

                if (ctx.post().hybrid() != null) {
                    String convID = ctx.post().hybrid().getText();
                    node.setID(convID);

                    Conversion conv;
                    if (convIDMap.containsKey(convID))
                        conv = convIDMap.get(convID);
                    else {
                        conv = new Conversion(parseConvID(convID));
                        convIDMap.put(convID, conv);
                    }

                    if (ctx.node().isEmpty()) {
                        for (ExtendedNewickParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
                            switch (attribCtx.attribKey.getText()) {
                                default:
                                	// Once we write more information in the Ext.Newick, we could parse it here.
                                    break;
                            }
                        }
                    }
                }

                for (ExtendedNewickParser.NodeContext childCtx : ctx.node())
                    node.addChild(visitNode(childCtx));

                if (ctx.post().label() != null) {
                    node.setID(ctx.post().label().getText());
                    String lbl = ctx.post().label().getText();
                    TaxonSet taxonSet = getTaxonset();
                    if ((getTaxonset() != null) && (taxonSet.asStringList().contains(lbl))) {
                        node.setNr(taxonSet.getTaxonIndex(lbl));
                    } else {
                        node.setNr(Integer.parseInt(ctx.post().label().getText())
                                   - nodeNumberoffset);
                    }
                }
                node.setHeight(Double.parseDouble(ctx.post().length.getText()));

                return node;
            }
        }.visit(parseTree);

        initAfterParsingFromNewick(root, m_nodes);

        cfEventList = new CFEventList(this);

        if (!dropNewickConvsInput.get()) {
            for (Conversion conv : convIDMap.values())
                addConversion(conv);
        }
    }

    public int _getNodeCount() {
    	return nodeCount;
    }

    public void initAfterParsingFromNewick(Node root, Node[] mnodes) {
        m_nodes = root.getAllChildNodesAndSelf().toArray(m_nodes);
        nodeCount = m_nodes.length;
        leafNodeCount = root.getAllLeafNodes().size();
        setRoot(root);
        initArrays();
    }

    public void initTreeArrays() {
    	initArrays();
    }

    protected int parseConvID(String sConvID) {
    	return Integer.parseInt(sConvID.substring(1));
    }

    /**
     * DEBUG CHECKS
     */

    @Override
    public int getChecksum() {
        Integer[] hashValues = new Integer[3];

        // If the ACG is the same, the following properties need to match:

        // Number of conversions should be the same
        hashValues[0] = getConvCount();

        // The root height should be the same (could be moved to Tree.getChecksum())
        hashValues[1] = Double.hashCode(getRoot().getHeight());

        // The mean conversion height should be the same
//        hashValues[2] = Double.hashCode(getMeanConvHeight());
        // (hashing the exact Double value can lead to false alarms due to rounding errors)
        hashValues[2] = (int) (100000 * getMeanConvHeight());

        return Arrays.deepHashCode(hashValues);
    }

    public double getMeanConvHeight() {
        if (getConvCount() == 0)
            return 0.0;

        double sumConvHeight = 0;
        for (Conversion conv : getConversions()) {
            sumConvHeight += conv.getHeight();
        }
        return sumConvHeight / getConvCount();
    }

}
