/**
 * 
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contactrees.CFEventList.Event;
import contactrees.util.parsers.ExtendedNewickBaseVisitor;
import contactrees.util.parsers.ExtendedNewickLexer;
import contactrees.util.parsers.ExtendedNewickParser;

/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("Conversion graph based around the clonal frame.")
public class ConversionGraph extends Tree {

    public Input<Integer> siteCountInput = new Input<>(
            "siteCount",
            "Number of sites contained in alignments associated with this graph.");
	
	public Input<Alignment> alignmentInput = new Input<>(
			"alignment",
			"Initialize graph using the alignment.",
			Input.Validate.XOR, siteCountInput);
	
    /**
     * List of conversion edges on graph (and a copy for restore).
     */
    protected List<Conversion> convs;
    protected List<Conversion> storedConvs;
      
    /**
     * Clonal frame event list and number of sites in alignment.
     */
    protected CFEventList cfEventList;
	protected int siteCount;

    @Override
    public void initAndValidate() {
    	// Initialise conversions lists
        convs = new ArrayList<Conversion>();
        storedConvs = new ArrayList<Conversion>();
        
        if (alignmentInput.get() != null) {
        	siteCount = alignmentInput.get().getSiteCount();
        } else {
        	siteCount = siteCountInput.get();
        }
        
        cfEventList = new CFEventList(this);
        
        super.initAndValidate();
    }

    /**
     * Add conversion event to graph, ensuring conversions list
     * remains sorted (by height).
     *
     * @param conv conversion events to add
     */
    public void addConversion(Conversion conv) {
        startEditing(null);
        
        conv.setConversionGraph(this);

        int i;
        for (i=0; i<convs.size(); i++)
            if (convs.get(i).height > conv.height)
                break;
        
        convs.add(i, conv);
    }
   
    /**
     * Remove recombination from graph.
     *
     * @param conv conversion to remove.
     */
    public void deleteConversion(Conversion conv) {
        startEditing(null);
        convs.remove(conv);
    }
    
    /**
     * Obtain total number of conversion events.
     *
     * @return Number of conversions.
     */
    public int getTotalConvCount() {
        return convs.size();
    }

    /**
     * Obtain index of conversion when conversions are listed in order
     * of height.
     *
     * @param conv conversion whose index is required
     * @return Conversion index
     */
    public int getConversionIndex(Conversion conv) {
        return convs.indexOf(conv);
    }
    
    /**
     * Obtain the list of conversions affecting each site.
     * 
     * @return A list of lists of conversions for every site.
     */
    public ArrayList<List<Conversion>> getConversionsPerSite() {
    	ArrayList<List<Conversion>> conversionsPerSite = new ArrayList<>();
    	
    	for (int s=0; s<siteCount; s++)
    		conversionsPerSite.add(new ArrayList<Conversion>());
    	
    	for (Conversion conv : convs) {
    		for (int s : conv.affectedSites)
    			conversionsPerSite.get(s).add(conv);
    	}
    	
    	return conversionsPerSite;
    }
    
    /**
     * Count the number of sites that are affected by a conversion.
     * 
     * @return Number of affected sites
     */
    public int countAffectedSites() {
    	boolean[] affected = new boolean[siteCount];
    	for (Conversion conv : convs) {
    		for (int s : conv.affectedSites) {
    			affected[s] = true;
    		}
    	}
    	
    	int affectedCount = 0;
    	for (boolean b : affected) if (b) affectedCount++;
    	return affectedCount;
    }
    

    /**
     * Count the number of conversion edges which do not affect any sites.
     * 
     * @return The number of edges without an affected site.
     */
    public int getUselessConvCount() {
        int count = 0;
        for (Conversion conv : convs) {
            if (conv.affectedSites.isEmpty())
                count += 1;
        }

        return count;
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
        	// TODO check whehter k and dt match! Should we use eEnd.lineages?
        }
        
        return length;
    }
    
    /**
     * Check validity of conversions.  Useful for probability densities
     * over the ACG to decide whether to return 0 based on an unphysical
     * state.
     * 
     * @return true if all conversions are valid w.r.t. clonal frame.
     */
    public boolean isInvalid() {
        for (Conversion conv : convs) {
            if (!conv.isValid()) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Produces an extended Newick representation of this ACG.  This
     * method is also used to serialize the state to a state file.
     *
     * @return an extended Newick representation of ACG.
     */
    @Override
    public String toString() {
        return getExtendedNewick();
    }

    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        fromExtendedNewick(node.getTextContent().replaceAll("&amp", "&"));
    }

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
        
        acg.convs = new ArrayList<>();
        acg.storedConvs = new ArrayList<>();

        for (Conversion conv : convs) {
            Conversion convCopy = conv.getCopy();
            convCopy.setConversionGraph(acg);
            convCopy.setNode1(acg.m_nodes[conv.getNode1().getNr()]);
            convCopy.setNode2(acg.m_nodes[conv.getNode2().getNr()]);
            acg.convs.add(convCopy);
        }
        for (Conversion conv : storedConvs) {
            Conversion convCopy = conv.getCopy();
            convCopy.setConversionGraph(acg);
            convCopy.setNode1(acg.m_nodes[conv.getNode1().getNr()]);
            convCopy.setNode2(acg.m_nodes[conv.getNode2().getNr()]);
            acg.storedConvs.add(convCopy);
        }
        
        acg.siteCount = siteCount;

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
                convs.add(convCopy);
            }
            
            // TODO what about copying storedConvs?

            if (cfEventList == null)
                cfEventList = new CFEventList(this);
            
            siteCount = acg.siteCount;
        }
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
    
    /**
     * Obtain extended Newick representation of ACG.  Includes Nexus metadata
     * on hybrid leaf nodes describing the alignment sites affected by the
     * conversion event.
     * 
     * @return Extended Newick string.
     */
    public String getExtendedNewick() {
        return getExtendedNewick(true);
    }

    /**
     * Obtain extended Newick representation of ACG.  Optionally Nexus metadata
     * on hybrid leaf nodes describing the alignment sites affected by the
     * conversion event.
     *
     * @param includeSiteStats if true, include count and ration of affected sites
     * @return Extended Newick string.
     */
    public String getExtendedNewick(boolean includeSiteStats) {
        return extendedNewickTraverse(root, includeSiteStats) + ";";

    }
    
    private String extendedNewickTraverse(Node node,
                                          boolean includeSiteStats) {
        StringBuilder sb = new StringBuilder();
        
        // Determine sequence of events along this node.
        class Event {
            boolean isArrival;
            double time;
            Conversion conv;
            
            public Event(boolean isArrival, double time, Conversion conv) {
                this.isArrival = isArrival;
                this.time = time;
                this.conv = conv;
            }
        }
        List<Event> events = new ArrayList<>();
        for (Conversion conv : convs) {
            if (conv.node1 == node)
                events.add(new Event(false, conv.getHeight(), conv));
            if (conv.node2 == node)
                events.add(new Event(true, conv.getHeight(), conv));
        }
        
        // Sort events from oldest to youngest.
        events.sort((Event e1, Event e2) -> {
            if (e1.time > e2.time)
                return -1;
            else
                return 1;
        });

        // Process events.
        
        int cursor = 0;
        
        double lastTime;
        if (node.isRoot())
            lastTime = Double.POSITIVE_INFINITY;
        else
            lastTime = node.getParent().getHeight();

        for (Event event : events) {

            double thisLength;
            if (Double.isInfinite(lastTime))
                thisLength = 0.0;
            else
            	// TODO Is this fine for multiple conversions along an edge? 
                thisLength = lastTime - event.time;

            if (event.isArrival) {
                String meta =  String.format(Locale.ENGLISH,
                        "[&conv=%d, relSize=%g",
                        convs.indexOf(event.conv),
                        event.conv.getSiteCount()/(double) siteCount
                );

                if (includeSiteStats) {
                	int affectedSiteCount = countAffectedSites();
        			double affectedSiteFraction = affectedSiteCount / (double) siteCount;
                    meta += String.format(Locale.ENGLISH,
                            ", affectedSites=%d, uselessSiteFraction=%g",
                            affectedSiteCount,
                            1.0-affectedSiteFraction);
                }

                if (event.conv.newickMetaDataMiddle != null)
                    meta += ", " + event.conv.newickMetaDataMiddle;

                meta += "]";

                String parentMeta;
                if (event.conv.newickMetaDataTop != null)
                    parentMeta = "[&" + event.conv.newickMetaDataTop + "]";
                else
                    parentMeta = "";

                sb.insert(cursor, "(,#" + getConversionIndex(event.conv)
                        + meta
                        + ":0.0" // TODO Is that allowed? Old version:  
//                        + (event.conv.height2-event.conv.height1) 
                        + ")"
                        + parentMeta
                        + ":" + thisLength);
                cursor += 1;
            } else {
                String meta;
                if (event.conv.newickMetaDataBottom != null)
                    meta = "[&" + event.conv.newickMetaDataBottom + "]";
                else
                    meta = "";

                sb.insert(cursor, "()#" + getConversionIndex(event.conv)
                        + meta
                        + ":" + thisLength);
                cursor += 1;
            }
            
            lastTime = event.time;
        }
        
        // Process this node and its children.

        if (!node.isLeaf()) {
            String subtree1 = extendedNewickTraverse(node.getChild(0), includeSiteStats);
            String subtree2 = extendedNewickTraverse(node.getChild(1), includeSiteStats);
            sb.insert(cursor, "(" + subtree1 + "," + subtree2 + ")");
            cursor += subtree1.length() + subtree2.length() + 3;
        }

        double thisLength;
        if (Double.isInfinite(lastTime))
            thisLength = 0.0;
        else
            thisLength = lastTime - node.getHeight();
        sb.insert(cursor, (node.getNr() + taxaTranslationOffset)
                + node.getNewickMetaData() + ":" + thisLength);
        
        return sb.toString();
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
                        conv.setHeight1(node.getHeight());
                        conv.setHeight2(node.getParent().getHeight());
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
                        conv = new Conversion();
                        convIDMap.put(convID, conv);
                    }

                    if (ctx.node().isEmpty()) {
                        String locusID;
                        for (ExtendedNewickParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
                            switch (attribCtx.attribKey.getText()) {
                                case "region":
                                    conv.setStartSite(Integer.parseInt(
                                            attribCtx.attribValue().vector().attribValue(0).getText()));
                                    conv.setEndSite(Integer.parseInt(
                                            attribCtx.attribValue().vector().attribValue(1).getText()));
                                    break;

                                case "locus":
                                    locusID = attribCtx.attribValue().getText();
                                    if (locusID.startsWith("\""))
                                        locusID = locusID.substring(1,locusID.length()-1);

                                    Locus locus = null;
                                    for (Locus thisLocus : getConvertibleLoci()) {
                                        if (thisLocus.getID().equals(locusID))
                                            locus = thisLocus;
                                    }

                                    if (locus == null)
                                        throw new IllegalArgumentException(
                                                "Locus with ID " + locusID + " not found.");

                                    conv.setLocus(locus);
                                    break;

                                default:
                                    break;
                            }
                        }
                    }
                }

                for (ExtendedNewickParser.NodeContext childCtx : ctx.node())
                    node.addChild(visitNode(childCtx));

                if (ctx.post().label() != null) {
                    node.setID(ctx.post().label().getText());
                    node.setNr(Integer.parseInt(ctx.post().label().getText())
                            - nodeNumberoffset);
                }

                node.setHeight(Double.parseDouble(ctx.post().length.getText()));

                return node;
            }
        }.visit(parseTree);

        m_nodes = root.getAllChildNodesAndSelf().toArray(m_nodes);
        nodeCount = m_nodes.length;
        leafNodeCount = root.getAllLeafNodes().size();

        setRoot(root);
        initArrays();

        for (Locus locus : getConvertibleLoci())
            convs.get(locus).clear();

        for (Conversion conv : convIDMap.values())
            addConversion(conv);

    }
    
    
    
    
    
    
    
    
    
    /*
    * StateNode implementation
    */
    
    @Override
    protected void store () {
        super.store();
        
        storedConvs.clear();

        for (Conversion conv : convs) {
        	Conversion convCopy = conv.getCopy();
        	
        	// TODO Why not do this in Conversion.getCopy() ?
        	convCopy.newickMetaDataBottom = conv.newickMetaDataBottom;
        	convCopy.newickMetaDataMiddle = conv.newickMetaDataMiddle;
        	convCopy.newickMetaDataTop = conv.newickMetaDataTop;

        	// Adapt the stored conversion to reference the stored nodes
            convCopy.setNode1(m_storedNodes[conv.getNode1().getNr()]);
            convCopy.setNode2(m_storedNodes[conv.getNode2().getNr()]);

            convCopy.setConversionGraph(this);

            storedConvs.add(convCopy);
        }
    }
    
    @Override
    public void restore() {
        super.restore();
        
        // Swap conversions with storedConversions
        List<Conversion> tmp = storedConvs;
        storedConvs = convs;
        convs = tmp;

        cfEventList.makeDirty();
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

    @Override
    public void log(long nSample, PrintStream out) {
        ConversionGraph arg = (ConversionGraph) getCurrent();
        
        out.print(String.format("tree STATE_%d = [&R] %s",
                nSample, arg.getExtendedNewick()));
    }
   
}
