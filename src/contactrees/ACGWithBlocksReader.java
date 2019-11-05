package contactrees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import beast.evolution.tree.Node;
import contactrees.util.parsers.ExtendedNewickBaseVisitor;
import contactrees.util.parsers.ExtendedNewickLexer;
import contactrees.util.parsers.ExtendedNewickParser;

public class ACGWithBlocksReader {

	public ConversionGraph acg;
	public BlockSet blockSet;

//	public ACGWithBlocksReader(ConversionGraph acg, BlockSet blockSet) {
	public ACGWithBlocksReader(int nBlocks) {
		acg = new ConversionGraph();
		acg.initAndValidate();
		
		blockSet = new BlockSet();
		ArrayList<Block> blocks = new ArrayList<>();
		for (int i=0; i<nBlocks; i++) {
			Block b = new Block();
			b.initAndValidate();
			blocks.add(b);
		}
		
		blockSet.initByName("network", acg, "block", blocks);
		
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

        Map<Integer, List<Conversion>> blockConvs = new HashMap<>();
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
                                case "affectedBlocks":
                                	List<Integer> blockIDs = parseIntList(attribCtx.attribValue().getText());
                                	for (int blockID : blockIDs) 
                                		blockSet.addBlockMove(conv, blockID);
//	                                	if (!blockConvs.containsKey(blockID))
//	                                		blockConvs.put(blockID, new ArrayList<Conversion>());
//	                                	blockConvs.get(blockID).add(conv);
//                                	}
                                	
                                    break;

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
                    node.setNr(Integer.parseInt(ctx.post().label().getText())
                            - nodeNumberoffset);
                }

                node.setHeight(Double.parseDouble(ctx.post().length.getText()));

                return node;
            }
        }.visit(parseTree);

        acg.initAfterParsingFromNewick(root, acg.getNodesAsArray());

        for (Conversion conv : convIDMap.values())
            acg.addConversion(conv);
        
//        // Initialize blockSet
//        for (int iBlock : blockConvs.keySet()) {
//        	for (Conversion conv : blockConvs.get(iBlock))
//        		blockSet.addBlockMove(conv, iBlock);
//        }
    }

    
    protected List<Integer> parseIntList(String listAsString) {
    	ArrayList<Integer> ints = new ArrayList<>();
    	
    	String listAsStringStripped = listAsString.substring(2, listAsString.length()-2);
    	if (listAsStringStripped.length() == 0)
    		return ints;

    	for (String s : listAsStringStripped.split(","))
    		ints.add(Integer.parseInt(s));

    	return ints;
    }
    
    protected int parseConvID(String sConvID) {
    	return Integer.parseInt(sConvID.substring(1));
    }
    
}