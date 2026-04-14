package contactrees.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("CF subtree exchange operator.")
public class CFSubtreeExchange extends CFOperator {

    public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
            "Whether to use narrow exchange. (Default true.)", true);

    int count = 0;

    @Override
    public double proposal() {
        double logHGF = 0.0;

        if (acg.getLeafNodeCount()<3)
            return Double.NEGATIVE_INFINITY;

        Node srcNode, srcNodeP, destNode, destNodeP;
        if (isNarrowInput.get()) {

            // Narrow exchange selection:
            do {
                srcNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
            } while (srcNode.getParent().isRoot());

            srcNodeP = srcNode.getParent();
            destNode = getSibling(srcNodeP);
            destNodeP = destNode.getParent();

        } else {

            // Wide exchange selection:
            srcNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
            srcNodeP = srcNode.getParent();

            do {
                destNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
            } while(destNode == srcNode || destNode.getParent() == srcNode.getParent()); 
            								// TODO We never allow wide to be narrow? 
            destNodeP = destNode.getParent();
        }
        
        assert srcNodeP != null;
        assert destNodeP != null;
        double t_srcNodeP = srcNodeP.getHeight();
        double t_destNodeP = destNodeP.getHeight();

        // Reject if substitution would result in negative branch lengths:
        // TODO include those in the while loop above?
        if (destNode == srcNodeP || srcNode == destNodeP
                || destNode.getHeight()>t_srcNodeP
                || srcNode.getHeight()>t_destNodeP)
            return Double.NEGATIVE_INFINITY;

        if (t_srcNodeP > t_destNodeP) {
            Node srcNodeS = getSibling(srcNode);

            logHGF += collapseConversions(srcNode, destNodeP, t_destNodeP);

            if (!srcNodeS.isRoot() && srcNodeS.getLength() == 0.0)
                srcNodeS = srcNodeS.getParent();

            logHGF -= expandConversions(destNode, srcNodeS, t_srcNodeP);

        } else {
            Node srcNodeS = getSibling(srcNode);

            logHGF -= expandConversions(srcNode, destNodeP, t_destNodeP);

            if (srcNodeP.isRoot())
                acg.setRoot(srcNodeP);

            logHGF += collapseConversions(destNode, srcNodeS, t_srcNodeP);
        }
        
        if (logHGF > Double.NEGATIVE_INFINITY) 
        	assert !acg.isInvalid() : "CFSTX proposed invalid state.";
    	
        return logHGF;
    }
}
