package contactrees.model.likelihood;

import beast.app.beauti.Beauti;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.likelihood.TreeLikelihood;
import contactrees.MarginalTree;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class CTreeLikelihood extends TreeLikelihood {

    public Input<MarginalTree> marginalTreeInput = new Input<>("marginalTree", "marginal tree based on actual tree");

    public BranchRateModel.Base marginalTreeClock;

    @Override
    public void initAndValidate() {

        if (!Beauti.isInBeauti()) {

            MarginalTree mTree;
            if (marginalTreeInput.get() != null) {
                mTree = marginalTreeInput.get();
                treeInput.set(mTree);
            } else {
                // If the marginalTreeInput is null, the treeInput has to be a marginal tree
                mTree = (MarginalTree) treeInput.get();
                marginalTreeInput.set(mTree);
            }

            marginalTreeClock =  branchRateModelInput.get();
            if (marginalTreeClock != null) {
                // Forward branchRateModel to the tree and use fixed strict clock model here instead
                mTree.setBranchRateModel(marginalTreeClock);
                branchRateModelInput.setValue(null, this);
            }
        }

        super.initAndValidate();

        if (Beauti.isInBeauti() && marginalTreeClock != null) {
            // Reset branchRateModelInput for XML file
            branchRateModelInput.setValue(marginalTreeClock, this);
        }

        // Update the marginal tree with the new model
        if (!Beauti.isInBeauti()) {
            MarginalTree mTree = marginalTreeInput.get();
            mTree.recalculate();
            mTree.makeOutdated();
            marginalTreeInput.setValue(null, this);
        }

    }


}
