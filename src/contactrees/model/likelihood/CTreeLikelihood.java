package contactrees.model.likelihood;

import beast.core.Description;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import contactrees.MarginalTree;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class CTreeLikelihood extends TreeLikelihood {

    public BranchRateModel.Base marginalTreeClock;

    public BranchRateModel.Base stealClockModel() {
        BranchRateModel.Base clockModel =  branchRateModelInput.get();
        branchRateModelInput.setValue(new StrictClockModel(), this);
        return clockModel;
    }

    @Override
    public void initAndValidate() {

        MarginalTree mTree = (MarginalTree) treeInput.get();
        marginalTreeClock =  branchRateModelInput.get();

        if (marginalTreeClock != null) {
            // Forward branchRateModel to the tree and use fixed strict clock model here instead
            mTree.setBranchRateModel(marginalTreeClock);
            branchRateModelInput.setValue(null, this);
        }

        super.initAndValidate();

        if (marginalTreeClock != null) {
            // Reset branchRateModelInput for XML file
            branchRateModelInput.setValue(marginalTreeClock, this);
        }

        // Update the marginal tree with the new model
        mTree.recalculate();
        mTree.makeOutdated();
    }


}
