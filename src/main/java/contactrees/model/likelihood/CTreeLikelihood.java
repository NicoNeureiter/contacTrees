package contactrees.model.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.evolution.branchratemodel.StrictClockModel;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import contactrees.MarginalTree;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class CTreeLikelihood extends TreeLikelihood {
	final public Input<MarginalTree> marginalTreeInput = new Input<>("marginalTree", "marginal tree based on actual tree");

    public Base marginalTreeClock;

//<<<<<<< Updated upstream
//
//    @Override
//    public void initAndValidate() {
//
//        if (!isInBeauti()) {
//
//            MarginalTree mTree;
//            if (marginalTreeInput.get() != null) {
//                mTree = marginalTreeInput.get();
//                treeInput.set(mTree);
//            } else {
//                // If the marginalTreeInput is null, the treeInput has to be a marginal tree
//                mTree = (MarginalTree) treeInput.get();
//                marginalTreeInput.set(mTree);
//            }
//
//            marginalTreeClock =  branchRateModelInput.get();
//            if (marginalTreeClock != null) {
//                // Forward branchRateModel to the tree and use fixed strict clock model here instead
//                mTree.setBranchRateModel(marginalTreeClock);
//                branchRateModelInput.setValue(null, this);
//            }
//=======
    boolean isInBeauti(){
        return ProgramStatus.name.equals("BEAUti");
    }

    @Override
    public void initAndValidate() {
    	if (!isInBeauti() && marginalTreeInput.get() != null) {
            MarginalTree mTree = marginalTreeInput.get();
            treeInput.set(mTree);
    	}
        marginalTreeClock =  branchRateModelInput.get();    		
        if (!isInBeauti() && marginalTreeClock != null) {
            // Forward branchRateModel to the tree and use fixed strict clock model here instead
            ((MarginalTree)treeInput.get()).setBranchRateModel(marginalTreeClock);
            branchRateModelInput.setValue(new StrictClockModel(), this);
        }

        super.initAndValidate();

        if (isInBeauti() && marginalTreeClock != null) {
            // Reset branchRateModelInput for XML file
            branchRateModelInput.setValue(marginalTreeClock, this);
        }

        // Update the marginal tree with the new model
        if (!isInBeauti()) {
            MarginalTree mTree = (MarginalTree) treeInput.get();
            mTree.recalculate();
            mTree.makeOutdated();
            marginalTreeInput.setValue(null, this);
        }
    }


}
