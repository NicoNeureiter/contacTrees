package contactrees.model.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import contactrees.MarginalTree;
import contactrees.model.MarginalClockModel;

@Description("Tree likelihood of a block's marginal tree. A clock model given as branchRateModel is "
        + "defined on the clonal frame and is wrapped in a MarginalClockModel.")
public class CTreeLikelihood extends TreeLikelihood {
	final public Input<MarginalTree> marginalTreeInput = new Input<>("marginalTree", "marginal tree based on actual tree");

    boolean isInBeauti(){
        return ProgramStatus.name.equals("BEAUti");
    }

    @Override
    public void initAndValidate() {
        if (isInBeauti()) {
            super.initAndValidate();
            return;
        }

        if (marginalTreeInput.get() != null)
            treeInput.set(marginalTreeInput.get());

        if (!(treeInput.get() instanceof MarginalTree mTree))
            throw new IllegalArgumentException("CTreeLikelihood " + getID() + " requires a MarginalTree as tree.");

        Base clock = branchRateModelInput.get();
        if (clock != null && !(clock instanceof MarginalClockModel)) {
            MarginalClockModel marginalClock = new MarginalClockModel();
            marginalClock.setID("MarginalClock." + mTree.getID());
            marginalClock.initByName("clock", clock, "tree", mTree);
            clock.getOutputs().remove(this);
            branchRateModelInput.setValue(marginalClock, this);
        }

        super.initAndValidate();
        marginalTreeInput.setValue(null, this);
    }

}
