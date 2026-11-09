package contactrees.model.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.ProgramStatus;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.substitutionmodel.BasicGeneralSubstitutionModel;
import contactrees.MarginalTree;
import contactrees.model.MarginalClockModel;

@Description("Tree likelihood of a block's marginal tree. A clock model given as branchRateModel is "
        + "defined on the clonal frame and is wrapped in a MarginalClockModel.")
public class CTreeLikelihood extends TreeLikelihood {
    final public Input<Boolean> shareMatricesInput = new Input<>("shareTransitionMatrices",
            "share transition probability matrices with all other likelihoods using the same substitution model "
            + "(only used by the Java likelihood core and for substitution models based on BasicGeneralSubstitutionModel)", true);

    SharedTransitionMatrixCache.CachingSubstitutionModel cachingModel;

    boolean isInBeauti(){
        return ProgramStatus.name.equals("BEAUti");
    }

    @Override
    public void initAndValidate() {
        if (isInBeauti()) {
            super.initAndValidate();
            return;
        }

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

        if (beagle == null && shareMatricesInput.get() && SharedTransitionMatrixCache.supports(substitutionModel)) {
            cachingModel = new SharedTransitionMatrixCache.CachingSubstitutionModel(
                    (BasicGeneralSubstitutionModel) substitutionModel);
            substitutionModel = cachingModel;
        }
    }

    @Override
    public double calculateLogP() {
        if (cachingModel != null)
            cachingModel.beginEvaluation();
        return super.calculateLogP();
    }

    public boolean usesSharedTransitionMatrices() {
        return cachingModel != null;
    }

}
