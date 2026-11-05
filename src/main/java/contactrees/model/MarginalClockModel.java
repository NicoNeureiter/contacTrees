package contactrees.model;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.spec.evolution.branchratemodel.Base;
import contactrees.MarginalTree;

/**
 * Branch rates for the branches of a {@link MarginalTree}, derived from a clock model defined
 * on the clonal frame. Marginal trees keep their true heights, so the rate of a marginal branch
 * is the time-weighted mean of the clonal frame rates along the path the branch follows.
 *
 * Currently only clocks whose rate does not depend on the branch are supported (strict clocks,
 * or a FreezableClock without frozen taxa around one); then every marginal branch has the clock
 * rate.
 *
 * @author Nico Neureiter
 */
@Description("Branch rates of a marginal tree, derived from a clock model on the clonal frame.")
public class MarginalClockModel extends Base {

    public Input<BranchRateModel> clockInput = new Input<>(
            "clock",
            "The clock model defined on the clonal frame.",
            Input.Validate.REQUIRED);
    public Input<MarginalTree> treeInput = new Input<>(
            "tree",
            "The marginal tree whose branch rates are provided.");

    protected BranchRateModel clock;

    @Override
    public void initAndValidate() {
        clock = clockInput.get();
        if (!isBranchIndependent(clock)) {
            throw new UnsupportedOperationException("MarginalClockModel " + getID() + ": clock model "
                    + clock.getClass().getSimpleName() + " has branch-dependent rates, which are not yet "
                    + "supported on marginal trees. Use a strict clock.");
        }
    }

    /**
     * @return true iff the clock returns the same rate for every branch.
     */
    public static boolean isBranchIndependent(BranchRateModel clock) {
        if (clock instanceof beast.base.spec.evolution.branchratemodel.StrictClockModel)
            return true;
        if (clock instanceof beast.base.evolution.branchratemodel.StrictClockModel)
            return true;
        if (clock instanceof FreezableClock freezable)
            return !freezable.anythingFrozen && isBranchIndependent(freezable.clock);
        return false;
    }

    @Override
    public double getRateForBranch(Node node) {
        // Branch independent: the rate of any branch is the clock rate.
        return clock.getRateForBranch(node);
    }

    @Override
    protected boolean requiresRecalculation() {
        return (clock instanceof CalculationNode cn) && cn.somethingIsDirty();
    }

}
