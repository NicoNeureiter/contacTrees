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
 * For clocks whose rate does not depend on the branch (strict clocks, or a FreezableClock
 * without frozen taxa around one) every marginal branch has the clock rate. For any other clock
 * defined on the clonal frame (relaxed clocks, frozen taxa, ...) the marginal tree records the
 * CF path of each marginal branch and computes the path-averaged rates
 * ({@link MarginalTree#setBranchRateClock}).
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
    protected MarginalTree tree;
    protected boolean branchIndependent;

    @Override
    public void initAndValidate() {
        clock = clockInput.get();
        tree = treeInput.get();
        branchIndependent = isBranchIndependent(clock);
        if (!branchIndependent) {
            if (tree == null)
                throw new IllegalArgumentException("MarginalClockModel " + getID() + ": the clock "
                        + clock.getClass().getSimpleName() + " has branch-dependent rates, which requires the tree input.");
            tree.setBranchRateClock(clock);
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
        if (branchIndependent)
            return clock.getRateForBranch(node);
        tree.ensureUpToDate();
        return tree.getBranchRate(node.getNr());
    }

    @Override
    protected boolean requiresRecalculation() {
        boolean clockDirty = (clock instanceof CalculationNode cn) && cn.somethingIsDirty();
        if (branchIndependent)
            return clockDirty;

        // Recompute the rates now; nodes whose rate changed are flagged dirty in the tree
        if (clockDirty)
            tree.invalidateRates();
        tree.ensureUpToDate();
        return clockDirty || tree.somethingIsDirty();
    }

    /*
     * After a proposal that only changed the clock, BEAST accepts or restores this wrapper but not
     * the marginal tree (which is not downstream of the clock), although the tree's rates changed.
     * Forward both; they are no-ops for a tree that BEAST handles itself.
     */

    @Override
    protected void store() {
        super.store();
        if (!branchIndependent)
            tree.store();
    }

    @Override
    protected void restore() {
        super.restore();
        if (!branchIndependent)
            tree.restore();
    }

    @Override
    protected void accept() {
        super.accept();
        if (!branchIndependent)
            tree.acceptProposal();
    }

}
