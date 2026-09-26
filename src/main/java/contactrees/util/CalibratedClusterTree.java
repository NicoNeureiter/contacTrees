package contactrees.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.evolution.tree.Node;
import beast.base.inference.distribution.ParametricDistribution;
import beastlabs.evolution.tree.ConstrainedClusterTree;

@Description("A ConstrainedClusterTree that starts the chain at the ages its calibrations imply. " +
        "The clustering produces node heights in substitutions, which ConstrainedClusterTree only " +
        "converts to time by dividing through its clock.rate input. It can also stretch a clade onto " +
        "an MRCAPrior, but its handlebounds() reads inverseCumulativeProbability(0) and (1), which are " +
        "infinite for the usual Normal or LogNormal calibration, so in practice the heights are left " +
        "in substitution units and the chain starts orders of magnitude away from its calibrations. " +
        "This subclass instead rescales each calibrated clade onto the median of its own MRCAPrior, " +
        "so no clock.rate has to be guessed and nothing has to be re-measured when the data change. " +
        "Set rootHeight to override this with a fixed root age, e.g. when there is no usable calibration.")
public class CalibratedClusterTree extends ConstrainedClusterTree {

    final public Input<Double> rootHeightInput = new Input<>("rootHeight",
            "scale the starting tree so that its root sits at this height, in the time units of the " +
            "tree. Optional: when omitted, each calibrated clade is scaled onto the median of its own " +
            "MRCAPrior instead, which needs no configuration and cannot fall out of step with the " +
            "calibration. Provide this only when there is no calibration to derive a scale from.");

    final public Input<Double> quantileInput = new Input<>("quantile",
            "quantile of each MRCAPrior to place its clade at; 0.5 (the median) by default. The median " +
            "is preferred over the mean because it stays inside the bulk of a skewed calibration such " +
            "as a LogNormal or Gamma.", 0.5);

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        double quantile = quantileInput.get();
        if (quantile <= 0.0 || quantile >= 1.0) {
            throw new IllegalArgumentException("quantile must lie strictly between 0 and 1; got " + quantile);
        }

        if (rootHeightInput.get() != null) {
            double target = rootHeightInput.get();
            if (target <= 0.0) {
                throw new IllegalArgumentException("rootHeight must be strictly positive; got " + target);
            }
            scaleClade(getRoot(), target);
        } else {
            scaleOntoCalibrations(quantile);
        }

        // super.initAndValidate() already pushed the unscaled tree into the `initial` tree;
        // repeat that now the heights have moved.
        initStateNodes();
    }

    /**
     * Rescale every calibrated clade so that its height is the requested quantile of its own
     * MRCAPrior.
     *
     * Calibrations are visited from the root down, matching ConstrainedClusterTree.handlebounds():
     * scaling an outer clade also moves the nested ones, so the nested calibration has to be
     * applied afterwards to have the last word.
     */
    private void scaleOntoCalibrations(double quantile) {
        // ConstrainedClusterTree keeps its calibrations in a local variable, so re-derive them
        // the same way. The constraint lists are only outputs here and are discarded.
        List<MRCAPrior> calibrations = collectCalibrations(taxaNames(), m_initial.get(),
                allConstraints.get(), calibrationsInput.get(),
                new ArrayList<boolean[]>(), new ArrayList<Integer>());

        Map<Node, MRCAPrior> nodeToBoundMap = new HashMap<>();
        findConstrainedNodes(calibrations, getRoot(), nodeToBoundMap);
        applyCalibrations(getRoot(), nodeToBoundMap, quantile);
    }

    private void applyCalibrations(Node node, Map<Node, MRCAPrior> nodeToBoundMap, double quantile) {
        if (node.isLeaf()) {
            return;
        }
        MRCAPrior calibration = nodeToBoundMap.get(node);
        if (calibration != null && calibration.distInput.get() != null) {
            ParametricDistribution distr = calibration.distInput.get();
            double target;
            try {
                target = distr.inverseCumulativeProbability(quantile) + distr.offsetInput.get();
            } catch (Exception e) {
                throw new IllegalArgumentException("Could not evaluate the " + quantile
                        + " quantile of calibration " + calibration.getID(), e);
            }
            if (Double.isFinite(target) && target > 0.0) {
                scaleClade(node, target);
            }
        }
        for (Node child : node.getChildren()) {
            applyCalibrations(child, nodeToBoundMap, quantile);
        }
    }

    /**
     * Scale the clade below {@code node} so that {@code node} ends up at {@code targetHeight},
     * then lift any ancestor that this would have left below its own child.
     */
    private void scaleClade(Node node, double targetHeight) {
        double height = node.getHeight();
        if (!(height > 0.0)) {
            // A degenerate clade has no scale to stretch; place it directly.
            node.setHeight(targetHeight);
        } else {
            scale(node, targetHeight / height);
        }
        liftAncestors(node);
    }

    /** ConstrainedClusterTree.taxaNames is package private, so resolve it the same way it does. */
    private List<String> taxaNames() {
        if (dataInput.get() != null) {
            return dataInput.get().getTaxaNames();
        }
        return m_taxonset.get().asStringList();
    }

    /** ConstrainedClusterTree.scale() is private, so repeat it here. */
    private void scale(Node node, double factor) {
        if (node.isLeaf() || node.isFake()) {
            return;
        }
        double height = node.getHeight() * factor;
        for (Node child : node.getChildren()) {
            scale(child, factor);
        }
        for (Node child : node.getChildren()) {
            height = Math.max(height, child.getHeight());
        }
        node.setHeight(height);
    }

    private void liftAncestors(Node node) {
        Node child = node;
        Node parent = child.getParent();
        while (parent != null && parent.getHeight() < child.getHeight()) {
            parent.setHeight(child.getHeight() + epsilonInput.get());
            child = parent;
            parent = child.getParent();
        }
    }
}
