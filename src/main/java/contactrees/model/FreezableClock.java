package contactrees.model;

import java.util.ArrayList;

import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;


public class FreezableClock extends Base {

    public Input<BranchRateModel> clockInput = new Input<>(
            "clock",
            "The clock model defining the rate along unfrozen branches.",
            Input.Validate.REQUIRED);

    public Input<ArrayList<Taxon>> frozenTaxaInput = new Input<>(
            "frozenTaxa",
            "Taxa for which the last branch should have a fixed branch rate of 0.",
            new ArrayList<>());

    BranchRateModel clock;
    ArrayList<String> frozenTaxa;
    boolean anythingFrozen;

    @Override
    public void initAndValidate() {
        clock = clockInput.get();

        frozenTaxa = new ArrayList<>();
        for (final Taxon taxon : frozenTaxaInput.get()) {
            frozenTaxa.add(taxon.getID());
        }

        anythingFrozen = (frozenTaxa.size() > 0);
    }

    @Override
    public double getRateForBranch(Node node) {
        if (anythingFrozen && frozenTaxa.contains(node.getID()))
            return 0.0;
        else
            return clock.getRateForBranch(node);
    }

    @Override
    public boolean requiresRecalculation() {
        return ((CalculationNode) clock).somethingIsDirty();
    }

    @Override
    protected void restore() {
        super.restore();
    }

    @Override
    protected void store() {
        super.store();
    }

}
