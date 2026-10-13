package contactrees.test;

import static org.junit.Assert.assertThrows;

import org.junit.Test;

import contactrees.Block;
import contactrees.MarginalNode;
import contactrees.MarginalTree;
import contactrees.operators.ACGOperator;

/**
 * {@code ACGOperator.proposal()} calls {@code MarginalTree.assertConsistent()} for every
 * validating marginal tree after every successful proposal -- but that check's own contract
 * (see {@link MarginalTreeValidatorTest}) is "the tree's CURRENT, not-yet-recalculated display
 * already matches a fresh rebuild". That is only true right after an operator whose doProposal()
 * proactively recalculates every marginal tree it affects (the Gibbs-choreography operators).
 * {@code CFUniform} (like every {@code CFOperator} subclass) changes a CF node's height directly
 * and relies on BEAST's own lazy requiresRecalculation() to catch up later -- by design, for
 * performance (see the plan's store/restore work). Calling assertConsistent() right after such an
 * operator's proposal() is therefore a false positive, not a real staleness bug.
 *
 * These tests use minimal synthetic {@link ACGOperator} subclasses that mirror the exact
 * mechanism (mutate a CF node's height, with or without recalculating), rather than a real
 * {@code CFUniform}, to isolate {@link ACGOperator#operatorUpdatesMarginalTree()}'s gating effect
 * from CFUniform's own random-selection logic (which hits an unrelated, pre-existing
 * getInternalNodeCount()==0 quirk on the acg2/acg fixtures that no prior test has exercised).
 *
 * @author Nico Neureiter
 */
public class ACGOperatorMarginalTreeSyncTest extends ContactreesTest {

    private MarginalTree newValidatingMarginalTree(Block block) {
        MarginalTree marginalTree = new MarginalTree();
        marginalTree.initByName("network", acg2, "block", block, "nodetype", MarginalNode.class.getName(),
                "validationLevel", 1);
        return marginalTree;
    }

    /** Mutates a CF node's height and returns, WITHOUT recalculating any marginal tree -- the
     *  default (operatorUpdatesMarginalTree() == false), matching CFUniform's own behavior. */
    private class NonSyncingHeightChanger extends ACGOperator {
        @Override
        protected double doProposal() {
            node2_6.setHeight(node2_6.getHeight() + 0.1);
            return 0.0;
        }
    }

    /** Same mutation, but overrides operatorUpdatesMarginalTree() to true AND actually keeps
     *  the given tree recalculated -- the class of operator the flag is meant for. */
    private class SyncingHeightChanger extends ACGOperator {
        private final MarginalTree marginalTree;
        SyncingHeightChanger(MarginalTree marginalTree) {
            this.marginalTree = marginalTree;
        }
        @Override
        protected boolean operatorUpdatesMarginalTree() {
            return true;
        }
        @Override
        protected double doProposal() {
            node2_6.setHeight(node2_6.getHeight() + 0.1);
            marginalTree.recalculate();
            return 0.0;
        }
    }

    /** Same mutation, claims (via the flag) to keep trees synced but does NOT -- a buggy
     *  operator. Proves the flag actually re-enables the check when true, not just always off. */
    private class FalselyClaimingSyncOperator extends ACGOperator {
        @Override
        protected boolean operatorUpdatesMarginalTree() {
            return true;
        }
        @Override
        protected double doProposal() {
            node2_6.setHeight(node2_6.getHeight() + 0.1);
            return 0.0;
        }
    }

    private ACGOperator initOperator(ACGOperator op) {
        op.initByName("acg", acg2, "weight", 1.0);
        return op;
    }

    @Test
    public void testNonSyncingOperatorDoesNotFalselyReportStaleTree() throws Exception {
        Block block = blockSet2.getBlocks().get(0); // untouched -- no moves
        newValidatingMarginalTree(block);

        ACGOperator op = initOperator(new NonSyncingHeightChanger());

        // Assertion is that this does not throw: op.proposal() must not report a false "stale
        // tree" just because this operator (like CFUniform) defers recalculation.
        op.proposal();
    }

    @Test
    public void testSyncingOperatorPassesAssertConsistent() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        MarginalTree marginalTree = newValidatingMarginalTree(block);

        ACGOperator op = initOperator(new SyncingHeightChanger(marginalTree));

        // The flag being true must still let a genuinely well-behaved operator pass.
        op.proposal();
    }

    @Test
    public void testFalseSyncClaimStillThrows() throws Exception {
        Block block = blockSet2.getBlocks().get(0);
        newValidatingMarginalTree(block);

        ACGOperator op = initOperator(new FalselyClaimingSyncOperator());

        // The flag must actually gate the check (not just be ignored/always-false): an operator
        // that claims to sync trees but doesn't must still be caught.
        assertThrows(IllegalStateException.class, op::proposal);
    }

}
