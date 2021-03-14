/**
 *
 */
package contactrees.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.special.Beta;
import org.apache.commons.math3.distribution.BetaDistribution;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.Conversion;
import contactrees.ConversionGraph;

/**
 *
 *
 * @author Nico Neureiter
 */
public class BetaBinomialMovePrior extends Distribution {

    final public Input<ConversionGraph> networkInput = new Input<>(
            "network",
            "The conversion graph containing the conversion edges.",
            Input.Validate.REQUIRED);

    final public Input<BlockSet> blockSetInput = new Input<>(
            "blockSet",
            "The blockSet containing all blocks that evolve along the clonal frame and "
            + "move over conversion edges in the network.",
            Input.Validate.REQUIRED);

    final public Input<RealParameter> alphaInput = new Input<>(
            "alpha",
            "Alpha parameter of the betaBinomial prior on the number of blocks moving over a conversion edge.",
            Input.Validate.REQUIRED);

    final public Input<RealParameter> betaInput = new Input<>(
            "beta",
            "Beta parameter of the betaBinomial prior on the number of blocks moving over a conversion edge.",
            Input.Validate.REQUIRED);

    ConversionGraph acg;
    BlockSet blockSet;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        acg = networkInput.get();
        blockSet = blockSetInput.get();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;
        double alpha = alphaInput.get().getValue();
        double beta = betaInput.get().getValue();
        int n = blockSet.getBlockCount();
        for (Conversion conv : acg.getConversions()) {
            int k = blockSet.countAffectedBlocks(conv);
            logP += Beta.logBeta(alpha + k, beta + n - k) - Beta.logBeta(alpha, beta);
        }
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        arguments.add(blockSetInput.get().getID());
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(networkInput.get().getID());
        conditions.add(alphaInput.get().getID());
        conditions.add(betaInput.get().getID());
        return conditions;
    }

    @Override
    public void sample(State state, Random random) {
        if (sampledFlag)
            return;
        sampledFlag = true;

        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

        acg = networkInput.get();
        blockSet = blockSetInput.get();
        double alpha = alphaInput.get().getValue();
        double beta = betaInput.get().getValue();

        for (Block block : blockSet) {
            // Remove old conversion moves
            block.removeAllMoves();
        }
        // Sample new conversion moves
        BetaDistribution pMovePrior = new BetaDistribution(alpha, beta);
        for (Conversion conv : acg.getConversions()) {
            double pMove = pMovePrior.sample();
            for (Block block : blockSet) {
                if (random.nextDouble() < pMove) {
                    block.addMove(conv);
                }
            }
        }
    }

}
