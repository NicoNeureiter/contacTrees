/**
 *
 */
package contactrees.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.numbers.gamma.LogBeta;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.statistics.distribution.BetaDistribution;
import org.apache.commons.statistics.distribution.ContinuousDistribution;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.type.RealScalar;
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

    final public Input<RealScalar<PositiveReal>> alphaInput = new Input<>(
            "alpha",
            "Alpha parameter of the betaBinomial prior on the number of blocks moving over a conversion edge.",
            Input.Validate.REQUIRED);

    final public Input<RealScalar<PositiveReal>> betaInput = new Input<>(
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
        double alpha = alphaInput.get().get();
        double beta = betaInput.get().get();
        int n = blockSet.getBlockCount();
        for (Conversion conv : acg.getConversions()) {
            int k = blockSet.countAffectedBlocks(conv);
            logP += LogBeta.value(alpha + k, beta + n - k) - LogBeta.value(alpha, beta);
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
        maybeAddInputToConditions(conditions, alphaInput);
        maybeAddInputToConditions(conditions, betaInput);
        return conditions;
    }

    /**
     * Add an input to the `conditions` list if it is a BEASTInterface.
     * Otherwise, we assume it is constant and doesn't need to be sampled.
     */
    void maybeAddInputToConditions(List<String> conditions, Input<?> input) {
        if (input.get() instanceof BEASTInterface beastInterface) {
            conditions.add(beastInterface.getID());
        }
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
        double alpha = alphaInput.get().get();
        double beta = betaInput.get().get();

        for (Block block : blockSet) {
            // Remove old conversion moves
            block.removeAllMoves();
        }
        // Sample new conversion moves (RNG seeded from the passed Random for reproducibility)
        UniformRandomProvider rng = RandomSource.XO_RO_SHI_RO_128_PP.create(random.nextLong());
        ContinuousDistribution.Sampler pMovePrior = BetaDistribution.of(alpha, beta).createSampler(rng);
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
