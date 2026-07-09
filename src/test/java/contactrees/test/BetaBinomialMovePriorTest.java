package contactrees.test;

import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.statistics.distribution.BetaDistribution;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.junit.Test;

import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import contactrees.model.BetaBinomialMovePrior;

/**
 * Guards the migration of BetaBinomialMovePrior off commons-math3 onto the
 * core-provided commons-numbers (LogBeta) and commons-statistics + commons-rng
 * (Beta sampling). Exercises both replaced code paths.
 */
public class BetaBinomialMovePriorTest extends ContactreesTest {

    /** calculateLogP() now uses org.apache.commons.numbers.gamma.LogBeta. */
    @Test
    public void logPIsFinite() {
        RealScalarParam<PositiveReal> alpha = new RealScalarParam<>(2.0, PositiveReal.INSTANCE);
        RealScalarParam<PositiveReal> beta = new RealScalarParam<>(16.0, PositiveReal.INSTANCE);

        BetaBinomialMovePrior prior = new BetaBinomialMovePrior();
        prior.initByName("network", acg, "blockSet", blockSet, "alpha", alpha, "beta", beta);

        double logP = prior.calculateLogP();
        assertTrue("logP should be finite, was " + logP, Double.isFinite(logP));
    }

    /** sample() now draws Beta variates via commons-statistics + commons-rng. */
    @Test
    public void betaSamplerProducesUnitIntervalValues() {
        Random random = new Random(42);
        UniformRandomProvider rng = RandomSource.XO_RO_SHI_RO_128_PP.create(random.nextLong());
        ContinuousDistribution.Sampler sampler = BetaDistribution.of(2.0, 16.0).createSampler(rng);
        for (int i = 0; i < 1000; i++) {
            double x = sampler.sample();
            assertTrue("Beta sample out of [0,1]: " + x, x >= 0.0 && x <= 1.0);
        }
    }
}
