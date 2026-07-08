open module contactrees {
    requires beast.pkgmgmt;
    requires beast.base;
    requires beast.fx;

    requires com.google.common;
    requires org.antlr.antlr4.runtime;
    requires org.apache.commons.statistics.distribution;

    requires java.desktop;
    requires java.management;
    requires java.naming;
    requires java.xml;

    requires commons.math3;

    exports contactrees;
    exports contactrees.acgannotator;
    exports contactrees.model;
    exports contactrees.model.likelihood;
    exports contactrees.operators;
    exports contactrees.util;
    exports contactrees.util.parsers;

    provides beast.base.core.BEASTInterface with
        contactrees.ACGStatsLogger,
        contactrees.ACGWithBlocks,
        contactrees.ACGWithMetaDataLogger,
        contactrees.Block,
        contactrees.BlockSet,
        contactrees.ConversionGraph,
        contactrees.MarginalNode,
        contactrees.MarginalNodeSlow,
        contactrees.MarginalTree,
        contactrees.MarginalTreeSlow,
        contactrees.model.ACGDistribution,
        contactrees.model.ACGSimulator,
        contactrees.model.BetaBinomialMovePrior,
        contactrees.model.ConversionMovePrior,
        contactrees.model.ConversionPrior,
        contactrees.model.ConversionRate,
        contactrees.model.FreezableClock,
        contactrees.model.likelihood.CTreeLikelihood,
        contactrees.model.SimulatedACGWithBlocks,
        contactrees.operators.ACGScaler,
        contactrees.operators.AddRemoveConversion,
        contactrees.operators.AddRemoveConversionGibbs,
        contactrees.operators.CFConversionSwap,
        contactrees.operators.CFSubtreeExchange,
        contactrees.operators.CFSubtreeSlide,
        contactrees.operators.CFUniform,
        contactrees.operators.CFWilsonBalding,
        contactrees.operators.ConversionSplit,
        contactrees.operators.ConvertedEdgeFlip,
        contactrees.operators.ConvertedEdgeHop,
        contactrees.operators.ConvertedEdgeHopGibbs,
        contactrees.operators.ConvertedEdgeSlide,
        contactrees.operators.GibbsSampleMovesPerBlock,
        contactrees.operators.GibbsSampleMovesPerConversion,
        contactrees.operators.ResampleBorrowings,
        contactrees.RandomACG,
        contactrees.test.SimulatorOperator,
        contactrees.test.StationarityTestSchedule;
}