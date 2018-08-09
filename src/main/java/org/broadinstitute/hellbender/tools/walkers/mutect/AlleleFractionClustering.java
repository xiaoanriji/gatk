package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

public class AlleleFractionClustering {
    private BetaDistributionShape highConfidenceDistribution;
    private BetaDistributionShape lowConfidenceDistribution;
    private double log10HighConfidencePrior;
    private double log10LowConfidencePrior;
    private double log10NothingPrior;

    public AlleleFractionClustering(final List<ImmutablePair<double[], double[]>> tumorLodsAndCounts,
                                    final int callableSites, final M2FiltersArgumentCollection MTFAC) {
        Utils.validateArg(MTFAC.highConfidenceLod >= MTFAC.lowConfidenceLod, "High confidence threshold can't be smaller than low-confidence threshold.");
        final List<double[]> highConfidenceCounts = tumorLodsAndCounts.stream()
                .filter(pair -> MathUtils.arrayMax(pair.getLeft()) > MTFAC.highConfidenceLod)
                .map(pair -> pair.getRight())
                .collect(Collectors.toList());
        final List<double[]> lowConfidenceCounts = tumorLodsAndCounts.stream()
                                                .filter(pair -> MathUtils.arrayMax(pair.getLeft()) > MTFAC.lowConfidenceLod && MathUtils.arrayMax(pair.getLeft()) <= MTFAC.highConfidenceLod)
                                                .map(pair -> pair.getRight())
                                                .collect(Collectors.toList());
        highConfidenceDistribution = fitShape(highConfidenceCounts);
        log10HighConfidencePrior = FastMath.log10((double) highConfidenceCounts.size() / callableSites);
        lowConfidenceDistribution = fitShape(lowConfidenceCounts);
        log10LowConfidencePrior = FastMath.log10((double) lowConfidenceCounts.size() / callableSites);
        log10NothingPrior = FastMath.log10(1 - log10HighConfidencePrior - log10LowConfidencePrior);

        // Now do an E step with fractional assignments to nothing, low-confidence, high confidence
        final List<double[]> responsibilities = tumorLodsAndCounts.stream().map(pair -> {
            final double tumorLog10Odds = pair.getLeft()[0];
            final double refCount = pair.getRight()[0];
            final double altCount = pair.getRight()[1];

            final double lowConfidenceLog10OddsCorrection = SomaticLikelihoodsEngine.log10OddsCorrection(
                    lowConfidenceDistribution, Dirichlet.flat(2), new double[] {altCount, refCount});
            final double highConfidenceLog10OddsCorrection = SomaticLikelihoodsEngine.log10OddsCorrection(
                    highConfidenceDistribution, Dirichlet.flat(2), new double[] {altCount, refCount});

            final double[] unweightedLog10Responsibilities = new double[] {log10NothingPrior,
                    log10LowConfidencePrior + tumorLog10Odds + lowConfidenceLog10OddsCorrection,
                    log10HighConfidencePrior + tumorLog10Odds + highConfidenceLog10OddsCorrection};

            return MathUtils.normalizeFromLog10ToLinearSpace(unweightedLog10Responsibilities);
            }).collect(Collectors.toList());


        log10LowConfidencePrior = FastMath.log10(responsibilities.stream().mapToDouble(r -> r[1]).sum() / callableSites);
        log10HighConfidencePrior = FastMath.log10(responsibilities.stream().mapToDouble(r -> r[2]).sum() / callableSites);
        log10NothingPrior = FastMath.log10(responsibilities.stream().mapToDouble(r -> r[0]).sum() / callableSites);

        //todo write M step for beta shapes by using fit shape w/ responsibilities
        // TODO wrap this in an iteration

    }

    private BetaDistributionShape fitShape(final List<double[]> counts) {
        //TODO fill this in
        return null;
    }
}
