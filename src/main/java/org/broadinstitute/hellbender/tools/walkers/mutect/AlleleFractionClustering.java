package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.*;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public class AlleleFractionClustering {
    private BetaDistributionShape highConfidenceDistribution = new BetaDistributionShape(1,1);
    private BetaDistributionShape lowConfidenceDistribution = new BetaDistributionShape(1,1);
    private final int callableSites;
    private double log10HighConfidencePrior;
    private double log10LowConfidencePrior;
    private double log10NothingPrior;

    public AlleleFractionClustering(final List<ImmutablePair<double[], double[]>> tumorLodsAndCounts,
                                    final int callableSites, final M2FiltersArgumentCollection MTFAC) {
        Utils.validateArg(MTFAC.highConfidenceLod >= MTFAC.lowConfidenceLod, "High confidence threshold can't be smaller than low-confidence threshold.");
        this.callableSites = callableSites;
        final List<double[]> allCounts = tumorLodsAndCounts.stream().map(pair -> pair.getRight()).collect(Collectors.toList());

        List<double[]> responsibilities = tumorLodsAndCounts.stream()
                .mapToDouble(pair -> MathUtils.arrayMax(pair.getLeft()))
                .mapToObj(lod -> lod < MTFAC.lowConfidenceLod ? new double[] {1, 0, 0} :
                        (lod < MTFAC.highConfidenceLod ? new double[] {0,1,0} : new double[] {0,0,1}))
                .collect(Collectors.toList());
        updatePriors(responsibilities);
        fitShape(allCounts, responsibilities);

        //TODO: real measure of convergence or at least extract magic constant!
        for (int n = 0; n < 3; n++) {
            // E step with fractional assignments to nothing, low-confidence, high confidence
            responsibilities = tumorLodsAndCounts.stream().map(pair -> {
                final double tumorLog10Odds = pair.getLeft()[0];
                final double refCount = pair.getRight()[0];
                final double altCount = pair.getRight()[1];

                final double lowConfidenceLog10OddsCorrection = SomaticLikelihoodsEngine.log10OddsCorrection(
                        lowConfidenceDistribution, Dirichlet.flat(2), new double[]{altCount, refCount});
                final double highConfidenceLog10OddsCorrection = SomaticLikelihoodsEngine.log10OddsCorrection(
                        highConfidenceDistribution, Dirichlet.flat(2), new double[]{altCount, refCount});

                final double[] unweightedLog10Responsibilities = new double[]{log10NothingPrior,
                        log10LowConfidencePrior + tumorLog10Odds + lowConfidenceLog10OddsCorrection,
                        log10HighConfidencePrior + tumorLog10Odds + highConfidenceLog10OddsCorrection};

                return MathUtils.normalizeFromLog10ToLinearSpace(unweightedLog10Responsibilities);
            }).collect(Collectors.toList());

            updatePriors(responsibilities);
            fitShape(allCounts, responsibilities);
        }
    }

    private void updatePriors(final List<double[]> responsibilities) {

        final double lowConfCount = responsibilities.stream().mapToDouble(r -> r[1]).sum();
        log10LowConfidencePrior = FastMath.log10(lowConfCount / callableSites);
        final double highConfCount = responsibilities.stream().mapToDouble(r -> r[2]).sum();
        log10HighConfidencePrior = FastMath.log10(highConfCount / callableSites);
        log10NothingPrior = FastMath.log10((callableSites - lowConfCount - highConfCount)/ callableSites);
    }

    private void fitShape(final List<double[]> counts, final List<double[]> responsibilities) {
        Utils.validateArg(counts.size() == responsibilities.size(), "must have one responsibility per count");

        final MutableDouble xLow = new MutableDouble(0);
        final MutableDouble yLow = new MutableDouble(0);
        final MutableDouble zLow = new MutableDouble(0);
        final MutableDouble xHigh = new MutableDouble(0);
        final MutableDouble yHigh = new MutableDouble(0);
        final MutableDouble zHigh = new MutableDouble(0);


        final MutableDouble altCountLow = new MutableDouble(0);
        final MutableDouble altCountHigh = new MutableDouble(0);
        final MutableDouble refCountLow = new MutableDouble(0);
        final MutableDouble refCountHigh = new MutableDouble(0);
        for (int n = 0; n < counts.size(); n++) {
            final double altCount = counts.get(n)[0];
            final double refCount = counts.get(n)[1];
            final double rLow = responsibilities.get(n)[1];
            final double rHigh = responsibilities.get(n)[1];

            altCountLow.add(rLow * altCount);
            altCountHigh.add(rHigh * altCount);
            refCountLow.add(rLow * refCount);
            refCountHigh.add(rHigh * refCount);

            final double weightedAlpha = rLow * lowConfidenceDistribution.getAlpha() +
                    rHigh * highConfidenceDistribution.getAlpha() + altCount;
            final double weightedBeta = rLow * lowConfidenceDistribution.getBeta() +
                    rHigh * highConfidenceDistribution.getBeta() + refCount;

            final BetaDistributionShape meanFieldPosterior = new BetaDistributionShape(weightedAlpha, weightedBeta);
            final double[] meanLogs = meanFieldPosterior.meanLogs();

            xLow.add(rLow * meanLogs[0]);
            yLow.add(rLow * meanLogs[1]);
            xHigh.add(rHigh * meanLogs[0]);
            yHigh.add(rHigh * meanLogs[1]);
            zLow.add(rLow);
            zHigh.add(rHigh);
        }

        final double lowMean = (0.5 + altCountLow.doubleValue()) / (altCountLow.doubleValue()+ refCountLow.doubleValue() + 1);
        final double highMean = (0.5 + altCountHigh.doubleValue()) / (altCountHigh.doubleValue()+ refCountHigh.doubleValue() + 1);

        // We get the mean a/(a+b) of our beta from the ratio of weighted alt counts to weighted total counts
        // now we optimize with respect to the overall scale, a.  That is, we optimize wrt a subject to b = (1-mean)*a/mean
        final Function<Double, Double> lowObjective = a -> objective(a, lowMean, xLow.doubleValue(), yLow.doubleValue(), zLow.doubleValue());
        final Function<Double, Double> highObjective = a -> objective(a, highMean, xHigh.doubleValue(), yHigh.doubleValue(), zHigh.doubleValue());

        final double alphaLow = OptimizationUtils.argmax(lowObjective, 0, 1000, lowConfidenceDistribution.getAlpha());
        final double alphaHigh = OptimizationUtils.argmax(highObjective, 0, 1000, highConfidenceDistribution.getAlpha());
        lowConfidenceDistribution = new BetaDistributionShape(alphaLow, getBeta(alphaLow, lowMean));
        highConfidenceDistribution = new BetaDistributionShape(alphaHigh, getBeta(alphaHigh, lowMean));
    }

    private double objective(final double alpha, final double mean, final double x, final double y, final double z) {
        final double beta = getBeta(alpha, mean);
        return new BetaDistributionShape(alpha,beta).log10Normalization() * z + x * alpha + y * beta;
    }

    private double getBeta(final double alpha, final double mean) {
        return (1 - mean) * alpha / mean;
    }
}
