package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Stores the results of the first pass of {@link FilterMutectCalls}, a purely online step in which each variant is
 * not "aware" of other variants, and learns various global properties necessary for a more refined second step.
 */
public class FilteringFirstPass {
    final List<FilterResult> filterResults = new ArrayList<>();
    final List<ImmutablePair<double[], double[]>> unfilteredTumorLodsAndCounts = new ArrayList<>();
    final Map<String, ImmutablePair<String, Integer>> filteredPhasedCalls = new HashMap<>();
    final Map<String, FilterStats> filterStats = new HashMap<>();
    AlleleFractionClustering afClustering = null;
    final String tumorSample;
    final int callableSites;    //TODO: emit this in M2 and grab from vcf just like tumor sample
    boolean readyForSecondPass = false;

    public FilteringFirstPass(final String tumorSample, final int callableSites) {
        this.tumorSample = tumorSample;
        this.callableSites = callableSites;
    }

    public boolean isReadyForSecondPass() { return readyForSecondPass; }

    public FilterStats getFilterStats(final String filterName){
        Utils.validateArg(filterStats.containsKey(filterName), "invalid filter name: " + filterName);
        return filterStats.get(filterName);
    }

    public boolean isOnFilteredHaplotype(final VariantContext vc, final int maxDistance) {

        final Genotype tumorGenotype = vc.getGenotype(tumorSample);

        if (!hasPhaseInfo(tumorGenotype)) {
            return false;
        }

        final String pgt = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "");
        final String pid = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "");
        final int position = vc.getStart();

        final Pair<String, Integer> filteredCall = filteredPhasedCalls.get(pid);
        if (filteredCall == null) {
            return false;
        }

        // Check that vc occurs on the filtered haplotype
        return filteredCall.getLeft().equals(pgt) && Math.abs(filteredCall.getRight() - position) <= maxDistance;
    }

    public void add(final FilterResult filterResult, final VariantContext vc) {
        filterResults.add(filterResult);
        final Genotype tumorGenotype = vc.getGenotype(tumorSample);

        final Set<String> appliedFilters = filterResult.getFilters();
        if (!appliedFilters.isEmpty() && hasPhaseInfo(tumorGenotype)) {
            final String pgt = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "");
            final String pid = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "");
            final int position = vc.getStart();
            filteredPhasedCalls.put(pid, new ImmutablePair<>(pgt, position));
        }

        // if a variant has no artifact filter applied (it could have a TLOD filter) we use it for the AF clustering model
        if (appliedFilters.isEmpty() || (appliedFilters.size() == 1 && appliedFilters.contains(GATKVCFConstants.TUMOR_LOD_FILTER_NAME))) {
            final double[] tumorLods = Mutect2FilteringEngine.getDoubleArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final double[] tumorADs = Arrays.stream(tumorGenotype.getAD()).mapToDouble(n->n).toArray();
            unfilteredTumorLodsAndCounts.add(new ImmutablePair<>(tumorLods, tumorADs));
        }
    }

    public void learnModelForSecondPass(final M2FiltersArgumentCollection MTFAC) {
        final double[] readOrientationPosteriors = getFilterResults().stream()
                .filter(r -> r.getFilters().isEmpty())
                .mapToDouble(r -> r.getReadOrientationPosterior())
                .toArray();

        final FilterStats readOrientationFilterStats = calculateThresholdForReadOrientationFilter(readOrientationPosteriors, MTFAC.maxFalsePositiveRate);
        filterStats.put(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, readOrientationFilterStats);

        afClustering = new AlleleFractionClustering(unfilteredTumorLodsAndCounts, callableSites, MTFAC);
        readyForSecondPass = true;
    }

    /**
     *
     * Compute the filtering threshold that ensures that the false positive rate among the resulting pass variants
     * will not exceed the requested false positive rate
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param requestedFPR We set the filtering threshold such that the FPR doesn't exceed this value
     * @return
     */
    public static FilterStats calculateThresholdForReadOrientationFilter(final double[] posteriors, final double requestedFPR){
        ParamUtils.isPositiveOrZero(requestedFPR, "requested FPR must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Arrays.sort(posteriors);

        final int numPassingVariants = posteriors.length;
        double cumulativeExpectedFPs = 0.0;

        for (int i = 0; i < numPassingVariants; i++){
            final double posterior = posteriors[i];

            // One can show that the cumulative error rate is monotonically increasing in i
            final double expectedFPR = (cumulativeExpectedFPs + posterior) / (i + 1);
            if (expectedFPR > requestedFPR){
                return i > 0 ?
                        new FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, posteriors[i-1],
                                cumulativeExpectedFPs, i-1, cumulativeExpectedFPs/i, requestedFPR) :
                        new FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, thresholdForFilteringAll,
                                0.0, 0, 0.0, requestedFPR);
            }

            cumulativeExpectedFPs += posterior;
        }

        // If the expected FP rate never exceeded the max tolerable value, then we can let everything pass
        return new FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, thresholdForFilteringNone,
                cumulativeExpectedFPs, numPassingVariants, cumulativeExpectedFPs/numPassingVariants, requestedFPR);
    }

    public double getSomaticProbability(final double tumorLog10Odds, final double refCount, final double altCount) {
        Utils.validateArg(readyForSecondPass, "somatic probability should only be called after learning from first pass.");
        return afClustering.getSomaticProbability(tumorLog10Odds, refCount, altCount);
    }

    public static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

    public List<FilterResult> getFilterResults() {
        return filterResults;
    }

    public static class FilterStats {
        private final String filterName;
        private final double threshold;
        private final double expectedNumFPs;
        private final int numPassingVariants;
        private final double expectedFPR;
        private final double requestedFPR;

        public FilterStats(final String filterName, final double threshold, final double expectedNumFPs,
                           final int numPassingVariants, final double expectedFPR, final double requestedFPR){
            this.filterName = filterName;
            this.threshold = threshold;
            this.expectedNumFPs = expectedNumFPs;
            this.numPassingVariants = numPassingVariants;
            this.expectedFPR = expectedFPR;
            this.requestedFPR = requestedFPR;
        }

        public String getFilterName() { return filterName; }

        public double getExpectedNumFPs() { return expectedNumFPs; }

        public int getNumPassingVariants() { return numPassingVariants; }

        public double getThreshold() { return threshold; }

        public double getExpectedFPR() { return expectedFPR; }

        public double getRequestedFPR() { return requestedFPR; }

    }

    private enum M2FilterStatsTableColumn {
        FILTER_NAME("filter_name"),
        THRESHOLD("threshold"),
        EXPECTED_FALSE_POSITIVES("expected_fps"),
        EXPECTED_FALSE_POSITIVE_RATE("expected_fpr"),
        REQUESTED_FALSE_POSITIVE_RATE("requested_fpr"),
        NUM_PASSING_VARIANTS("num_passing_variants");

        private String columnName;

        M2FilterStatsTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static class Mutect2FilterStatsWriter extends TableWriter<FilterStats> {
        private Mutect2FilterStatsWriter(final File output) throws IOException {
            super(output, M2FilterStatsTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final FilterStats stats, final DataLine dataLine) {
            dataLine.set(M2FilterStatsTableColumn.FILTER_NAME.toString(), stats.getFilterName())
                    .set(M2FilterStatsTableColumn.THRESHOLD.toString(), stats.getThreshold())
                    .set(M2FilterStatsTableColumn.EXPECTED_FALSE_POSITIVES.toString(), stats.getExpectedNumFPs())
                    .set(M2FilterStatsTableColumn.EXPECTED_FALSE_POSITIVE_RATE.toString(), stats.getExpectedFPR())
                    .set(M2FilterStatsTableColumn.REQUESTED_FALSE_POSITIVE_RATE.toString(), stats.getRequestedFPR())
                    .set(M2FilterStatsTableColumn.NUM_PASSING_VARIANTS.toString(), stats.getNumPassingVariants());
        }
    }

    public void writeM2FilterSummary(final File outputTable) {
        try (Mutect2FilterStatsWriter writer = new Mutect2FilterStatsWriter(outputTable)) {
            writer.writeAllRecords(filterStats.values());
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable), e);
        }
    }
}
