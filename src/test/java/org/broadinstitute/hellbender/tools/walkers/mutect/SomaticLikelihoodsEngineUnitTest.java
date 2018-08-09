package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by David Benjamin on 3/9/17.
 */
public class SomaticLikelihoodsEngineUnitTest extends GATKBaseTest {
    @Test
    public void testAlleleFractionsPosterior() {

        //likelihoods completely favor allele 0 over allele 1 for every read, so
        // we should get no counts for allele 1
        final Dirichlet prior1 = Dirichlet.of(1, 1);
        final RealMatrix mat1 = new Array2DRowRealMatrix(new double[][] {{0, 0, 0, 0}, {-10, -10, -10, -10}});
        final Dirichlet posterior1 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat1, prior1);
        final double[] expectedCounts1 = new double[] {4, 0};

        final Dirichlet expectedPosterior1 = prior1.addCounts(expectedCounts1);
        Assert.assertEquals(posterior1.distance1(expectedPosterior1),0, 1.0e-6);

        //prior is extremely strong and outweighs ambiguous likelihoods
        final Dirichlet prior2 = Dirichlet.of(1e8, 1);
        final RealMatrix mat2 = new Array2DRowRealMatrix(new double[][] {{0, 0, 0, 0}, {0, 0, 0, 0}});
        final Dirichlet posterior2 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat2, prior2);
        final double[] expectedCounts2 = new double[] {4, 0};

        final Dirichlet expectedPosterior2 = prior2.addCounts(expectedCounts2);
        Assert.assertEquals(posterior2.distance1(expectedPosterior2),0, 1.0e-6);

        //prior is extremely weak and likelihoods speak for themselves
        final Dirichlet prior3 = Dirichlet.of(1e-6, 1e-6);
        final RealMatrix mat3 = new Array2DRowRealMatrix(new double[][] {{0, 0, 0, -10}, {-10, -10, -10, 0}});
        final Dirichlet posterior3 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat3, prior3);
        final double[] expectedCounts3 = new double[] {3, 1};

        final Dirichlet expectedPosterior3 = prior3.addCounts(expectedCounts3);
        Assert.assertEquals(posterior3.distance1(expectedPosterior3),0, 1.0e-6);

        // test convergence i.e. posterior = prior + effective counts
        final Dirichlet prior4 = Dirichlet.of(0.2, 1.7);
        final RealMatrix mat4 = new Array2DRowRealMatrix(new double[][] {{0.1, 5.2, 0.5, 0.2}, {2.6, 0.6, 0.5, 0.4}});
        final Dirichlet posterior4 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat4, prior4);
        final double[] effectiveCounts = SomaticLikelihoodsEngine.getEffectiveCounts(mat4, posterior4);
        Assert.assertEquals(prior4.addCounts(effectiveCounts).distance1(posterior4), 0, 1.0e-3);
    }

    @Test
    public void testEvidence() {
        // one exact limit for the evidence is when the likelihoods of each read are so peaked (i.e. the most likely allele
        // of each read is much likelier than all other alleles) that the sum over latent read-to-allele assignments
        // (that is, over the indicator z in the notes) is dominated by the max-likelihood allele configuration
        // and thus the evidence reduces to exactly integrating out the Dirichlet allele fractions

        final Dirichlet prior = Dirichlet.of(1,2);
        final RealMatrix log10Likelihoods = new Array2DRowRealMatrix(new double[][] {{0.1, 4.0, 3.0, -10}, {-12, -9, -5.0, 0.5}});
        final double calculatedLog10Evidence = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods, prior);
        final double[] maxLikelihoodCounts = new double[] {3, 1};
        final double expectedLog10Evidence = prior.log10Normalization() - prior.addCounts(maxLikelihoodCounts).log10Normalization()
                + new IndexRange(0,log10Likelihoods.getColumnDimension()).sum(read -> log10Likelihoods.getColumnVector(read).getMaxValue());
        Assert.assertEquals(calculatedLog10Evidence, expectedLog10Evidence, 1e-5);

        // when there's just one read we can calculate the likelihood exactly

        final Dirichlet prior2 = Dirichlet.of(1,2);
        final RealMatrix log10Likelihoods2 = new Array2DRowRealMatrix(new double[][] {{0.1}, {0.5}});
        final double calculatedLog10Evidence2 = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods2, prior2);
        final double[] delta0 = new double[] {1, 0};
        final double[] delta1 = new double[] {0, 1};
        final double expectedLog10Evidence2 = MathUtils.log10SumLog10(log10Likelihoods2.getEntry(0,0) +
                prior2.log10Normalization() - prior2.addCounts(delta0).log10Normalization(),
                log10Likelihoods2.getEntry(1,0) +
                prior2.log10Normalization() - prior2.addCounts(delta1).log10Normalization());
        Assert.assertEquals(calculatedLog10Evidence2, expectedLog10Evidence2, 0.05);


    }

}