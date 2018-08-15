package org.broadinstitute.hellbender.tools.funcotator.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class FuncotationEngine {

    private static final Logger logger = LogManager.getLogger(FuncotationEngine.class);

    private FuncotationMetadata inputMetadata;
    private List<DataSourceFuncotationFactory> dataSourceFactories;

    public FuncotationEngine(final FuncotationMetadata metadata, final List<DataSourceFuncotationFactory> funcotationFactories) {
        inputMetadata = metadata;
        dataSourceFactories = funcotationFactories;
        dataSourceFactories.sort(DataSourceUtils::datasourceComparator);
    }

    public List<DataSourceFuncotationFactory> getDataSourceFactories() {
        return dataSourceFactories;
    }

    /**
     * Creates a {@link FuncotationMap} for the given {@code variantContext}.
     *
     * @param variantContext   {@link VariantContext} to annotate.  Never {@code null}.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variantContext}.  Never {@code null}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variantContext}.  Never {@code null}.
     * @return an instance of FuncotationMap that maps transcript IDs to lists of funcotations for the given variantContext context.
     */
    public FuncotationMap createFuncotationMapForVariant(final VariantContext variantContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        Utils.nonNull(variantContext);
        Utils.nonNull(referenceContext);
        Utils.nonNull(featureContext);

        //==============================================================================================================
        // First create only the transcript (Gencode) funcotations:

        if (retrieveGencodeFuncotationFactoryStream().count() > 1) {
            logger.warn("Attempting to annotate with more than one GENCODE datasource.  If these have overlapping transcript IDs, errors may occur.");
        }

        final List<GencodeFuncotation> transcriptFuncotations = retrieveGencodeFuncotationFactoryStream()
                .map(gf -> gf.createFuncotations(variantContext, referenceContext, featureContext))
                .flatMap(List::stream)
                .map(gf -> (GencodeFuncotation) gf).collect(Collectors.toList());

        //==============================================================================================================
        // Create the funcotations for non-Gencode data sources:

        // Create a place to keep our funcotations:
        final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(transcriptFuncotations);

        // Perform the rest of the annotation.  Note that this code manually excludes the Gencode Funcotations.
        for (final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

            // Note that this guarantees that we do not add GencodeFuncotations a second time.
            if (!funcotationFactory.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE)) {
                final List<String> txIds = funcotationMap.getTranscriptList();

                for (final String txId: txIds) {
                    funcotationMap.add(txId, funcotationFactory.createFuncotations(variantContext, referenceContext,
                            featureContext, funcotationMap.getGencodeFuncotations(txId)));
                }
            }
        }

        //==============================================================================================================
        // Create the funcotations for the input and add to all txID mappings.

        final List<String> txIds = funcotationMap.getTranscriptList();

        for (final String txId: txIds) {
            funcotationMap.add(txId, FuncotatorUtils.createFuncotations(variantContext, inputMetadata, FuncotatorConstants.DATASOURCE_NAME_FOR_INPUT_VCFS));
        }

        return funcotationMap;
    }

    /**
     * Shutdown the engine.  Closes all datasource factories.
     */
    public void close() {
        for ( final DataSourceFuncotationFactory factory : dataSourceFactories ) {
            if ( factory != null ) {
                factory.close();
            }
        }
    }

    private Stream<DataSourceFuncotationFactory> retrieveGencodeFuncotationFactoryStream() {
        return dataSourceFactories.stream()
                .filter(f -> f.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE));
    }
}
