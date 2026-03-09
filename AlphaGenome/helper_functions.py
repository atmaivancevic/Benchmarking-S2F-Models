# Standard library imports
import io
import itertools
import os

# Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import clear_output

# AlphaGenome imports
import alphagenome
from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import track_data
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components

###############################################################################################################

def plot_baseline(variant_name, start=None, end=None):
    """
    Plot baseline (reference) genome state for all selected assays.
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        start: Custom start coordinate (default: variant_start - 5kb)
        end: Custom end coordinate (default: variant_end + 5kb)
    """
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    variant_start = variant_row['POS']
    variant_end = variant_row['POS'] + len(variant_row['REF'])
    
    # Use custom coordinates if provided, otherwise default to ±5kb
    if start is None:
        start = variant_start - 5000
    if end is None:
        end = variant_end + 5000
    
    print(f"Creating baseline view for {variant_name}")
    
    # Define interval
    interval = genome.Interval(chromosome=chrom, start=start, end=end, strand='+')
    
    # Convert PRIMARY_ASSAYS to OutputType set
    requested_outputs = {getattr(dna_client.OutputType, assay) for assay in PRIMARY_ASSAYS}
    
    # Get predictions
    print(f"Making baseline predictions for {PRIMARY_ASSAYS}...")
    baseline = dna_model.predict_interval(
        interval=interval.resize(2**20),
        requested_outputs=requested_outputs,
        ontology_terms=[PRIMARY_CELL_LINE_ID],
    )
    
    # Get transcripts and create markers
    transcripts = longest_transcript_extractor.extract(interval)
    markers = [genome.Variant(chromosome=chrom, position=pos, reference_bases='N', 
                              alternate_bases='N', name=name) 
               for pos, name in [(variant_start, 'Variant_Start'), (variant_end, 'Variant_End')]]
    
    # Map assay types to baseline attributes
    assay_map = {
        'RNA_SEQ': 'rna_seq', 'CHIP_HISTONE': 'chip_histone', 'CHIP_TF': 'chip_tf',
        'DNASE': 'dnase', 'ATAC': 'atac', 'CAGE': 'cage', 'SPLICE_SITES': 'splice_sites',
        'SPLICE_SITE_USAGE': 'splice_site_usage', 'SPLICE_JUNCTIONS': 'splice_junctions',
        'CONTACT_MAPS': 'contact_maps', 'PROCAP': 'procap',
    }
    
    # Build tracks
    tracks = [plot_components.TranscriptAnnotation(transcripts)]
    for assay in PRIMARY_ASSAYS:
        if assay in assay_map:
            attr_name = assay_map[assay]
            if hasattr(baseline, attr_name):
                try:
                    data = getattr(baseline, attr_name)
                    ylabel_template = '{biosample_name}\n{name}' if hasattr(data, 'metadata') and 'biosample_name' in data.metadata.columns else '{name}'
                    tracks.append(plot_components.Tracks(tdata=data, ylabel_template=ylabel_template, filled=True))
                except (AttributeError, ValueError, KeyError):
                    print(f" Skipping {assay}: not compatible with Tracks")
    
    plot_components.plot(
        tracks,
        annotations=[plot_components.VariantAnnotation(markers)],
        interval=interval,
        title=f'Baseline Genome State at {variant_name} in {PRIMARY_CELL_LINE}',
    )
    
    # Rasterize fill collections to keep PDF small and Illustrator-compatible
    for ax in plt.gcf().get_axes():
        for collection in ax.collections:
            collection.set_rasterized(True)
    
    # Save to current directory
    plt.gcf().set_size_inches(18, 12)
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_baseline.pdf'
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    
###############################################################################################################

def plot_deletion(variant_name, start=None, end=None):
    """
    Plot genome state WITH deletion (mimics CRISPR knockout) for all selected assays.
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        start: Custom start coordinate (default: variant_start - 5kb)
        end: Custom end coordinate (default: variant_end + 5kb)
    """
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    variant_start = variant_row['POS']
    variant_end = variant_row['POS'] + len(variant_row['REF'])
    ref_seq = variant_row['REF']
    
    # Use custom coordinates if provided, otherwise default to ±5kb
    if start is None:
        start = variant_start - 5000
    if end is None:
        end = variant_end + 5000
    
    print(f"Predicting genome state WITH deletion for {variant_name}")
    
    # Create deletion variant
    deletion_variant = genome.Variant(
        chromosome=chrom,
        position=variant_start,
        reference_bases=ref_seq,
        alternate_bases='',  # Empty = deletion
        name=f'{variant_name}_Deletion'
    )
    
    # Define interval
    interval = genome.Interval(chromosome=chrom, start=start, end=end, strand='+')
    
    # Convert PRIMARY_ASSAYS to OutputType set
    requested_outputs = {getattr(dna_client.OutputType, assay) for assay in PRIMARY_ASSAYS}
    
    # Make predictions with the deletion
    print(f"Making predictions with deletion for {PRIMARY_ASSAYS}...")
    deletion_output = dna_model.predict_variant(
        interval=interval.resize(2**20),
        variant=deletion_variant,
        requested_outputs=requested_outputs,
        ontology_terms=[PRIMARY_CELL_LINE_ID],
    )
    
    # Get transcripts
    transcripts = longest_transcript_extractor.extract(interval)
    
    # Map assay types to baseline attributes
    assay_map = {
        'RNA_SEQ': 'rna_seq', 'CHIP_HISTONE': 'chip_histone', 'CHIP_TF': 'chip_tf',
        'DNASE': 'dnase', 'ATAC': 'atac', 'CAGE': 'cage', 'SPLICE_SITES': 'splice_sites',
        'SPLICE_SITE_USAGE': 'splice_site_usage', 'SPLICE_JUNCTIONS': 'splice_junctions',
        'CONTACT_MAPS': 'contact_maps', 'PROCAP': 'procap',
    }
    
    # Build tracks (showing ALTERNATE allele, i.e the deletion state)
    tracks = [plot_components.TranscriptAnnotation(transcripts)]
    for assay in PRIMARY_ASSAYS:
        if assay in assay_map:
            attr_name = assay_map[assay]
            if hasattr(deletion_output.alternate, attr_name):
                try:
                    # Use the alternate (deletion) state
                    alt_data = getattr(deletion_output.alternate, attr_name)
                    
                    ylabel_template = '{biosample_name}\n{name}' if hasattr(alt_data, 'metadata') and 'biosample_name' in alt_data.metadata.columns else '{name}'
                    tracks.append(plot_components.Tracks(tdata=alt_data, ylabel_template=ylabel_template, filled=True))
                except (AttributeError, ValueError, KeyError):
                    print(f" Skipping {assay}: not compatible with Tracks")
    
    plot_components.plot(
        tracks,
        annotations=[plot_components.VariantAnnotation([deletion_variant])],
        interval=interval,
        title=f'Genome State WITH {variant_name} Deletion in {PRIMARY_CELL_LINE}',
    )
    
    # Save to current directory
    plt.gcf().set_size_inches(18, 14)
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_with_deletion.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f" Saved: {filename}")
    plt.show()
    
###############################################################################################################    

def plot_difference(variant_name, start=None, end=None):
    """
    Plot the difference from baseline for all selected assays, e.g. the effect of the deletion.
    
    Positive values (i.e. tracks facing upward) indicate GAIN of signal.
    Negative values (i.e. tracks facing downward) indicate LOSS of signal.
    
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        start: Custom start coordinate (default: variant_start - 5kb)
        end: Custom end coordinate (default: variant_end + 5kb)
    """
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    variant_start = variant_row['POS']
    variant_end = variant_row['POS'] + len(variant_row['REF'])
    ref_seq = variant_row['REF']
    
    # Use custom coordinates if provided, otherwise default to ±5kb
    if start is None:
        start = variant_start - 5000
    if end is None:
        end = variant_end + 5000
    
    print(f"Predicting deletion effects for {variant_name}")
    
    # Create deletion variant
    deletion_variant = genome.Variant(
        chromosome=chrom,
        position=variant_start,
        reference_bases=ref_seq,
        alternate_bases='',  # Empty = deletion
        name=f'{variant_name}_Deletion'
    )
    
    # Define interval
    interval = genome.Interval(chromosome=chrom, start=start, end=end, strand='+')
    
    # Convert PRIMARY_ASSAYS to OutputType set
    requested_outputs = {getattr(dna_client.OutputType, assay) for assay in PRIMARY_ASSAYS}
    
    # Make predictions with the deletion
    print(f"Making deletion predictions for {PRIMARY_ASSAYS}...")
    deletion_output = dna_model.predict_variant(
        interval=interval.resize(2**20),
        variant=deletion_variant,
        requested_outputs=requested_outputs,
        ontology_terms=[PRIMARY_CELL_LINE_ID],
    )
    
    # Get transcripts
    transcripts = longest_transcript_extractor.extract(interval)
    
    # Map assay types to baseline attributes
    assay_map = {
        'RNA_SEQ': 'rna_seq', 'CHIP_HISTONE': 'chip_histone', 'CHIP_TF': 'chip_tf',
        'DNASE': 'dnase', 'ATAC': 'atac', 'CAGE': 'cage', 'SPLICE_SITES': 'splice_sites',
        'SPLICE_SITE_USAGE': 'splice_site_usage', 'SPLICE_JUNCTIONS': 'splice_junctions',
        'CONTACT_MAPS': 'contact_maps', 'PROCAP': 'procap',
    }
    
    # Build tracks (showing DIFFERENCE: alternate - reference)
    tracks = [plot_components.TranscriptAnnotation(transcripts)]
    for assay in PRIMARY_ASSAYS:
        if assay in assay_map:
            attr_name = assay_map[assay]
            if hasattr(deletion_output.alternate, attr_name) and hasattr(deletion_output.reference, attr_name):
                try:
                    # Calculate difference: alternate - reference
                    alt_data = getattr(deletion_output.alternate, attr_name)
                    ref_data = getattr(deletion_output.reference, attr_name)
                    diff_data = alt_data - ref_data
                    
                    ylabel_template = '{biosample_name}\n{name}' if hasattr(diff_data, 'metadata') and 'biosample_name' in diff_data.metadata.columns else '{name}'
                    tracks.append(plot_components.Tracks(tdata=diff_data, ylabel_template=ylabel_template, filled=True))
                except (AttributeError, ValueError, KeyError):
                    print(f"  Skipping {assay}: not compatible with Tracks")
    
    plot_components.plot(
        tracks,
        annotations=[plot_components.VariantAnnotation([deletion_variant])],
        interval=interval,
        title=f'Effect of {variant_name} Deletion in {PRIMARY_CELL_LINE} (Change from Baseline)',
    )
    
    # Save to current directory
    plt.gcf().set_size_inches(18, 14)
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_deletion_effect.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f" Saved: {filename}")
    plt.show()    
    
###############################################################################################################    
    
def plot_difference_corrected_rna(variant_name, start=None, end=None):
    """
    Plot RNA-seq deletion effects with frameshift correction.
    
    Deletions cause a coordinate shift in predictions downstream of the variant.
    This creates visualization errors when the variant is not at the right edge of the plot.
    This function uses AlphaGenome's recommended padding approach to correct the alignment,
    but this method only really works for high-resolution tracks like RNA-seq.
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        start: Custom start coordinate (default: variant_start - 5kb)
        end: Custom end coordinate (default: variant_end + 5kb)
    """
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    variant_start = variant_row['POS']
    variant_end = variant_row['POS'] + len(variant_row['REF'])
    ref_seq = variant_row['REF']
    
    # Use custom coordinates if provided, otherwise default to ±5kb
    if start is None:
        start = variant_start - 5000
    if end is None:
        end = variant_end + 5000
    
    print(f"Creating RNA-seq difference plot for {variant_name} (with frameshift correction)")
    
    # Create deletion variant
    deletion_variant = genome.Variant(
        chromosome=chrom,
        position=variant_start,
        reference_bases=ref_seq,
        alternate_bases='',
        name=f'{variant_name}_Deletion'
    )
    
    # Define interval
    interval = genome.Interval(chromosome=chrom, start=start, end=end, strand='+')
    
    # Make predictions for RNA-seq only
    print(f"Making RNA-seq predictions...")
    predictions = dna_model.predict_variant(
        interval=interval.resize(2**20),
        variant=deletion_variant,
        requested_outputs={dna_client.OutputType.RNA_SEQ},
        ontology_terms=[PRIMARY_CELL_LINE_ID],
    )
    
    # Get transcripts
    transcripts = longest_transcript_extractor.extract(interval)
    
    # Apply AlphaGenome's padding correction for RNA-seq
    deletion_length = len(deletion_variant.reference_bases)
    
    aligned_alternate = track_data.TrackData(
        np.pad(
            predictions.alternate.rna_seq.values,
            ((deletion_length, 0), (0, 0)),
        )[:-deletion_length, :],
        predictions.alternate.rna_seq.metadata,
        interval=predictions.alternate.rna_seq.interval,
    )
    
    # Calculate the DIFFERENCE: aligned alternate - reference
    diff_data = aligned_alternate - predictions.reference.rna_seq
    
    # Create tracks with the difference
    tracks = [
        plot_components.TranscriptAnnotation(transcripts),
        plot_components.Tracks(
            tdata=diff_data,
            ylabel_template='{biosample_name}\n{name}',
            filled=True
        )
    ]
    
    plot_components.plot(
        tracks,
        annotations=[plot_components.VariantAnnotation([deletion_variant], alpha=0.8)],
        interval=interval,
        title=f'RNA-seq: Effect of {variant_name} Deletion in {PRIMARY_CELL_LINE} (Change from Baseline)',
    )
    
    # Save
    plt.gcf().set_size_inches(18, 8)
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_rna_seq_difference.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f" Saved: {filename}")
    plt.show()

###############################################################################################################    
    
def plot_overlaid_corrected_rna(variant_name, start=None, end=None):
    """
    Plot reference vs deletion RNA-seq with frameshift correction (overlaid view).
    
    Same as plot_difference_corrected_rna_only, but with an overlaid view instead. 
    Shows reference (gray) overlaid with deletion (red).
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        start: Custom start coordinate (default: variant_start - 5kb)
        end: Custom end coordinate (default: variant_end + 5kb)
    """
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    variant_start = variant_row['POS']
    variant_end = variant_row['POS'] + len(variant_row['REF'])
    ref_seq = variant_row['REF']
    
    # Use custom coordinates if provided, otherwise default to ±5kb
    if start is None:
        start = variant_start - 5000
    if end is None:
        end = variant_end + 5000
    
    print(f"Creating RNA-seq overlay for {variant_name} (with frameshift correction)")
    
    # Create deletion variant
    deletion_variant = genome.Variant(
        chromosome=chrom,
        position=variant_start,
        reference_bases=ref_seq,
        alternate_bases='',
        name=f'{variant_name}_Deletion'
    )
    
    # Define interval
    interval = genome.Interval(chromosome=chrom, start=start, end=end, strand='+')
    
    # Make predictions for RNA-seq only
    print(f"Making RNA-seq predictions...")
    predictions = dna_model.predict_variant(
        interval=interval.resize(2**20),
        variant=deletion_variant,
        requested_outputs={dna_client.OutputType.RNA_SEQ},
        ontology_terms=[PRIMARY_CELL_LINE_ID],
    )
    
    # Get transcripts
    transcripts = longest_transcript_extractor.extract(interval)
    
    # Apply AlphaGenome's padding correction for RNA-seq
    deletion_length = len(deletion_variant.reference_bases)
    
    aligned_alternate = track_data.TrackData(
        np.pad(
            predictions.alternate.rna_seq.values,
            ((deletion_length, 0), (0, 0)),
        )[:-deletion_length, :],
        predictions.alternate.rna_seq.metadata,
        interval=predictions.alternate.rna_seq.interval,
    )
    
    # Create overlaid tracks
    tracks = [
        plot_components.TranscriptAnnotation(transcripts),
        plot_components.OverlaidTracks(
            {
                'Reference': predictions.reference.rna_seq,
                'Deletion': aligned_alternate,
            },
            colors={'Reference': 'dimgrey', 'Deletion': 'red'},
        )
    ]
    
    plot_components.plot(
        tracks,
        annotations=[plot_components.VariantAnnotation([deletion_variant], alpha=0.8)],
        interval=interval,
        title=f'RNA-seq: Reference vs Deletion for {variant_name} in {PRIMARY_CELL_LINE}',
    )
    
    # Save
    plt.gcf().set_size_inches(18, 8)
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_rna_seq_overlay.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f" Saved: {filename}")
    plt.show()   
    
###############################################################################################################    

def score_deletion(variant_name, gene_name=None, assay=None):
    """
    Score deletion effects for a variant using AlphaGenome variant scoring.
    Set the alternate allele as an empty string.
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        gene_name: Optional gene(s) to focus on (e.g., 'ATG12' or 'ATG12,AP3S1')
        assay: Optional specific assay to score (e.g., 'RNA_SEQ', 'ATAC', 'CHIP_HISTONE')
    
    Returns:
        DataFrame with deletion effects
    """
    import pandas as pd
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return None
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    pos = variant_row['POS']
    ref = variant_row['REF']
    
    print(f"Scoring deletion effects for {variant_name}")
    print(f"Location: {chrom}:{pos:,}-{pos + len(ref):,}")
    if assay:
        print(f"Assay: {assay}")
    if gene_name:
        print(f"Genes: {gene_name}")
    
    # Create deletion variant
    deletion_variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases='',  # Empty = deletion
        name=f'{variant_name}_Deletion'
    )
    
    # Map assay names to scorer configuration
    # Based on recommended scorers: https://www.alphagenomedocs.com/api/models.html#variant-scorers
    assay_config = {
        'RNA_SEQ': ('GeneMaskLFCScorer', None),
        'CAGE': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'PROCAP': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'ATAC': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'DNASE': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'CHIP_TF': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'CHIP_HISTONE': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'CONTACT_MAPS': ('ContactMapScorer', None),
    }
    
    # Select scorers
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    
    if assay and assay in assay_config:
        scorer_type, aggregation = assay_config[assay]
        selected_scorers = []
        
        for s in all_scorers.values():
            # Match scorer type and output
            if (type(s).__name__ == scorer_type and 
                hasattr(s, 'requested_output') and 
                str(s.requested_output).split('.')[-1] == assay):
                
                # Match aggregation if specified
                if aggregation is None:
                    selected_scorers.append(s)
                elif hasattr(s, 'aggregation_type') and str(s.aggregation_type).split('.')[-1] == aggregation:
                    selected_scorers.append(s)
        
        if not selected_scorers:
            selected_scorers = list(all_scorers.values())
    else:
        selected_scorers = list(all_scorers.values())
    
    # Score variant
    print("Running variant scoring...")
    variant_scores = dna_model.score_variant(
        interval=deletion_variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=deletion_variant,
        variant_scorers=selected_scorers,
        organism=dna_client.Organism.HOMO_SAPIENS
    )
    
    # Convert to DataFrame and filter
    effects_df = variant_scorers.tidy_scores(variant_scores, match_gene_strand=True)
    filtered_df = effects_df[effects_df['biosample_name'] == PRIMARY_CELL_LINE].copy()
    
    # Parse gene names (support both single gene and comma-separated list)
    if gene_name:
        gene_list = [g.strip() for g in gene_name.split(',')]
        filtered_df = filtered_df[filtered_df['gene_name'].isin(gene_list)]
    
    if assay:
        filtered_df = filtered_df[filtered_df['output_type'] == assay]
    
    # Remove variant_id and scored_interval columns to keep output clean
    cols_to_drop = ['variant_id', 'scored_interval']
    filtered_df = filtered_df.drop(columns=[col for col in cols_to_drop if col in filtered_df.columns])
    
    # Sort by quantile_score (most negative first = biggest decrease)
    filtered_df = filtered_df.sort_values('quantile_score', ascending=True)
    
    # Display results
    print(f"\n{'='*70}")
    print(f"RESULTS: {len(filtered_df)} effects in {PRIMARY_CELL_LINE}")
    print(f"{'='*70}")
    
    if gene_name and assay and not filtered_df.empty:
        for gene in gene_list:
            gene_df = filtered_df[filtered_df['gene_name'] == gene]
            if not gene_df.empty:
                print(f"\n{gene} - {assay}:")
                for _, row in gene_df.iterrows():
                    direction = "DECREASE" if row['quantile_score'] < -0.5 else ("INCREASE" if row['quantile_score'] > 0.5 else "Minimal effect")
                    print(f"  Raw: {row['raw_score']:.4f} | Quantile: {row['quantile_score']:.4f} | {direction}")
    elif not filtered_df.empty:
        summary = filtered_df.groupby('gene_name' if gene_name else 'output_type')['quantile_score'].agg(['count', 'mean']).round(3)
        print(summary)
    
    # Save results to CSV
    if not filtered_df.empty:
        # Build filename
        filename_parts = [variant_name, PRIMARY_CELL_LINE, 'deletion_scores']
        if assay:
            filename_parts.append(assay)
        if gene_name:
            # Use all gene names in filename (replace commas with underscores)
            filename_parts.append(gene_name.replace(',', '_').replace(' ', ''))
        filename = '_'.join(filename_parts) + '.csv'
        
        filtered_df.to_csv(filename, index=False)
        print(f"\n Saved results to: {filename}")
    
    return filtered_df

###############################################################################################################    

def score_deletion_alt(variant_name, gene_name=None, assay=None):
    """
    Score deletion effects using a left-aligned alternate allele instead of an empty string.
    Use the first nucleotide of the reference sequence as the alternate allele.
    This may provide better predictions for some variants.
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        gene_name: Optional gene(s) to focus on (e.g., 'ATG12' or 'ATG12,AP3S1')
        assay: Optional specific assay to score (e.g., 'RNA_SEQ', 'ATAC', 'CHIP_HISTONE')
    
    Returns:
        DataFrame with deletion effects
    """
    import pandas as pd
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return None
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    pos = variant_row['POS']
    ref = variant_row['REF']
    
    # Use first nucleotide of reference as alternate (left-aligned)
    alt = ref[0]
    
    print(f"Scoring deletion effects for {variant_name} (left-aligned, alt={alt})")
    print(f"Location: {chrom}:{pos:,}-{pos + len(ref):,}")
    if assay:
        print(f"Assay: {assay}")
    if gene_name:
        print(f"Genes: {gene_name}")
    
    # Create deletion variant (left-aligned with first base as alternate)
    deletion_variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,  # First nucleotide of reference
        name=f'{variant_name}_Deletion_LeftAlign'
    )
    
    # Map assay names to scorer configuration
    assay_config = {
        'RNA_SEQ': ('GeneMaskLFCScorer', None),
        'CAGE': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'PROCAP': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'ATAC': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'DNASE': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'CHIP_TF': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'CHIP_HISTONE': ('CenterMaskScorer', 'DIFF_LOG2_SUM'),
        'CONTACT_MAPS': ('ContactMapScorer', None),
    }
    
    # Select scorers
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    
    if assay and assay in assay_config:
        scorer_type, aggregation = assay_config[assay]
        selected_scorers = []
        
        for s in all_scorers.values():
            if (type(s).__name__ == scorer_type and 
                hasattr(s, 'requested_output') and 
                str(s.requested_output).split('.')[-1] == assay):
                
                if aggregation is None:
                    selected_scorers.append(s)
                elif hasattr(s, 'aggregation_type') and str(s.aggregation_type).split('.')[-1] == aggregation:
                    selected_scorers.append(s)
        
        if not selected_scorers:
            selected_scorers = list(all_scorers.values())
    else:
        selected_scorers = list(all_scorers.values())
    
    # Score variant
    print("Running variant scoring...")
    variant_scores = dna_model.score_variant(
        interval=deletion_variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=deletion_variant,
        variant_scorers=selected_scorers,
        organism=dna_client.Organism.HOMO_SAPIENS
    )
    
    # Convert to DataFrame and filter
    effects_df = variant_scorers.tidy_scores(variant_scores, match_gene_strand=True)
    filtered_df = effects_df[effects_df['biosample_name'] == PRIMARY_CELL_LINE].copy()
    
    # Parse gene names
    if gene_name:
        gene_list = [g.strip() for g in gene_name.split(',')]
        filtered_df = filtered_df[filtered_df['gene_name'].isin(gene_list)]
    
    if assay:
        filtered_df = filtered_df[filtered_df['output_type'] == assay]
    
    # Remove variant_id and scored_interval columns
    cols_to_drop = ['variant_id', 'scored_interval']
    filtered_df = filtered_df.drop(columns=[col for col in cols_to_drop if col in filtered_df.columns])
    
    # Sort by quantile_score
    filtered_df = filtered_df.sort_values('quantile_score', ascending=True)
    
    # Display results
    print(f"\n{'='*70}")
    print(f"RESULTS: {len(filtered_df)} effects in {PRIMARY_CELL_LINE}")
    print(f"{'='*70}")
    
    if gene_name and assay and not filtered_df.empty:
        for gene in gene_list:
            gene_df = filtered_df[filtered_df['gene_name'] == gene]
            if not gene_df.empty:
                print(f"\n{gene} - {assay}:")
                for _, row in gene_df.iterrows():
                    direction = "DECREASE" if row['quantile_score'] < -0.5 else ("INCREASE" if row['quantile_score'] > 0.5 else "Minimal effect")
                    print(f"  Raw: {row['raw_score']:.4f} | Quantile: {row['quantile_score']:.4f} | {direction}")
    elif not filtered_df.empty:
        summary = filtered_df.groupby('gene_name' if gene_name else 'output_type')['quantile_score'].agg(['count', 'mean']).round(3)
        print(summary)
    
    # Save results to CSV
    if not filtered_df.empty:
        filename_parts = [variant_name, PRIMARY_CELL_LINE, 'deletion_alt_scores']
        if assay:
            filename_parts.append(assay)
        if gene_name:
            filename_parts.append(gene_name.replace(',', '_').replace(' ', ''))
        filename = '_'.join(filename_parts) + '.csv'
        
        filtered_df.to_csv(filename, index=False)
        print(f"\n Saved results to: {filename}")
    
    return filtered_df

###############################################################################################################

def plot_variant_effects_ma(variant_name, assay='RNA_SEQ', track_filter=None, color_threshold=0.5):
    """
    Create an MA-style plot showing predicted deletion effects for ALL genes in the region.
    
    X-axis: Transformed |Quantile score| (0-0.999 mapped to 0-0.5, 0.999-1.0 mapped to 0.5-1.0)
    Y-axis: Raw score (magnitude and direction of effect)
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        assay: Assay type (default: 'RNA_SEQ')
        track_filter: How to handle multiple tracks per gene (REQUIRED):
                     'totalRNA' - Use only total RNA-seq tracks
                     'polyA' - Use only polyA RNA-seq tracks
                     'average' - Average all tracks for each gene
                     'max_effect' - Use track with maximum absolute quantile score per gene
        color_threshold: Quantile score threshold for coloring (default: 0.5)
                        Genes with |quantile| > threshold are colored red/blue and labeled
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    
    if track_filter is None:
        raise ValueError("track_filter is required. Choose: 'totalRNA', 'polyA', 'average', or 'max_effect'")
    
    # Get scores for ALL genes (don't filter by gene_name)
    df = score_deletion(variant_name, assay=assay)
    
    if df is None or df.empty:
        print("No data to plot")
        return
    
    # Handle different track filtering options
    if track_filter == 'totalRNA':
        df_filtered = df[df['track_name'].str.contains('total', case=False, na=False)]
        print(f"Filtered to total RNA-seq tracks: {len(df_filtered)} predictions")
        gene_summary = df_filtered.groupby('gene_name').agg({
            'raw_score': 'mean',
            'quantile_score': 'mean'
        }).reset_index()
        
    elif track_filter == 'polyA':
        df_filtered = df[df['track_name'].str.contains('polyA', case=False, na=False)]
        print(f"Filtered to polyA RNA-seq tracks: {len(df_filtered)} predictions")
        gene_summary = df_filtered.groupby('gene_name').agg({
            'raw_score': 'mean',
            'quantile_score': 'mean'
        }).reset_index()
        
    elif track_filter == 'average':
        print(f"Averaging all tracks: {len(df)} predictions")
        gene_summary = df.groupby('gene_name').agg({
            'raw_score': 'mean',
            'quantile_score': 'mean'
        }).reset_index()
        
    elif track_filter == 'max_effect':
        print(f"Using max absolute effect per gene: {len(df)} predictions")
        # For each gene, take the prediction with largest absolute quantile score
        idx = df.groupby('gene_name')['quantile_score'].apply(lambda x: x.abs().idxmax())
        gene_summary = df.loc[idx][['gene_name', 'raw_score', 'quantile_score']].reset_index(drop=True)
        
    else:
        raise ValueError(f"Invalid track_filter: {track_filter}. Choose: 'totalRNA', 'polyA', 'average', or 'max_effect'")
    
    print(f"Plotting {len(gene_summary)} genes")
    
    # Add absolute quantile score column
    gene_summary['abs_quantile'] = gene_summary['quantile_score'].abs()
    
    # Piecewise linear transformation: 
    # 0-0.999 maps to 0-0.5, 0.999-1.0 maps to 0.5-1.0
    def transform_quantile(q):
        if q <= 0.999:
            return q * (0.5 / 0.999)  # Scale 0-0.999 to 0-0.5
        else:
            return 0.5 + (q - 0.999) * (0.5 / 0.001)  # Scale 0.999-1.0 to 0.5-1.0
    
    gene_summary['transformed_quantile'] = gene_summary['abs_quantile'].apply(transform_quantile)
    
    # Categorize genes by effect using the specified threshold
    def categorize_effect(row):
        if row['quantile_score'] < -color_threshold:
            return 'decrease'
        elif row['quantile_score'] > color_threshold:
            return 'increase'
        else:
            return 'minimal'
    
    gene_summary['effect_category'] = gene_summary.apply(categorize_effect, axis=1)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Define colors
    color_map = {
        'decrease': '#FF0000',  # Standard R red
        'minimal': '#808080',   # Standard gray
        'increase': '#0000FF'   # Standard R blue
    }
    colors = gene_summary['effect_category'].map(color_map)
    
    # Size by extremity (larger = more extreme)
    sizes = gene_summary['abs_quantile'].apply(lambda x: 30 + x * 120)
    
    # Scatter plot with transformed x-axis
    scatter = ax.scatter(
        gene_summary['transformed_quantile'],
        gene_summary['raw_score'],
        c=colors,
        s=sizes,
        alpha=0.75,
        edgecolors='none',
        linewidth=0.5
    )
    
    # Add gene labels for all colored points (red and blue)
    for _, row in gene_summary.iterrows():
        if row['effect_category'] in ['decrease', 'increase']:
            ax.annotate(
                row['gene_name'],
                (row['transformed_quantile'], row['raw_score']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=9,
                alpha=0.8,
                fontweight='bold'
            )
    
    # Add threshold lines
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    ax.axvline(x=0.5, color='black', linestyle='--', linewidth=1.5, alpha=0.5, label='|quantile| = 0.999')
    
    # Labels
    ax.set_xlabel('|Quantile Score| (Effect Extremity)', fontsize=12)
    ax.set_ylabel('Raw Variant Score (Effect Magnitude)', fontsize=12)
    
    # Add custom x-axis ticks showing original quantile values
    # Left half: 0 to 0.999
    left_ticks = [0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999]
    left_positions = [transform_quantile(x) for x in left_ticks]
    # Right half: 0.999 to 1.0
    right_ticks = [0.9992, 0.9994, 0.9996, 0.9998, 1.0]
    right_positions = [transform_quantile(x) for x in right_ticks]
    
    all_positions = left_positions + right_positions
    all_labels = [f'{x:.2f}' if x < 0.99 else f'{x:.3f}' for x in left_ticks] + [f'{x:.4f}' for x in right_ticks]
    
    ax.set_xticks(all_positions)
    ax.set_xticklabels(all_labels, fontsize=9)
    
    # Title based on track filter
    track_label = {
        'totalRNA': 'total RNA-seq',
        'polyA': 'polyA RNA-seq',
        'average': 'average across tracks',
        'max_effect': 'max absolute effect'
    }
    title = f'MA Plot: Predicted {assay} Effects ({track_label[track_filter]})\n{variant_name} Deletion in {PRIMARY_CELL_LINE}'
    ax.set_title(title, fontsize=14)
    
    # Set x-axis limits
    ax.set_xlim(-0.05, 1.05)
    
    # Legend with counts
    from matplotlib.patches import Patch
    n_decrease = (gene_summary['effect_category'] == 'decrease').sum()
    n_minimal = (gene_summary['effect_category'] == 'minimal').sum()
    n_increase = (gene_summary['effect_category'] == 'increase').sum()
    
    legend_elements = [
        Patch(facecolor='red', label=f'Predicted decrease (|q| > {color_threshold}, n={n_decrease})'),
        Patch(facecolor='lightgray', label=f'Minimal effect (n={n_minimal})'),
        Patch(facecolor='blue', label=f'Predicted increase (|q| > {color_threshold}, n={n_increase})')
    ]
    ax.legend(handles=legend_elements, loc='best')
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle=':')
    
    plt.tight_layout()
    
    # Save
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_{assay}_{track_filter}_MA_plot.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f" Saved: {filename}")
    plt.show()

###############################################################################################################

def plot_variant_effects_coord_scatterplot(variant_name, assay='RNA_SEQ', track_filter=None, 
                                           color_threshold=0.9, gene_types='all', 
                                           ylim=None, xlim=None):
    """
    Plot predicted deletion effects along genomic coordinates (like CRISPR plot).
    
    X-axis: Genomic coordinates (TSS position)
    Y-axis: Raw score (effect magnitude and direction)
    Colors: Based on quantile score threshold
    
    Args:
        variant_name: Variant ID (e.g., 'LTR10.ATG12')
        assay: Assay type (default: 'RNA_SEQ')
        track_filter: How to handle multiple tracks per gene (REQUIRED):
                     'totalRNA' - Use only total RNA-seq tracks
                     'polyA' - Use only polyA RNA-seq tracks
                     'average' - Average all tracks for each gene
                     'max_effect' - Use track with maximum absolute quantile score per gene
        color_threshold: Quantile score threshold for coloring (default: 0.9)
        gene_types: Which gene types to include (default: 'all'):
                   'all' - All gene types (protein_coding, lncRNA, pseudogenes, etc.)
                   'protein_coding' - Only protein-coding genes
                   Or a list of gene types: ['protein_coding', 'lncRNA']
        ylim: Optional manual y-axis limits as [ymin, ymax] (default: None = auto symmetric)
              Example: ylim=[-0.2, 0.2]
        xlim: Optional manual x-axis limits as [xmin, xmax] (default: None = auto 1.5 MB window)
              Example: xlim=[112_500_000, 120_000_000]
    """
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'
    import pandas as pd
    import numpy as np
    from matplotlib.patches import Rectangle, Patch
    
    if track_filter is None:
        raise ValueError("track_filter is required. Choose: 'totalRNA', 'polyA', 'average', or 'max_effect'")
    
    # Find variant in DataFrame
    variant_row = variant_df[variant_df['ID'] == variant_name]
    if variant_row.empty:
        print(f"Error: {variant_name} not found!")
        return
    
    variant_row = variant_row.iloc[0]
    chrom = variant_row['CHROM']
    variant_start = variant_row['POS']
    variant_end = variant_row['POS'] + len(variant_row['REF'])
    variant_center = (variant_start + variant_end) // 2
    
    # Get scores for ALL genes
    df = score_deletion(variant_name, assay=assay)
    
    if df is None or df.empty:
        print("No data to plot")
        return
    
    # Handle different track filtering options
    if track_filter == 'totalRNA':
        df_filtered = df[df['track_name'].str.contains('total', case=False, na=False)]
        print(f"Filtered to total RNA-seq tracks: {len(df_filtered)} predictions")
        gene_summary = df_filtered.groupby('gene_name').agg({
            'raw_score': 'mean',
            'quantile_score': 'mean',
            'gene_id': 'first'
        }).reset_index()
        
    elif track_filter == 'polyA':
        df_filtered = df[df['track_name'].str.contains('polyA', case=False, na=False)]
        print(f"Filtered to polyA RNA-seq tracks: {len(df_filtered)} predictions")
        gene_summary = df_filtered.groupby('gene_name').agg({
            'raw_score': 'mean',
            'quantile_score': 'mean',
            'gene_id': 'first'
        }).reset_index()
        
    elif track_filter == 'average':
        print(f"Averaging all tracks: {len(df)} predictions")
        gene_summary = df.groupby('gene_name').agg({
            'raw_score': 'mean',
            'quantile_score': 'mean',
            'gene_id': 'first'
        }).reset_index()
        
    elif track_filter == 'max_effect':
        print(f"Using max absolute effect per gene: {len(df)} predictions")
        idx = df.groupby('gene_name')['quantile_score'].apply(lambda x: x.abs().idxmax())
        gene_summary = df.loc[idx][['gene_name', 'raw_score', 'quantile_score', 'gene_id']].reset_index(drop=True)
    
    print(f"\n=== Got {len(gene_summary)} unique genes from scoring ===")
    gene_summary = gene_summary.copy()
    
    # Use the GTF directly to get TSS positions for ALL gene types
    print(f"\n=== Fetching TSS coordinates from GTF ===")
    
    # Filter GTF to genes on the same chromosome
    gtf_chr = gtf[gtf['Chromosome'] == chrom].copy()
    print(f"Found {len(gtf_chr)} GTF entries on {chrom}")
    
    # Get gene-level information (TSS = transcription start site) AND gene_type
    gene_tss_map = {}
    gene_type_map = {}
    
    # Group by gene_id to get TSS and gene_type
    for gene_id, gene_data in gtf_chr.groupby('gene_id'):
        gene_id_base = gene_id.split('.')[0]
        
        # Get strand and gene_type
        strand = gene_data['Strand'].iloc[0]
        gene_type = gene_data['gene_type'].iloc[0] if 'gene_type' in gene_data.columns else None
        
        # Get all transcript starts and ends
        if strand == '+':
            tss = gene_data['Start'].min()
        else:
            tss = gene_data['End'].max()
        
        gene_tss_map[gene_id_base] = tss
        gene_tss_map[gene_id] = tss
        gene_type_map[gene_id_base] = gene_type
        gene_type_map[gene_id] = gene_type
        
        # Also store by gene_name if available
        if 'gene_name' in gene_data.columns and pd.notna(gene_data['gene_name'].iloc[0]):
            gene_name = gene_data['gene_name'].iloc[0]
            gene_tss_map[gene_name] = tss
            gene_type_map[gene_name] = gene_type
    
    print(f"Found TSS positions for {len(set(gene_tss_map.values()))} unique genes in GTF")
    
    # Map TSS and gene_type to gene_summary
    for idx, row in gene_summary.iterrows():
        gene_name = row['gene_name']
        gene_id_base = row['gene_id'].split('.')[0] if pd.notna(row['gene_id']) else None
        
        # Try gene_name first
        if gene_name in gene_tss_map:
            gene_summary.loc[idx, 'position'] = gene_tss_map[gene_name]
            gene_summary.loc[idx, 'gene_type'] = gene_type_map[gene_name]
        # Try gene_id
        elif gene_id_base and gene_id_base in gene_tss_map:
            gene_summary.loc[idx, 'position'] = gene_tss_map[gene_id_base]
            gene_summary.loc[idx, 'gene_type'] = gene_type_map[gene_id_base]
        # Try full gene_id with version
        elif pd.notna(row['gene_id']) and row['gene_id'] in gene_tss_map:
            gene_summary.loc[idx, 'position'] = gene_tss_map[row['gene_id']]
            gene_summary.loc[idx, 'gene_type'] = gene_type_map[row['gene_id']]
        # Try without version number from gene_name (in case gene_name has ENSG ID)
        elif gene_name.startswith('ENSG'):
            gene_name_base = gene_name.split('.')[0]
            if gene_name_base in gene_tss_map:
                gene_summary.loc[idx, 'position'] = gene_tss_map[gene_name_base]
                gene_summary.loc[idx, 'gene_type'] = gene_type_map[gene_name_base]
    
    # Filter by gene_types
    genes_before_filter = len(gene_summary)
    if gene_types != 'all':
        if gene_types == 'protein_coding':
            gene_summary = gene_summary[gene_summary['gene_type'] == 'protein_coding'].copy()
            print(f"Filtered to protein_coding genes: {len(gene_summary)}/{genes_before_filter} genes")
        elif isinstance(gene_types, list):
            gene_summary = gene_summary[gene_summary['gene_type'].isin(gene_types)].copy()
            print(f"Filtered to {gene_types} genes: {len(gene_summary)}/{genes_before_filter} genes")
        else:
            print(f"Warning: Unknown gene_types value '{gene_types}', using all genes")
    else:
        print(f"Using all gene types: {len(gene_summary)} genes")
    
    # Show statistics
    genes_with_pos = gene_summary['position'].notna().sum()
    genes_without_pos = gene_summary['position'].isna().sum()
    
    print(f"\n=== Position matching results ===")
    print(f"Matched: {genes_with_pos}/{len(gene_summary)} genes ({100*genes_with_pos/len(gene_summary):.1f}%)")
    
    if genes_without_pos > 0:
        missing = gene_summary[gene_summary['position'].isna()][['gene_name', 'gene_id']]
        print(f"Still missing {genes_without_pos} genes:")
        for _, row in missing.iterrows():
            print(f"  - {row['gene_name']} ({row['gene_id']})")
    
    # Remove genes without positions
    gene_summary_with_pos = gene_summary.dropna(subset=['position']).copy()
    
    if gene_summary_with_pos.empty:
        print("\nError: Could not map any genes to positions!")
        return gene_summary
    
    print(f"\n=== Creating plot with {len(gene_summary_with_pos)} genes ===")
    
    # Categorize genes by effect using quantile threshold
    def categorize_effect(row):
        if abs(row['quantile_score']) > color_threshold:
            if row['quantile_score'] < 0:
                return 'decrease'
            else:
                return 'increase'
        else:
            return 'minimal'
    
    gene_summary_with_pos['effect_category'] = gene_summary_with_pos.apply(categorize_effect, axis=1)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(5.5, 4))
    
    # Define colors
    color_map = {
        'decrease': '#FF0000',  # Standard R red
        'minimal': '#808080',   # Standard gray
        'increase': '#0000FF'   # Standard R blue
    }
    
    # Different sizes for colored vs gray dots
    def get_size(category):
        if category == 'minimal':
            return 20  # Smaller for gray
        else:
            return 40  # Much larger for colored

    colors = gene_summary_with_pos['effect_category'].map(color_map)
    sizes = gene_summary_with_pos['effect_category'].apply(get_size)
    
    # Scatter plot - Y-AXIS IS RAW SCORE, NO OUTLINES
    scatter = ax.scatter(
        gene_summary_with_pos['position'],
        gene_summary_with_pos['raw_score'],
        c=colors,
        s=sizes,
        alpha=0.75,
        edgecolors='none',
        linewidth=0
    )
    
    # Add gene labels for significant effects
    for _, row in gene_summary_with_pos.iterrows():
        if row['effect_category'] in ['decrease', 'increase']:
            ax.annotate(
                row['gene_name'],
                (row['position'], row['raw_score']),
                xytext=(0, -10 if row['raw_score'] < 0 else 5),
                textcoords='offset points',
                fontsize=11,
                alpha=0.9,
                ha='center',
            )
    
    # Set x-axis limits
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
        print(f"Using manual x-axis limits: [{xlim[0]:,}, {xlim[1]:,}]")
    else:
        plot_window = 1_500_000  # ~1 MB + 0.5 MB padding
        plot_start = variant_center - plot_window // 2
        plot_end = variant_center + plot_window // 2
        ax.set_xlim(plot_start, plot_end)
    
    # Get initial y-axis limits (needed for rectangle sizing)
    y_min_initial, y_max_initial = ax.get_ylim()
    y_range_initial = y_max_initial - y_min_initial
    
    # Draw the enhancer as a filled black rectangle at y=0
    rect_height = y_range_initial * 0.045
    min_display_width = 450_000  # Minimum visible width at this scale
    rect_width = max(variant_end - variant_start, min_display_width)
    rect = Rectangle(
        (variant_start, -rect_height/2),
        rect_width,
        rect_height,
        facecolor='black',
        edgecolor='black',
        linewidth=0,
        zorder=3
    )
    ax.add_patch(rect)
    # Add horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    # Labels
    ax.set_xlabel('Genome coordinates (bp)', fontsize=10)
    ax.set_ylabel('Raw Variant Score (Effect Magnitude)', fontsize=10)
    ax.set_title(f'{variant_name} Predicted Effects in {PRIMARY_CELL_LINE}', fontsize=12, fontweight='bold')
    
    # Format x-axis with comma separators
    ax.ticklabel_format(style='plain', axis='x')
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
    plt.xticks(rotation=0)
    
    # Add chromosome label
    ax.text(-0.13, -0.08, chrom, transform=ax.transAxes,
            fontsize=15, fontweight='bold', va='top', clip_on=False)
    
    # Grid
    # ax.grid(True, alpha=0.3, linestyle=':', axis='y')
    
    # Set y-axis limits: manual if provided, otherwise auto symmetric
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
        print(f"Using manual y-axis limits: [{ylim[0]}, {ylim[1]}]")
    else:
        # Make y-axis symmetric around zero (centered on y=0)
        max_abs_y = max(abs(y_min_initial), abs(y_max_initial))
        ax.set_ylim(-max_abs_y * 1.05, max_abs_y * 1.05)  # 5% padding
        print(f"Using auto symmetric y-axis limits: [{-max_abs_y * 1.05:.3f}, {max_abs_y * 1.05:.3f}]")
    
    # 1.5 Mb scale bar (bottom left, like CRISPRi plot)
    scale_len = 1_500_000
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    bar_x_start = xlims[0] + (xlims[1] - xlims[0]) * 0.03
    bar_y = ylims[0] + (ylims[1] - ylims[0]) * 0.06
    ax.plot([bar_x_start, bar_x_start + scale_len], [bar_y, bar_y],
            color='black', linewidth=1.5, solid_capstyle='butt')
    ax.text(bar_x_start, bar_y + (ylims[1] - ylims[0]) * 0.06,
            '1.5Mb', fontsize=10, va='top')
    
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, labelsize=8)
    
    import matplotlib.ticker as ticker
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    
    plt.tight_layout()
    
    # Save
    filename = f'{variant_name}_{PRIMARY_CELL_LINE}_{assay}_{track_filter}_genomic_plot.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {filename}")
    plt.show()
    
    #return gene_summary

###############################################################################################################