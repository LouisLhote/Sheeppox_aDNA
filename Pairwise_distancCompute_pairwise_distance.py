#!/usr/bin/env python3
"""
Script to compute pairwise similarity (ANI-like) for global alignment and sliding windows.
- Computes global pairwise similarity for all samples (0-100 scale, closer = higher)
- Uses sliding windows of 500bp with 100bp step
- For each window, computes pairwise similarity for all samples
- Groups by species and computes mean per species pair
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
from itertools import combinations
import seaborn as sns
from typing import Tuple, List, Dict
import os

def parse_msa(msa_file):
    """Parse MSA file and return sequences as dictionary."""
    sequences = {}
    print(f"Reading MSA file: {msa_file}...")
    
    with open(msa_file, 'r') as f:
        current_seq = None
        current_seq_data = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_seq is not None:
                    sequences[current_seq] = ''.join(current_seq_data).upper()
                # Start new sequence
                current_seq = line[1:].split()[0]  # Get accession ID
                current_seq_data = []
            else:
                current_seq_data.append(line)
        
        # Don't forget the last sequence
        if current_seq is not None:
            sequences[current_seq] = ''.join(current_seq_data).upper()
    
    print(f"  Loaded {len(sequences)} sequences")
    if sequences:
        print(f"  Alignment length: {len(list(sequences.values())[0])} bp")
    return sequences

def trim_sequences(sequences, trim_start=0, trim_end=0):
    """
    Trim sequences by removing specified number of bp from start and end.
    
    Args:
        sequences: Dictionary of sequence_name -> sequence_string
        trim_start: Number of bp to remove from start (default: 0)
        trim_end: Number of bp to remove from end (default: 0)
    
    Returns:
        Dictionary of trimmed sequences
    """
    if trim_start == 0 and trim_end == 0:
        return sequences
    
    trimmed = {}
    original_length = len(list(sequences.values())[0]) if sequences else 0
    
    for name, seq in sequences.items():
        if len(seq) < trim_start + trim_end:
            print(f"  Warning: Sequence {name} is too short to trim ({len(seq)} bp)")
            trimmed[name] = seq
        else:
            trimmed[name] = seq[trim_start:len(seq)-trim_end] if trim_end > 0 else seq[trim_start:]
    
    new_length = len(list(trimmed.values())[0]) if trimmed else 0
    print(f"  Trimmed sequences: removed {trim_start}bp from start, {trim_end}bp from end")
    print(f"  Original length: {original_length} bp")
    print(f"  New length: {new_length} bp")
    print(f"  Removed total: {original_length - new_length} bp ({100*(original_length - new_length)/original_length:.1f}%)")
    
    return trimmed

def parse_merged_regions(regions_file=None, regions_list=None):
    """
    Parse merged regions from a file or list.
    
    Args:
        regions_file: Path to file with merged regions (one per line, format: "start - end" or "start - end (length: ...)")
        regions_list: List of tuples [(start, end), ...] or list of strings ["start - end", ...]
    
    Returns:
        List of tuples [(start, end), ...] with 0-based indexing
    """
    regions = []
    
    if regions_file and os.path.exists(regions_file):
        print(f"  Reading merged regions from file: {regions_file}")
        with open(regions_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('-'):
                    continue
                # Parse format: "start - end" or "start - end (length: ...)"
                parts = line.split('-')
                if len(parts) >= 2:
                    try:
                        start = int(parts[0].strip().replace(',', ''))
                        end_part = parts[1].split()[0].strip().replace(',', '')
                        end = int(end_part)
                        # Convert to 0-based indexing
                        regions.append((start - 1, end))
                    except ValueError:
                        continue
    elif regions_list:
        print(f"  Using provided merged regions list")
        for region in regions_list:
            if isinstance(region, tuple):
                # Convert to 0-based if needed (assuming 1-based input)
                start, end = region
                regions.append((start - 1, end))
            elif isinstance(region, str):
                # Parse string format "start - end"
                parts = region.split('-')
                if len(parts) >= 2:
                    try:
                        start = int(parts[0].strip().replace(',', ''))
                        end = int(parts[1].split()[0].strip().replace(',', ''))
                        regions.append((start - 1, end))
                    except ValueError:
                        continue
    
    # Sort and merge overlapping regions
    if regions:
        regions = sorted(regions)
        merged = []
        for start, end in regions:
            if merged and start <= merged[-1][1]:
                # Merge with previous region
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        regions = merged
    
    return regions

def extract_regions(sequences, regions):
    """
    Extract specified regions from sequences.
    
    Args:
        sequences: Dictionary of sequence_name -> sequence_string
        regions: List of tuples [(start, end), ...] with 0-based indexing
    
    Returns:
        Dictionary of extracted sequences (concatenated regions)
    """
    if not regions:
        return {}
    
    extracted = {}
    for name, seq in sequences.items():
        extracted_parts = []
        for start, end in regions:
            # Ensure indices are within bounds
            start = max(0, min(start, len(seq)))
            end = max(start, min(end, len(seq)))
            if start < end:
                extracted_parts.append(seq[start:end])
        # Concatenate all extracted regions
        extracted[name] = ''.join(extracted_parts)
    
    total_length = sum(end - start for start, end in regions)
    print(f"  Extracted {len(regions)} region(s), total length: {total_length} bp")
    for start, end in regions:
        print(f"    Position {start+1}-{end} ({end-start} bp)")
    
    return extracted

def mask_regions(sequences, regions, mask_char='N'):
    """
    Mask specified regions in sequences by replacing with mask_char.
    
    Args:
        sequences: Dictionary of sequence_name -> sequence_string
        regions: List of tuples [(start, end), ...] with 0-based indexing
        mask_char: Character to use for masking (default: 'N')
    
    Returns:
        Dictionary of masked sequences
    """
    if not regions:
        return sequences
    
    masked = {}
    seq_length = len(list(sequences.values())[0]) if sequences else 0
    total_masked = 0
    
    for name, seq in sequences.items():
        masked_seq = list(seq)
        for start, end in regions:
            # Ensure indices are within bounds
            start = max(0, min(start, len(masked_seq)))
            end = max(start, min(end, len(masked_seq)))
            if start < end:
                masked_seq[start:end] = [mask_char] * (end - start)
                total_masked += (end - start)
        masked[name] = ''.join(masked_seq)
    
    print(f"  Masked {len(regions)} region(s)")
    for start, end in regions:
        print(f"    Position {start+1}-{end} ({end-start} bp)")
    print(f"  Total masked: {total_masked} bp ({100*total_masked/seq_length:.2f}% of alignment)")
    
    return masked

def adjust_regions_after_trim(regions, trim_start, trim_end, original_length):
    """
    Adjust region coordinates after trimming.
    
    Args:
        regions: List of tuples [(start, end), ...] with 0-based indexing (original coordinates)
        trim_start: Number of bp removed from start
        trim_end: Number of bp removed from end
        original_length: Original sequence length
    
    Returns:
        List of adjusted regions (regions that fall within trimmed sequence)
    """
    if not regions:
        return []
    
    adjusted = []
    new_length = original_length - trim_start - trim_end
    
    for start, end in regions:
        # Adjust for trimming
        new_start = start - trim_start
        new_end = end - trim_start
        
        # Check if region is still within bounds after trimming
        if new_end <= 0:
            # Region was completely removed by trimming
            continue
        elif new_start < 0:
            # Region partially removed by start trimming
            new_start = 0
        
        if new_start >= new_length:
            # Region was completely removed by trimming
            continue
        elif new_end > new_length:
            # Region partially removed by end trimming
            new_end = new_length
        
        if new_start < new_end:
            adjusted.append((new_start, new_end))
    
    return adjusted

def calculate_pairwise_similarity(seq1, seq2, start_pos=None, end_pos=None):
    """
    Calculate pairwise similarity (ANI-like) between two sequences.
    Returns similarity as percentage (0-100), where closer samples are closer to 100.
    Excludes gaps and ambiguous positions.
    If start_pos and end_pos are provided, only compares that region.
    """
    if start_pos is not None and end_pos is not None:
        chunk1 = seq1[start_pos:end_pos]
        chunk2 = seq2[start_pos:end_pos]
    else:
        chunk1 = seq1
        chunk2 = seq2
    
    # Find positions where both sequences have valid nucleotides (not gaps, not N)
    valid_positions = []
    for i in range(len(chunk1)):
        if (chunk1[i] in 'ATCG' and chunk2[i] in 'ATCG' and 
            chunk1[i] != 'N' and chunk2[i] != 'N'):
            valid_positions.append(i)
    
    if len(valid_positions) == 0:
        return np.nan  # No valid positions to compare
    
    # Count matches
    matches = sum(1 for i in valid_positions if chunk1[i] == chunk2[i])
    similarity = matches / len(valid_positions) * 100  # Percentage (0-100)
    
    return similarity

def get_species_mapping(metadata_file):
    """Create mapping from accession to species."""
    metadata = pd.read_csv(metadata_file, sep='\t')
    metadata.columns = ['accession', 'relative_to_2024', 'Virus']
    
    # Create mapping
    accession_to_species = dict(zip(metadata['accession'], metadata['Virus']))
    
    # Also create base name mapping (handle _rev suffix)
    base_name_to_species = {}
    for meta_acc, sp in accession_to_species.items():
        base_name = meta_acc.replace('_rev', '').replace('rev', '')
        if base_name not in base_name_to_species:
            base_name_to_species[base_name] = sp
    
    return accession_to_species, base_name_to_species

def get_species_for_accession(acc, accession_to_species, base_name_to_species):
    """Get species for an accession, handling various naming conventions."""
    # Handle reverse sequences: _R_SAMPLE means the sample is SAMPLE
    if acc.startswith('_R_'):
        base_acc = acc.replace('_R_', '')
        if base_acc in base_name_to_species:
            return base_name_to_species[base_acc]
        elif base_acc in accession_to_species:
            return accession_to_species[base_acc]
        # Try prefix matching
        for base_name, sp in base_name_to_species.items():
            if base_acc == base_name or base_acc.startswith(base_name) or base_name.startswith(base_acc):
                return sp
    
    # Handle TCLL002 -> TCL002
    if acc == 'TCLL002':
        if 'TCL002' in base_name_to_species:
            return base_name_to_species['TCL002']
        elif 'TCL002' in accession_to_species:
            return accession_to_species['TCL002']
    
    # Try direct match
    if acc in accession_to_species:
        return accession_to_species[acc]
    
    # Try base name match
    if acc in base_name_to_species:
        return base_name_to_species[acc]
    
    # Try prefix matching
    for meta_acc, sp in accession_to_species.items():
        if acc.startswith(meta_acc) or meta_acc.startswith(acc):
            return sp
    
    for base_name, sp in base_name_to_species.items():
        if acc.startswith(base_name) or base_name.startswith(acc):
            return sp
    
    return None

def compute_global_pairwise_similarities(sequences, accession_to_species, base_name_to_species):
    """Compute global pairwise similarities (ANI-like) for all samples."""
    print("\nComputing global pairwise similarities...")
    
    sample_names = list(sequences.keys())
    n_samples = len(sample_names)
    total_comparisons = n_samples * (n_samples - 1) // 2
    print(f"  Total comparisons to perform: {total_comparisons}")
    
    # Compute pairwise similarities
    similarities = []
    species_pairs = []
    pairwise_records: List[Dict] = []
    comparison_count = 0
    
    for i, acc1 in enumerate(sample_names):
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  Processing sample {i+1}/{n_samples}...")
        for j, acc2 in enumerate(sample_names):
            if i < j:  # Only compute once per pair
                comparison_count += 1
                if comparison_count % 100 == 0:
                    print(f"    Completed {comparison_count}/{total_comparisons} comparisons ({100*comparison_count/total_comparisons:.1f}%)...")
                sim = calculate_pairwise_similarity(sequences[acc1], sequences[acc2])
                if not np.isnan(sim):
                    sp1 = get_species_for_accession(acc1, accession_to_species, base_name_to_species)
                    sp2 = get_species_for_accession(acc2, accession_to_species, base_name_to_species)
                    if sp1 and sp2:
                        similarities.append(sim)
                        # Sort species names to ensure consistent pair naming
                        pair = tuple(sorted([sp1, sp2]))
                        species_pairs.append(pair)
                        pairwise_records.append({
                            'sample1': acc1,
                            'sample2': acc2,
                            'species1': sp1,
                            'species2': sp2,
                            'species_pair': f"{pair[0]}_vs_{pair[1]}",
                            'similarity': sim
                        })
    
    # Create DataFrame
    if len(similarities) > 0:
        global_df = pd.DataFrame({
            'similarity': similarities,
            'species_pair': [f"{p[0]}_vs_{p[1]}" for p in species_pairs]
        })
        pairwise_df = pd.DataFrame(pairwise_records)
    else:
        global_df = pd.DataFrame(columns=['similarity', 'species_pair'])
        pairwise_df = pd.DataFrame(columns=['sample1', 'sample2', 'species1', 'species2', 'species_pair', 'similarity'])
    
    print(f"\nGlobal pairwise similarities computed:")
    print(f"  Total comparisons: {len(similarities)}")
    print(f"\nMean similarities by species pair:")
    for pair_name in global_df['species_pair'].unique():
        pair_similarities = global_df[global_df['species_pair'] == pair_name]['similarity']
        print(f"  {pair_name}: {pair_similarities.mean():.2f}% ± {pair_similarities.std():.2f}%")
    
    return global_df, pairwise_df

def compute_sliding_window_similarities(sequences, accession_to_species, base_name_to_species, 
                                        window_size=500, step_size=100):
    """Compute pairwise similarities (ANI-like) for sliding windows."""
    if not sequences:
        return pd.DataFrame()
    
    alignment_length = len(list(sequences.values())[0])
    sample_names = list(sequences.keys())
    n_samples = len(sample_names)
    
    print(f"\nComputing sliding window similarities...")
    print(f"  Window size: {window_size}bp")
    print(f"  Step size: {step_size}bp")
    print(f"  Alignment length: {alignment_length}bp")
    
    results = []
    windows = []
    
    # Generate windows
    for start in range(0, alignment_length, step_size):
        end = min(start + window_size, alignment_length)
        if end - start < window_size * 0.5:  # Skip windows that are too small
            break
        windows.append((start, end))
    
    print(f"  Total windows: {len(windows)}")
    
    # Process each window
    for window_idx, (start, end) in enumerate(windows):
        if (window_idx + 1) % 50 == 0:
            print(f"  Processing window {window_idx + 1}/{len(windows)}...")
        
        window_center = (start + end) / 2
        chunk_number = window_idx + 1  # 1-indexed chunk number
        
        # Compute pairwise similarities for this window
        window_similarities = defaultdict(list)
        
        for i, acc1 in enumerate(sample_names):
            for j, acc2 in enumerate(sample_names):
                if i < j:  # Only compute once per pair
                    sim = calculate_pairwise_similarity(sequences[acc1], sequences[acc2], start, end)
                    if not np.isnan(sim):
                        sp1 = get_species_for_accession(acc1, accession_to_species, base_name_to_species)
                        sp2 = get_species_for_accession(acc2, accession_to_species, base_name_to_species)
                        if sp1 and sp2:
                            # Sort species names to ensure consistent pair naming
                            pair = tuple(sorted([sp1, sp2]))
                            pair_name = f"{pair[0]}_vs_{pair[1]}"
                            window_similarities[pair_name].append(sim)
        
        # Compute mean for each species pair
        for pair_name, sims in window_similarities.items():
            if sims:
                results.append({
                    'chunk_number': chunk_number,
                    'window_start': start,
                    'window_end': end,
                    'window_center': window_center,
                    'window_size': end - start,
                    'chunk_position': f"{start}-{end}",  # Human-readable chunk position
                    'species_pair': pair_name,
                    'mean_similarity': np.mean(sims),
                    'std_similarity': np.std(sims),
                    'n_comparisons': len(sims)
                })
    
    return pd.DataFrame(results)

def plot_global_similarities(global_df, output_file='global_pairwise_similarities.png'):
    """Plot global pairwise similarities by species pair."""
    if global_df.empty:
        print("No data to plot for global similarities.")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Box plot
    ax1 = axes[0]
    species_pairs = sorted(global_df['species_pair'].unique())
    data_to_plot = [global_df[global_df['species_pair'] == pair]['similarity'].values 
                     for pair in species_pairs]
    
    bp = ax1.boxplot(data_to_plot, tick_labels=species_pairs, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)
    
    ax1.set_ylabel('Pairwise Similarity (%)', fontsize=12)
    ax1.set_title('Global Pairwise Similarities (ANI-like) by Species Pair', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    # Zoom in on relevant range (80-100)
    min_val = max(0, global_df['similarity'].min() - 2)
    max_val = min(105, global_df['similarity'].max() + 2)
    ax1.set_ylim([min_val, max_val])
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Violin plot
    ax2 = axes[1]
    sns.violinplot(data=global_df, x='species_pair', y='similarity', ax=ax2)
    ax2.set_ylabel('Pairwise Similarity (%)', fontsize=12)
    ax2.set_title('Global Pairwise Similarities Distribution', fontsize=14, fontweight='bold')
    # Zoom in on relevant range
    ax2.set_ylim([min_val, max_val])
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGlobal similarities plot saved: {output_file}")
    plt.close()

def plot_species_distance_bars(global_df, output_file='global_pairwise_distances.png'):
    """Plot pairwise distance (100 - similarity) by species pair using boxplot/violin plot."""
    if global_df.empty:
        print("No data to plot for global distances.")
        return
    df = global_df.copy()
    df['distance'] = 100 - df['similarity']
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Box plot
    ax1 = axes[0]
    species_pairs = sorted(df['species_pair'].unique())
    data_to_plot = [df[df['species_pair'] == pair]['distance'].values 
                     for pair in species_pairs]
    
    bp = ax1.boxplot(data_to_plot, tick_labels=species_pairs, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('salmon')
        patch.set_alpha(0.7)
    
    ax1.set_ylabel('Pairwise Distance (100 - similarity)', fontsize=12)
    ax1.set_title('Global Pairwise Distances by Species Pair', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    # Zoom in on relevant range
    min_val = max(0, df['distance'].min() - 2)
    max_val = min(25, df['distance'].max() + 2)  # Distances are typically 0-25
    ax1.set_ylim([min_val, max_val])
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Violin plot
    ax2 = axes[1]
    sns.violinplot(data=df, x='species_pair', y='distance', ax=ax2, color='salmon', alpha=0.7)
    ax2.set_ylabel('Pairwise Distance (100 - similarity)', fontsize=12)
    ax2.set_title('Global Pairwise Distances Distribution', fontsize=14, fontweight='bold')
    # Zoom in on relevant range
    ax2.set_ylim([min_val, max_val])
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGlobal distance plot saved: {output_file}")
    plt.close()

def plot_global_distance_distribution(global_df, output_file='global_pairwise_distance_distribution.png'):
    """Plot distribution of global pairwise distances."""
    if global_df.empty:
        print("No data to plot for global distance distribution.")
        return
    df = global_df.copy()
    df['distance'] = 100 - df['similarity']
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Histogram
    ax1 = axes[0]
    ax1.hist(df['distance'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Pairwise Distance (100 - similarity)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Global Pairwise Distance Distribution (Histogram)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # KDE plot by species pair
    ax2 = axes[1]
    species_pairs = sorted(df['species_pair'].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(species_pairs)))
    for idx, pair in enumerate(species_pairs):
        pair_data = df[df['species_pair'] == pair]['distance']
        sns.kdeplot(data=pair_data, label=pair, ax=ax2, color=colors[idx], linewidth=2)
    ax2.set_xlabel('Pairwise Distance (100 - similarity)', fontsize=12)
    ax2.set_ylabel('Density', fontsize=12)
    ax2.set_title('Global Pairwise Distance Distribution by Species Pair (KDE)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGlobal distance distribution plot saved: {output_file}")
    plt.close()

def save_global_similarities_summary(global_df, output_file):
    """Save summary statistics for global pairwise similarities."""
    if global_df.empty:
        return
    
    summary_data = []
    for pair in sorted(global_df['species_pair'].unique()):
        pair_data = global_df[global_df['species_pair'] == pair]['similarity']
        summary_data.append({
            'species_pair': pair,
            'mean_similarity': pair_data.mean(),
            'std_similarity': pair_data.std(),
            'min_similarity': pair_data.min(),
            'max_similarity': pair_data.max(),
            'median_similarity': pair_data.median(),
            'q25_similarity': pair_data.quantile(0.25),
            'q75_similarity': pair_data.quantile(0.75),
            'n_comparisons': len(pair_data)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_file, index=False)
    print(f"  Summary statistics saved to: {output_file}")

def save_window_similarities_summary(window_df, output_file):
    """Save summary statistics for sliding window pairwise similarities."""
    if window_df.empty:
        return
    
    summary_data = []
    for pair in sorted(window_df['species_pair'].unique()):
        pair_data = window_df[window_df['species_pair'] == pair]
        summary_data.append({
            'species_pair': pair,
            'mean_similarity_across_windows': pair_data['mean_similarity'].mean(),
            'std_similarity_across_windows': pair_data['mean_similarity'].std(),
            'min_similarity': pair_data['mean_similarity'].min(),
            'max_similarity': pair_data['mean_similarity'].max(),
            'median_similarity': pair_data['mean_similarity'].median(),
            'q25_similarity': pair_data['mean_similarity'].quantile(0.25),
            'q75_similarity': pair_data['mean_similarity'].quantile(0.75),
            'n_windows': len(pair_data)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_file, index=False)
    print(f"  Summary statistics saved to: {output_file}")

def plot_sliding_window_similarities(window_df, output_file='sliding_window_pairwise_similarities.png', masked_regions=None):
    """Plot sliding window pairwise similarities."""
    if window_df.empty:
        print("No data to plot for sliding windows.")
        return
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot each species pair
    species_pairs = sorted(window_df['species_pair'].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(species_pairs)))
    
    for idx, pair in enumerate(species_pairs):
        pair_data = window_df[window_df['species_pair'] == pair].sort_values('window_center')
        ax.plot(pair_data['window_center'], pair_data['mean_similarity'], 
                'o-', label=pair, color=colors[idx], linewidth=2, markersize=4, alpha=0.7)
    
    # Add masked regions as shaded areas
    if masked_regions:
        ylim = ax.get_ylim()
        for start, end in masked_regions:
            ax.axvspan(start, end, alpha=0.2, color='red', label='Masked region' if start == masked_regions[0][0] else '')
    
    ax.set_xlabel('Window Position in Alignment (bp)', fontsize=12)
    ax.set_ylabel('Mean Pairwise Similarity (%)', fontsize=12)
    title = f'Pairwise Similarities (ANI-like) Across Sliding Windows (500bp windows, 100bp step)'
    if masked_regions:
        title += '\n(Red shaded areas = masked regions)'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, linestyle='--')
    # Zoom in on relevant range (80-100)
    min_val = max(0, window_df['mean_similarity'].min() - 2)
    max_val = min(105, window_df['mean_similarity'].max() + 2)
    ax.set_ylim([min_val, max_val])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSliding window similarities plot saved: {output_file}")
    plt.close()

def plot_sliding_window_heatmap(window_df, output_file='sliding_window_similarities_heatmap.png', masked_regions=None):
    """Plot sliding window similarities as a heatmap across the genome."""
    if window_df.empty:
        print("No data to plot for sliding window heatmap.")
        return
    
    # Pivot the data for heatmap
    pivot_data = window_df.pivot_table(
        values='mean_similarity',
        index='window_center',
        columns='species_pair',
        aggfunc='mean'
    )
    
    # Sort by window center
    pivot_data = pivot_data.sort_index()
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create heatmap with zoomed range (80-100 for better visibility)
    vmin = max(0, pivot_data.min().min() - 2)
    vmax = min(100, pivot_data.max().max() + 2)
    im = ax.imshow(pivot_data.T, aspect='auto', cmap='RdYlGn', vmin=vmin, vmax=vmax, interpolation='bilinear')
    
    # Add masked regions as vertical lines
    if masked_regions:
        n_windows = len(pivot_data)
        window_positions = pivot_data.index.values
        for start, end in masked_regions:
            # Find closest window indices
            start_idx = np.argmin(np.abs(window_positions - start))
            end_idx = np.argmin(np.abs(window_positions - end))
            ax.axvline(start_idx, color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax.axvline(end_idx, color='red', linestyle='--', linewidth=2, alpha=0.7)
            # Add label
            if start == masked_regions[0][0]:
                ax.text(start_idx, -0.5, 'Masked', color='red', fontsize=9, ha='center', rotation=90)
    
    # Set labels
    ax.set_xlabel('Window Position in Alignment (bp)', fontsize=12)
    ax.set_ylabel('Species Pair', fontsize=12)
    title = 'Pairwise Similarities (ANI-like) Heatmap Across Genome\n(500bp windows, 100bp step)'
    if masked_regions:
        title += '\n(Red dashed lines = masked regions)'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Set y-axis labels
    ax.set_yticks(range(len(pivot_data.columns)))
    ax.set_yticklabels(pivot_data.columns)
    
    # Set x-axis labels (show every nth window to avoid crowding)
    n_ticks = 10
    tick_positions = np.linspace(0, len(pivot_data) - 1, n_ticks, dtype=int)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([f"{int(pivot_data.index[pos])}" for pos in tick_positions], rotation=45, ha='right')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Mean Pairwise Similarity (%)', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSliding window heatmap saved: {output_file}")
    plt.close()

def _kmeans_numpy(data: np.ndarray, k: int, random_state: int = 42, max_iter: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Lightweight k-means implementation to avoid external dependencies."""
    rng = np.random.default_rng(random_state)
    n_samples = data.shape[0]
    # Initialize centers with random unique samples
    initial_idx = rng.choice(n_samples, size=min(k, n_samples), replace=False)
    centers = data[initial_idx]
    labels = np.zeros(n_samples, dtype=int)
    
    for _ in range(max_iter):
        # Compute distances to centers
        dists = np.linalg.norm(data[:, None, :] - centers[None, :, :], axis=2)
        new_labels = np.argmin(dists, axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        # Update centers; if a cluster is empty, reinitialize it to a random sample
        for c in range(centers.shape[0]):
            members = data[labels == c]
            if len(members) == 0:
                centers[c] = data[rng.integers(0, n_samples)]
            else:
                centers[c] = members.mean(axis=0)
    return labels, centers

def cluster_sliding_windows(window_df, n_clusters=4, random_state=42,
                            output_csv='sliding_window_clusters.csv',
                            output_plot='sliding_window_clusters.png'):
    """
    Cluster genomic windows based on their pairwise similarity profiles.
    Returns a DataFrame with cluster assignments and saves CSV/plot.
    """
    if window_df.empty:
        print("No window data available for clustering.")
        return pd.DataFrame()
    
    # Pivot to matrix: rows = windows, cols = species pairs
    pivot = window_df.pivot_table(
        index=['chunk_number', 'window_start', 'window_end', 'window_center'],
        columns='species_pair',
        values='mean_similarity',
        aggfunc='mean'
    )
    pivot = pivot.sort_index()
    
    # Drop columns that are entirely NaN
    pivot = pivot.dropna(axis=1, how='all')
    if pivot.empty:
        print("Pivot table is empty after dropping NaN columns; skipping clustering.")
        return pd.DataFrame()
    
    # Fill remaining NaN with column means
    col_means = pivot.mean(axis=0)
    pivot_filled = pivot.fillna(col_means)
    
    # Standardize features
    data = pivot_filled.to_numpy()
    col_std = data.std(axis=0)
    col_std[col_std == 0] = 1.0
    data_std = (data - data.mean(axis=0)) / col_std
    
    # Run k-means
    k = min(n_clusters, data_std.shape[0])  # do not exceed number of windows
    labels, centers = _kmeans_numpy(data_std, k=k, random_state=random_state)
    
    # Build output
    result = pivot_filled.reset_index()
    result['cluster'] = labels
    # Add cluster size
    cluster_sizes = result['cluster'].value_counts().to_dict()
    result['cluster_size'] = result['cluster'].map(cluster_sizes)
    
    # Save CSV
    result.to_csv(output_csv, index=False)
    print(f"\nSliding window clusters saved to: {output_csv}")
    
    # Save cluster summary
    cluster_summary = []
    for cluster_id in sorted(result['cluster'].unique()):
        cluster_data = result[result['cluster'] == cluster_id]
        cluster_summary.append({
            'cluster_id': cluster_id,
            'n_windows': len(cluster_data),
            'mean_window_center': cluster_data['window_center'].mean(),
            'min_window_center': cluster_data['window_center'].min(),
            'max_window_center': cluster_data['window_center'].max(),
            'mean_window_start': cluster_data['window_start'].mean(),
            'mean_window_end': cluster_data['window_end'].mean()
        })
    
    summary_df = pd.DataFrame(cluster_summary)
    summary_csv = output_csv.replace('.csv', '_summary.csv')
    summary_df.to_csv(summary_csv, index=False)
    print(f"  Cluster summary saved to: {summary_csv}")
    
    # Plot clusters along genome (window center vs cluster id)
    fig, ax = plt.subplots(figsize=(14, 4))
    unique_clusters = sorted(result['cluster'].unique())
    colors_map = plt.cm.tab10(np.linspace(0, 1, len(unique_clusters)))
    cluster_colors = {c: colors_map[i] for i, c in enumerate(unique_clusters)}
    
    for cluster_id in unique_clusters:
        cluster_data = result[result['cluster'] == cluster_id]
        ax.scatter(cluster_data['window_center'], cluster_data['cluster'], 
                  c=[cluster_colors[cluster_id]], label=f'Cluster {cluster_id} (n={len(cluster_data)})',
                  s=30, alpha=0.8)
    
    ax.set_xlabel('Window Position in Alignment (bp)', fontsize=12)
    ax.set_ylabel('Cluster ID', fontsize=12)
    ax.set_title(f'Window Clusters (k={k}) based on pairwise similarity profiles', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10, loc='best', ncol=min(len(unique_clusters), 4))
    ax.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    print(f"Sliding window cluster plot saved: {output_plot}")
    plt.close()
    
    return result

def plot_window_regions_pca(window_df, output_plot='sliding_window_regions_pca.png'):
    """
    Perform PCA on window regions based on their pairwise similarity profiles.
    Each window is represented by its similarity values across species pairs.
    """
    if window_df.empty:
        print("No window data available for PCA.")
        return
    
    # Pivot to matrix: rows = windows, cols = species pairs
    pivot = window_df.pivot_table(
        index=['chunk_number', 'window_start', 'window_end', 'window_center'],
        columns='species_pair',
        values='mean_similarity',
        aggfunc='mean'
    )
    pivot = pivot.sort_index()
    
    # Drop columns that are entirely NaN
    pivot = pivot.dropna(axis=1, how='all')
    if pivot.empty:
        print("Pivot table is empty; skipping window regions PCA.")
        return
    
    # Fill remaining NaN with column means
    col_means = pivot.mean(axis=0)
    pivot_filled = pivot.fillna(col_means)
    
    # Standardize features
    data = pivot_filled.to_numpy()
    col_std = data.std(axis=0)
    col_std[col_std == 0] = 1.0
    data_std = (data - data.mean(axis=0)) / col_std
    
    # Center the data
    data_centered = data_std - data_std.mean(axis=0, keepdims=True)
    
    # SVD for PCA
    U, S, Vt = np.linalg.svd(data_centered, full_matrices=False)
    components = U[:, :2] * S[:2]
    
    # Calculate variance explained
    variance_explained = (S[:2] ** 2) / (S ** 2).sum() * 100
    
    # Get window information for annotation
    pivot_reset = pivot_filled.reset_index()
    window_centers = pivot_reset['window_center'].values
    window_starts = pivot_reset['window_start'].values
    window_ends = pivot_reset['window_end'].values
    
    # Detect outliers using IQR method for both PC1 and PC2
    pc1 = components[:, 0]
    pc2 = components[:, 1]
    
    # Calculate IQR for PC1
    q1_pc1, q3_pc1 = np.percentile(pc1, [25, 75])
    iqr_pc1 = q3_pc1 - q1_pc1
    lower_bound_pc1 = q1_pc1 - 1.5 * iqr_pc1
    upper_bound_pc1 = q3_pc1 + 1.5 * iqr_pc1
    
    # Calculate IQR for PC2
    q1_pc2, q3_pc2 = np.percentile(pc2, [25, 75])
    iqr_pc2 = q3_pc2 - q1_pc2
    lower_bound_pc2 = q1_pc2 - 1.5 * iqr_pc2
    upper_bound_pc2 = q3_pc2 + 1.5 * iqr_pc2
    
    # Identify outliers (points outside IQR bounds in either dimension)
    outlier_mask = (pc1 < lower_bound_pc1) | (pc1 > upper_bound_pc1) | \
                   (pc2 < lower_bound_pc2) | (pc2 > upper_bound_pc2)
    outlier_indices = np.where(outlier_mask)[0]
    
    # Plot PCA
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Plot all points
    scatter = ax.scatter(components[:, 0], components[:, 1], 
                        c=window_centers, cmap='viridis', 
                        s=30, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Highlight and annotate outliers
    if len(outlier_indices) > 0:
        ax.scatter(components[outlier_indices, 0], components[outlier_indices, 1],
                  c='red', s=100, alpha=0.8, edgecolors='darkred', linewidth=2,
                  marker='X', label=f'Outliers (n={len(outlier_indices)})', zorder=5)
        
        # Annotate outliers with window positions
        for idx in outlier_indices:
            # Format window position as "start-end" or just center if preferred
            pos_text = f"{int(window_starts[idx])}-{int(window_ends[idx])}"
            # Alternative: just show center: pos_text = f"{int(window_centers[idx])}"
            ax.annotate(pos_text,
                       (components[idx, 0], components[idx, 1]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, color='darkred', fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7, edgecolor='darkred'),
                       arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))
    
    ax.set_xlabel(f'PC1 ({variance_explained[0]:.1f}% variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({variance_explained[1]:.1f}% variance)', fontsize=12)
    ax.set_title('PCA of Genomic Window Regions\n(based on pairwise similarity profiles)', 
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add legend if outliers exist
    if len(outlier_indices) > 0:
        ax.legend(loc='best', fontsize=10)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Window Position (bp)', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    print(f"\nWindow regions PCA plot saved: {output_plot}")
    if len(outlier_indices) > 0:
        print(f"  Annotated {len(outlier_indices)} outlier windows with positions")
    plt.close()

def build_sample_distance_outputs(pairwise_df: pd.DataFrame,
                                  sim_matrix_csv='global_pairwise_similarity_matrix.csv',
                                  dist_matrix_csv='global_pairwise_distance_matrix.csv',
                                  pca_plot='global_pairwise_pca.png'):
    """Build sample-level similarity/distance matrices and a PCA plot."""
    if pairwise_df.empty:
        print("No pairwise sample data to build matrices.")
        return
    
    samples = sorted(set(pairwise_df['sample1']).union(set(pairwise_df['sample2'])))
    idx_map = {s: i for i, s in enumerate(samples)}
    n = len(samples)
    
    sim_matrix = np.full((n, n), np.nan)
    np.fill_diagonal(sim_matrix, 100.0)
    
    for _, row in pairwise_df.iterrows():
        i = idx_map[row['sample1']]
        j = idx_map[row['sample2']]
        sim_matrix[i, j] = row['similarity']
        sim_matrix[j, i] = row['similarity']
    
    sim_df = pd.DataFrame(sim_matrix, index=samples, columns=samples)
    sim_df.to_csv(sim_matrix_csv)
    print(f"\nGlobal pairwise similarity matrix saved: {sim_matrix_csv}")
    
    dist_matrix = 100 - sim_matrix
    np.fill_diagonal(dist_matrix, 0.0)
    dist_df = pd.DataFrame(dist_matrix, index=samples, columns=samples)
    dist_df.to_csv(dist_matrix_csv)
    print(f"Global pairwise distance matrix saved: {dist_matrix_csv}")
    
    # PCA on similarity matrix (fill NaN with column means)
    sim_filled = sim_df.copy()
    for col in sim_filled.columns:
        col_mean = sim_filled[col].mean()
        sim_filled[col] = sim_filled[col].fillna(col_mean)
    data = sim_filled.to_numpy()
    # Center
    data = data - data.mean(axis=0, keepdims=True)
    # SVD for PCA
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    components = U[:, :2] * S[:2]
    
    plt.figure(figsize=(8, 6))
    plt.scatter(components[:, 0], components[:, 1], c='steelblue', alpha=0.8)
    for i, sample in enumerate(samples):
        plt.text(components[i, 0], components[i, 1], sample, fontsize=8)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of Pairwise Similarities')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(pca_plot, dpi=300, bbox_inches='tight')
    print(f"PCA plot saved: {pca_plot}")
    plt.close()

def main():
    # File paths
    msa_file = 'CAPRI_aln.fa'
    metadata_file = 'METADATAS.tsv'
    
    # Create output directory
    output_dir = 'PAIRWISE'
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}/")
    
    # Samples to exclude from all analyses (matches by raw or base name)
    exclude_samples = {'BA27', 'YG16', '79RR'}
    def _normalize_name(name):
        base = name
        if base.startswith('_R_'):
            base = base.replace('_R_', '')
        return base.replace('_rev', '').replace('rev', '')
    def should_exclude(name):
        norm = _normalize_name(name)
        # Also try stripping common rev markers on already normalized strings
        alt_norm = norm.replace('_rev', '').replace('rev', '')
        return norm in exclude_samples or alt_norm in exclude_samples
    
    # Parse MSA
    print("="*80)
    print("PARSING MSA FILE")
    print("="*80)
    sequences = parse_msa(msa_file)
    # Drop excluded samples
    filtered_sequences = {acc: seq for acc, seq in sequences.items() if not should_exclude(acc)}
    removed = set(sequences) - set(filtered_sequences)
    if removed:
        print(f"\nExcluding {len(removed)} sample(s): {', '.join(sorted(removed))}")
    sequences = filtered_sequences
    original_length = len(list(sequences.values())[0]) if sequences else 0
    
    # Get species mapping (needed for all analyses)
    print("\n" + "="*80)
    print("LOADING SPECIES METADATA")
    print("="*80)
    accession_to_species, base_name_to_species = get_species_mapping(metadata_file)
    print(f"  Loaded metadata for {len(accession_to_species)} accessions")
    
    # Define merged regions (1-based coordinates)
    # Merged regions from recombinant region summary:
    mask_regions_list = [
        (1, 543),           # Events: 1, 2
        (111315, 113565),   # Events: 6
        (145378, 145488),   # Events: 4
        (151208, 152461)    # Events: 5, 1, 3
    ]
    # Set to None to disable masking, or provide a file path: mask_regions_file = 'merged_regions.txt'
    mask_regions_file = None
    
    # Parse regions (converts to 0-based)
    original_regions = []
    if mask_regions_list or mask_regions_file:
        print("\n" + "="*80)
        print("PARSING MERGED REGIONS")
        print("="*80)
        original_regions = parse_merged_regions(regions_file=mask_regions_file, regions_list=mask_regions_list)
    
    # Extract each masked region individually BEFORE masking/trimming for separate analysis
    if original_regions:
        print("\n" + "="*80)
        print("EXTRACTING MERGED REGIONS INDIVIDUALLY FOR SEPARATE ANALYSIS")
        print("="*80)
        
        # Create folder for extracted regions
        extracted_regions_dir = os.path.join(output_dir, 'extracted_regions')
        os.makedirs(extracted_regions_dir, exist_ok=True)
        
        # Process each region individually
        for region_idx, (start, end) in enumerate(original_regions, 1):
            region_name = f"region_{region_idx}_{start+1}_{end}"
            region_folder = os.path.join(extracted_regions_dir, region_name)
            os.makedirs(region_folder, exist_ok=True)
            
            print(f"\n  Processing {region_name} (position {start+1}-{end}, {end-start} bp)")
            
            # Extract this single region
            single_region = [(start, end)]
            extracted_region_sequences = extract_regions(sequences, single_region)
            
            if extracted_region_sequences:
                # Compute pairwise similarities on this extracted region
                print(f"\n  Computing pairwise similarities for {region_name}...")
                region_global_df, region_pairwise_df = compute_global_pairwise_similarities(
                    extracted_region_sequences, accession_to_species, base_name_to_species)
                
                # Save this region's analysis
                region_global_df.to_csv(os.path.join(region_folder, 'global_pairwise_similarities.csv'), index=False)
                print(f"    Saved: {region_folder}/global_pairwise_similarities.csv")
                region_pairwise_df.to_csv(os.path.join(region_folder, 'global_pairwise_sample_pairs.csv'), index=False)
                print(f"    Saved: {region_folder}/global_pairwise_sample_pairs.csv")
                
                # Save summary statistics
                save_global_similarities_summary(region_global_df, 
                                               os.path.join(region_folder, 'summary_statistics.csv'))
                
                # Plot this region's similarities
                plot_global_similarities(region_global_df, 
                                       output_file=os.path.join(region_folder, 'global_pairwise_similarities.png'))
                plot_species_distance_bars(region_global_df, 
                                         output_file=os.path.join(region_folder, 'global_pairwise_distances.png'))
                plot_global_distance_distribution(region_global_df, 
                                                 output_file=os.path.join(region_folder, 'global_pairwise_distance_distribution.png'))
                
                print(f"    Plots saved in: {region_folder}/")
    
    # Mask merged regions FIRST (before trimming)
    if original_regions:
        print("\n" + "="*80)
        print("MASKING MERGED REGIONS")
        print("="*80)
        sequences = mask_regions(sequences, original_regions, mask_char='N')
    
    # Trim sequences: remove first 3000bp and last 3000bp (AFTER masking)
    print("\n" + "="*80)
    print("TRIMMING SEQUENCES")
    print("="*80)
    trim_start = 3000
    trim_end = 3000
    sequences = trim_sequences(sequences, trim_start=trim_start, trim_end=trim_end)
    
    # Adjust masked regions coordinates for plotting (after trimming)
    masked_regions_for_plotting = []
    if original_regions:
        masked_regions_for_plotting = adjust_regions_after_trim(original_regions, trim_start, trim_end, original_length)
    
    # Compute global pairwise similarities
    print("\n" + "="*80)
    print("COMPUTING GLOBAL PAIRWISE SIMILARITIES")
    print("="*80)
    global_df, pairwise_df = compute_global_pairwise_similarities(sequences, accession_to_species, base_name_to_species)
    
    # Save global similarities
    global_df.to_csv(os.path.join(output_dir, 'global_pairwise_similarities.csv'), index=False)
    print(f"\nGlobal similarities saved to: {output_dir}/global_pairwise_similarities.csv")
    # Save pairwise (sample-level) table
    pairwise_df.to_csv(os.path.join(output_dir, 'global_pairwise_sample_pairs.csv'), index=False)
    print(f"Global pairwise sample table saved to: {output_dir}/global_pairwise_sample_pairs.csv")
    
    # Save summary statistics
    save_global_similarities_summary(global_df, os.path.join(output_dir, 'global_pairwise_similarities_summary.csv'))
    
    # Plot global similarities
    plot_global_similarities(global_df, output_file=os.path.join(output_dir, 'global_pairwise_similarities.png'))
    plot_species_distance_bars(global_df, output_file=os.path.join(output_dir, 'global_pairwise_distances.png'))
    plot_global_distance_distribution(global_df, output_file=os.path.join(output_dir, 'global_pairwise_distance_distribution.png'))
    
    # Compute sliding window similarities
    print("\n" + "="*80)
    print("COMPUTING SLIDING WINDOW SIMILARITIES")
    print("="*80)
    window_df = compute_sliding_window_similarities(sequences, accession_to_species, base_name_to_species,
                                                   window_size=500, step_size=100)
    
    # Save sliding window similarities (includes window positions)
    # Reorder columns to put chunk information first
    column_order = ['chunk_number', 'window_start', 'window_end', 'window_center', 'window_size', 
                   'chunk_position', 'species_pair', 'mean_similarity', 'std_similarity', 'n_comparisons']
    # Only include columns that exist
    column_order = [col for col in column_order if col in window_df.columns]
    window_df = window_df[column_order]
    
    window_df.to_csv(os.path.join(output_dir, 'sliding_window_pairwise_similarities.csv'), index=False)
    print(f"\nSliding window similarities saved to: {output_dir}/sliding_window_pairwise_similarities.csv")
    print("  CSV includes: chunk_number, window_start, window_end, window_center, window_size, chunk_position, species_pair, mean_similarity, std_similarity, n_comparisons")
    
    # Save summary statistics
    save_window_similarities_summary(window_df, os.path.join(output_dir, 'sliding_window_pairwise_similarities_summary.csv'))
    
    # Plot sliding window similarities (with masked regions indicated)
    plot_sliding_window_similarities(window_df, 
                                    output_file=os.path.join(output_dir, 'sliding_window_pairwise_similarities.png'),
                                    masked_regions=masked_regions_for_plotting)
    
    # Plot sliding window heatmap (with masked regions indicated)
    plot_sliding_window_heatmap(window_df, 
                                output_file=os.path.join(output_dir, 'sliding_window_similarities_heatmap.png'),
                                masked_regions=masked_regions_for_plotting)
    
    # Cluster windows based on similarity profiles
    cluster_sliding_windows(window_df, n_clusters=4, random_state=42,
                            output_csv=os.path.join(output_dir, 'sliding_window_clusters.csv'),
                            output_plot=os.path.join(output_dir, 'sliding_window_clusters.png'))
    
    # PCA of window regions
    plot_window_regions_pca(window_df, output_plot=os.path.join(output_dir, 'sliding_window_regions_pca.png'))
    
    # Build sample similarity/distance matrices and PCA
    build_sample_distance_outputs(pairwise_df,
                                  sim_matrix_csv=os.path.join(output_dir, 'global_pairwise_similarity_matrix.csv'),
                                  dist_matrix_csv=os.path.join(output_dir, 'global_pairwise_distance_matrix.csv'),
                                  pca_plot=os.path.join(output_dir, 'global_pairwise_pca.png'))
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print("\nSliding Window Similarities by Species Pair:")
    for pair in sorted(window_df['species_pair'].unique()):
        pair_data = window_df[window_df['species_pair'] == pair]
        print(f"\n  {pair}:")
        print(f"    Mean similarity across windows: {pair_data['mean_similarity'].mean():.2f}%")
        print(f"    Std deviation: {pair_data['mean_similarity'].std():.2f}%")
        print(f"    Min: {pair_data['mean_similarity'].min():.2f}%")
        print(f"    Max: {pair_data['mean_similarity'].max():.2f}%")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)

if __name__ == "__main__":
    main()

