#!/usr/bin/env python3
"""
alpha_variant_scan.py — Configurable & Stable AlphaGenome Variant Expression Scanner

need: alphagenome, pandas, numpy, matplotlib.
"""

from __future__ import annotations
import argparse, os, sys, warnings
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    P = argparse.ArgumentParser(
        prog="alpha_variant_scan",
        description="Scan predicted expression effects of variants across organs using AlphaGenome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    P.add_argument("--variants", required=True, help="TSV/CSV/VCF of variants.")
    P.add_argument("--organs", nargs="*", default=None, help="UBERON terms. Default=demo list.")
    P.add_argument("--threshold", type=float, default=0.5, help="|ALT/REF−1| score threshold.")
    P.add_argument("--min-length", dest="min_length", type=int, default=1000, help="Min bp length of merged region.")
    P.add_argument("--merge-distance", dest="merge_distance", type=int, default=300, help="Max bp gap to merge regions.")
    P.add_argument("--window-size", dest="window_size", type=int, default=100, help="Sliding window size (bp).")
    P.add_argument("--scan-span", dest="scan_span", type=int, default=50000, help="Bp each side of variant center to scan.")
    P.add_argument("--plot-non-sig", dest="plot_non_sig", action="store_true", help="Plot even when not significant.")
    P.add_argument("--scan-all-tracks", dest="scan_all_tracks", action="store_true", help="Scan all tracks (recommended).")
    P.add_argument("--epsilon", type=float, default=1e-8, help="Add to REF to avoid divide-by-zero.")
    P.add_argument("--output-table", default="alphagenome_scan_results.csv", help="Output summary table path (csv/tsv/xlsx).")
    P.add_argument("--output-dir", default="alphagenome_scan_plots", help="Dir for plot images.")
    P.add_argument("--api-key", dest="api_key", default=None, help="AlphaGenome API key; else from env.")
    P.add_argument("--gtf", default='https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather', help="GENCODE transcript annotation feather path/URL.")
    P.add_argument("--chrom-col", dest="col_chrom", default="CHROM", help="Column name: chromosome.")
    P.add_argument("--pos-col", dest="col_pos", default="POS", help="Column name: 1-based position.")
    P.add_argument("--ref-col", dest="col_ref", default="REF", help="Column name: reference bases.")
    P.add_argument("--alt-col", dest="col_alt", default="ALT", help="Column name: alternate bases.")
    return P

# ------------------------------------------------------------------
# I/O
# ------------------------------------------------------------------

def load_variants_table(path: str, col_chrom: str, col_pos: str, col_ref: str, col_alt: str) -> pd.DataFrame:
    fp = Path(path)
    if not fp.exists():
        raise FileNotFoundError(path)
    suf = fp.suffix.lower()
    if suf in {".tsv", ".txt"}:
        df = pd.read_csv(fp, sep="\t")
    elif suf == ".csv":
        df = pd.read_csv(fp)
    elif suf == ".vcf":
        df = _load_vcf(fp)
    else:
        # fallback: try tab then comma
        try:
            df = pd.read_csv(fp, sep="\t")
        except Exception:
            df = pd.read_csv(fp)
    mapping = {col_chrom: "CHROM", col_pos: "POS", col_ref: "REF", col_alt: "ALT"}
    miss = [c for c in mapping if c not in df.columns]
    if miss:
        raise ValueError(f"Missing required columns: {miss}")
    df = df.rename(columns=mapping)[["CHROM", "POS", "REF", "ALT"]]
    return df


def _load_vcf(fp: Path) -> pd.DataFrame:
    chrom = []
    pos = []
    ref = []
    alt = []
    with fp.open() as fh:
        for ln in fh:
            if ln.startswith('#'):
                continue
            f = ln.rstrip('\n').split('\t')
            chrom.append(f[0])
            pos.append(int(f[1]))
            ref.append(f[3])
            alt.append(f[4].split(',')[0])  # first ALT only
    return pd.DataFrame({"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt})


def ensure_output_dir(path: str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

# ------------------------------------------------------------------
# Model / Annotation
# ------------------------------------------------------------------

def load_transcript_extractor(gtf_path: str):
    gtf = pd.read_feather(gtf_path)
    gtf_t = gene_annotation.filter_protein_coding(gtf)
    gtf_t = gene_annotation.filter_to_longest_transcript(gtf_t)
    return transcript_utils.TranscriptExtractor(gtf_t)


def get_dna_model(api_key: Optional[str] = None):
    if api_key is None:
        try:
          api_key = colab_utils.get_api_key()
        except Exception as e:
        # If not Colab
          print(f"Error getting API key from Colab: {e}")
          api_key = ""  # Insert your API key here if not working
    return dna_client.create(api_key)

# ------------------------------------------------------------------
# Core computations
# ------------------------------------------------------------------

def align_reference_for_indel(variant, interval, vout, length_alter: int):
    """Shift REF track in-place to align with ALT for indels (original logic)."""
    if length_alter > 0:  # deletion
        vout.reference.rna_seq.values[
            (variant.position - interval.start):(interval.end - interval.start - length_alter)
        ] = vout.reference.rna_seq.values[
            (variant.position - interval.start + length_alter):(interval.end - interval.start)
        ]
        vout.reference.rna_seq.values[
            (interval.end - interval.start - length_alter):(interval.end - interval.start)
        ] = np.nan
    elif length_alter < 0:  # insertion
        vout.reference.rna_seq.values[
            (variant.position - interval.start - length_alter):(interval.end - interval.start)
        ] = vout.reference.rna_seq.values[
            (variant.position - interval.start):(interval.end - interval.start + length_alter)
        ]
        vout.reference.rna_seq.values[
            (variant.position - interval.start):(variant.position - interval.start - length_alter)
        ] = np.nan
    # SNV => no shift


def compute_window_scores(alt_vals: np.ndarray, ref_vals: np.ndarray, start_idx: int, end_idx: int,
                          window_size: int, epsilon: float) -> np.ndarray:
    """Return (n_windows, n_tracks) of ALT/REF−1 window means."""
    n_bases, n_tracks = alt_vals.shape
    if start_idx < 0: start_idx = 0
    if end_idx > n_bases: end_idx = n_bases
    n_windows = end_idx - start_idx - window_size + 1
    if n_windows <= 0:
        return np.empty((0, n_tracks))
    # cumulative sums (prepend 0 row)
    alt_cs = np.cumsum(np.vstack([np.zeros((1, n_tracks)), alt_vals]), axis=0)
    ref_cs = np.cumsum(np.vstack([np.zeros((1, n_tracks)), ref_vals]), axis=0)
    # slice rolling windows
    a = alt_cs[start_idx + window_size:start_idx + window_size + n_windows] - alt_cs[start_idx:start_idx + n_windows]
    r = ref_cs[start_idx + window_size:start_idx + window_size + n_windows] - ref_cs[start_idx:start_idx + n_windows]
    a /= float(window_size)
    r /= float(window_size)
    return (a / (r + epsilon)) - 1.0


def call_regions(scores: np.ndarray, threshold: float, min_length: int, merge_distance: int) -> List[Tuple[int,int]]:
    """Identify high-score regions (inclusive window indices)."""
    idx = np.where(np.abs(scores) > threshold)[0]
    if idx.size == 0:
        return []
    # compress consecutive indices
    regions = []
    start = prev = idx[0]
    for i in idx[1:]:
        if i == prev + 1:
            prev = i
        else:
            regions.append((start, prev))
            start = prev = i
    regions.append((start, prev))  # tail
    # merge
    merged = []
    cur_s, cur_e = regions[0]
    for s, e in regions[1:]:
        if s - cur_e <= merge_distance:
            cur_e = e
        else:
            if cur_e - cur_s + 1 >= min_length:
                merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    if cur_e - cur_s + 1 >= min_length:
        merged.append((cur_s, cur_e))
    return merged

# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------

def plot_variant_tracks(variant, interval, vout, transcript_extractor, plot_size: int, outpath: Path):
    longest_transcripts = transcript_extractor.extract(interval)
    plot_components.plot([
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.OverlaidTracks(
            tdata={'REF': vout.reference.rna_seq, 'ALT': vout.alternate.rna_seq},
            colors={'REF': 'blue', 'ALT': 'red'},
        ),
    ], interval=vout.reference.rna_seq.interval.resize(plot_size),
       annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)])
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()


def plot_scores(scores: np.ndarray, regions: List[Tuple[int,int]], threshold: float, outpath: Path, title: str):
    plt.figure(figsize=(8,3))
    plt.plot(scores, label='score')
    plt.axhline(threshold, color='grey', lw=0.5, ls='--')
    plt.axhline(-threshold, color='grey', lw=0.5, ls='--')
    for s,e in regions:
        plt.axvspan(s, e, color='red', alpha=0.3)
    plt.xlabel('window index'); plt.ylabel('ALT/REF−1'); plt.title(title)
    plt.legend(loc='upper right', fontsize='x-small')
    plt.tight_layout(); plt.savefig(outpath, dpi=300); plt.close()

# ------------------------------------------------------------------
# Table writer
# ------------------------------------------------------------------

def write_table(df: pd.DataFrame, path: str):
    ext = Path(path).suffix.lower()
    if ext in {'.tsv', '.txt'}:
        df.to_csv(path, sep='\t', index=False)
    elif ext in {'.xlsx', '.xls'}:
        df.to_excel(path, index=False)
    else:
        # default csv
        df.to_csv(path, index=False)

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    # load inputs
    variants_df = load_variants_table(args.variants, args.col_chrom, args.col_pos, args.col_ref, args.col_alt)
    transcript_extractor = load_transcript_extractor(args.gtf)
    dna_model = get_dna_model(args.api_key)

    organs = args.organs or [
        'UBERON:0000992', 'UBERON:0002371', 'UBERON:0000948',
        'UBERON:0000955', 'UBERON:0001134', 'UBERON:0001264']

    out_dir = ensure_output_dir(args.output_dir)

    results_rows = []
    
    for ontology in organs:
        number_rank = 0
        for _, row in variants_df.iterrows():
            number_rank += 1
            variant = genome.Variant(
                chromosome=row['CHROM'], position=row['POS'],
                reference_bases=row['REF'], alternate_bases=row['ALT'])
        
            # 1Mb centered interval
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        
            vout = dna_model.predict_variant(
                interval=interval,
                variant=variant,
                requested_outputs=[dna_client.OutputType.RNA_SEQ],
                ontology_terms=[ontology],
            )
        
            # track count
            alt_vals = vout.alternate.rna_seq.values
            ref_vals = vout.reference.rna_seq.values
            n_tracks = alt_vals.shape[1]
            if n_tracks == 0:
                warnings.warn(f"No tracks available for {ontology}; skipping variant {variant}.")
                continue
        
            # indel alignment
            length_alter = len(variant.reference_bases) - len(variant.alternate_bases)
            if length_alter != 0:
                align_reference_for_indel(variant, interval, vout, length_alter)
                ref_vals = vout.reference.rna_seq.values  # updated
        
            # scan window bounds (0-based indexes into interval arrays)
            center_idx = variant.position - interval.start
            start_idx = center_idx - args.scan_span
            end_idx = center_idx + args.scan_span  # exclusive after adjust in compute
        
            # compute window scores (vectorized)
            win_scores = compute_window_scores(alt_vals, ref_vals, start_idx, end_idx, args.window_size, args.epsilon)
            if win_scores.size == 0:
                warnings.warn(f"Scan span/window too large/small near edges for {variant}.")
                continue
        
            # track names
            try:
                track_names = vout.reference.rna_seq.metadata.name
            except Exception:
                track_names = [f"track_{i}" for i in range(n_tracks)]
        
            # scan all tracks always; early-stop if user did NOT set scan_all_tracks? We'll still scan for table,
            # but we skip plotting after first significant if flag not set (speed).
            first_sig_found = False
            for ti in range(n_tracks):
                scores = win_scores[:, ti]
                regions = call_regions(scores, args.threshold, args.min_length, args.merge_distance)
                is_sig = len(regions) > 0
        
                # gather rows (one per region; if none, n_regions=0 + region_index=-1)
                if is_sig:
                    for ri, (rs, re) in enumerate(regions):
                        # convert window index -> bp relative to variant
                        rel_start = (start_idx + rs) - center_idx
                        rel_end = (start_idx + re + args.window_size - 1) - center_idx
                        abs_start = variant.position + rel_start
                        abs_end = variant.position + rel_end
                        seg_scores = scores[rs:re+1]
                        mean_s = float(np.nanmean(seg_scores))
                        max_s = float(np.nanmax(seg_scores))
                        min_s = float(np.nanmin(seg_scores))
                        if mean_s > 0 and min_s > 0:
                            direction = 'up'
                        elif mean_s < 0 and max_s < 0:
                            direction = 'down'
                        else:
                            direction = 'mixed'
                        results_rows.append({
                            'chrom': variant.chromosome,
                            'pos': variant.position,
                            'ref': variant.reference_bases,
                            'alt': variant.alternate_bases,
                            'ontology': ontology,
                            'track_name': track_names[ti] if ti < len(track_names) else f'track_{ti}',
                            'is_significant': True,
                            'n_regions': len(regions),
                            'region_index': ri,
                            'rel_start_bp': int(rel_start),
                            'rel_end_bp': int(rel_end),
                            'abs_start_bp': int(abs_start),
                            'abs_end_bp': int(abs_end),
                            'mean_score': mean_s,
                            'max_score': max_s,
                            'min_score': min_s,
                            'direction': direction,
                            'plot_file': None,  # filled after plotting variant-level below
                        })
                else:
                    results_rows.append({
                        'chrom': variant.chromosome,
                        'pos': variant.position,
                        'ref': variant.reference_bases,
                        'alt': variant.alternate_bases,
                        'ontology': ontology,
                        'track_name': track_names[ti] if ti < len(track_names) else f'track_{ti}',
                        'is_significant': False,
                        'n_regions': 0,
                        'region_index': -1,
                        'rel_start_bp': np.nan,
                        'rel_end_bp': np.nan,
                        'abs_start_bp': np.nan,
                        'abs_end_bp': np.nan,
                        'mean_score': float(np.nanmean(scores)),
                        'max_score': float(np.nanmax(scores)),
                        'min_score': float(np.nanmin(scores)),
                        'direction': 'none',
                        'plot_file': None,
                    })
        
                if is_sig and (not args.scan_all_tracks) and (not first_sig_found):
                    first_sig_found = True
                    # we will still scan others for table (already done!), but can skip extra heavy plotting of scores if desired
        
            # decide plot size (reuse original heuristic)
            if abs(length_alter) >= 2**14:
                plot_size = abs(length_alter) * 4
            else:
                plot_size = 2**15
        
            # plot organ-level REF/ALT overlay once per variant-organ if any sig OR user asked plot_non_sig
            sig_any = any((r['is_significant'] and r['chrom']==variant.chromosome and r['pos']==variant.position and r['ontology']==ontology) for r in results_rows[-n_tracks:])
            if sig_any or args.plot_non_sig:
                plot_path = out_dir / (f"{ontology}_{number_rank}_{variant.chromosome}_{variant.position}.png")
                plot_variant_tracks(variant, interval, vout, transcript_extractor, plot_size, plot_path)
                # back-fill plot_file for recent rows
                for r in results_rows[-n_tracks:]:
                    if r['chrom']==variant.chromosome and r['pos']==variant.position and r['ontology']==ontology:
                        r['plot_file'] = str(plot_path)

    # --------------------------
    # build DataFrame & write
    # --------------------------
    if not results_rows:
        print("No results generated.")
        return 0

    res_df = pd.DataFrame(results_rows)
    write_table(res_df, args.output_table)

    # variant×organ summary
    agg = (res_df.groupby(['chrom','pos','ref','alt','ontology'])['is_significant']
                 .any().reset_index().rename(columns={'is_significant':'is_significant_any'}))
    # track list for each variant×organ
    sig_tracks = (res_df[res_df.is_significant]
                    .groupby(['chrom','pos','ref','alt','ontology'])['track_name']
                    .apply(lambda s: ','.join(sorted(set(s))))
                    .reset_index())
    agg = agg.merge(sig_tracks, how='left', on=['chrom','pos','ref','alt','ontology'])
    agg['track_name'] = agg['track_name'].fillna('')
    # write
    agg_path = Path(args.output_table).with_name(Path(args.output_table).stem + '_variant_organ_summary' + Path(args.output_table).suffix)
    write_table(agg, str(agg_path))

    print(f"Wrote detailed results to {args.output_table}")
    print(f"Wrote variant×organ summary to {agg_path}")
    print("Done.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
