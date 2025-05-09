# ctcf_boundaries/cli.py

import argparse
import pandas as pd
import bioframe
import cooler
from cooltools import insulation
from bioframe.ops import _verify_columns
from bioframe import ops

def make_chromarms(chromsizes, midpoints, cols_chroms=("chrom", "length"), cols_mids=("chrom", "mid"), suffixes=("_p", "_q")):
    columns_to_drop = ["index", "sub_index_"]
    if len(cols_chroms) == 2:
        ck1, sk1 = cols_chroms
    elif len(cols_chroms) == 3:
        ck1, sk1, ek1 = cols_chroms

    if isinstance(chromsizes, (pd.Series, dict)):
        chromsizes = dict(chromsizes)
        df_chroms = pd.DataFrame({ck1: list(chromsizes.keys()), "length": list(chromsizes.values())})
    elif isinstance(chromsizes, pd.DataFrame):
        df_chroms = chromsizes.copy()
    else:
        raise ValueError("unknown input type for chromsizes")

    if len(cols_chroms) == 2:
        _verify_columns(df_chroms, [ck1, sk1])
        columns_to_drop += [sk1]
        df_chroms["end"] = df_chroms[sk1].values
        df_chroms["start"] = 0
        sk1, ek1 = "start", "end"
    elif len(cols_chroms) == 3:
        ck1, sk1, ek1 = cols_chroms
        _verify_columns(df_chroms, [ck1, sk1, ek1], unique_cols=True)
        if any(df_chroms[sk1].values != 0):
            raise ValueError("all values in starts column must be zero")

    ck2, sk2 = cols_mids
    if isinstance(midpoints, (pd.Series, dict)):
        midpoints = dict(midpoints)
        df_mids = pd.DataFrame.from_dict(midpoints, orient="index", columns=[sk2])
        df_mids.reset_index(inplace=True)
        df_mids.rename(columns={"index": ck2}, inplace=True)
    elif isinstance(midpoints, pd.DataFrame):
        df_mids = midpoints.copy()
    else:
        raise ValueError("unknown input type for midpoints")
    _verify_columns(df_mids, [ck2, sk2])
    df_mids["start"] = df_mids[sk2]
    df_mids["end"] = df_mids[sk2]

    df_chromarms = ops.subtract(
        df_chroms, df_mids, cols1=(ck1, sk1, ek1), cols2=(ck2, "start", "end"), return_index=True
    )
    if df_chromarms["sub_index_"].max() > 1:
        raise ValueError("chromosome split into more than two arms, double-check midpoints")
    df_chromarms["name"] = df_chromarms[ck1] + [
        suffixes[int(i)] for i in df_chromarms["sub_index_"].values
    ]
    return df_chromarms[[ck1, sk1, ek1, "name"]]


def init_files(cool_file, peaks_file=None, resolution=10000, genome='hg38'):
    import csv

    try:
        clr = cooler.Cooler(f'{cool_file}::/resolutions/{resolution}')
    except:
        clr = cooler.Cooler(cool_file)

    # Fetch genome arms
    genome_chromsizes = bioframe.fetch_chromsizes(genome)
    genome_cens = bioframe.fetch_centromeres(genome)
    genome_arms = make_chromarms(genome_chromsizes, genome_cens)

    #    Reorder to match the cooler’s chromosomes exactly
    genome_arms = genome_arms.set_index("chrom").loc[clr.chromnames].reset_index()

    # Validate as viewframe
    genome_arms = bioframe.make_viewframe(genome_arms)

    peaks = None
    if peaks_file:
        with open(peaks_file, 'r') as f:
            header = next(csv.reader(f, delimiter="\t"))

        required_cols = {'chrom', 'start', 'end'}
        if not required_cols.issubset(set(header)):
            raise ValueError(f"BED file must contain at least: {required_cols}. Found: {header}")

        peaks = bioframe.read_table(peaks_file, schema=None, names=header, header=0)
        peaks = peaks.query("chrom in @clr.chromnames")
        peaks['mid'] = (peaks.end + peaks.start) // 2

    return clr, genome_arms, peaks


def call_boundaries(clr, view, resolution=10000, nproc=4):
    window = resolution * 10
    return insulation(clr, window_bp=[window], view_df=view, nproc=nproc)


def bin_overlaps(df, overlapper, binsize=10000, colnames=["chrom", "start", "end"], colnames_overlapper=["chrom", "start", "end"], added_column="is_transcribed"):
    df = df.copy()
    df["bin"] = df[colnames[1]] // binsize

    overlapper_bins = (
        overlapper[[colnames_overlapper[0], colnames_overlapper[1]]]
        .assign(bin=lambda x: x[colnames_overlapper[1]] // binsize)
        .rename(columns={colnames_overlapper[0]: "chrom"})
    )
    overlapper_bins_end = (
        overlapper[[colnames_overlapper[0], colnames_overlapper[2]]]
        .assign(bin=lambda x: x[colnames_overlapper[2]] // binsize)
        .rename(columns={colnames_overlapper[0]: "chrom", colnames_overlapper[2]: "start"})
    )
    overlapper_bins = pd.concat([overlapper_bins, overlapper_bins_end], ignore_index=True).drop_duplicates()

    df["key"] = list(zip(df["chrom"], df["bin"]))
    overlapper_bins["key"] = list(zip(overlapper_bins["chrom"], overlapper_bins["bin"]))

    df[added_column] = df["key"].isin(set(overlapper_bins["key"]))
    df.drop(columns=["bin", "key"], inplace=True)

    return df


def generate_ctcf_boundaries(insulation_df, peaks, out_path, resolution=10000):
    window = resolution * 10
    ctcf_boundaries = bin_overlaps(insulation_df, peaks, binsize=resolution, colnames=["chrom", "start", "end"], added_column="has_ctcf")
    ctcf_boundaries = ctcf_boundaries[ctcf_boundaries["has_ctcf"]]
    ctcf_boundaries = ctcf_boundaries[ctcf_boundaries[f"is_boundary_{window}"]]
    ctcf_boundaries.to_csv(out_path, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Call CTCF boundaries from a .cool file and peak annotations.")
    parser.add_argument("cool_file", help="Path to the .cool file")
    parser.add_argument("ctcf_peaks", help="Path to BED file with CTCF peaks")
    parser.add_argument("out_path", help="Path to output file for boundaries")
    parser.add_argument("--resolution", type=int, default=10000, help="Resolution of .cool file (default: 10000)")
    parser.add_argument("--nproc", type=int, default=4, help="Number of processors to use (default: 4)")
    parser.add_argument("--genome", default="hg38", help="Genome assembly to use (default: hg38)")

    args = parser.parse_args()

    clr, view, peaks = init_files(args.cool_file, args.ctcf_peaks, resolution=args.resolution, genome=args.genome)
    insulation_df = call_boundaries(clr, view, resolution=args.resolution, nproc=args.nproc)
    generate_ctcf_boundaries(insulation_df, peaks, args.out_path, resolution=args.resolution)


if __name__ == "__main__":
    main()