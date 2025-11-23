#!/usr/bin/env python3


"""
Functions to parse FastQC output reports and produce tables.
"""


from pathlib import Path
from typing import TypeAlias, Dict, List, Union

import math
import os
import pandas as pd
import zipfile


def extract_module(lines, module_name: str):
    """
    Extracts and parses a single block in fastqc_data.txt
    Blocks begin with ">>module_name" and end with "<<END_MODULE"
    """
    # Defines start and end of the module
    start_tag = f">>{module_name}"
    end_tag = ">>END_MODULE"
    # Boolean to check if lines are inside the module
    inside = False
    # List of lines corresponding to the module
    block = []
    # Go through lines
    for line in lines:
        # Remove beggining and end spaces/end-of-line characters
        line = line.strip()
        # Check if line starts with the start tag
        if line.startswith(start_tag):
            # Check if inside is already True
            if inside:
                print('WARNING: start_tag found twice in a row.')
            else:
                # Mark that we are inside the module
                inside = True
            continue
        # Stop the loop if we find the end tag after seeing the start tag
        if inside and line.startswith(end_tag):
            break
        # While inside the module, save lines to block
        if inside:
            block.append(line)

    # remove header lines starting with '#'
    block_no_comments = [l for l in block if not l.startswith("#")]
    return block_no_comments


def parse_fastqc_directory(
    directory: Path|str,
    file_suffix: str = "_fastqc.zip"
) -> pd.DataFrame:
    """
    Parses all FastQC .zip reports in a directory and generates 
    a table with them.

    Output: DataFrame with columns:
      sample, mean_quality, median_quality, GC_mean, read_length_min,
      read_length_max, read_length_mean, read_length_median
    """
    # Initialize the output list
    rows = []
    # Go through elements in directory
    for f in os.listdir(directory):
        if f.endswith(file_suffix):
            zip_path = os.path.join(directory, f)
            try:
                parsed = parse_fastqc_zip(zip_path)

                # Summary metrics
                pbq = parsed["per_base_quality"]
                psq = parsed["per_sequence_quality"]
                gc = parsed["per_sequence_gc"]
                rl = parsed["read_length_summary"]

                sample = os.path.splitext(os.path.basename(f))[0]

                # Calculate summary values
                mean_quality = (pbq["mean"] * 1).mean()
                median_quality = pbq["median"].median()
                gc_mean = (gc["gc_percent"] * gc["count"]).sum() / \
                    gc["count"].sum()

                rl_min = rl['min']
                rl_max = rl['max']
                rl_mean = rl['mean']
                rl_median = rl['median']

                rows.append([sample, mean_quality, median_quality, 
                             gc_mean, rl_min, rl_max, rl_mean, rl_median])

            except Exception as e:
                raise ValueError(f"ERROR parsing {f}: {e}")

    return pd.DataFrame(
        rows,
        columns=["sample", "mean_quality", "median_quality", "GC_mean", 
                 "read_length_min", "read_length_max", "read_length_mean", 
                 "read_length_median"],
    )


def parse_fastqc_zip(
    zip_path: str|Path,
    report_name: str = "fastqc_data.txt"
) -> Dict[str, Union[pd.DataFrame, Dict[str, Union[int, float]]]]:
    """
    Parses fastqc_data.txt inside a FastQC .zip report.

    Returns a dict with:
      - per_base_quality : DataFrame
      - per_sequence_quality : DataFrame
      - per_sequence_gc : DataFrame
      - read_length : Dict[str : int|float]
    """
    # Open the fastqc report .zip file
    with zipfile.ZipFile(zip_path, "r") as z:
        # Find report_name inside the zip
        fastqc_files = [f for f in z.namelist() if f.endswith(report_name)]
        # Raise error if report_name is not found
        if not fastqc_files:
            raise ValueError(f"No {report_name} found in: {zip_path}")
        # Open the first file found with report_name
        with z.open(fastqc_files[0]) as f:
            lines = f.read().decode("utf-8").splitlines()

    try:
        # Per base quality
        block = extract_module(lines, "Per base sequence quality")
        cols = ["base", "mean", "median", "lower_quartile", "upper_quartile", "percentile_10", "percentile_90"]

        clean_rows = []
        for i, line in enumerate(block):
            if not line or line.strip().startswith("#"):
                # skip empty or commented lines
                continue
            parts = line.split("\t")
            if len(parts) != len(cols):
                # skip lines that don't have exactly the expected number of columns
                continue
            base = parts[0].strip()
            numeric_parts = parts[1:]
            # Validate numeric columns explicitly so we fail only on the offending line
            try:
                numeric_vals = [float(x) for x in numeric_parts]
            except ValueError:
                # skip this line if any numeric conversion fails
                continue
            # keep the row: base as string, numeric columns as floats
            clean_rows.append([base] + numeric_vals)

        # Build DataFrame from cleaned rows (preserves your column names)
        per_base_quality = pd.DataFrame(clean_rows, columns=cols)
    except Exception as e:
        raise ValueError(f"Error parsing Per base sequence quality: {e}")

    try:
        # Per sequence quality scores
        block = extract_module(lines, "Per sequence quality scores")
        clean_gc = []
        for l in block:
            parts = l.split("\t")
            if len(parts) != 2:
                continue
            a, b = parts
            # Try to parse both columns as integers
            if a.strip().isdigit() and b.strip().isdigit():
                clean_gc.append([int(a), int(b)])
        per_sequence_quality = pd.DataFrame(
            clean_gc,
            columns=["quality", "count"],
        )
    except Exception as e:
        raise ValueError(f"Error parsing Per sequence quality scores: {e}")

    try:
        # Per sequence GC content
        block = extract_module(lines, "Per sequence GC content")
        clean_gc = []
        for l in block:
            parts = l.split("\t")
            if len(parts) != 2:
                continue
            a, b = parts
            # Try to parse both columns as integers
            try:
                clean_gc.append([int(float(a)), int(float(b))])
            except Exception as e:
                raise ValueError(f"Error parsing GC content line '{l}': {e}")
        per_sequence_gc = pd.DataFrame(
            clean_gc,
            columns=["gc_percent", "count"],
        )
    except Exception as e:
        raise ValueError(f"Error parsing Per sequence GC content: {e}")

    try:
        # Sequence length distribution
        block = extract_module(lines, "Sequence Length Distribution")
        if not block:
            # Module missing or empty: keep the same output variables but empty
            length_df = pd.DataFrame(columns=["length", "count", "len_min", "len_max", "cum"])
            read_length_summary = {"min": None, "max": None, "mean": None, "median": None}
        else:
            # parse rows: expects "length<TAB>count"
            length_rows = [l.split("\t") for l in block if l.strip()]
            length_df = pd.DataFrame(length_rows, columns=["length", "count"])
            # Apply parse_range to the length column
            parsed = length_df["length"].apply(lambda s: pd.Series(parse_range(s)))
            parsed.columns = ["len_min", "len_max"]
            # Add parsed columns and set their types
            length_df = pd.concat([length_df, parsed], axis=1)
            length_df["count"] = length_df["count"].astype(float).astype(int)
            length_df["len_min"] = length_df["len_min"].astype(int)
            length_df["len_max"] = length_df["len_max"].astype(int)
            # Summary min/max
            try:
                read_length_min = int(length_df["len_min"].min()) if not length_df.empty else None
                read_length_max = int(length_df["len_max"].max()) if not length_df.empty else None
            except Exception as e:
                raise ValueError(f"Error calculating min/max read lengths: {e}")

            try:
                # weighted mean: assumes uniform distribution within each bin
                total_count = int(length_df["count"].sum())
            except Exception as e:
                raise ValueError(f"Error calculating total count for mean read length: {e}")
            # If there are 0 counts, there is no mean
            if total_count == 0:
                mean_read_length = None
            else:
                # Calculate the mean value per bin
                mean_per_bin = (length_df["len_min"] + length_df["len_max"]) / 2.0
                mean_read_length = \
                    float((mean_per_bin * length_df["count"]).sum() / total_count)

            # median: find the bin containing the cumulative half and interpolate inside it
            if total_count == 0:
                # If there are 0 counts, there is no median
                median_read_length = None
            else:
                # Ensure sorted and cumulative column
                length_df = length_df.sort_values(
                    ["len_min", "len_max"]
                ).reset_index(drop=True)
                length_df["cum"] = length_df["count"].cumsum()
                half = total_count / 2.0

                # If the last cumulative is less than half, that's inconsistent data
                if length_df["cum"].iloc[-1] < half:
                    raise ValueError(
                        f"Inconsistent cumulative counts: total_count={total_count} "
                        f"but final cumulative={length_df['cum'].iloc[-1]}"
                    )

                # Find the first bin whose cumulative count reaches/exceeds the median position
                try:
                    idx = int(length_df["cum"].searchsorted(half))
                except Exception as e:
                    raise ValueError(f"Failed to locate median bin: {e}")

                bin_row = length_df.iloc[idx]
                bin_count = int(bin_row["count"])
                if bin_count == 0:
                    raise ValueError(
                        f"Median lies in a bin with zero count at index {idx}, "
                        f"bin={bin_row.to_dict()}"
                    )

                len_min = int(bin_row["len_min"])
                len_max = int(bin_row["len_max"])
                n_ints = len_max - len_min + 1
                if n_ints <= 0:
                    raise ValueError(
                        f"Invalid bin range len_min={len_min}, len_max={len_max}"
                    )

                cum_prev = int(length_df["cum"].iloc[idx-1]) if idx > 0 else 0
                # how many counts into the bin are needed
                need = half - cum_prev

                # counts per integer length inside the bin (uniform assumption)
                counts_per_int = float(bin_count) / float(n_ints)
                if counts_per_int <= 0:
                    raise ValueError(
                        f"Non-positive counts_per_int for bin at index {idx}, "
                        f"bin={bin_row.to_dict()}"
                    )

                integers_needed = math.ceil(need / counts_per_int)
                # If this exceeds bin size, something is inconsistent
                if integers_needed < 1 or integers_needed > n_ints:
                    raise ValueError(
                        f"Median integer index out of bin range: need={need}, "
                        f"integers_needed={integers_needed}, n_ints={n_ints}, "
                        f"bin={bin_row.to_dict()}"
                    )
                # Compute the actual median integer length
                median_read_length = len_min + (integers_needed - 1)
                median_read_length = int(median_read_length)
            read_length_summary: Dict[str, Union[int, float]] = {
                "min": read_length_min,
                "max": read_length_max,
                "mean": mean_read_length,
                "median": median_read_length,
            }
    except Exception as e:
        raise ValueError(f"Error parsing Sequence Length Distribution: {e}")

    return {
        "per_base_quality": per_base_quality,
        "per_sequence_quality": per_sequence_quality,
        "per_sequence_gc": per_sequence_gc,
        "read_length_summary": read_length_summary
    }


def parse_range(s: str, sep: str='-') -> list[int, int]:
    """
    Parses a string in range format (val1-val2).
    Handles int format (val).
    """
    try:
        # Remove leading spaces and end of line characters
        s = s.strip()
        if sep in s:
            parts = [p.strip() for p in s.split(sep)]
            return int(parts[0]), int(parts[1])
        else:
            v = int(s)
            return v, v
    except Exception as e:
        raise ValueError(f"Error parsing range from string '{s}': {e}")
