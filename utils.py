#!/usr/bin/env python3


"""
Utility scripts to be used in various parts of the project.
"""


from pathlib import Path


def sample_id_from_folder(
    dir_path: Path,
    cut_string: str = "_R",
    file_extension: str = ".fastq.gz"
) -> dict[str: list[str]]:
    """
    Extracts sample IDs from FASTQ file names in the given directory.

    Args:
        dir_path (Path): Path to the directory containing FASTQ files.
        cut_string (str): String to cut the sample ID from the file name.
        file_extension (str): Extension of the FASTQ files to consider.
    Returns:
        dict: A dictionary with sample IDs as keys and a list of fastq 
              file names as values.
    
    """
    # First pass: extract unique sample IDs
    sample_ids = set()
    # Get all files within dir_path
    for file in dir_path.iterdir():
        if file.is_file() and file.name.endswith(file_extension):
            sample_id = file.name.split(cut_string)[0]
            sample_ids.add(sample_id)
    # Second pass: associate files with sample IDs
    sample_dict = {}
    for sample_id in sample_ids:
        associated_files = []
        for file in dir_path.iterdir():
            if file.is_file() and file.name.endswith(file_extension):
                if file.name.startswith(sample_id+str(cut_string)):
                    associated_files.append(file.name)
        sample_dict[sample_id] = associated_files
    return sample_dict

