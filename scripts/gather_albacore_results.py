#!/usr/bin/env python
import logging
from glob import glob
import os
import shutil
from collections import defaultdict
from typing import List, Tuple, Dict
import pandas as pd


def get_fastq_files(input_dirs: List[str], search_dir: str) -> List[str]:
    """Traverse through input directories and make a list of all fastq files.

    Traverse though `<input_dir><search_dir>` and get all passing fastq file
    paths. All file paths are placed in `fastq_files`, which is globally
    accessible.

    :param input_dirs:
    :param search_dir:
    :return:
    """
    fastq_files = []
    for input_dir in input_dirs:
        search_path = os.path.join(input_dir, search_dir, "**", '*.fastq')
        fastq_files += glob(search_path, recursive=True)
    return fastq_files


def cat_fastq_files(fastq_files: List[str], output_dir: str, subdir: str) -> None:
    """Concatenate all fastq files with the same run ID and barcode.

    Firstly, create sets of all barcodes and run IDs that the files in
    `fastq_files` are a part of. Then iterate through each combination of
    barcode and run ID, and write all files that match into a fastq file,
    effectively combining all fastq files that have the same barcode and run ID.

    :param fastq_files:
    :param output_dir:
    :param subdir:
    :return:
    """
    full_output_dir = os.path.join(output_dir, subdir)
    # create directory if does not already exist
    os.makedirs(full_output_dir, exist_ok=True)

    for f in fastq_files:
        *_, barcode, fastq_filename = f.split(os.path.sep)
        runid = fastq_filename.split("_")[2]

        output_fastq_path = os.path.join(full_output_dir, f'{barcode}_{runid}.fastq')
        with open(output_fastq_path, 'a') as fh_out, open(f, "r") as fh_in:
            fh_out.write(fh_in.read())


def cat_sequence_summaries(input_dirs: List[str], sequencing_summary: str) -> None:
    """Combine all sequencing summaries into one file

    Get sequence summaries from the directories in `input_dirs` and concatenate
    them all into one string, being written to the master sequencing summary
    file located `sequencing_summary` (specified by user).

    :param input_dirs:
    :param sequencing_summary:
    :return:
    """
    sequencing_summary_files = [txt_path for input_dir in input_dirs for txt_path in
                                glob(os.path.join(input_dir, 'sequencing_summary.txt'))]
    df: pd.DataFrame = pd.concat([pd.read_table(f) for f in sequencing_summary_files])
    df.to_csv(sequencing_summary, sep='\t', index=None)


def get_num_reads_processed(sequencing_summary: str) -> int:
    """Get the total number of reads processed by Albacore.

    :param sequencing_summary: Sequencing summary file path
    :return: number of lines in the sequencing summary file, (minus the header line).
    """
    with open(sequencing_summary, "r") as fp:
        num_reads_processed = sum(1 for l in fp) - 1
    return num_reads_processed


def get_num_reads_stats(sequencing_summary: str) -> Tuple[int, int, int]:
    """Get the number of failed and passed reads.

    Go through the sequencing summary file and for each line record whether
    the `passes_filtering` field is true or false.

    :param sequencing_summary: Sequencing summary file path
    :return: Tuple of reads passed, failed and with no match
    """
    df = pd.read_table(sequencing_summary)

    num_reads_passed = df.passes_filtering.sum()
    num_reads_failed = df.passes_filtering.size - num_reads_passed
    num_no_match = (df.calibration_strand_genome_template == 'no_match').sum()

    return num_reads_passed, num_reads_failed, num_no_match


def get_barcode_stats(output_dir: str) -> Dict[str, Dict[str, int]]:
    """Generate and return a dict containing barcode run statistics.

    Read the `fastq` files in the given output directory and use these to
    generate a dict with a key for each barcode. The value for each barcode is a
    a dict with `pass`, `fail`, and `total` entries, containing the number of
    reads that passed, the number of reads failed, and total number of reads,
    respectively.

    :param output_dir:
    :return:
    """
    barcode_stats = defaultdict(lambda: defaultdict(int))
    for status in ['pass', 'fail']:
        for fastq_file in glob(os.path.join(output_dir, status, '*.fastq')):
            barcode_name, *_ = os.path.basename(fastq_file).split("_")
            barcode_stats[barcode_name][status] += count_reads(fastq_file)

    return barcode_stats


def count_reads(fastq_file: str) -> int:
    """Count number of reads in a FASTQ file.

    Assuming correctly formatted FASTQ file where number of reads is the number of lines divided by 4.

    :param fastq_file: FASTQ file path
    :return: Number of reads in FASTQ file
    """
    with open(fastq_file, "r") as fp:
        numlines = sum(1 for l in fp)
    return int(numlines / 4)



def generate_run_health(sequencing_summary, run_health_path) -> None:
    """Write run health statistics file.

    Generate a report including reads processed by Albacore, number of passing
    reads and failing reads, and more. Save this data in nice tabular format to
    `output_path`.

    :param sequencing_summary:
    :param run_health_path:
    """
    num_reads_processed = get_num_reads_processed(sequencing_summary)
    (num_reads_passed, num_reads_failed, num_reads_no_match) = get_num_reads_stats(sequencing_summary)

    output_dir = os.path.dirname(run_health_path)
    barcode_stats = get_barcode_stats(output_dir)

    with open(run_health_path, "w") as fp:
        # save info in comments
        fp.write(f'#Total reads processed by Albacore:\t{num_reads_processed}\n'
                 f'#Reads passed:\t{num_reads_passed}\n'
                 f'#Reads failed:\t{num_reads_failed}\n'
                 f'#Reads that failed with "no_match":\t{num_reads_no_match}\n'
                 f'Barcode\tPassing Reads\tFailing Reads\tTotal Reads\n')
        total_pass = 0
        total_fail = 0
        for barcode, stats in barcode_stats.items():
            n_pass = stats['pass']
            n_fail = stats['fail']
            total_pass += n_pass
            total_fail += n_fail
            fp.write(f'{barcode}\t{n_pass}\t{n_fail}\t{n_pass + n_fail}\n')
        fp.write(f'TOTAL\t{total_pass}\t{total_fail}\t{total_pass + total_fail}\n')


def get_albacore_results(input_dirs: List[str],
                         sequencing_summary: str,
                         run_health: str,
                         deletefailed: bool = False) -> None:
    """

    :param input_dirs:
    :param sequencing_summary:
    :param run_health:
    :param deletefailed:
    :return:
    """
    # get directory of ouput `sequencing_summary` so we can output the reads
    # there as well.
    output_dir = os.path.dirname(sequencing_summary)

    # Ensure input_dirs is always a list, even when there's only one entry.
    if not isinstance(input_dirs, (list, tuple)):
        input_dirs = [input_dirs]

    fastq_pass_files = get_fastq_files(input_dirs, "workspace/pass")

    # combine pass fastq files into result directory
    cat_fastq_files(fastq_pass_files, output_dir, "pass")

    fastq_failed_files = get_fastq_files(input_dirs, "workspace/fail")

    # combine fail fastq files into result directory
    cat_fastq_files(fastq_failed_files, output_dir, "fail")

    cat_sequence_summaries(input_dirs, sequencing_summary)
    generate_run_health(sequencing_summary, run_health)

    # check to see if we delete the failed fastqs reads
    if deletefailed:
        shutil.rmtree(os.path.join(output_dir, 'fail'))
