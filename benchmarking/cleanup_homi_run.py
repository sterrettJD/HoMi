#!/usr/bin/env python3
"""
Deletes intermediate files from a HoMi run. Snakemake does not run HoMi steps
again because of these deletions. You can use it on a run that is not complete.
The benchmarking pipeline runs it after each local HoMi run.

Deletes:
    * the adapter trimmed reads ({PROJ}.noadpt/) and the trimmed reads
      ({PROJ}.f*.r*/), if host read removal is done for all samples and both
      MultiQC reports exist.
    * the Kraken2 classified and unclassified reads and the Kraken2 per read
      outputs (*.kraken.txt). No HoMi step reads these files.
    * the files in each HUMAnN sample temp folder, if the aggregated HUMAnN outputs
      exist. The MetaPhlAn bugs lists and the logs are kept, because HoMi declares
      the bugs lists as outputs.
    * the nonhost reads, if the reads breakdown report, the aggregated Kraken2
      outputs and the aggregated HUMAnN outputs exist. The host read removal
      reports are kept, because the reads breakdown report uses them. If HUMAnN
      did not run, the nonhost reads are kept, so that HUMAnN can run later.

The benchmarking plots use only {PROJ}_reads_breakdown.csv, Combined-taxonomy.tsv
and all_bugs_list.tsv. This script keeps these files.

Usage, from benchmarking/:
    python cleanup_homi_run.py semi/dna_semi_HoMi_config.yaml           # list only
    python cleanup_homi_run.py semi/dna_semi_HoMi_config.yaml --delete  # delete
"""

import argparse
import glob
import os
import shutil

import yaml

# run_humann_nonhost and aggregate_humann_outs_nonhost declare the bugs lists as outputs
HUMANN_TEMP_FILES_TO_KEEP = ("_metaphlan_bugs_list.tsv", "_metaphlan_bugs_list_v3.tsv", ".log")

def get_args():
    parser = argparse.ArgumentParser(description="Delete intermediate files from a HoMi run.")
    parser.add_argument("config", help="HoMi config file")
    parser.add_argument("--delete", action="store_true",
                        help="delete the files. Without this option, the script only lists them.")
    return parser.parse_args()


def read_config(config_path):
    with open(config_path) as config_file:
        return yaml.safe_load(config_file)


def get_trim_trunc_path(config):
    # the same path that the HoMi snakefile uses
    return (f"{config['PROJ'].strip()}.f{config['trim_fwd']}.{config['trunc_fwd']}"
            f".r{config['trim_rev']}.{config['trunc_rev']}")


def get_samples(config):
    # HoMi links the raw reads of each sample to {PROJ}/{sample}.R1.fq.gz
    raw_reads = glob.glob(os.path.join(config["PROJ"].strip(), "*.R1.fq.gz"))
    return [os.path.basename(path).removesuffix(".R1.fq.gz") for path in raw_reads]


def get_missing(paths):
    return [path for path in paths if not os.path.exists(path)]


def get_trimmed_reads_paths(config):
    """
    Returns the trimmed read folders, or an empty list if HoMi still needs them.
    Host read removal and FastQC read these folders.
    """
    proj = config["PROJ"].strip()
    trim_trunc_path = get_trim_trunc_path(config)
    samples = get_samples(config)
    required = ([os.path.join(f"{trim_trunc_path}.nonhost", f"{sample}.report") for sample in samples]
                + [os.path.join(f"{proj}.noadpt.fastqc", "multiqc_report"),
                   os.path.join(f"{trim_trunc_path}.fastqc", "multiqc_report")])
    missing = get_missing(required)
    if not samples or missing:
        print("Keeping the trimmed reads, because host read removal or MultiQC is not done")
        return []
    return [f"{proj}.noadpt", trim_trunc_path]


def get_kraken_paths(config):
    kraken_path = f"{get_trim_trunc_path(config)}.nonhost.kraken"
    return (glob.glob(os.path.join(kraken_path, "*classified_*.fq"))
            + glob.glob(os.path.join(kraken_path, "*.kraken.txt")))


def is_kept_humann_temp_file(path):
    return path.endswith(HUMANN_TEMP_FILES_TO_KEEP)


def get_humann_temp_paths(config):
    """
    Returns the files in the HUMAnN sample temp folders, except the bugs lists and
    the logs, or an empty list if the aggregated HUMAnN outputs do not exist.
    """
    humann_path = f"{get_trim_trunc_path(config)}.nonhost.humann"
    if not os.path.exists(os.path.join(humann_path, "all_genefamilies.tsv")):
        print("Keeping the HUMAnN temp files, because the aggregated HUMAnN outputs do not exist")
        return []
    temp_files = glob.glob(os.path.join(humann_path, "*", "*_humann_temp", "*"))
    return [path for path in temp_files if not is_kept_humann_temp_file(path)]


def get_nonhost_reads_paths(config):
    """
    Returns the nonhost reads and the concatenated nonhost reads, or an empty list
    if Kraken2, HUMAnN or the reads breakdown report is not done.
    """
    trim_trunc_path = get_trim_trunc_path(config)
    required = [f"{config['PROJ'].strip()}_reads_breakdown.csv",
                os.path.join(f"{trim_trunc_path}.nonhost.kraken", "Combined-taxonomy.tsv"),
                os.path.join(f"{trim_trunc_path}.nonhost.humann", "all_genefamilies.tsv")]
    if get_missing(required):
        print("Keeping the nonhost reads, because Kraken2, HUMAnN or the reads breakdown is not done")
        return []
    return (glob.glob(os.path.join(f"{trim_trunc_path}.nonhost", "*.fq.gz"))
            + glob.glob(os.path.join(f"{trim_trunc_path}.nonhost.concat", "*.fq.gz")))


def get_size_gb(path):
    if os.path.isfile(path):
        return os.path.getsize(path) / 1e9
    total = 0
    for dirpath, _, filenames in os.walk(path):
        total += sum(os.path.getsize(os.path.join(dirpath, name)) for name in filenames)
    return total / 1e9


def delete_path(path):
    # a folder can be a symlink to another disk, so delete the target and the link
    if os.path.islink(path) and os.path.isdir(path):
        shutil.rmtree(os.path.realpath(path))
        os.remove(path)
    elif os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


def main():
    args = get_args()
    config = read_config(args.config)
    paths = [path for path in (get_trimmed_reads_paths(config)
                               + get_kraken_paths(config)
                               + get_humann_temp_paths(config)
                               + get_nonhost_reads_paths(config))
             if os.path.lexists(path)]

    total_gb = 0
    for path in paths:
        size_gb = get_size_gb(path)
        total_gb += size_gb
        print(f"{size_gb:8.1f} GB  {path}")
        if args.delete:
            delete_path(path)

    action = "Deleted" if args.delete else "Can delete (use --delete)"
    print(f"{action}: {total_gb:.1f} GB in {len(paths)} paths")


if __name__ == "__main__":
    main()
