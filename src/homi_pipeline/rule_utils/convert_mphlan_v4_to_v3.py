#!/usr/bin/env python3

import argparse
from aggregate_metaphlan_bugslists import get_filepaths


def get_args():
    """
    handles arg parsing for this script

    returns the parsed args
    """
    parser = argparse.ArgumentParser(
        prog="Convert MetaPhlan v4 reports to v3 format",
        description="Convert MetaPhlan v4 reports to v3-like format, or at least close enough to trick Pavian. It really just modifies the first line"
    )
    # input directory
    parser.add_argument("-i", "--indir",
                        help="Input directory that contains all of the sample directories",
                        required=True)
    # a "padder" string that exists between the sample id and the metaphlan output name type
    parser.add_argument("-p", "--padder",
                        help="a 'padder' string that exists between the sample id and the humann output name type",
                        default="")

    parsed_args = parser.parse_args()
    return parsed_args


def fix_files(filepaths):
    for file in filepaths:
        print(f"Fixing: {file}")
        with open(file, "r") as f:
            contents = f.readlines()

        # replace the first line to pass the Pavian checks for Metaphlan v3
        contents[0] = "#mpa_v3\n"

        # create new filepath
        # should take "dir/subdir/bugs_list.tsv" -> "dir/subdir/bugs_list_v3.tsv"
        split_by_dot = file.split(".")
        split_by_dot[-2] = split_by_dot[-2] + "_v3"
        v3_path = ".".join(split_by_dot)

        # write it to a new file
        print(f"Writing to {v3_path}")
        with open(v3_path, "w") as outfile:
            outfile.writelines(contents)


if __name__ == "__main__":
    args = get_args()
    filepaths, subdirs = get_filepaths(args.indir, args.padder)
    fix_files(filepaths)
