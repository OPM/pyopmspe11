# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Compare the content of two folders"""

import filecmp
import argparse
import os


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main script to compare the content of two folders",
    )
    parser.add_argument(
        "-a",
        "--afolder",
        default="result1",
        type=str.strip,
        help="Name of the first folder",
    )
    parser.add_argument(
        "-b",
        "--bfolder",
        default="result2",
        type=str.strip,
        help="Name of the second folder",
    )
    parser.add_argument(
        "-c",
        "--extensions",
        default="INC,csv,DATA",
        type=str.strip,
        help="Extension of files to compare",
    )
    args = vars(parser.parse_known_args()[0])
    afolder = args["afolder"]
    bfolder = args["bfolder"]
    extensions = (args["extensions"].strip()).split(",")
    for root, _, files in os.walk(afolder):
        for file in files:
            for ext in extensions:
                if file.endswith(ext):
                    afile = os.path.join(root, file)
                    bfile = os.path.join(
                        bfolder + "/" + "/".join((root.split("/")[1:])), file
                    )
                    assert os.path.join(
                        bfolder + "/" + "/".join((root.split("/")[1:])), file
                    ), f"File {bfile} does not exist"
                    is_equal = filecmp.cmp(afile, bfile, shallow=False)
                    if not is_equal:
                        print(f"Config mismatch for {afile} and {bfile}")


if __name__ == "__main__":
    main()
