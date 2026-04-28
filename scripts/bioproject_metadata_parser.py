#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import argparse
import subprocess

__copyright__ = ""
__license__ = ""
__version__ = ""
__authors__ = ["cgcole"]
__requires__ = [""]
__provides__ = [""]
__description__ = "This program generates the sample file that will be used to run the snakemake assembly program."


def fetch_bioproject_metadata(bioproject):
    os.makedirs("bioproject_metadata", exist_ok=True)
    output_path = f"bioproject_metadata/{bioproject}_SraRunInfo.csv"

    esearch = subprocess.Popen(
        ["esearch", "-db", "sra", "-query", bioproject], stdout=subprocess.PIPE
    )

    with open(output_path, "w") as f:
        subprocess.run(["efetch", "-format", "runinfo"], stdin=esearch.stdout, stdout=f)

    esearch.stdout.close()

    return pd.read_csv(output_path)


def parse_bioproject_metadata(bioproject):

    bioproject_metadata_table = fetch_bioproject_metadata(bioproject)

    bioproject_metadata_table = (
        bioproject_metadata_table[
            [
                "bases",
                "BioProject",
                "BioSample",
                "SampleName",
                "LibraryStrategy",
                "LibrarySource",
                "LibraryLayout",
                "Run",
                "ReleaseDate",
            ]
        ]
        .query('LibraryStrategy == "WGS" and LibrarySource == "GENOMIC"')
        .assign(
            short_R1_filename=lambda x: np.where(
                x["LibraryLayout"] == "PAIRED", x["Run"] + "_1.fastq.gz", np.nan
            ),
            short_R2_filename=lambda x: np.where(
                x["LibraryLayout"] == "PAIRED", x["Run"] + "_2.fastq.gz", np.nan
            ),
            long_filename=lambda x: np.where(
                x["LibraryLayout"] == "SINGLE", x["Run"] + ".fastq.gz", np.nan
            ),
        )
        .groupby(["BioSample", "LibraryLayout"])
        .apply(lambda x: x.nlargest(1, "bases"))
        .reset_index()
        .assign(
            sample=lambda x: x["SampleName"].str.replace(r"[-\s.]", "_", regex=True),
            library_ID=lambda x: pd.to_datetime(x["ReleaseDate"]).dt.strftime("%Y%m%d")
            + "-"
            + x["BioProject"],
        )
        .groupby("BioSample")
        .agg(
            BioProject=("BioProject", "first"),
            sample=("sample", "first"),
            library_ID=("library_ID", "first"),
            short_R1_filename=("short_R1_filename", "first"),
            short_R2_filename=("short_R2_filename", "first"),
            long_filename=("long_filename", "first"),
        )
        .reset_index()
    )

    bioproject_hybrid_samples = bioproject_metadata_table[
        bioproject_metadata_table[
            ["short_R1_filename", "short_R2_filename", "long_filename"]
        ]
        .notna()
        .all(axis=1)
    ]

    bioproject_long_samples = bioproject_metadata_table[
        bioproject_metadata_table[["short_R1_filename", "short_R2_filename"]]
        .isna()
        .any(axis=1)
        & bioproject_metadata_table["long_filename"].notna()
    ]

    bioproject_short_samples = bioproject_metadata_table[
        bioproject_metadata_table[["short_R1_filename", "short_R2_filename"]]
        .notna()
        .all(axis=1)
        & bioproject_metadata_table["long_filename"].isna()
    ]

    if os.path.exists("samples/hybrid_assembly/bioproject_hybrid_samples.csv"):
        existing = pd.read_csv("samples/hybrid_assembly/bioproject_hybrid_samples.csv")
        combined = pd.concat([existing, bioproject_hybrid_samples]).drop_duplicates(
            subset=["BioSample"]
        )
        combined.to_csv(
            "samples/hybrid_assembly/bioproject_hybrid_samples.csv", index=False
        )
    else:
        os.makedirs(
            os.path.dirname("samples/hybrid_assembly/bioproject_hybrid_samples.csv"),
            exist_ok=True,
        )
        bioproject_hybrid_samples.to_csv(
            "samples/hybrid_assembly/bioproject_hybrid_samples.csv", index=False
        )

    if os.path.exists("samples/long_assembly/bioproject_long_samples.csv"):
        existing = pd.read_csv("samples/long_assembly/bioproject_long_samples.csv")
        combined = pd.concat([existing, bioproject_long_samples]).drop_duplicates(
            subset=["BioSample"]
        )
        combined.to_csv(
            "samples/long_assembly/bioproject_long_samples.csv", index=False
        )
    else:
        os.makedirs(
            os.path.dirname("samples/long_assembly/bioproject_long_samples.csv"),
            exist_ok=True,
        )
        bioproject_long_samples.to_csv(
            "samples/long_assembly/bioproject_long_samples.csv", index=False
        )

    if os.path.exists("samples/short_assembly/bioproject_short_samples.csv"):
        existing = pd.read_csv("samples/short_assembly/bioproject_short_samples.csv")
        combined = pd.concat([existing, bioproject_short_samples]).drop_duplicates(
            subset=["BioSample"]
        )
        combined.to_csv(
            "samples/short_assembly/bioproject_short_samples.csv", index=False
        )
    else:
        os.makedirs(
            os.path.dirname("samples/short_assembly/bioproject_short_samples.csv"),
            exist_ok=True,
        )
        bioproject_short_samples.to_csv(
            "samples/short_assembly/bioproject_short_samples.csv", index=False
        )

    print("Finished downloading and parsing ", bioproject)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("bioproject", help="BioProject accesion")
    args = parser.parse_args()
    parse_bioproject_metadata(args.bioproject)
