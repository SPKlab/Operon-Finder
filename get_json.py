from functools import cache
from gzip import decompress
from shutil import move, rmtree
from os import makedirs
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve
from glob import glob
from subprocess import run, check_output
from concurrent.futures import ThreadPoolExecutor

import streamlit as st
from JsonToCoordinates import parse_string_scores
from test import main

@cache
def operon_clusters(genome_id: str, pegs: frozenset) -> list[set[int]]:
    progress_bar = st.progress(0.05)
    genome_data_changed = False
    gene_figure_name = {f"fig|{genome_id}.peg.{i}" for i in pegs}

    json_folder = f".json_files/compare_region/{genome_id}"
    makedirs(json_folder, exist_ok=True)
    if sum(1 for _ in Path(json_folder).iterdir()) < len(gene_figure_name):
        with ThreadPoolExecutor(30, "JSONFetcher") as executor:
            for gene in gene_figure_name:
                json_path = f"{json_folder}/{gene}.json"
                if not Path(json_path).exists():
                    args = [
                        "curl",
                        "--fail",
                        "--max-time",
                        "300",
                        "--data-binary",
                        '{"method": "SEED.compare_regions_for_peg", "params": ["'
                        + gene
                        + '", 5000, 20, "pgfam", "representative+reference"], "id": 1}',
                        "https://p3.theseed.org/services/compare_region",
                        "--compressed",
                        "-o",
                        json_path,
                    ]
                    executor.submit(run, args)
                    genome_data_changed = True

    progress_bar.progress(0.10)

    genome_file_path = Path("genomes/" + genome_id + ".PATRIC.gff")
    if not genome_file_path.exists():
        run(
            [
                "wget",
                f"ftp://ftp.patricbrc.org/genomes/{genome_id}/{genome_id}.PATRIC.gff",
                "-O",
                genome_file_path,
            ]
        )
        genome_data_changed = True

    organism = genome_id.split(".")[0]

    if not glob("strings/" + organism + ".protein.links.v11.*.txt"):
        raw_path = f"{organism}.protein.links.v11.5.txt.gz"
        urlretrieve(
            f"https://stringdb-static.org/download/protein.links.v11.5/{raw_path}",
            raw_path,
        )
        path = Path(raw_path)

        target_file_path = path.parent.joinpath("strings").joinpath(
            path.name.removesuffix(".gz")
        )
        with open(path, "rb") as compressed_file, open(
            target_file_path, "w", encoding="utf8"
        ) as decompressed_file:
            decom_str = decompress(compressed_file.read()).decode("utf-8")
            decompressed_file.write(decom_str)

        path.unlink()
        genome_data_changed = True

    progress_bar.progress(0.13)

    from JsonToCoordinates import to_coordinates

    test_operons_path = f"images_custom/test_operons/{genome_id}"
    if genome_data_changed:
        Path(test_operons_path).unlink(missing_ok=True)

    if not Path(test_operons_path).exists():
        coords_filename = to_coordinates(json_folder, genome_id)
        print("Coordinates created")

        progress_bar.progress(0.15)
        makedirs(test_operons_path, exist_ok=True)
        run(["java", "CoordsToJpg.java", coords_filename, test_operons_path])
        Path(coords_filename).unlink()
        progress_bar.progress(0.20)

    operons = main(genome_id, progress_bar)
    peg_next  = {}
    prev = -1
    for peg in sorted(pegs):
        peg_next[prev] = peg
        prev = peg
    peg_nums = sorted(
        [[peg_num, peg_next[peg_num]] for peg_num in operons]
    )
    # pair = [op, op[:op.rfind('.')+1] + str(peg_num + 1)]

    clusters: list[set[int]] = []
    for peg_num in peg_nums:
        if clusters and peg_num[0] in clusters[-1]:
            clusters[-1].add(peg_num[1])
        else:
            clusters.append(set(peg_num))

    return clusters
