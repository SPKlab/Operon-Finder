from functools import cache
from gzip import decompress
from shutil import move, rmtree
from os import makedirs
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve
from glob import glob
from subprocess import run, check_output
from concurrent.futures import ThreadPoolExecutor, as_completed
from json import loads, dump
from helpers import data

import streamlit as st
from JsonToCoordinates import parse_string_scores

def get_operons(genome_id:str, pegs: frozenset) -> dict[str, float]:
    placeholder = st.empty()
    placeholder.info("Please wait while we fetch the data and predict operons. It might take upto 15 minutes.")

    progress_bar = st.progress(0.05)
    genome_data_changed = False
    gene_figure_name = {f"fig|{genome_id}.peg.{i}" for i in pegs}

    json_folder = f".json_files/{genome_id}/compare_region"
    makedirs(json_folder, exist_ok=True)
    if len(list(Path(json_folder).iterdir())) < len(gene_figure_name)-50:
        with ThreadPoolExecutor(60, "JSONFetcher") as executor:
            tasks =[]
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
                    tasks.append(executor.submit(run, args))
                    genome_data_changed = True
            for i, future in enumerate(as_completed(tasks)):
                progress_bar.progress(i/len(tasks)*0.50)

    progress_bar.progress(0.50)

    from JsonToCoordinates import to_coordinates

    test_operons_path = f"images_custom/test_operons/{genome_id}"
    if genome_data_changed:
        data.updated()
        Path(test_operons_path).unlink(missing_ok=True)

    if not Path(test_operons_path).exists() or len(list(Path(test_operons_path).glob('*.jpg'))) < len(pegs) - 50:
        coords_filename = to_coordinates(json_folder, genome_id)
        print("Coordinates created")

        progress_bar.progress(0.55)
        makedirs(test_operons_path, exist_ok=True)
        run(["java", "CoordsToJpg.java", coords_filename, test_operons_path])
        Path(coords_filename).unlink()
        progress_bar.progress(0.65)

    placeholder.empty()

    from test import main
    return main(genome_id, progress_bar)

@cache
def operon_probs(genome_id: str, pegs: frozenset) -> dict[str, float]:
    makedirs(f'.json_files/{genome_id}', exist_ok=True)
    predict_json = Path(f'.json_files/{genome_id}/operons.json')
    if predict_json.exists():
        operons = loads(predict_json.read_bytes())
        1+1
    else:
        operons = get_operons(genome_id, pegs)
        with open(predict_json, 'w') as f:
            dump(operons, f)
        data.updated()
    # JSON keys can only be strings
    return {int(gene_id): prob for gene_id, prob in operons.items()}

def operon_clusters(genome_id: str, pegs: frozenset[int], min_prob: float) -> tuple[list[set[int]], dict[int, float]]:
    peg_next  = {}
    prev = -1
    for peg in sorted(pegs):
        peg_next[prev] = peg
        prev = peg
    probs = operon_probs(genome_id, pegs)
    peg_nums = sorted(
        [[peg_num, peg_next[peg_num]] for peg_num, prob in probs.items() if prob >= min_prob]
    )
    # pair = [op, op[:op.rfind('.')+1] + str(peg_num + 1)]

    clusters: list[set[int]] = []
    for peg_num, next_peg_num in peg_nums:
        if clusters and peg_num in clusters[-1]:
            clusters[-1].add(next_peg_num)
        else:
            clusters.append({peg_num, next_peg_num})

    return clusters, probs
