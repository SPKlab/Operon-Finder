from collections import namedtuple
from os import makedirs, mkdir
from pathlib import Path
from gzip import decompress
from urllib.request import urlretrieve
from glob import glob
from functools import cache
from subprocess import run, check_output
from json import dump, dumps, load, loads
from typing import Iterable, Optional
from time import time

from attr import dataclass


@cache
def curl_output(*args: str)->bytes:
    return check_output(('curl', '--compressed') + args)

PatricMeta = namedtuple('PatricMeta', ['refseq', 'desc', 'protein_id'])
LocInfo = namedtuple('LocInfo', ['start', 'end'])

def to_pid( genome_id: str) -> tuple[dict[int, PatricMeta], int, str, dict[str, LocInfo]]:
    genome_data = get_genome_data(genome_id)
    feature_data = genome_data["docs"]

    gene_locations = {}
    full_data = {}
    assert feature_data
    for feature in feature_data:
        desc: str = feature["product"]
        desc_col_loc = desc.find(': ')
        if desc_col_loc != -1:
            desc = desc[desc_col_loc + 2:]

        refseq = feature.get("refseq_locus_tag", "None")
        protein_id = feature.get("protein_id", "None")
        patric_id = int(feature["patric_id"].split(".")[-1])
        full_data[patric_id] = PatricMeta(
            desc=desc, refseq=refseq, protein_id=protein_id
        )
        
        gene_locations[patric_id] = LocInfo(start=feature['start'], end=feature['end'])

    sequence_accession_id = feature["sequence_id"]
    gene_count: int = genome_data["numFound"]
    return full_data, gene_count, sequence_accession_id, gene_locations

def query_keywords(query: str) -> set[str]:
    return {qs.lower() for qs in query.split(' ') if qs}

def fetch_string_scores(genome_id: str) -> None:
    genome_data_changed = False

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

class _Data:
    def __init__(self):
        self.last_update = time()
        self.changed = False
    def updated(self, changed = True):
        self.last_update = time()
        self.changed = changed
data = _Data()

@cache
def get_genome_data(genome_id: str):
    genome_data_dir = f'.json_files/{genome_id}'
    genome_data_path = Path(f'{genome_data_dir}/genome.json')
    if genome_data_path.exists():
        genome_data = loads(genome_data_path.read_bytes())
    else:
        genome_data = loads(curl_output(
                "https://patricbrc.org/api/genome_feature/?&http_accept=application/solr+json&http_download=true",
                "-H",
                "Content-Type: application/x-www-form-urlencoded",
                "--data-raw",
                f"rql=eq%28genome_id%252C{genome_id}%29%2526and%28eq%28feature_type%252C%252522CDS%252522%29%252Ceq%28annotation%252C%252522PATRIC%252522%29%29%2526sort%28%252Bfeature_id%29%2526limit%2825000%29",
                "--compressed",
        ))["response"]
        makedirs(genome_data_dir, exist_ok=True)
        with open(genome_data_path, 'w') as f:
            dump(genome_data, f)
        data.updated()
    return genome_data
