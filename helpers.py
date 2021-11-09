from collections import namedtuple
from os import makedirs, mkdir
from pathlib import Path
from functools import cache
from subprocess import check_output
from json import dump, dumps, load, loads
from typing import Iterable, Optional

from attr import dataclass


@cache
def curl_output(*args: str)->bytes:
    return check_output(('curl', '--compressed') + args)

# @dataclass
# class PatricMeta:
#     refseq: str
#     desc: str
#     protein_id: str

PatricMeta = namedtuple('PatricMeta', ['refseq', 'desc', 'protein_id'])

def to_pid(
    genome_id: str, refseqs: Optional[frozenset[str]], query: Optional[str], protein_ids: Optional[frozenset[str]]
) -> tuple[dict[int, PatricMeta], int]:
    if refseqs:
        refseqs = frozenset({refseq.lower() for refseq in refseqs})
    query_strings = query_keywords(query) if query and query.strip(' ') else None
    if protein_ids:
        protein_ids = frozenset({pid.lower() for pid in protein_ids})

    genome_data = get_genome_data(genome_id)
    feature_data = genome_data["docs"]

    full_data = {}
    for feature in feature_data:
        desc: str = feature["product"]
        desc_col_loc = desc.find(': ')
        if desc_col_loc != -1:
            desc = desc[desc_col_loc + 2:]

        refseq = feature.get("refseq_locus_tag", "None")
        protein_id = feature.get("protein_id", "None")
        if (
            (query_strings and not all(s in desc.lower() for s in query_strings))
            or (refseqs and (refseq.lower() not in refseqs))
            or (protein_ids and protein_id.lower() not in protein_ids)
        ):
            continue
        patric_id = int(feature["patric_id"].split(".")[-1])
        full_data[patric_id] = PatricMeta(
            desc=desc, refseq=refseq, protein_id=protein_id
        )

    gene_count: int = genome_data["numFound"]
    return full_data, gene_count

def query_keywords(query: str) -> set[str]:
    return {qs.lower() for qs in query.split(' ') if qs}


@cache
def get_genome_data(genome_id: str):
    genome_data_dir = '.json_files/genome'
    genome_data_path = Path(f'{genome_data_dir}/{genome_id}.json')
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
    return genome_data
