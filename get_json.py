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
from json import loads

@cache
def operon_clusters(genome_id: str, pegs: Optional[frozenset[int]])->list[set[int]]:
    if pegs is None:
        x = loads(check_output([
            "curl", 'https://patricbrc.org/api/genome_feature/',
              "-H", 'accept: application/solr+json',
              "-H", 'Content-Type: application/rqlquery+x-www-form-urlencoded',
              "--data-raw", 'and(eq(genome_id,83332.12),eq(annotation,%22PATRIC%22))&limit(1)&facet((field,feature_type),(mincount,1))',
              "--compressed" # TODO: Add it to other commands as well
        ]))
        max_pegs = x["facet_counts"]["facet_fields"]["feature_type"][1]
        pegs = frozenset(range(1, max_pegs + 1))

    genes = [f"fig|{genome_id}.peg.{i}" for i in pegs]

    json_folder = ".json_files"
    makedirs(json_folder, exist_ok=True)
    with ThreadPoolExecutor(5, "JSONFetcher") as executor:
        for gene in genes:
            json_path = f"{json_folder}/{gene}.json"
            if not Path(json_path).exists():
                args = [
                        "curl",
                        "--fail",
                        "--max-time", "300",
                        "--data-binary",
                        '{"method": "SEED.compare_regions_for_peg", "params": ["'
                        + gene
                        + '", 5000, 20, "pgfam", "representative+reference"], "id": 1}',
                        "https://p3.theseed.org/services/compare_region",
                        "--compressed",
                        "-o", json_path,
                    ]
                executor.submit(run, args)

    organism_genomes = {".".join(gene.split("|")[1].split(".")[:2]) for gene in genes}

    for organism_genome in organism_genomes:
        genome_file_path = Path("genomes/" + organism_genome + ".PATRIC.gff")
        if not genome_file_path.exists():
            run(
            [
                "wget",
                f"ftp://ftp.patricbrc.org/genomes/{organism_genome}/{organism_genome}.PATRIC.gff",
                "-O",
                genome_file_path,
            ]
        )
            Path('strings/parsed.json').unlink(missing_ok=True)

    organisms = {organism_genome.split(".")[0] for organism_genome in organism_genomes}

    for organism in organisms:
        if not glob("strings/" + organism + ".protein.links.v11.+.txt"):
            raw_path = f"{organism}.protein.links.v11.5.txt.gz"
            urlretrieve( f"https://stringdb-static.org/download/protein.links.v11.5/{raw_path}", raw_path)
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
            Path('strings/parsed.json').unlink(missing_ok=True)

    from JsonToCoordinates import to_coordinates
    to_coordinates(json_folder)

    makedirs("images", exist_ok=True)
    run(["java", "CoordsToJpg.java"])

    makedirs("images_custom", exist_ok=True)
    test_operons_path = 'images_custom/test_operons'
    rmtree(test_operons_path)
    move('images', test_operons_path)

    from test import main
    main()

    operons = Path('operon_pegs.txt').read_text().splitlines()
    peg_nums = sorted([[(_peg_num:=int(op.split('.')[-1])), _peg_num+1] for op in operons])
    # pair = [op, op[:op.rfind('.')+1] + str(peg_num + 1)]

    clusters: list[set[int]] = []
    for peg_num in peg_nums:
        if clusters and peg_num[0] in clusters[-1]:
            clusters[-1].add(peg_num[1])
        else:
            clusters.append(set(peg_num))
    
    return clusters
    