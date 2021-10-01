from gzip import decompress
from shutil import move, rmtree
from os import makedirs
from pathlib import Path
from urllib.request import urlretrieve
from glob import glob
from subprocess import run
from concurrent.futures import ThreadPoolExecutor

def operon_clusters(genome_id: str, pegs: list[int])->list[list[int]]:
    genes = [f"fig|{genome_id}.peg.{i}" for i in pegs]

    json_folder = ".json_files"
    makedirs(json_folder, exist_ok=True)
    with ThreadPoolExecutor(5, "JSONFetcher") as executor:
        for gene in genes:
            json_path = f"{json_folder}/{gene}.json"
            if not Path(json_path).exists():
                args = [
                        "curl",
                        "--max-time",
                        "300",
                        "--data-binary",
                        '{"method": "SEED.compare_regions_for_peg", "params": ["'
                        + gene
                        + '", 5000, 20, "pgfam", "representative+reference"], "id": 1}',
                        "https://p3.theseed.org/services/compare_region",
                        "-o",
                        json_path,
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

    clusters: list[list[int]] = []
    for peg_num in peg_nums:
        if clusters and peg_num[0] == clusters[-1][1]:
            clusters[-1].append(peg_num[1])
        else:
            clusters.append(peg_num)
    
    return clusters
    