from time import sleep
from base64 import b64decode
from tempfile import NamedTemporaryFile
from multiprocessing import Process
import requests
import traceback
from threading import Thread
import sys
import asyncio
from functools import lru_cache
from gzip import decompress, compress
from shutil import move, rmtree
from os import makedirs, environ
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve
from glob import glob
from subprocess import run, check_output, Popen
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from json import loads, dump, dumps
from helpers import get_session, get_compare_region_data, get_compare_region_json_path, send_alert, source_email, send_alert_background, logged_background, to_pid, logger
from pid import PidFile, PidFileError

def listen_pdb(port_hint):
    import pdb_attach
    for p in range(port_hint, port_hint+100):
        try:
            pdb_attach.listen(p)
            break
        except OSError:
            pass
    else:
        raise Exception(f"Coudn't find a port")

def get_email_alerts_dir_path(genome_id):
    return Path(f'.json_files/{genome_id}/email_alerts')

def get_operon_progress_path(genome_id):
    return Path(f'.json_files/{genome_id}/operons_progress')

def get_operon_path(genome_id):
    return Path(f'.json_files/{genome_id}/operons.json')

def operons_in_progress(genome_id):
    try:
        with PidFile('.lock_'+genome_id):
            return False
    except PidFileError:
        return True

def get_operons_background_process(genome_id:str):
    Popen([sys.executable, "-c", f"from get_json import get_operons; get_operons({repr(genome_id)})"])


def get_operons(genome_id:str) -> dict[str, float]:
    logger.info("Enter get operons")
    for _ in range(10):
        try:
            with PidFile('.lock_'+genome_id):
                listen_pdb(40000)
                logger.info("Start get operons")
                predict_json = get_operon_path(genome_id)
                if predict_json.exists():
                    return loads(predict_json.read_bytes())

                progress_file = get_operon_progress_path(genome_id)
                progress_writer = lambda progress: progress_file.write_text(str(progress))

                progress_writer(0.0)

                compare_region_json_path = get_compare_region_json_path(genome_id)
                test_operons_path = Path(f"images_custom/test_operons/{genome_id}")
                if not compare_region_json_path.exists():
                    rmtree(test_operons_path, ignore_errors=True)
                pegs = list(to_pid(genome_id).full_data.keys())
                compare_region_data = get_compare_region_data(genome_id, pegs, progress_writer)
                logger.info("compare region")

                progress_writer(0.50)


                if not test_operons_path.exists() or len(list(test_operons_path.glob('*.jpg'))) < len(pegs) - 50:
                    from JsonToCoordinates import to_coordinates
                    logger.info("Enter tocoord")
                    coords_filename = to_coordinates(compare_region_data, genome_id, progress_writer)
                    logger.info("Coordinates created")

                    makedirs(test_operons_path, exist_ok=True)
                    run(["java", "CoordsToJpg.java", coords_filename, test_operons_path.as_posix()])
                    Path(coords_filename).unlink()
                    progress_writer(0.65)


                from test import main
                logger.info("before predict")
                for _ in range(10*60):
                    try:
                        with PidFile('.main_predictor_lock'):
                            operons = main(genome_id, progress_writer)
                            break
                    except PidFileError:
                        sleep(1)
                logger.info("predicted")

                progress_writer(1.0)

                predict_json.write_text(dumps(operons))

                logger.info("bye")
        except PidFileError:
            sleep(0.1)
        else:
            email_dir = get_email_alerts_dir_path(genome_id)
            if email_dir.exists():
                logger.info("Sending email")
                for b64_validated_email_path in email_dir.iterdir():
                    recipient = b64decode(b64_validated_email_path.name.encode()).decode()
                    send_alert_background(recipient, genome_id, None)
                    b64_validated_email_path.unlink()
                    logger.info(f"Sent to email {recipient}")
            return operons
    raise PidFileError


def operon_clusters(pegs: frozenset[int], min_prob: float, probs: dict[int, float]) -> list[set[int]]:
    peg_next  = {}
    prev = -1
    for peg in sorted(pegs):
        peg_next[prev] = peg
        prev = peg
    peg_nums = sorted(
        [[peg_num, peg_next[peg_num]] for peg_num, prob in probs.items() if peg_num in peg_next and prob >= min_prob]
    )
    # pair = [op, op[:op.rfind('.')+1] + str(peg_num + 1)]

    clusters: list[set[int]] = []
    for peg_num, next_peg_num in peg_nums:
        if clusters and peg_num in clusters[-1]:
            clusters[-1].add(next_peg_num)
        else:
            clusters.append({peg_num, next_peg_num})

    return clusters
