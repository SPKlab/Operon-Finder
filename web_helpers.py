from os import makedirs, environ
from base64 import b64encode
import sys
from time import sleep
from helpers import source_email, send_alert_background, logger
import streamlit as st
from get_json import get_operons_background_process, get_operon_path, get_operon_progress_path, operons_in_progress, get_email_alerts_dir_path
from json import loads
from email_validator import validate_email, EmailNotValidError

class LruDict:
    def __init__(self, size):
        self.size = size
        self.d = {}
    def __getitem__(self, key):
        self.d[key] = self.d.pop(key)
        return self.d[key]
    def __setitem__(self, key, val):
        if len(self.d) >= self.size:
            del self.d[next(iter(self.d))]
        self.d[key] = val
    def __contains__(self, key):
        return key in self.d
    def __iter__(self):
        return iter(self.d)



operon_probs_cache = LruDict(128)
def operon_probs(genome_id: str, peg_count) -> dict[str, float]:
    if genome_id not in operon_probs_cache:
        makedirs(f'.json_files/{genome_id}', exist_ok=True)
        predict_json = get_operon_path(genome_id)
        if not predict_json.exists():
            placeholders = []
            stpl = lambda: placeholders.append(st.empty()) or placeholders[-1]
            stpl().info(f"Please wait while we fetch the data and predict operons. It might take upto {round(peg_count/5/60)} minutes.")
            email = stpl().text_input("Get email alert on completion", value=st.session_state.get("email", ""), placeholder='Email address', help="Make sure to check spam/junk folder. Email will be recieved from spklab.iitg@gmail.com").strip()
            validated_email = None
            if email:
                try:
                    validated_email = validate_email(email).email
                except EmailNotValidError as e:
                    stpl().error(f"Invalid email address\n\n{e}")
                else:
                    alert_dir = get_email_alerts_dir_path(genome_id)
                    alert_dir.mkdir(parents=True, exist_ok=True)
                    email_filename_bytes = b64encode(validated_email.encode())
                    if len(email_filename_bytes) > 255: # Filename length limit on ext, ntfs, etc. Email length limit 64@255
                        stpl().error("Please use another shorter email address.")
                    else:
                        alert_dir.joinpath(email_filename_bytes.decode()).touch(exist_ok=True)
                        stpl().success(f"On completion the email alert will be sent to {validated_email}. Make sure to check the junk/spam folder.")
            progress_file = get_operon_progress_path(genome_id)
            if not progress_file.exists():
                progress_file.write_text(str(0.01))
            progress_bar = stpl().progress(0.01)

            if not operons_in_progress(genome_id):
                get_operons_background_process(genome_id)
                while not operons_in_progress(genome_id):
                    sleep(0.1)
            logger.info("Starting loop")
            while not predict_json.exists() and operons_in_progress(genome_id):
                for _ in range(10):
                    try:
                        progress = float(progress_file.read_text())
                        break
                    except ValueError:
                        sleep(0.1)
                else:
                    continue
                progress_bar.progress(progress)
                sleep(1)

            if predict_json.exists():
                logger.info("Json exists")
            else:
                st.error(f"Some error occured, please retry and report the genome id to {source_email}")
                raise Exception(f"Error with {genome_id=}")
            for p in placeholders:
                p.empty()

        operons = loads(predict_json.read_bytes())
        operon_probs_cache[genome_id] = {int(gene_id): prob for gene_id, prob in operons.items()}
    probs = operon_probs_cache[genome_id]
    return probs


