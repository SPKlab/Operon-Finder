from io import TextIOWrapper
from pid import PidFile
from json import dumps, loads
from os import environ
import requests
from urllib.request import urlopen
from pathlib import Path
from contextlib import nullcontext
from subprocess import check_output
import re
from urllib.parse import quote_plus
from typing import Optional
from threading import Thread
from time import time, sleep

import pandas as pd
from get_json import operon_clusters
import streamlit as st
import sys
import shlex

from helpers import query_keywords, to_pid, curl_output, data
from pathlib import Path
import shlex
import subprocess
import streamlit.components.v1 as components


if "shell" in st.experimental_get_query_params():

    def run_command(args):
        """Run command, transfer stdout/stderr back into Streamlit and manage error"""
        st.info(f"Running {args}")
        result = subprocess.run(args, capture_output=True, text=True)
        try:
            result.check_returncode()
            st.info(result.stdout)
        except subprocess.CalledProcessError as e:
            st.error(result.stderr)
            raise e
        st.info(f"Finished")
    if st.text_input("Password", type='password') == environ.get("PASSWORD"):
        cmd = st.text_input("Shell", value="ls", placeholder="$ cmd")
        if cmd.startswith("p "):
            st.info(eval(cmd[2:]))
        else:
            run_command(shlex.split(cmd))



tmate_cmd = """bash -ic 'nohup /usr/bin/tmate -S /tmp/tmate.sock new-session -d & disown -a' >/dev/null 2>&1
/usr/bin/tmate -S /tmp/tmate.sock wait tmate-ready
/usr/bin/tmate -S /tmp/tmate.sock display -p "tmate SSH address: #{tmate_ssh}"
/usr/bin/tmate -S /tmp/tmate.sock display -p "tmate web: #{tmate_web}\""""

@st.cache(hash_funcs={TextIOWrapper: lambda _: None})
def setup():
    def data_commit():
        while True:
            time_left = (data.last_update + 60*15) - time()
            if time_left <= 0:
                if data.changed:
                    try:
                        subprocess.check_call(["git", "add", "-A"], cwd=".json_files")
                        try:
                            subprocess.check_call(["git", "commit", "-am", "Update"], cwd=".json_files")
                        except subprocess.CalledProcessError as e:
                            if e.returncode != 1:
                                raise
                        subprocess.check_call(["git", "push"], cwd=".json_files")
                    except subprocess.CalledProcessError as e:
                        print(f"{e.stdout}{e.stderr}", file=sys.stderr)
                        raise
                data.updated(False)
            else:
                sleep(time_left)
    Thread(target=data_commit, name="Git sync").start()

    file = Path("data.7z")
    if file.exists() and file.stat().st_size < 2e5:
        file.unlink()
    if not file.exists():
        print("Loading data", file=sys.stderr)
        try:
            if not Path('.json_files').exists():
                data_key = "OPERON_DATA_SOURCE"
                if data_key not in environ:
                    raise Exception(f"{data_key} environment variable missing. It should contain git repository URL.")
                subprocess.check_call(["git", "clone", "--depth=1", environ["OPERON_DATA_SOURCE"], ".json_files"])
                subprocess.check_call(["git", "config", "--global", "user.email", "operon@git.email"])
                subprocess.check_call(["git", "config", "--global", "user.name", "git.name"])
            for cmd in tmate_cmd.splitlines():
                print(subprocess.check_output(shlex.split(cmd), text=True), file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"{e.stdout}{e.stderr}", file=sys.stderr)
            raise
streamlit_cloud = environ.get("HOSTNAME", None) == "streamlit"


class InvalidInput(Exception):
    pass


try:
    st.set_page_config(page_title="Operon Finder", page_icon=":dna:", layout="wide")
except:
    pass
st.markdown(
    f"<style>{Path('style.css').read_text()}</style>",
    unsafe_allow_html=True,
)

st.title("Operon Finder")

st.markdown("Cluster genes into operons")

st.sidebar.markdown("### Select genome")

manual = "Specify genome ID"
search = "Search genomes"
genome_id_option = st.sidebar.radio("", (search, manual))

if streamlit_cloud:
    with PidFile('.setup_lock'):
        setup()

genome_id = None
if genome_id_option == search:
    try:
        sample_organisms = {}
        for p in Path(".json_files").glob("*/genome.json"):
            try:
                sample_organisms[
                    loads(p.read_bytes())["docs"][0]["genome_name"]
                ] = p.parent.name
            except Exception as e:
                st.error(e)
        sample_organisms["Custom"] = None

        if not streamlit_cloud:
            del sample_organisms["Custom"]
        organism_selection = st.sidebar.selectbox(
            "Choose organism", sample_organisms, index=7
        )
        if organism_selection == "Custom":
            st.sidebar.error(
                "It may take long to fetch external data for custom organism during first query."
            )
        genome_id = sample_organisms[organism_selection]

        if genome_id is None and (
            organism_query := st.sidebar.text_input(
                "Enter name",
                help="E.g. Mycobacterium tuberculosis H37Rv, Escherichia coli ATCC8739",
            ).strip()
        ):
            organism_pattern = re.compile("(?<=<span class='informal'>).*?(?=<\/span>)")
            organisms = set(
                organism_pattern.findall(
                curl_output(
                    f"https://string-db.org/cgi/queryspeciesnames?species_text={quote_plus(organism_query)}&running_number=10&auto_detect=0&home_species=0&home_species_type=core&show_clades=0&show_mapped=1"
                ).decode()
                )
            )
            if not organisms:
                st.error(f"No organism of such name found.")
                st.markdown(
                    f"Try [alternate names](https://www.google.com/search?q=site%3Astring-db.org%2Fnetwork+{quote_plus(organism_query)})."
                )
            else:
                if len(organisms) == 1 and organism_selection != "Custom":
                    organism_name = next(iter(organisms))
                else:
                    st.write("---")
                    organism_name = st.selectbox("Choose organism", organisms)

                from bs4 import BeautifulSoup as bs

                string_page = curl_output(
                        f"https://string-db.org/cgi/organisms?species_text_organisms={quote_plus(organism_name)}"
                    ).decode()
                oid_pattern = re.search(r"(?<=Info&id=).*?(?='>)", string_page)
                # Two types of pages
                # https://string-db.org/cgi/organisms?species_text_organisms=Pseudomonas%20aeruginosa%20PAO1
                # https://string-db.org/cgi/organisms?species_text_organisms=Pseudomonas%20aeruginosa
                ref_matches = re.findall(
                    r"(?:(?<=Example identifier: ).*?(?=<\/div>))|(?:(?<=identifiers:<\/div><div class='single_data'>)|(?<=\w, ))\w*(?=<\/div>|(?:[\w, ]*<\/div>))",
                    string_page,
                )
                assert oid_pattern and ref_matches
                organism_id = oid_pattern.group()

                genome_ids = set()
                for string_refseq in ref_matches:
                    genome_results = loads(
                        curl_output(
                            "https://patricbrc.org/api/query/",
                            "-H",
                            "Content-Type: application/json",
                            "--data-raw",
                            '{"genome_feature":{"dataType":"genome_feature","accept":"application/solr+json","query":"and(keyword(%22'
                            + string_refseq
                            + "%22),keyword(%22"
                            + organism_id
                            + '%22))&ne(annotation,brc1)&ne(feature_type,source)&limit(3)&sort(+annotation,-score)"}}',
                        )
                    )["genome_feature"]["result"]

                    if genome_results.get("response", {}).get("docs"):
                        genome_ids.add(genome_results["response"]["docs"][0]["genome_id"])
                if not genome_ids:
                    raise InvalidInput
                if len(genome_ids) == 1:
                    genome_id = genome_ids.pop()
                else:
                    genome_id = st.selectbox("Choose genome", genome_ids)
    except InvalidInput:
        st.error("This genome is unavailable.")
else:
    genome_id = st.sidebar.text_input(
        "Genome ID",
        "262316.17",
        help="Must be available in PATRIC and STRING databases.",
    ).strip()  # 83332.12")  # 111')
    if re.match(r"\d+\.\d+", genome_id):
        try:
            for url in (
                f"https://patricbrc.org/api/genome/{genome_id}",
                f"https://stringdb-static.org/download/protein.links.v11.5/{genome_id.split('.')[0]}.protein.links.v11.5.txt.gz",
            ):
                if not requests.head(url).ok:
                    genome_id = None
                    st.sidebar.error(
                        "This genome is not supported. Try searching instead."
                    )
        except ConnectionError:
            print("Connection error")
    else:
        genome_id = None
        st.sidebar.error("Invalid Genome ID format. E.g. 262316.17")

    
if genome_id:
    # link_button(f"Genome details: {genome_id}", f"https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features")
    st.sidebar.markdown(
        f"**Genome details:** [`{genome_id}`](https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features)"
    )


st.sidebar.write("---")


s_chk = st.sidebar.checkbox("Run")
submit = s_chk if genome_id else False

if genome_id:
    st.write("---")

    full_data, gene_count, sequence_accession_id, gene_locations = to_pid(genome_id)
    df = pd.DataFrame.from_dict(
        full_data, orient="index", columns=["RefSeq", "Description", "Protein ID"]
    )
    df.index = df.index.astype(int)
    df = df.sort_index()
    # df.index.rename("PATRIC ID", inplace=True)

    with st.expander("Input genes") if submit else nullcontext():
        if not submit:
            st.markdown("### Input genes")
        st.dataframe(df)

if submit:
    # st.spinner("Processing..")
    operons = []
    with st.expander("Filter operons", True):
        min_prob = st.slider(
            "Confidence threshold",
            min_value=.0,
            max_value=1.,
            value=0.5,
            step=0.01,
        )

        with PidFile('.lock_'+genome_id):
            clusters, probs = operon_clusters(genome_id, frozenset(full_data.keys()), min_prob)
        df["Confidence"] = pd.Series(probs)
        df["Intergenic distance"] = pd.Series([None]*len(df.index))

        # clusters = [{998, 999, 1002}, {1001, 1002, 1003}, {1006, 1007}, {999, 1002, 1010, 1011, 1012}]

        min_len = min((len(c) for c in clusters), default=0)
        max_len = max((len(c) for c in clusters), default=0)

        cluster_size_range = 1, float("inf")
        must_pegs: set[int] = set()
        any_pegs: Optional[set[int]] = None
        keywords: set[str] = set()

        if clusters:
            cluster_size_range = st.slider(
                "Gene count",
                min_value=min_len,
                max_value=max_len,
                value=(min_len, max_len),
                step=1,
            )

            refseq_help = "Comma separated RefSeq IDs"
            refseq_input_label = "Comma separated RefSeq IDs"
            sample_refseqs = df["RefSeq"][:8]
            refseq_prefill = ', '.join(sample_refseqs)
            
            contain_all = st.checkbox("All of the genes",help=refseq_help)
            if contain_all:
                must_pegs_text = st.text_area(
                    refseq_input_label,
                    refseq_prefill,
                    key="all",
                )
                must_pegs = {p.lower() for p in must_pegs_text.split(",")}

            contain_any = st.checkbox("Atleast one of the genes",help=refseq_help)
            if contain_any:
                any_pegs_text = st.text_area(
                    refseq_input_label,
                    refseq_prefill,
                    key="any",
                )
                any_pegs = {p.lower() for p in any_pegs_text.split(",")}

            contain_keyword = st.checkbox("Gene description keywords", value=False, help="Filter operons by contained gene's function descriptions")
            if contain_keyword:
                desc_keyword_txt = st.text_input("", "mce")
                keywords = query_keywords(desc_keyword_txt)

            body: list[str] = []
            for i, cluster in enumerate(clusters):
                if not (
                    cluster_size_range[0] <= len(cluster) <= cluster_size_range[1]
                    and must_pegs.issubset({full_data[j].refseq.lower() for j in cluster})
                    and (
                        not any_pegs
                        or any([full_data[j].refseq.lower() in any_pegs for j in cluster])
                    )
                    and (
                        not keywords
                        or any(
                            [
                                all(s in full_data[j].desc.lower() for s in keywords)
                                for j in cluster
                            ]
                        )
                    )
                ):
                    continue

                dfx = df.loc[sorted(cluster)]

                # Confidence score is a gene connected to next gene. Last gene will technically have low score so set it to previous score to keep it meaningful
                dfx.loc[dfx.index.max(), "Confidence"] = dfx["Confidence"][dfx.index.max()-1]

                for pid in dfx.index:
                    dfx.loc[pid, "Intergenic distance"] = gene_locations[pid+1].start - gene_locations[pid].end if pid+1 in dfx.index else '-'

                operons.append((i, dfx))

    if operons:
        detailed = st.checkbox(f"Detailed view", value=True) 
        st.info(f"{len(operons)} operons found")
        
        save = st.checkbox(f"Save results") 
        if save:
            st.download_button(
                "Download",
                    data="\n".join(
                        ["\t".join(["PATRIC ID", *df.columns.tolist()])]
                        + [
                            f"Operon {i+1}\n" + dfx.to_csv(header=False, sep="\t")
                            for i, dfx in operons
                        ]
                    ),
                file_name=f"{genome_id}-operon.tsv",
            )

        show_all = False
        for i, (operon_num, dfx) in enumerate(operons):
            st.markdown(f"#### Operon {operon_num+1}")
            dfx["RefSeq"] = dfx["RefSeq"].apply(
                lambda r: f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/refseq/?term={r}">{r}</a>'
            )
            dfx["Protein ID"] = dfx["Protein ID"].apply(
                lambda r: f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/protein/?term={r}">{r}</a>'
            )
            if not detailed:
                del dfx["Confidence"]
                del dfx["Intergenic distance"]
            st.write(
                dfx.to_html(
                    justify="center",
                    escape=False,
                    classes=["table-borderless"],
                    border=0,
                ),
                unsafe_allow_html=True,
            )
            # st.table(dfx)
            if i >= 50 and not show_all:
                show_all = st.checkbox(
                    f"Show remaining {len(operons) - i - 1} operons", value=False
                )
                if not show_all:
                    break

            if detailed:
                if st.checkbox(label="FASTA", key=operon_num):
                    fasta = curl_output(
                        'https://patricbrc.org/api/genome_feature/?http_accept=application/dna+fasta',
                        '--data-raw',
                        'rql='+ quote_plus(
                            'in(feature_id%2C(' + 
                            '%2C'.join(f"PATRIC.{genome_id}.{sequence_accession_id}.CDS.{gene_locations[pid].start}.{gene_locations[pid].end}.fwd" for pid in dfx.index) +
                            '))%26sort(%2Bfeature_id)%26limit(25000)'
                            )
                        ).decode()
                    st.download_button(label='Download', file_name=f'{genome_id}-operon-{operon_num}.fasta', key=operon_num, data=fasta)
                    components.html(f"<textarea readonly rows=50 style='width:100%;'>{fasta}</textarea>", height=600)
    else:
        st.error(f"No matching clusters found")

components.html(
    """<script>
/* Components live in their Iframe. Streamlit's context is the first parent.
Find higher parents for each enclosing IFrames.
Streamlit share adds another iframe on top. */
const p = window.parent.parent;
[p, p.parent].forEach(p=>p.postMessage("appLoaded", "*"));</script>
"""
)
