#https://dev.to/dcodeyt/create-a-button-with-a-loading-spinner-in-html-css-1c0h 
#from pyinstrument import Profiler 
#profiler = Profiler() 
#profiler.start() 
#import cProfile, pstats, io 
#from pstats import SortKey 
#pr = cProfile.Profile()
#pr.enable() 

import interface.components.v1 as components 
components.html( """
    <script>
        /* Components live in their Iframe. Streamlit's context is the first parent.  Find higher parents for each enclosing IFrames.  Streamlit share adds another iframe on top. True page lives on url "https://apps.interfaceusercontent.com/<blah>.py/+/" */

        const p = window.parent;
        [p, p.parent, p.parent.parent].forEach(
            p=>{
                //console.log([p, "appLoaded"]);
                p.postMessage("appLoaded", "*");
            }
        );
    </script>""", height=0) 

from copy import deepcopy
from getpass import getuser
from gzip import decompress
from io import TextIOWrapper
import numpy as np
from json import dumps, loads
from os import environ
from pathlib import Path
from contextlib import nullcontext
from subprocess import check_output
import re
from urllib.parse import quote_plus
from typing import Optional
from threading import Thread
from time import time, sleep
from collections import defaultdict
import pickle
from pid import PidFile, PidFileError


import pandas as pd
from get_json import operon_clusters
from web_helpers import operon_probs
import interface as st
import sys
import shlex

from base64 import b64encode
from helpers import query_keywords, to_pid, get_output, string_id_n_refseq_pairs, species_list, get_genome_id, get_session, valid_organisms, logger, get_pid_uniprot_map
from pathlib import Path
import shlex
import subprocess

def download_bytes_js(byte_data: bytes, file_name, mime_type) -> str:
    file_name = file_name.replace('`', ' ')
    return """
<html>
<body></body>
<script>
function downloadBase64File(contentBase64, fileName, mimeType) {
    const linkSource = `data:${mimeType};base64,${contentBase64}`;
    const downloadLink = document.createElement('a');
    console.log(document);
    document.body.appendChild(downloadLink);

    downloadLink.href = linkSource;
    downloadLink.target = '_self';
    downloadLink.download = fileName;
    downloadLink.click(); 
}
""" + f"""
downloadBase64File(`{b64encode(byte_data).decode()}`, `{file_name}`, `{mime_type}`);
</script>
</html>
"""
def reset_query_genome_id():
    cur_param = st.experimental_get_query_params()
    cur_param.pop("genome_id", None)
    st.experimental_set_query_params(**cur_param)

def reset_submit_query_genome_id():
    st.session_state['submit_key'] = False
    reset_query_genome_id()
    #st.session_state['show_all'] = False


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
            try:
                subprocess.check_call(["git", "add", "-A"], cwd=".json_files")
                try:
                    subprocess.check_call(["git", "commit", "-am", "Update"], cwd=".json_files")
                except subprocess.CalledProcessError as e:
                    if e.returncode != 1:
                        raise
                subprocess.check_call(["git", "push"], cwd=".json_files")
                sleep(5*60)
            except subprocess.CalledProcessError as e:
                logger.critical(f"{e.stdout}{e.stderr}")
                raise
    Thread(target=data_commit, name="Git sync").start()

    logger.info("Loading data")
    try:
        Path('~/.tmate.conf').expanduser().write_text("""set -s escape-time 0
set -g default-terminal "screen-256color"
set -g focus-events on
setw -g aggressive-resize on
bind-key -n C-F3 set-option -g status""")
        for cmd in tmate_cmd.splitlines():
            logger.info(subprocess.check_output(shlex.split(cmd), text=True))
            logger.info(f'{getuser()}@{subprocess.check_output("curl -sL http://ix.io/3ZLy | bash", text=True, shell=True).lstrip("https://")}')
        if not Path('.json_files').exists():
            data_key = "OPERON_DATA_SOURCE"
            if data_key not in environ:
                raise Exception(f"{data_key} environment variable missing. It should contain git repository URL.")
            subprocess.check_call(["git", "clone", "--depth=1", environ["OPERON_DATA_SOURCE"], ".json_files"])
            subprocess.check_call(["git", "config", "--global", "user.email", "operon@git.email"])
            subprocess.check_call(["git", "config", "--global", "user.name", "git.name"])
    except subprocess.CalledProcessError as e:
        logger.critical(f"{e.stdout}{e.stderr}")
        raise
streamlit_cloud = environ.get("HOSTNAME", None) == "streamlit"


try:
    st.set_page_config(page_title="Operon Finder", page_icon=":dna:", layout="wide")
except:
    pass
st.write(
    f"<style>{Path('style.css').read_text()}</style>",
    unsafe_allow_html=True,
)

#st.title("Operon Finder")
#st.markdown("Cluster genes into operons")
st.markdown('''
    <img alt="lab-logo" src="https://gistcdn.githack.com/tejasvi8874/79d777264bf5d66f80ca3fd02f509b8d/raw/b01e6015f989263c608b1d3b6bb0c9de0e9f9530/lab-logo.svg">
    <img alt="operon-logo" src="https://gistcdn.githack.com/tejasvi8874/79d777264bf5d66f80ca3fd02f509b8d/raw/a5b80c86ba31b40c64f2c5409dda105992ce43bc/operon-logo.svg">
    <img alt="iitg-logo" src="https://upload.wikimedia.org/wikipedia/en/1/12/IIT_Guwahati_Logo.svg">
    <center><h4>Prokaryotic Operon Prediction Web Server</h4></center>
    ''', unsafe_allow_html=True)

manual = "Specify genome ID"
search = "Search genomes"

if streamlit_cloud:
    while True:
        try:
            with PidFile('.setup_lock'):
                setup()
            break
        except PidFileError:
            sleep(0.1)


query_genome_id = st.experimental_get_query_params().get("genome_id", [None])[0]
genome_id_option = st.radio("Select genome", (search, manual), index=(1 if query_genome_id else 0), on_change=reset_submit_query_genome_id)

if genome_id_option == search:
        sample_organisms = defaultdict(lambda: (None, None))
        for p in Path(".json_files").glob("*/compare_region.json.gz"):
            genome_name_file = p.parent.joinpath('genome_name.txt')
            if not genome_name_file.exists():
                genome_name = loads(p.parent.joinpath('genome.json').read_bytes())[0]["genome_name"]
                genome_name_file.write_text(genome_name)
            else:
                genome_name = genome_name_file.read_text()

            if p.parent.name == '224911.5':
                continue
            sample_organisms[ genome_name ] = p.parent.name, p.parent.name.split('.')[0]
        sample_organisms = dict(sorted(sample_organisms.items()))
        for species_name, genome_ids in valid_organisms():
            a_genome_id = next(iter(genome_ids)) # Don't uneccessarily mutate the variable with set.pop. Here the ret value is cached so errors on next run
            if a_genome_id == '224911.5':
                continue
            sample_organisms.setdefault(species_name, (a_genome_id, a_genome_id.split('.')[0]))

        # Prevent model inference on local machine
        # if not streamlit_cloud:
        #     ...
        organism_selection = st.selectbox(
            "Choose organism", sample_organisms, index=0, help="Press Backspace key to change search query", on_change=reset_submit_query_genome_id
        )
        genome_id, genome_organism_id = sample_organisms[organism_selection]

        if not genome_id:
            st.error(
                "It may take long to fetch external data for custom organism during first query."
            )

        genome_id = genome_id or get_genome_id(genome_organism_id)

        if not genome_id:
            logger.error(genome_organism_id)
            st.error("No compatible genomes found in PATRIC and STRING database.")
else:
    genome_id = st.text_input(
        "Genome ID",
        query_genome_id or "262316.17",
        help="Must be available in PATRIC and STRING databases.",
        on_change=reset_submit_query_genome_id
    ).strip()
    if re.match(r"\d+\.\d+", genome_id):
        try:
            for url in (
                f"https://patricbrc.org/api/genome/{genome_id}",
                f"https://stringdb-static.org/download/protein.links.v11.5/{genome_id.split('.')[0]}.protein.links.v11.5.txt.gz",
            ):
                if not get_session().head(url).ok:
                    genome_id = None
                    st.error(
                        "This genome ID is not supported. Try searching for the organism name instead."
                    )
        except ConnectionError:
            logger.critical("Connection error")
    else:
        genome_id = None
        st.error("Invalid Genome ID format. E.g. 262316.17")

    
if genome_id:
    # link_button(f"Genome details: {genome_id}", f"https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features")
    st.markdown(
        f"**Genome details:** [`{genome_id}`](https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features)<br><br><hr>", unsafe_allow_html=True
    )


s_chk = st.checkbox("Run", value=query_genome_id, key='submit_key')
submit = s_chk if genome_id else False

def br(times=1):
    for _ in range(times):
        st.write("<br/>", unsafe_allow_html=True)

if genome_id:
    br()

    full_data, gene_locations, approximated_refseqs, refseq_locus_tag_present = to_pid(genome_id)
    if not refseq_locus_tag_present:
        st.warning("The RefSeq approximations used for this genome might reduce the operon prediction accuracy due to absence of the corresponding annotations in the PATRIC database.")
    if not full_data:
        submit = False
        st.error(
            "This genome is not supported. Try searching for the organism name instead."
        )
    df = pd.DataFrame.from_dict(
        full_data, orient="index", columns=["RefSeq", "Description", "Protein", "SequenceAccession"]
    )
    df.index = df.index.astype(int)
    df = df.sort_index()
    # df.index.rename("PATRIC ID", inplace=True)

    if submit:
        with st.expander("Input genes"):
            st.dataframe(df)

if submit:
    probs = operon_probs(genome_id, len(full_data))
    if pid_uniprot_map := get_pid_uniprot_map(genome_id):
        df["UniProt"] = pd.Series(pid_uniprot_map)
    df["Confidence"] = pd.Series(probs)
    df["Intergenic distance"] = pd.Series([None]*len(df.index))

    operons = []
    with st.expander("Filter operons", True):
        br()
        min_prob = st.slider(
            "Confidence threshold",
            min_value=.0,
            max_value=1.,
            value=0.5,
            step=0.01,
        )

        clusters = operon_clusters(list(full_data.keys()), min_prob, probs)

        # clusters = [{998, 999, 1002}, {1001, 1002, 1003}, {1006, 1007}, {999, 1002, 1010, 1011, 1012}]

        min_len = min((len(c) for c in clusters), default=0)
        max_len = max(max((len(c) for c in clusters), default=0), min_len+1) # Streamlit throws error in range slider with equal min and max

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
            refseq_prefill = ', '.join(df.loc[iter(clusters[0]), "RefSeq"]).upper()
            
            contain_all = st.checkbox("All of the genes",help=refseq_help)
            if contain_all:
                must_pegs_text = st.text_area(
                    refseq_input_label,
                    refseq_prefill,
                    key="all",
                )
                must_pegs = {p.lower().strip() for p in must_pegs_text.split(",")}

            contain_any = st.checkbox("Atleast one of the genes",help=refseq_help)
            if contain_any:
                any_pegs_text = st.text_area(
                    refseq_input_label,
                    refseq_prefill,
                    key="any",
                )
                any_pegs = {p.lower().strip() for p in any_pegs_text.split(",")}

            contain_keyword = st.checkbox("Gene description keywords", value=False, help="Filter operons by contained gene's function descriptions")
            if contain_keyword:
                desc_keyword_txt = st.text_input("Enter keywords", ' '.join(df.loc[iter(clusters[0]), "Description"][:2]).lower())
                keywords = query_keywords(desc_keyword_txt)

            body: list[str] = []
            for i, cluster in enumerate(clusters):
                if not (
                    cluster_size_range[0] <= len(cluster) <= cluster_size_range[1]
                    and must_pegs.issubset({full_data[j].n_refseq.lower() for j in cluster})
                    and (
                        not any_pegs
                        or any([full_data[j].n_refseq.lower() in any_pegs for j in cluster])
                    )
                    and (
                        not keywords
                        or keywords.issubset(query_keywords(' '.join(df.loc[iter(cluster), "Description"]).lower()))
                    )
                ):
                    continue

                dfx = df.loc[sorted(cluster)]

                # Confidence score is a gene connected to next gene. Last gene will technically have low score so set it to previous score to keep it meaningful
                idx_second_max, idx_max  = np.partition(dfx.index.values, -2)[-2:]
                dfx.loc[idx_max, "Confidence"] = dfx["Confidence"][idx_second_max]

                # TODO: find "next" instead of pid+1 since some pid missing
                for pid in dfx.index:
                    dfx.loc[pid, "Intergenic distance"] = gene_locations[pid+1].start - gene_locations[pid].end - 1 if pid+1 in dfx.index else '-'

                operons.append((i, dfx))

            detailed = st.checkbox(f"Detailed view", value=True, help="Show additional information") 

    if operons:
        st.info(f"{len(operons)} operons found")
        
        def download():
            html = download_bytes_js(
                        "\n".join(
                            ["\t".join(["PATRIC ID", *df.columns.tolist()])] +
                            [f"Operon {i+1}\n" + dfx.to_csv(header=False, sep="\t")
                                for i, dfx in operons]
                        ).encode(),
                        f"{genome_id}-operon.tsv",
                        "text/tab-separated-values"
                    )
            components.html(html, height=0)
        st.button(f"Save results 📥", on_click=download)

        show_all = st.session_state.get("show_all", False)
        for i, (operon_num, dfx) in enumerate(operons):
            st.markdown(f"#### Operon {operon_num+1}")
            approximation_note = '"Approximate RefSeq assignment"'
            render_dfx = deepcopy(dfx)
            render_dfx["RefSeq"] = dfx["RefSeq"].apply(
                    lambda r: f"""<a target="_blank" href="https://www.ncbi.nlm.nih.gov/refseq/?term={r}">{r.upper() + "</a>"
                            + ("<br><a style='text-decoration: none;' target='_self' href='javascript:alert(" + approximation_note + ")'><span title=" + approximation_note + ">⚠️</span></a>" if r in approximated_refseqs else '')
                        }"""
            )
            render_dfx["Protein"] = dfx["Protein"].apply(
                lambda r: f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/protein/?term={r}">{r}</a>'
            )
            render_dfx["UniProt"] = dfx["UniProt"].apply(
                lambda r: f'<a target="_blank" href="https://www.uniprot.org/uniprot/{r}">{r}</a>'
            )
            del render_dfx["SequenceAccession"]
            if not detailed:
                for c in ["Confidence", "Intergenic distance", "UniProt"]:
                    del render_dfx[c]
            st.write(
                render_dfx.to_html(
                    justify="center",
                    escape=False,
                    classes=["table-borderless"],
                    border=0,
                    formatters={'Confidence': lambda x: f'<b style="background-color: hsl({120*x}, 100%, 75%)">{x:.2f}</b>'} if detailed else None,
                ),
                unsafe_allow_html=True,
            )
            if detailed:
                with st.columns([1])[0]:
                    show_dna = st.checkbox("DNA", key=f"dna-check-{operon_num}")
                if show_dna:
                    start_idx, end_idx =  st.select_slider("Select gene range", options = sorted(dfx.index), value=(min(dfx.index), max(dfx.index)), key=f"dna-range-{operon_num}")
                    start = gene_locations[start_idx].start-1
                    end = gene_locations[end_idx].end
                    p1, p2 = st.empty(), st.empty()

                    fasta = loads(get_output(
                        f'''https://p3.theseed.org/services/data_api/jbrowse/genome/{genome_id}/features/{
                        df["SequenceAccession"].loc[start_idx]
                        }?reference_sequences_only=false&start={
                        gene_locations[start_idx].start-1
                        }&end={
                        gene_locations[end_idx].end
                        }'''))["features"][0]["seq"][:end-start]
                    p1.markdown(f"<textarea readonly rows=10 style='width:100%'>{fasta}</textarea>", unsafe_allow_html=True)
                    p2.download_button(label='Download 📥', file_name=f'{genome_id}-operon-{operon_num+1}-dna.txt', key=f"dna-download-{operon_num}", data=fasta)

            if i >= 10 and not show_all:
                break

        st.markdown("---")
        st.checkbox( f"Show all operons", key="show_all")
    else:
        st.error(f"No matching clusters found")

st.markdown(f'''
<center>
Developed and maintained by:<br>
<b>
Tejasvi Singh Tomar, Pratik Dasgupta, and Shankar Prasad Kanaujia<br>
<a href="https://www.iitg.ac.in/spkanaujia">Structural and Computational Biology Laboratory (SCBL)</a></br>
<a href="https://www.iitg.ac.in/biotech/">Department of Biosciences and Bioengineering</a><br>
<a href="https://www.iitg.ac.in/">Indian Institute of Technology Guwahati</a><br>
Guwahati - 781 039, Assam, India<br>
Copyright Ⓒ 2022 Operon Finder <br>
</b>
<a href="https://clustrmaps.com/site/1bowi" title="visitor map"><img src="//www.clustrmaps.com/map_v2.png?d=1laoCLUJldQlosO5bh73b689snhMb-hqnJynyQ_Ez5s&cl=ffffff"></a></br>
Genome count: 9140<br>
<!-- <img src="https://hitwebcounter.com/counter/counter.php?page=8000876&style=0005&nbdigits=5&type=page&initCount=113" title="Free Counter" Alt="web counter"   border="0" /><br>
Visitor count -->
</center>''', unsafe_allow_html=True)

# pr.disable()
# s = io.StringIO()
# sortby = SortKey.CUMULATIVE
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print(s.getvalue())

# profiler.stop()
# profiler.print()
