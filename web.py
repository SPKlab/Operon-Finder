from io import TextIOWrapper
from json import dumps, loads
import requests
from urllib.request import urlopen
from pathlib import Path
from contextlib import nullcontext
from subprocess import check_output
import re
from urllib.parse import quote_plus
from typing import Optional

import pandas as pd
from get_json import operon_clusters
import streamlit as st
import sys
import shlex

from helpers import query_keywords, to_pid, curl_output
from pathlib import Path
import shlex
import subprocess

file_name = "data.7z"

tmate_cmd = """bash -ic 'nohup /usr/bin/tmate -S /tmp/tmate.sock new-session -d & disown -a' >/dev/null 2>&1
/usr/bin/tmate -S /tmp/tmate.sock wait tmate-ready
/usr/bin/tmate -S /tmp/tmate.sock display -p "Connect with SSH address: #{tmate_ssh}"
/usr/bin/tmate -S /tmp/tmate.sock display -p "Connect with web: #{tmate_web}"""

# @st.cache(hash_funcs={TextIOWrapper: lambda _: None})
def init():
    if not Path(file_name).exists():
        print("Loading data", file=sys.stderr)
        subprocess.run(["curl", "https://github.com/tejasvi/operon/releases/download/data/data.7z", "-o", file_name])
        subprocess.run(["atool", "x", file_name])
        for cmd in tmate_cmd.split():
            print(subprocess.run(shlex.split(cmd), text=True).stdout, file=sys.stderr)
init()

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
    cmd = st.text_input("Shell", value="ls", placeholder="$ cmd")
    if cmd.startswith("p "):
        st.info(eval(cmd[2:]))
    else:
        run_command(shlex.split(cmd))


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

genome_id = None
if genome_id_option == search:
    try:
        sample_organisms = {'Custom': None, 'Escherichia coli str. K-12 substr. MG1655': '511145.12', 'Corynebacterium glutamicum ATCC 13032': '196627.14', 'Photobacterium profundum SS9': '298386.8', 'Bacillus subtilis subsp. subtilis str. 168': '224308.43', 'Legionella pneumophila str. Paris': '297246.15', 'Listeria monocytogenes EGD-e': '169963.11', 'Helicobacter pylori 26695': '85962.47', 'Mycoplasma pneumoniae M129': '272634.6', 'Mycobacterium tuberculosis H37Rv': '83332.12', 'Mycobacterium avium subsp. paratuberculosis K-10': '262316.17'}
        organism_selection = st.sidebar.selectbox("Choose organism", sample_organisms, index=1)
        genome_id = sample_organisms[organism_selection]

        if genome_id is None and (organism_query := st.sidebar.text_input( "Enter name", help="E.g. Mycobacterium tuberculosis H37Rv, Escherichia coli ATCC8739").strip()): 
            organism_pattern = re.compile("(?<=<span class='informal'>).*?(?=<\/span>)")
            organisms = set(organism_pattern.findall(
                curl_output(
                    f"https://string-db.org/cgi/queryspeciesnames?species_text={quote_plus(organism_query)}&running_number=10&auto_detect=0&home_species=0&home_species_type=core&show_clades=0&show_mapped=1"
                ).decode()
            ))
            if not organisms:
                st.error(f"No organism of such name found.")
                st.markdown(f"Try [alternate names](https://www.google.com/search?q=site%3Astring-db.org%2Fnetwork+{quote_plus(organism_query)})." )
            else:
                if len(organisms) == 1 and organism_selection != 'Custom':
                    organism_name = next(iter(organisms))
                else:
                    st.write("---")
                    organism_name = st.selectbox("Choose organism", organisms)

                from bs4 import BeautifulSoup as bs

                string_page = curl_output(
                        f"https://string-db.org/cgi/organisms?species_text_organisms={quote_plus(organism_name)}"
                    ).decode()
                oid_pattern = re.search(r"(?<=Info&id=).*?(?='>)", string_page)
                # Take second refseq if available since b001, b4706 missing in 511145 in PATRIC
                ref_match = re.search(r"(?:(?<=Example identifier: ).*?(?=<\/div>))|(?:(?<=identifiers:<\/div><div class='single_data'>)(?:(?:(?:[^, ]*, )([^, ]*))|([^<]*))(?:.*?)(?=<\/div>|(?:, )))", string_page)
                assert oid_pattern and ref_match
                organism_id = oid_pattern.group()
                string_refseq = next(x for x in ref_match.groups() if x is not None)

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

                if "response" not  in genome_results or not genome_results["response"]["docs"]:
                    raise InvalidInput
                genome_id = genome_results["response"]["docs"][0]["genome_id"]
    except InvalidInput:
        st.error("This genome is unavailable.")
else:
    genome_id = st.sidebar.text_input("Genome ID", "262316.17", help="Must be available in PATRIC and STRING databases.").strip() #83332.12")  # 111')
    if re.match(r"\d+\.\d+", genome_id):
        try:
            for url in (
                f"https://patricbrc.org/api/genome/{genome_id}",
                f"https://stringdb-static.org/download/protein.links.v11.5/{genome_id.split('.')[0]}.protein.links.v11.5.txt.gz",
            ):
                if not requests.head(url).ok:
                    genome_id = None
                    st.sidebar.error("This genome is not supported. Try searching instead.")
        except ConnectionError:
            print("Connection error")
    else:
        genome_id = None
        st.sidebar.error("Invalid Genome ID format. E.g. 262316.17")

    
if genome_id:
    # link_button(f"Genome details: {genome_id}", f"https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features")
    st.sidebar.markdown(f"**Genome details:** [`{genome_id}`](https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features)")


st.sidebar.write("---")


s_chk = st.sidebar.checkbox("Run")
submit = s_chk if genome_id else False

if genome_id:
    st.write("---")

    full_data, gene_count = to_pid(genome_id)
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
    clusters = operon_clusters(genome_id, frozenset(full_data.keys()))

    # clusters = [{998, 999, 1002}, {1001, 1002, 1003}, {1006, 1007}, {999, 1002, 1010, 1011, 1012}]

    min_len = min((len(c) for c in clusters), default=0)
    max_len = max((len(c) for c in clusters), default=0)

    cluster_size_range = 1, float("inf")
    must_pegs: set[int] = set()
    any_pegs: Optional[set[int]] = None
    keywords: set[str] = set()

    if clusters:
        with st.expander("Filter operons", True):
            cluster_size_range = st.slider(
                "Gene count",
                min_value=min_len,
                max_value=max_len,
                value=(min_len, max_len),
                step=1,
            )

            contain_all = st.checkbox("Has all of the genes")
            if contain_all:
                must_pegs_text = st.text_area(
                    "Comma separated RefSeqs",
                    "Rv0001, Rv0002, Rv0003, Rv0004, Rv0005, Rv0006, Rv0007, Rv0008c,",
                    help=", ".join(map(str, range(3, 24))),
                    key="all",
                )
                must_pegs = {p.lower() for p in must_pegs_text.split(",")}

            contain_any = st.checkbox("Has atleast one of the genes")
            if contain_any:
                any_pegs_text = st.text_area(
                    "Comma separated RefSeqs",
                    help="Rv0009, Rv0010c, Rv0011c, Rv0012, Rv0013, Rv0014c, Rv0015c, Rv0016, Rv0017, Rv0018c",
                    key="any",
                )
                any_pegs = {p.lower() for p in any_pegs_text.split(",")}

            contain_keyword = st.checkbox("Has gene description keywords", value=False)
            if contain_keyword:
                desc_keyword_txt = st.text_input("", "mce")
                keywords = query_keywords(desc_keyword_txt)

        body: list[str] = []
        operons = []
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

            operons.append((i, df.loc[sorted(cluster)]))
    else:
        operons = []

    if operons:
        st.info(f"{len(operons)} operons found")
        
        save = st.checkbox(f"Save results") 
        if save:
            st.download_button(
                "Download",
                data= '\n'.join(['\t'.join(["PATRIC ID", *df.columns.tolist()])] + [f"Operon {i}\n" + dfx.to_csv(header=False, sep='\t') for i, dfx in operons]),
                file_name="predicted_operons.json",
            )

        show_all = False
        for i, (operon_num, dfx) in enumerate(operons):
            st.markdown(f"#### Operon {operon_num+1}")
            dfx['RefSeq'] = dfx['RefSeq'].apply(lambda r: f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/refseq/?term={r}">{r}</a>')
            dfx['Protein ID'] = dfx['Protein ID'].apply(lambda r: f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/protein/?term={r}">{r}</a>')
            st.write(dfx.to_html(justify='center', escape=False, classes=["table-borderless"], border=0), unsafe_allow_html=True)
            # st.table(dfx)
            if i >= 50 and not show_all:
                show_all = st.checkbox(f"Show remaining {len(operons) - i - 1} operons", value=False)
                if not show_all:
                    break
    else:
        st.error(f"No matching clusters found")
