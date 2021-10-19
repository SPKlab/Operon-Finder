from typing import Optional
from get_json import operon_clusters
import streamlit as st

st.title('Operon Finder')

st.write("Cluster genes into operons")

genome_id = st.text_input('Genome ID', '83332.12')

manual_label = 'Manual IDs'
specify_range = 'Specify Range'
method = st.radio("Protein Encoding Gene ids",(manual_label, specify_range, "Entire genome"))

start_gene_id = 998
end_gene_id = 1015

peg_nums: Optional[frozenset[int]]

if method == manual_label:
    peg_num_text = st.text_area('Comma separated Protein Encoding Gene numbers', ', '.join(map(str, range(start_gene_id, end_gene_id+1))))
    peg_nums = frozenset({int(p) for p in peg_num_text.split(',') if p.strip().isdigit()})
elif method == specify_range:
    start = st.number_input('Start Protein Encoding Gene number', value=start_gene_id, min_value=1)
    end = st.number_input('End Protein Encoding Gene number', value=end_gene_id, min_value=1)
    peg_nums = frozenset(range(int(start), int(end)+1))
else:
    peg_nums = None
    
submit = st.checkbox("Submit")
st.write("---")

st.spinner("Processing..")
    
if submit:
    clusters = operon_clusters(genome_id, peg_nums)

    # clusters = [{998, 999, 1002}, {1001, 1002, 1003}, {1006, 1007}, {999, 1002, 1010, 1011, 1012}]

    min_len = min(len(c) for c in clusters)
    max_len = max(len(c) for c in clusters)

    cluster_size_range = 1, float('inf')
    must_pegs: set[int] = set()
    any_pegs: Optional[set[int]] = None

    st.write("## Clustered genes")
    apply_filter = st.checkbox("Filter clusters")
    if apply_filter:
        st.write("---")
        cluster_size_range = st.slider("Cluster Size Range", min_value=min_len, max_value=max_len,value=(min_len, max_len), step=1)

        contain_all = st.checkbox("Contains all of the specified genes")
        if contain_all:
            must_pegs_text = st.text_area('Comma separated Protein Encoding Gene numbers',
            '193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 653, 654, 655, 656, 657, 658, 659, 1526, 1527, 2192, 2196, 2197, 2198, 2199, 2200, 2201, 2202, 2203, 2672, 3899, 3900, 3901, 3902, 3903, 3904, 3905, 3906 ',
             help=', '.join(map(str, range(3, 24))), key="all")
            must_pegs = {int(p) for p in must_pegs_text.split(',') if p.strip().isdigit()}

        contain_any = st.checkbox("Contains atleast one of the specified genes")
        if contain_any:
            any_pegs_text = st.text_area('Comma separated Protein Encoding Gene numbers', help=', '.join(map(str, range(3, 24))),key="any")
            any_pegs = {int(p) for p in any_pegs_text.split(',') if p.strip().isdigit()}
        st.write("---")

    body: list[str] = []
    for i, cluster in enumerate(clusters):
        if not (cluster_size_range[0] <= len(cluster) <= cluster_size_range[1]
            and must_pegs.issubset(cluster)
            and (not any_pegs or (any_pegs.intersection(cluster)))
        ):
            continue
            
        
        body.append(f"#### Cluster {i+1}")
        body.append(', '.join([f"***{j}***" if (j in must_pegs or (any_pegs and j in any_pegs)) else str(j) for j in cluster]))
    st.markdown('\n'.join(body))