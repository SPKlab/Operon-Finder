from get_json import operon_clusters
import streamlit as st

st.title('Operon Finder')

st.write("Cluster genes into operons")
genome_id = st.text_input('Genome ID', '83332.12')

manual_label = 'Manual IDs'
method = st.radio("Protein Encoding Gene ids",(manual_label, 'Specify Range'))

start_gene_id = 998
end_gene_id = 1007

with st.form("my_form"):
    if method == manual_label:
        peg_num_text = st.text_area('Comma separated Protein Encoding Gene numbers', ', '.join(map(str, range(start_gene_id, end_gene_id+1))))
        peg_nums = [int(p) for p in peg_num_text.split(',') if p.strip().isdigit()]
    else:
        start = st.number_input('Start Protein Encoding Gene number', value=start_gene_id)
        end = st.number_input('End Protein Encoding Gene number', value=end_gene_id)
        peg_nums = list(range(int(start), int(end)+1))
        
    submitted = st.form_submit_button("Submit")
    if submitted:
        clusters = operon_clusters(genome_id, peg_nums)
        body = ["## Clustered genes"]
        for i, cluster in enumerate(clusters):
            body.append(f"#### Cluster {i+1}")
            body.append(', '.join([str(j) for j in cluster]))
        st.markdown('\n'.join(body))



