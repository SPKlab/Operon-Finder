# from argparse import ArgumentParser


# from gooey import Gooey

# @Gooey()
# def main():
#     parser = ArgumentParser()
#     parser.add_argument(
#         "-i",
#         "--input_file", 
#         required=True,
#         help="File to be converted"
#     )
#     parser.add_argument(
#         "-o",
#         "--output_file",
#         required=True,
#         help="Path for the converted file"
#     )
#     args = parser.parse_args()
#     convert(args.input_file, args.output_file)


# def convert(infile, outfile):
#     pass
#     # (
#     #     ffmpeg
#     #     .input(infile)
#     #     .output(outfile, vcodec="copy", acodec="copy")
#     #     .run()
#     # )
    

# if __name__ == "__main__":
#     main()

from argparse import ArgumentParser
from get_json import operon_clusters
from gooey import Gooey
import sys, os

import warnings
warnings.filterwarnings('ignore')

@Gooey(program_name='Operon Finder')
def main():
    parser = ArgumentParser(description="Cluster genes into operons") 
    parser.add_argument('-g', '--Genome_ID', help='Genome ID obtained from PATRIC', default='83332.12')
  
    start_gene_id = 998
    end_gene_id = 1007
    parser.add_argument('-p', '--PEG_IDs', help='Comma separated Protein Encoding Gene numbers', default=', '.join(map(str, range(start_gene_id, end_gene_id+1))))
    args = parser.parse_args()

    peg_nums = [int(p) for p in args.PEG_IDs.split(',') if p.strip().isdigit()]
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    clusters = operon_clusters(args.Genome_ID, peg_nums)

    body = ["CLUSTERED GENES:\n"]
    for i, cluster in enumerate(clusters):
        body.append(f"Cluster {i+1}:")
        body.append(', '.join([str(j) for j in cluster]) + '\n')
    sys.stdout = sys.__stdout__
    print('\n'.join(body))

main()