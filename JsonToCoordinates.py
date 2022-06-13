from concurrent.futures import ProcessPoolExecutor, as_completed
from string import ascii_letters
from gzip import decompress
from helpers import get_output, normalize_refseq, string_id_n_refseq_pairs, to_pid, get_prefix_counter, logger, StringRefseq
import multiprocessing
import json
from tempfile import NamedTemporaryFile
from pathlib import Path
import os
import random

def valid(f):
        if f is not None and  len(f) > 1 and f[1] is not None and len(f[1]) > 1 and f[1][1] is not None and not isinstance(f[1][1], int) and 'blast' not in f[1][1].lower() and 'PLF' in f[1][1]:
                return True
        return False
def PGFvalid(f):
        if f is not None and  len(f) > 0 and f[0] is not None and len(f[0]) > 1 and f[0][1] is not None and not isinstance(f[0][1], int) and 'blast' not in f[0][1].lower() and 'PGF' in f[0][1]:
                return True
        return False

colors = ["vlred", "vlblue", "lyellow", "orange", "gold", "lgrey", "lbrown", "purple", "maroon",  "salmon", "magneta"]

random.seed(7)
random.shuffle(colors)

import re

string_pat = re.compile(r'^\d+?\.(.+?) \d+?\.(.+?) (\d+)', re.MULTILINE)
def parse_string_scores(genome_id: str)->dict[str,float]:
    string_refseq_obj = StringRefseq(genome_id)
    refseq_order_pid = string_refseq_obj.refseq_order_pid
    organism = genome_id.split('.')[0]

    string = {} 
    for match in string_pat.finditer(decompress(get_output(f"https://stringdb-static.org/download/protein.links.v11.5/{organism}.protein.links.v11.5.txt.gz")).decode()):
        g1, g2, score = match.groups()
        r1, r2 = [string_refseq_obj.get_refseq(g) for g in (g1, g2)]
        if r1 and r2 and (refseq_order_pid[r1][0] + 1 == refseq_order_pid[r2][0]):
            string[f"fig|{genome_id}.peg.{refseq_order_pid[r1][1]}"] = float(score)/1000

    assert string
    logger.info("Parsed STRING scores")
    return string

def to_coordinates(compare_region_data: list, genome_id: str, progress_clb=None) -> str:
    str_json_file = Path(f'.json_files/string/{genome_id}.json')
    if str_json_file.exists():
        with open(str_json_file) as f:
           string = json.load(f)
    else:
        string = parse_string_scores(genome_id)
        str_json_file.parent.mkdir(parents=True, exist_ok=True)
        with open(str_json_file, 'w') as f:
           json.dump(string, f)

###
    directory = f'coords/{genome_id}'
    os.makedirs(directory, exist_ok=True)
    out = NamedTemporaryFile(mode='w+t', prefix=genome_id, delete=False, dir=directory)
    out.close()
    out_file_name = out.name
    out_lock = multiprocessing.Manager().Lock()
    

    with ProcessPoolExecutor(os.cpu_count()*4) as executor:
        for i, r in enumerate(as_completed([executor.submit(writer, figdata, out_lock, out_file_name, string) for figdata in compare_region_data])):
            r.result()
            if progress_clb:
                progress_clb(0.5 + (i+1)/len(compare_region_data)/10)

    return out_file_name

def writer(figdata: str, lock, out_file_name: str, string: dict[str, float]):
    data = figdata
    if data['result'] is None:
            if data['error'] is not None:
                    raise Exception("JSON not fetched correctly")
            return
    num_pegs = len(data['result'][0][0]['features']) - 1
    if num_pegs == 0:
            return

### Assign color #########
    count = 0
    colormap = {}
    strands = {}
    is_query = 0
    query2 = False
    for r in data['result'][0]: # for all genomes
            is_query += 1
            if is_query == 1:
                    real_pin = r['pinned_peg']
                    real_pin_strand = r['pinned_peg_strand']

            r_pin = r['pinned_peg']
            for p in range(len(r['features'])):
                    if r['features'][p]['fid'] == r_pin:
                            break
            p += 1 # the one after the pin 
            if (p < len(r['features']) -1) and real_pin_strand == '-' and r['pinned_peg_strand'] == '-':
                            r['pinned_peg'] = r['features'][p]['fid']

            features = r['features'][:-1]
            if real_pin_strand == '-' and r['pinned_peg_strand'] == '-':
                    features = features[::-1]
            if real_pin_strand == '-' and r['pinned_peg_strand'] == '+':
                    if p >=2:
                            r['pinned_peg'] = r['features'][p-2]['fid']
        
            for f in features:
                    if f['type'] != 'peg':
                            continue
                    if is_query == 1:
                            if f['fid'] == r['pinned_peg']:
                                    query2 = True
                                    try:
                                        pgf = f['attributes'][0][1]
                                    except IndexError as e:
                                            continue
                                    colormap[pgf] = 'blue'
                            elif query2:
                                    query2 = False
                                    if PGFvalid(f['attributes']):
                                            pgf = f['attributes'][0][1]
                                            colormap[pgf] = 'red'
                    if PGFvalid(f['attributes']):
                        pgf = f['attributes'][0][1] #  actually 
                        if pgf not in colormap:
                            if count < 11:
                                colormap[pgf] = colors[count]
                                count += 1
                            else:
                                colormap[pgf] = 'black'

                    # now strand
                        if is_query == 1:
                                if pgf not in strands: # just for the query
                                    strands[pgf] = f['strand']

###################################################
    r = data['result'][0][0] #only the query genome
    out_lines = []
# for every new figure:
    out_lines.append('new\tfigure\n')
# write focus peg names on line 1, then start,end,strand,color,name for each arrow on line 2
    name1 = real_pin
    name2 = '.'.join(name1.split('.')[:-1]) + '.' + str(int(name1.split('.')[-1])+1)



    features = r['features'][:-1]
#get offset for query peg 
    feat = [] 
    for f in features:
            if f['type'] == 'peg':
                    feat.append(f)
    features = feat

    if len(features) ==0:
            return
    for f in features:
            if f['fid'] == name1:
                    q_size = int(f['size'] -1)
                    q_strand = f['strand']
                    if q_strand == '+':
                            q_start,q_end = f['beg'], f['end']
                    else:
                            q_start, q_end = f['end'], f['beg']

                    if features[0]['strand'] == '+':
                            f_start, f_end = features[0]['beg'], features[0]['end']
                    else:
                            f_start, f_end = features[0]['end'], features[0]['beg']
                    if features[-1]['strand'] == '+':
                            l_end = features[-1]['end']
                    else:
                            l_end = features[-1]['beg']

                    f_start = int(f_start)
                    f_end = int(f_end)
                    l_end = int(l_end)
                    total = abs(l_end - f_start)

                
                    q_start = int(q_start)
                    q_end = int(q_end)
                    if q_strand == '+':
                            q_offset = int(abs(q_start - f_start)/float(total) * 300)
                    else:
                            q_offset = int(abs(q_end- l_end)/float(total) * 300)
 
    for r in data['result'][0]:
            score = 0.1
            if name1 in string:
                    score = string[name1]
            bname1 = real_pin
            var = 1
            bname2 = '.'.join(bname1.split('.')[:-1]) + '.' + str(int(bname1.split('.')[-1])+var) 
            features = r['features'][:-1] 

        # get offset for current query
            cname = r['pinned_peg']
            feat = []
            for f in features:
                    if f['type'] == 'peg':
                            feat.append(f)
            features = feat
            if len(features) == 0:
                    continue
            out_lines.append(name1 + '\t' + name2 + '\t' + str(score) + '\n')
            if r['pinned_peg_strand'] == '-':
                    features = features[::-1]

        # no point in having those in a loop?
            if r['pinned_peg_strand'] == '-' and len(features) > 0:
                    if features[0]['strand'] == '-':# or features[0]['strand'] == '-':
                            f_start, f_end = features[0]['beg'], features[0]['end']
                    else:
                            f_start, f_end = features[0]['end'], features[0]['beg']
            elif len(features) > 0:
                    if features[0]['strand'] == '+':# or features[0]['strand'] == '-':
                            f_start, f_end = features[0]['beg'], features[0]['end']
                    else:
                            f_start, f_end = features[0]['end'], features[0]['beg']
            if r['pinned_peg_strand'] == '-':
                    if features[-1]['strand'] == '-':# or features[-1]['strand'] == '-':
                            l_end = int(features[-1]['end'])
                    else:
                            l_end = int(features[-1]['beg'])
            else:
                    if len(features) > 0 and features[-1]['strand'] == '+':# or features[-1]['strand'] == '-':
                            l_end = int(features[-1]['end'])
                    elif len(features) > 0:
                            l_end = int(features[-1]['beg'])
                
            fh_start = int(f_start)
            fh_end = 'x'
            sh_start = 'x'
            f_end = int(f_end)
            sh_end = int(l_end)
            total = abs(sh_end - fh_start)
            if total == 0:
                    break

            before = True
            focus = False
            after = False
            query2 = False
            c_start = 'x'
            c_end = 'x'
            total_after = 0
            c_strand = r['pinned_peg_strand']
            for f in features:
                    if f['fid'] == r['pinned_peg']:
                            before = False
                            focus = True
                            after = False
                    if before:
                            if r['pinned_peg_strand'] == '-':
                                    if f['strand'] == '-':# or f['strand'] == '-':
                                            c_start,c_end = f['beg'], f['end']
                                    else:
                                            c_start, c_end = f['end'], f['beg']
                            else:
                                    if f['strand'] == '+':# or f['strand'] == '-':
                                            c_start,c_end = f['beg'], f['end']
                                    else:
                                            c_start, c_end = f['end'], f['beg']
                            c_start = int(c_start)
                            c_end = int(c_end)
                    elif focus:
                            focus = False
                            query2 = True
                            after = True
                        # use the end of before to establish total_before 
                            fh_end = c_end
                            total_before = 0
                            if fh_end != 'x': # means there is a before
                                    total_before = abs(c_end - fh_start)        
                            if r['pinned_peg_strand'] == '-':
                                    if f['strand'] == '-':
                                            c_start,c_end = f['beg'], f['end']
                                    else:
                                            c_start, c_end = f['end'], f['beg']
                            else:
                                    if f['strand'] == '+':
                                            c_start,c_end = f['beg'], f['end']
                                    else:
                                            c_start, c_end = f['end'], f['beg']

                            cq_start = int(c_start)
                            cq_end = int(c_end)
                    elif after:
                            if query2:
                                    query2 = False
                                # use current end with c_start to establish total_query
                                    if r['pinned_peg_strand'] == '-':
                                            if f['strand'] == '-':
                                                    c_end = f['end']
                                            else:
                                                    c_end = f['beg']
                                    else:
                                            if f['strand'] == '+':
                                                    c_end = f['end']
                                            else:
                                                    c_end = f['beg']
                                    c_end = int(c_end)
                                    total_query = abs(c_end - cq_start) # c_start is the focus peg start
                            else:
                                # use current start with l_end to establish total_after
                                    if r['pinned_peg_strand'] == '-':
                                            if f['strand'] == '-':
                                                    c_start = f['beg']
                                            else:
                                                    c_start = f['end']
                                    else:
                                            if f['strand'] == '+':
                                                    c_start = f['beg']
                                            else:
                                                    c_start = f['end']
                                    c_start = int(c_start)

                                    if sh_start == 'x':
                                            sh_start = c_start
                                    total_after = abs(l_end - sh_start)

            before = True
            focus = False
            after = False
            query2 = False
            focus_300_offset = 0
            focus_150_size = 0
            for f in features:
                # check if peg has known pgf 
                    if f['type'] != 'peg':
                            continue
                    if PGFvalid(f['attributes']): #valid
                            pgf = f['attributes'][0][1] # PGF, 0 1 for plf
                    else:
                            continue
                                        
                # assign color to peg 
                    c = ''
                    if pgf in colormap:
                            c += colormap[pgf]
                    else:
                            c += 'black'
                
                    strand = f['strand']

                    size = int(f['size'] -1)
                    if r['pinned_peg_strand'] == '-':
                            if f['strand'] == '-':
                                    start,end = f['beg'], f['end']
                            else:
                                    start, end = f['end'], f['beg']
                    else:
                            if f['strand'] == '+':
                                    start,end = f['beg'], f['end']
                            else:
                                    start, end = f['end'], f['beg']

                    if total == 0:
                            logger.info(bname1)
                            continue
                    start = int(start)
                    end = int(end)

                    name = f['fid']
                    if pgf in strands:
                            strand = strands[pgf]
                    out_size = 0

                # for genes before blast peg, scale offset and size to 75 
                    if f['fid'] == r['pinned_peg']:
                            before = False
                            focus = True
                            after = False
                    if before and total_before != 0:
                            if c_strand == '+' or c_strand == '-':
                                    offset = int(abs(start - fh_start)/float(total_before) * 50) 
                            else:
                                    offset = int(abs(end- fh_end)/float(total_before) * 50)# + 225
                            out_size = float(size)/total_before * 50
                # for focus peg and the one after, scale offset and size to 150, add 75 to offset 
                    elif focus:
                            focus = False
                            offset = 50
                            out_size = float(size)/total_query * 200
                            focus_150_size = out_size
                            query2 = True
                            after = True
                # for after the focus peg, scale offset and size to 75, add 225 to offset 
                    elif after:
                            if query2:
                                    query2 = False
                                    if c_strand == '+' or c_strand == '-':
                                            offset = int(abs(start - cq_start)/float(total_query) * 200) 
                                    else:
                                            offset = int(abs(end- cq_end)/float(total_query) * 200) 
                                    offset += 50
                                    out_size = float(size)/total_query * 200
                            else:
                                    if total_after == 0:
                                            continue
                                    if c_strand == '+' or c_strand == '-':
                                            offset = int(abs(start - sh_start)/float(total_after) * 50) +250
                                    else:
                                            offset = int(abs(end- sh_end)/float(total_after) * 50)# - 225
                                    out_size = float(size)/total_after * 50
                    if offset < 0 or out_size < 0 or strand is None or c is None or name is None:
                            continue
                    out_lines.append(str(int(offset)) + '\t' + str(int(out_size)) + '\t' + strand + '\t' + c + '\t' + name + '\t')
            out_lines.append('\n')
    with lock:
            with  open(out_file_name, 'a') as f:
                    f.writelines(out_lines)
