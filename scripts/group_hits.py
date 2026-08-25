#!/usr/bin/env python3
import argparse
import csv

def group_hits(blast_file, out_ltr, out_genic, context_range, ltr_range):
    # Read blast hits
    # Columns expected (from Snakefile): 6 qseqid sseqid sstart send sstrand pident qcovs
    hits = []
    try:
        with open(blast_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) >= 5:
                    qseqid = row[0]
                    sseqid = row[1]
                    sstart = int(row[2])
                    send = int(row[3])
                    sstrand = row[4]
                    hits.append({
                        'sseqid': sseqid,
                        'sstart': sstart,
                        'send': send,
                        'sstrand': sstrand
                    })
    except FileNotFoundError:
        pass
        
    # Sort hits by subject sequence and start coordinate
    hits.sort(key=lambda x: (x['sseqid'], min(x['sstart'], x['send'])))

    loc_ltr_lines = []
    loc_genic_lines = []
    
    anterior = 'nada'
    inicio = 0
    fim = 0
    fita = 'nada'
    possivelLTR = ''
    
    for hit in hits:
        curr_sseqid = hit['sseqid']
        curr_start = hit['sstart']
        curr_end = hit['send']
        curr_strand = hit['sstrand']
        
        # logic from minERVa_023.py
        if anterior == curr_sseqid and (inicio - context_range) < curr_start < (fim + context_range):
            if curr_strand == 'plus' or curr_start < curr_end:
                if inicio > curr_start:
                    inicio = max(1, curr_start - ltr_range)
                if fim < curr_end:
                    fim = curr_end + ltr_range
            else:
                if inicio > curr_end:
                    inicio = max(1, curr_end - ltr_range)
                if fim < curr_start:
                    fim = curr_start + ltr_range
            possivelLTR = ''
        else:
            if possivelLTR == '' and anterior != 'nada':
                loc_genic_lines.append(f"{anterior} {inicio}-{fim} {fita}\n")
            elif anterior != 'nada':
                loc_ltr_lines.append(possivelLTR)
                
            anterior = curr_sseqid
            if curr_strand == 'plus' or curr_start < curr_end:
                inicio = max(1, curr_start - ltr_range)
                fim = curr_end + ltr_range
                fita = curr_strand
            else:
                inicio = max(1, min(curr_start, curr_end) - ltr_range)
                fim = max(curr_start, curr_end) + ltr_range
                fita = curr_strand
            possivelLTR = f"{anterior} {inicio}-{fim} {fita}\n"
            
    if possivelLTR == '' and anterior != 'nada':
        loc_genic_lines.append(f"{anterior} {inicio}-{fim} {fita}\n")
    elif anterior != 'nada':
        loc_ltr_lines.append(possivelLTR)

    with open(out_ltr, 'w') as f_ltr:
        f_ltr.writelines(loc_ltr_lines)
        
    with open(out_genic, 'w') as f_genic:
        f_genic.writelines(loc_genic_lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True, help="BLAST output TSV")
    parser.add_argument("--out_ltr", required=True, help="Output LTR loc file")
    parser.add_argument("--out_genic", required=True, help="Output genic loc file")
    parser.add_argument("--context_range", type=int, default=20000, help="Context range")
    parser.add_argument("--ltr_range", type=int, default=500, help="LTR range")
    args = parser.parse_args()
    
    group_hits(args.blast, args.out_ltr, args.out_genic, args.context_range, args.ltr_range)
