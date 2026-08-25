#!/usr/bin/env python3
import argparse
import re

def delimit_sequences(input_fasta, output_fasta, ltr_range):
    with open(input_fasta, 'r') as f:
        lines = f.readlines()

    PATTERNS = [
        (re.compile(r"[ :\.-]"), "_"),              # Substitui caracteres especiais
        (re.compile(r",.+$|isolate_"), ""),         # Remove texto após vírgula, número de acesso e isolate
        (re.compile(r"chromosome_"), "c"),          # Abrevia "chromosome"
    ]

    ervnovo = []
    for line in lines:
        if line.startswith('>'):
            for encontrar, substituir in PATTERNS:
                line = encontrar.sub(substituir, line)
        ervnovo.append(line)

    count = 0
    with open(output_fasta, 'w') as out:
        header = ''
        for line in ervnovo:
            if line.startswith('>'):
                header = line.rstrip()
            else:
                seqatual = line.rstrip() # Note: line originally might have \n, keep rstrip to handle exact length
                if len(seqatual) > ltr_range * 2:
                    # Original minERVa used hardcoded indices [496:500] and [494:500] assuming LTR_range=500
                    # Let's adapt them to use ltr_range dynamically
                    sig1_start = ltr_range - 4
                    sig1_end = ltr_range
                    sig1_rev_start = -(ltr_range + 1)
                    sig1_rev_end = -(ltr_range - 3)
                    
                    sig2_start = ltr_range - 6
                    sig2_end = ltr_range
                    sig2_rev_start = -(ltr_range + 1)
                    sig2_rev_end = -(ltr_range - 5)
                    
                    if seqatual[sig1_start:sig1_end] == seqatual[sig1_rev_start:sig1_rev_end] and seqatual[sig1_start:sig1_end] != '':
                        count += 1
                        seqtotal = f'{header}_{seqatual[sig1_start:sig1_end]}\n{seqatual[sig1_start:-(ltr_range - 3)]}\n'
                    elif seqatual[sig2_start:sig2_end] == seqatual[sig2_rev_start:sig2_rev_end] and seqatual[sig2_start:sig2_end] != '':
                        count += 1
                        seqtotal = f'{header}_{seqatual[sig2_start:sig2_end]}\n{seqatual[sig2_start:-(ltr_range - 5)]}\n'
                    else:
                        seqtotal = f'{header}_NS\n{seqatual}\n'
                else:
                    seqtotal = f'{header}_NS\n{seqatual}\n'
                    
                out.write(seqtotal)

    print(f'Signatures: {count} found in sequences.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output trimmed FASTA")
    parser.add_argument("--ltr_range", type=int, default=500, help="LTR range margin")
    
    args = parser.parse_args()
    delimit_sequences(args.input, args.output, args.ltr_range)
