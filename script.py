#!/usr/bin/env python3

import re
import csv
import sys

# Config
length_min = 35
loop_min = 3
loop_max = 10
bulge_max = 2
input_file = '../gene_lists/RNAsubopt_first_3utr.txt'
output_file = '../gene_lists/3utr_hairpin_trimmed_subopt.tsv'

log = open('parse_rnafold_hairpin.log', 'w') # for debugging
sys.stdout = log

def hairpin_check(candidate_hairpin, length_min, loop_min, bulge_max):
    middle_match = re.search((r'\((\.+)\)'), candidate_hairpin)
    if middle_match:
        middle_start = middle_match.start(1)  # refers only to the enclosed dots (group(1)), not the bracket
        middle_end = middle_match.end(1)
    if not middle_match or len(middle_match.group(1)) < loop_min\
            or len(middle_match.group(1)) > loop_max:
        print(f'Rejected!')
        return 'N', candidate_hairpin

    num_dots = (bulge_max + 1) * '\.'
    too_many_dots_upstream = re.compile(rf'\(({num_dots}+)\(')
    too_many_dots_downstream = re.compile(rf'\)({num_dots}+)\)')
    reversed_upstream = candidate_hairpin[:middle_start][::-1]  # reverse upstream for right-to-left search
    match_upstream = too_many_dots_upstream.search(reversed_upstream)
    match_downstream = too_many_dots_downstream.search(candidate_hairpin[middle_end:])
    if match_upstream:
        # Calculates original index in the non-reversed upstream part
        candidate_hairpin = candidate_hairpin[len(reversed_upstream) - match_upstream.start() - 1:]
    if match_downstream:
        middle_match = re.search((r'\((\.+)\)'), candidate_hairpin)  # recalculate, as positions changed
        if middle_match:
            middle_end = middle_match.end(1)
        candidate_hairpin = candidate_hairpin[:middle_end + match_downstream.start() + 1]

    open_count_pre = candidate_hairpin.count('(')
    close_count_pre = candidate_hairpin.count(')')
    diff = abs(open_count_pre - close_count_pre)

    if open_count_pre > close_count_pre:  # trimming to balance brackets
        list_hairpin = list(candidate_hairpin)  # converting to list to easily skip dots
        removed = 0
        for i in range(len(list_hairpin)):
            if list_hairpin[i] == '(':
                list_hairpin[i] = ''
                removed += 1
                if removed == diff:
                    break
        candidate_hairpin = ''.join(list_hairpin)  # converting back

    if open_count_pre < close_count_pre:
        list_hairpin = list(candidate_hairpin[::-1])  # same as above but in reverse
        removed = 0
        for i in range(len(list_hairpin)):
            if list_hairpin[i] == ')':
                list_hairpin[i] = ''
                removed += 1
                if removed == diff:
                    break
        list_hairpin.reverse()  # re-reverse
        candidate_hairpin = ''.join(list_hairpin)  # converting back

    candidate_hairpin = candidate_hairpin.strip('.')  # stripping possible trailing '.' after balancing

    print(f'Post-trimming: {candidate_hairpin}')

    hairpin_length = len(candidate_hairpin)
    if hairpin_length < length_min:
        print(f'Hairpin too short!')
        return 'N', candidate_hairpin
    if hairpin_length >= length_min:
        print(f'Hairpin found!')
        return 'Y', candidate_hairpin
    else:
        print(f'Rejected!')
        return 'N', candidate_hairpin

def parse_rnafold_hairpin(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    data = []
    gene_id = None
    rna = None

    for i, line in enumerate(lines):
        if line.startswith('>'):
            gene_id_match = re.search(r'gene_id=([^;]+);', line)
            if gene_id_match:
                gene_id = gene_id_match.group(1)
        elif line.startswith(('A', 'U', 'G', 'C')):
            rna = line.strip()
        elif line.startswith(('.', '(', ')')):
            # Extract MFE value
            mfe_match = re.search(r'\s+([-\d.]+)\s+', line)
            mfe = mfe_match.group(1) if mfe_match else None

            # Process structure
            temp_structure = line.strip()
            structure_parts = temp_structure.split()
            if structure_parts:
                trimmed_structure = structure_parts[0] + '('  # Add '(' at the end

            print(f"Processing: {gene_id}, {trimmed_structure}")  # Debug print

            if trimmed_structure and gene_id:
                start_position = None
                final_position = None
                candidate_hairpin = None
                result = 'N'
                for position, char in enumerate(trimmed_structure):
                    if char == '(' and start_position is None:
                        start_position = position
                    elif char == ')':
                        final_position = position
                    elif char == '(' and final_position is not None:
                        candidate_hairpin = trimmed_structure[start_position:final_position + 1]
                        print(f"Candidate hairpin: {candidate_hairpin}")
                        result, candidate_hairpin = hairpin_check(candidate_hairpin, length_min, loop_min, bulge_max)
                        print(f'Data appended.')
                        if result == 'Y':
                            # Hairpin position updated
                            pattern = re.escape(candidate_hairpin)
                            match = re.search(pattern, trimmed_structure[:-1])
                            if match:
                                exact_start = match.start()
                                exact_end = match.end()
                                rna_segment = rna[exact_start:exact_end] if rna else None
                                data.append([gene_id, result, trimmed_structure[:-1],
                                             candidate_hairpin, exact_start, exact_end,
                                             rna_segment, mfe])
                                break
                        else:
                            start_position = position
                            final_position = None
                if result == 'N':
                    data.append([gene_id, result, trimmed_structure[:-1], candidate_hairpin,  # Remove the added '(' when storing
                               start_position, final_position, rna, mfe])

    # Write output
    headers = ['GeneID', 'Hairpin', 'Trimmed_Str', 'Candidate_Hairpin_Str', 'StartPos', 'EndPos', 'RNA', 'MFE']

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
        writer.writerow(headers)
        writer.writerows(data)
        print(f"Written to TSV: {output_file}")

    return data

parse_rnafold_hairpin(input_file, output_file)
log.close()
```
