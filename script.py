#!/usr/bin/env python3


import re
import csv
import sys

# Config
paren_min = 30
middle_min = 3
bulge_max = 2
INPUT_FILE = 'your_input.tsv'
OUTPUT_TYPE = 'tsv'
OUTPUT_FILE = 'your_output.tsv'

# logging
log = open('parse_rnafold_hairpin.log', 'w')
sys.stdout = log


# identifies hairpins, trims accordingly, configurable
def hairpin_check(candidate_hairpin, paren_min, middle_min, bulge_max):
    middle_match = re.search((r'\((\.+)\)'), candidate_hairpin)
    if middle_match:
        middle_start = middle_match.start(1)  # refers only to the enclosed dots (group(1)), not the bracket
        middle_end = middle_match.end(1)
    if not middle_match or len(middle_match.group(1)) < middle_min:
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

    open_count_post = candidate_hairpin.count('(')
    close_count_post = candidate_hairpin.count(')')
    if open_count_post + close_count_post < paren_min:
        print(f'Hairpin too short!')
        return 'N', candidate_hairpin
    if open_count_post + close_count_post >= paren_min and open_count_post == close_count_post:
        print(f'Hairpin found!')
        return 'Y', candidate_hairpin
    else:
        print(f'Rejected!')
        return 'N', candidate_hairpin


def parse_rnafold_hairpin(input_file, output_type, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    data = []
    gene_id = None

    for line in lines:
        if line.startswith('>'):
            gene_id_match = re.search(r'geneID=([^;]+);', line)
            if gene_id_match:
                gene_id = gene_id_match.group(1)
        elif line.startswith(('A', 'U', 'G', 'C')):
            rna = line.strip()[:126]
        elif line.startswith(('.', '(', ')')) and gene_id:
            temp_structure = line.strip()[:126]  # trimming the first 125
            last_paren = temp_structure.rfind('(')  # if not trimmed, remove free energy calc from line's end
            temp_structure = temp_structure[:last_paren].strip()
            trimmed_structure = temp_structure + '('  # append a (, useful for the candidate_hairpin logic

            print(f"Processing: {gene_id}, {trimmed_structure}")  # Debug print

            if trimmed_structure and gene_id:
                start_position = None  # init
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
                        result, candidate_hairpin = hairpin_check(candidate_hairpin, paren_min, middle_min, bulge_max)
                        print(f'Data appended.')
                        if result == 'Y':
                            rna = rna[start_position:final_position + 1]
                            data.append([gene_id, result, trimmed_structure, candidate_hairpin,
                                         start_position, final_position, rna])
                            break
                        else:
                            start_position = position  # if set to 'None', we would 'lose' one pos in the enumeration
                            final_position = None
                if result == 'N':
                    data.append([gene_id, result, trimmed_structure, candidate_hairpin,
                                 start_position, final_position, rna])
    if output_type.lower() == 'csv':
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['GeneID', 'Hairpin', 'Trimmed Str', 'Candidate Hairpin Str', 'StartPos', 'EndPos', 'RNA'])
            writer.writerows(data)
            print(f"Written to CSV: {output_file}")
    if output_type.lower() == 'tsv':
        with open(output_file, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(['GeneID', 'Hairpin', 'Trimmed Str', 'Candidate Hairpin Str', 'StartPos', 'EndPos', 'RNA'])
            writer.writerows(data)
            print(f"Written to TSV: {output_file}")


parse_rnafold_hairpin(INPUT_FILE, OUTPUT_TYPE, OUTPUT_FILE)
log.close()
