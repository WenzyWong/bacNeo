#!/usr/bin/env python
#################################################################################
# Modified from original source: https://github.com/jenniferlu717/KrakenTools/extract_kraken_reads.py
# Main modifications:
## 1. Remove the -t/--taxon-id parameter 
## 2. Collect all unique taxon IDs before processing kraken files
## 3. Use collected taxon IDs to replace the IDs originally entered with the-t parameter
## 4. Remove the unused -o2 parameter when parsing paired-end fastq files

# Modifications made by Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2024-12-18
#################################################################################
import os, sys, argparse
import gzip
from time import gmtime
from time import strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#################################################################################
#Tree Class 
class Tree(object):
    'Tree node.'
    def __init__(self, taxid, level_num, level_id, children=None, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)

def process_kraken_output(kraken_line):
    l_vals = kraken_line.split('\t')
    if len(l_vals) < 5:
        return [-1, '']
    if "taxid" in l_vals[2]:
        temp = l_vals[2].split("taxid ")[-1]
        tax_id = temp[:-1]
    else:
        tax_id = l_vals[2]

    read_id = l_vals[1]
    if (tax_id == 'A'):
        tax_id = 81077
    else:
        tax_id = int(tax_id)
    return [tax_id, read_id]

def process_kraken_report(report_line):
    l_vals = report_line.strip().split('\t')
    if len(l_vals) < 5:
        return []
    try:
        int(l_vals[1])
    except ValueError:
        return []
    try:
        taxid = int(l_vals[-3]) 
        level_type = l_vals[-2]
        map_kuniq = {'species':'S', 'genus':'G','family':'F',
            'order':'O','class':'C','phylum':'P','superkingdom':'D',
            'kingdom':'K'}
        if level_type not in map_kuniq:
            level_type = '-'
        else:
            level_type = map_kuniq[level_type]
    except ValueError:
        taxid = int(l_vals[-2])
        level_type = l_vals[-3]
    spaces = 0
    for char in l_vals[-1]:
        if char == ' ':
            spaces += 1
        else:
            break
    level_num = int(spaces/2)
    return[taxid, level_num, level_type]

def collect_taxids(kraken_file):
    """Collect all unique taxonomy IDs from kraken output file"""
    taxids = {}
    with open(kraken_file, 'r') as k_file:
        for line in k_file:
            tax_id, _ = process_kraken_output(line)
            if tax_id != -1:
                if tax_id not in taxids:
                    taxids[tax_id] = 0
    return taxids

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', dest='kraken_file', required=True,
        help='Kraken output file to parse')
    parser.add_argument('-s','-s1', '-1', '-U', dest='seq_file1', required=True,
        help='FASTA/FASTQ File containing the raw sequence letters.')
    parser.add_argument('-s2', '-2', dest='seq_file2', default="",
        help='2nd FASTA/FASTQ File containing the raw sequence letters (paired).')
    parser.add_argument('-o', "--output",dest='output_file', required=True,
        help='Output FASTA/Q file containing the reads and sample IDs')
    parser.add_argument('-o2',"--output2", dest='output_file2', required=False, default='',
        help='Output FASTA/Q file containig the second pair of reads [optional]') 
    parser.add_argument('--append', dest='append', action='store_true',
        help='Append the sequences to the end of the output FASTA file specified.')
    parser.add_argument('--noappend', dest='append', action='store_false',
        help='Create a new FASTA file containing sample sequences and IDs (rewrite if existing) [default].')
    parser.add_argument('--max', dest='max_reads', required=False, 
        default=100000000, type=int,
        help='Maximum number of reads to save [default: 100,000,000]')
    parser.add_argument('-r','--report',dest='report_file', required=False,
        default="",
        help='Kraken report file. [required only if --include-parents/children is specified]')
    parser.add_argument('--include-parents',dest="parents", required=False, 
        action='store_true',default=False,
        help='Include reads classified at parent levels of the specified taxids')
    parser.add_argument('--include-children',dest='children', required=False,
        action='store_true',default=False,
        help='Include reads classified more specifically than the specified taxids')
    parser.add_argument('--exclude', dest='exclude', required=False,
        action='store_true',default=False,
        help='Instead of finding reads matching specified taxids, finds all reads NOT matching specified taxids') 
    parser.add_argument('--fastq-output', dest='fastq_out', required=False,
        action='store_true',default=False,
        help='Print output FASTQ reads [requires input FASTQ, default: output is FASTA]')
    parser.set_defaults(append=False)

    args=parser.parse_args()
    
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')
    
    # Removed the check for requiring -o2 with paired-end input

    # Collect unique taxids from kraken output
    sys.stdout.write(">> STEP 0: COLLECTING TAXONOMY IDs FROM KRAKEN OUTPUT\n")
    save_taxids = collect_taxids(args.kraken_file)
    sys.stdout.write(f"\t{len(save_taxids)} unique taxonomy IDs found\n")
    
    main_lvls = ['R','K','D','P','C','O','F','G','S']

    if args.parents or args.children:
        if args.report_file == "": 
            sys.stderr.write(">> ERROR: --report not specified.")
            sys.exit(1)
        sys.stdout.write(">> STEP 1: PARSING REPORT FILE %s\n" % args.report_file)
        base_nodes = {} 
        r_file = open(args.report_file,'r')
        prev_node = -1
        for line in r_file:
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            [taxid, level_num, level_id] = report_vals
            if taxid == 0:
                continue 
            if taxid == 1:
                level_id = 'R'
                root_node = Tree(taxid, level_num, level_id)
                prev_node = root_node
                if taxid in save_taxids:
                    base_nodes[taxid] = root_node
                continue
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent 
            if level_id == '-' or len(level_id) > 1:
                if prev_node.level_id in main_lvls:
                    level_id = prev_node.level_id + '1'
                else:
                    num = int(prev_node.level_id[-1]) + 1
                    level_id = prev_node.level_id[:-1] + str(num)
            curr_node = Tree(taxid, level_num, level_id, None, prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node
            if taxid in save_taxids:
                base_nodes[taxid] = curr_node 
        r_file.close()
        
        if args.parents:
            for tid in base_nodes:
                curr_node = base_nodes[tid]
                while curr_node.parent != None:
                    curr_node = curr_node.parent
                    save_taxids[curr_node.taxid] = 0
                    
        if args.children:
            for tid in base_nodes:
                curr_nodes = base_nodes[tid].children
                while len(curr_nodes) > 0:
                    curr_n = curr_nodes.pop()
                    if curr_n.taxid not in save_taxids:
                        save_taxids[curr_n.taxid] = 0
                    if curr_n.children != None:
                        for child in curr_n.children:
                            curr_nodes.append(child)
                    
    sys.stdout.write(">> STEP 2: PARSING KRAKEN FILE FOR READIDS %s\n" % args.kraken_file)
    count_kraken = 0
    read_line = -1
    exclude_taxids = {} 
    if args.exclude:
        exclude_taxids = save_taxids 
        save_taxids = {} 

    k_file = open(args.kraken_file, 'r')
    sys.stdout.write('\t0 reads processed')
    sys.stdout.flush()
    save_readids = {}
    save_readids2 = {} 
    for line in k_file:
        count_kraken += 1
        if (count_kraken % 10000 == 0):
            sys.stdout.write('\r\t%0.2f million reads processed' % float(count_kraken/1000000.))
            sys.stdout.flush()
        [tax_id, read_id] = process_kraken_output(line)
        if tax_id == -1:
            continue
        if (tax_id in save_taxids) and not args.exclude:
            save_taxids[tax_id] += 1
            save_readids2[read_id] = 0
            save_readids[read_id] = 0 
        elif (tax_id not in exclude_taxids) and args.exclude:
            if tax_id not in save_taxids:
                save_taxids[tax_id] = 1
            else:
                save_taxids[tax_id] += 1
            save_readids2[read_id] = 0
            save_readids[read_id] = 0 
        if len(save_readids) >= args.max_reads:
            break 

    k_file.close()
    sys.stdout.write('\r\t%0.2f million reads processed\n' % float(count_kraken/1000000.))
    sys.stdout.write('\t%i read IDs saved\n' % len(save_readids))

    seq_file1 = args.seq_file1
    seq_file2 = args.seq_file2
    if(seq_file1[-3:] == '.gz'):
        s_file1 = gzip.open(seq_file1,'rt')
    else:
        s_file1 = open(seq_file1,'rt')
    first = s_file1.readline()
    if len(first) == 0:
        sys.stderr.write("ERROR: sequence file's first line is blank\n")
        sys.exit(1)
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        sys.exit(1)
    s_file1.close()
    if filetype != 'fastq' and args.fastq_out:
        sys.stderr.write('ERROR: for FASTQ output, input file must be FASTQ\n')
        sys.exit(1)

    if(seq_file1[-3:] == '.gz'):
        s_file1 = gzip.open(seq_file1,'rt')
        if len(seq_file2) > 0:
            s_file2 = gzip.open(seq_file2,'rt')
    else:
        s_file1 = open(seq_file1, 'r')
        if len(seq_file2) > 0:
            s_file2 = open(seq_file2, 'r')

    sys.stdout.write(">> STEP 3: READING SEQUENCE FILES AND WRITING READS\n")
    sys.stdout.write('\t0 read IDs found (0 mill reads processed)')
    sys.stdout.flush()

    if (args.append):
        o_file = open(args.output_file, 'a')
        if args.output_file2 != '':
            o_file2 = open(args.output_file2, 'a')
    else:
        o_file = open(args.output_file, 'w')
        if args.output_file2 != '':
            o_file2 = open(args.output_file2, 'w')

    count_seqs = 0
    count_output = 0
    for record in SeqIO.parse(s_file1,filetype):
        count_seqs += 1
        if (count_seqs % 1000 == 0):
            sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
            sys.stdout.flush()
        test_id = str(record.id)
        test_id2 = test_id
        if ("/1" in test_id) or ("/2" in test_id):
            test_id2 = test_id[:-2]
        if test_id in save_readids or test_id2 in save_readids:
            count_output += 1
            sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
            sys.stdout.flush()
            if args.fastq_out:
                SeqIO.write(record, o_file, "fastq")
            else:
                SeqIO.write(record, o_file, "fasta")
        if len(save_readids) == count_output:
            break
    s_file1.close()
    o_file.close()
    sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)\n' % (count_output, float(count_seqs/1000000.)))
    sys.stdout.flush()
    
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    sys.exit(0)

if __name__ == "__main__":
    main()

#################################END OF PROGRAM##################################
#################################################################################