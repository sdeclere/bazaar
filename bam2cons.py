#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" fasta consensus from sorted indexed bam for a particular taxid 
    usage : arg1 -> bam 
            arg2 -> taxid 
"""

import sys
import pysam
import re 

from collections import Counter, OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ac_threshold = 3
af_threshold = 0.9

def do_consensus(bam, consensus):
    allele_counter = Counter()
    for pileup_column in bam.pileup(stepper='nofilter',max_depth=500000,truncate=False,min_base_quality=0):
        chrom = pileup_column.reference_name
        pos = pileup_column.reference_pos
        
        # assert current pos is uncalled 
        assert consensus[chrom][pos] == "N"
        # reinit counter
        allele_counter.clear()

        for pileup_read in pileup_column.pileups:
            if pileup_read.is_del:
                allele = "-"
            else:
                allele = pileup_read.alignment.query_sequence[pileup_read.query_position]
            allele_counter[allele] += 1

            max_allele = "N"
            max_count, total_count = 0, 0
            for allele, count in allele_counter.items():
                if count > max_count:
                    max_count = count
                    max_allele = allele
                total_count += count

            assert max_allele in "ACGTN-"
            if (max_count >= ac_threshold and
                max_count / total_count >= af_threshold):
                consensus[chrom][pos] = max_allele
    return consensus
        

def main(args):
    consensus = OrderedDict()
    for record in SeqIO.parse(sys.argv[2], "fasta"):
        consensus[record.id] = ["N"] * len(record)
    print('%d refs loaded from codex' % len(consensus.keys()))

    bam_file = args[1]
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        consensus =  do_consensus(bam, consensus)
    
    records = []
    for chrom, seq in consensus.items() :
        printed_seq = "".join(seq)   
        records.append( SeqRecord(Seq(printed_seq), id=bam_file, description="from : " + chrom) ) 

    out = bam_file.strip().split('.')[0]
    SeqIO.write(records, "%s_cons.fasta" % out, "fasta")

if __name__ == '__main__':
    main(sys.argv)
