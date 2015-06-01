#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = os.path.splitext(args.input)[0] + '.short.fa'

    ofile = open(args.output, "w")
    for record in SeqIO.parse(args.input, "fasta"):
        print record.id
        newseq    = Seq(str(record.seq))
        if len(str(record.seq)) > 1600:
            newseq    = Seq(str(record.seq)[:800] + 'N' + str(record.seq)[-800:])
        newrecord = SeqRecord(newseq, id=record.id, description="")
        SeqIO.write(newrecord, ofile, "fasta")
    ofile.close()

 

if __name__ == '__main__':
    main()

