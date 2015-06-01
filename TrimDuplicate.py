#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import networkx as nx

def usage():
    test="name"
    message='''
python TrimDuplicate.py --input test.table --fasta Rice.TE.short.fa

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = record
    return fastaid

#mping   430     1       430     OS_mPIFHarbinger_11     430     430     1       1       --      0       430     852     0.0     --      --
#mping   430     1       430     mping   430     1       430     1       --      0       430     852     0.0     --      --
#mping   430     1       253     ping    1601    1       253     0.99    --      0       253     494     1e-140  --      --
#mping   430     253     430     ping    1601    1424    1601    1       --      0       178     353     7e-98   --      --
#mping   430     1       47      pong    1601    1       47      0.89    --      0       47      54.0    9e-08   --      --
def readtable(infile, G):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Query_id'):
                unit = re.split(r'\t',line)
                if unit[0] != unit[4] and unit[1] == unit[11] and unit[5] == unit[11]:
                    G.add_edge(unit[0], unit[4])
                elif unit[0] != unit[4]:
                    G.add_nodes_from([unit[0], unit[4]])
    return G


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    #try:
    #    len(args.input) > 0
    #except:
    #    usage()
    #    sys.exit(2)

    G=nx.Graph()
    #G.add_node("spam")
    #G.add_edge(1,2)
    #G.add_edge(1,3)
    #print G.nodes()
    #print G.edges()
    #print nx.connected_components(G)
   
    readtable(args.input, G)
    #print G.nodes()
    #print 'number of node: %s' %(len(G.nodes()))
    #print 'number of edges: %s' %(len(G.edges()))   
    #print nx.connected_components(G)
    
    fasta_record = fasta_id(args.fasta)
    ofile = open ('%s.unique.fa' %(os.path.splitext(args.fasta)[0]), 'w') 
    for repeats in nx.connected_components(G):
        if len(repeats) > 1:
            top = 'N'*100
            for rep in repeats:
                top = rep if len(rep) < len(top) else top
            SeqIO.write(fasta_record[top], ofile, 'fasta')
            print '%s\tDuplicate\t%s' %(top, ','.join(repeats))
        else:
            SeqIO.write(fasta_record[repeats[0]], ofile, 'fasta')
    ofile.close()

if __name__ == '__main__':
    main()

