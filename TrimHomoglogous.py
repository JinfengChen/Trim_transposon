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
    pairs = defaultdict(lambda: defaultdict(lambda : int()))
    pairs_end = defaultdict(lambda: defaultdict(lambda : defaultdict(lambda : int())))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Query_id'):
                unit = re.split(r'\t',line)
                unit[2], unit[3] = sorted([unit[2], unit[3]], key=int)
                unit[6], unit[7] = sorted([unit[6], unit[7]], key=int)
                #if already compared this pairs before
                if pairs[unit[4]].has_key(unit[0]):
                    continue
                else:
                    pairs[unit[0]][unit[4]] = 1
                #if already have start or end compared
                #if pairs_end[unit[0]].has_key(unit[4]):
                #    if pairs_end[unit[0]][unit[4]][1] == 1:
                #        continue
                #    elif pairs_end[unit[0]][unit[4]][2] == 1:
                #        continue
                if unit[0] != unit[4] and unit[1] == unit[11] and unit[5] == unit[11]:
                    #duplicate
                    if G.has_edge(unit[0], unit[4]):
                         G[unit[0]][unit[4]]['weight'] += 1
                    else:
                         #G[unit[0]][unit[4]]['weight']  = 0.5
                         G.add_edge(unit[0], unit[4], weight=1)
                    
                elif unit[0] != unit[4] and (int(unit[2]) <= 3 and int(unit[6]) <= 3) and int(unit[11]) >= 100 and float(unit[8]) >= 0.99:
                    #start match start
                    #G.add_edge(unit[0], unit[4])
                    #if G.has_edge(unit[0], unit[4]):
                    #    G[unit[0]][unit[4]]['weight'] += 0.5
                    #else:
                    #    G.add_edge(unit[0], unit[4], weight=0.5)
                    #check if already have end to end compare
                    if 1:
                        if pairs_end[unit[0]][unit[4]][1] == 1 or pairs_end[unit[0]][unit[4]][3] == 1 or pairs_end[unit[0]][unit[4]][4] == 1:
                            continue
                        else:
                            if G.has_edge(unit[0], unit[4]):
                                G[unit[0]][unit[4]]['weight'] += 0.5
                            else:
                                G.add_edge(unit[0], unit[4], weight=0.5)    
                            pairs_end[unit[0]][unit[4]][1] = 1
                elif unit[0] != unit[4] and (int(unit[3]) >= int(unit[1]) - 3 and int(unit[7]) >= int(unit[5]) - 3) and int(unit[11]) >= 100 and float(unit[8]) >= 0.99:
                    #end match end
                    #if G.has_edge(unit[0], unit[4]):
                    #    G[unit[0]][unit[4]]['weight'] += 0.5
                    #else:
                    #    #G[unit[0]][unit[4]]['weight']  = 0.5
                    #    G.add_edge(unit[0], unit[4], weight=0.5)
                    #check if already have end to end compare
                    if 1:
                        if pairs_end[unit[0]][unit[4]][2] == 1 or pairs_end[unit[0]][unit[4]][3] == 1 or pairs_end[unit[0]][unit[4]][4] == 1:
                            continue
                        else:
                            if G.has_edge(unit[0], unit[4]):
                                G[unit[0]][unit[4]]['weight'] += 0.5
                            else:
                                G.add_edge(unit[0], unit[4], weight=0.5)    
                            pairs_end[unit[0]][unit[4]][2] = 1
                elif unit[0] != unit[4] and (int(unit[2]) <= 3 and int(unit[7]) >= int(unit[5]) - 3) and int(unit[11]) >= 100 and float(unit[8]) >= 0.99:
                    #start match end
                    #if G.has_edge(unit[0], unit[4]):
                    #    G[unit[0]][unit[4]]['weight'] += 0.5
                    #else:
                    #    #G[unit[0]][unit[4]]['weight']  = 0.5
                    #    G.add_edge(unit[0], unit[4], weight=0.5) 
                    #check if already have end to end compare
                    if 1:
                        if pairs_end[unit[0]][unit[4]][3] == 1 or pairs_end[unit[0]][unit[4]][1] == 1 or pairs_end[unit[0]][unit[4]][2] == 1:
                            continue
                        else:
                            if G.has_edge(unit[0], unit[4]):
                                G[unit[0]][unit[4]]['weight'] += 0.5
                            else:
                                G.add_edge(unit[0], unit[4], weight=0.5)    
                            pairs_end[unit[0]][unit[4]][3] = 1
                elif unit[0] != unit[4] and (int(unit[3]) >= int(unit[1]) - 3 and int(unit[6]) <= 3) and int(unit[11]) >= 100 and float(unit[8]) >= 0.99:
                    #end match start
                    #if G.has_edge(unit[0], unit[4]):
                    #    G[unit[0]][unit[4]]['weight'] += 0.5
                    #else:
                    #    #G[unit[0]][unit[4]]['weight']  = 0.5
                    #    G.add_edge(unit[0], unit[4], weight=0.5)
                    #check if already have end to end compare
                    if 1:
                        if pairs_end[unit[0]][unit[4]][4] == 1 or pairs_end[unit[0]][unit[4]][1] == 1 or pairs_end[unit[0]][unit[4]][2] == 1:
                            continue
                        else:
                            if G.has_edge(unit[0], unit[4]):
                                G[unit[0]][unit[4]]['weight'] += 0.5
                            else:
                                G.add_edge(unit[0], unit[4], weight=0.5)                   
                            pairs_end[unit[0]][unit[4]][4] = 1
                #elif unit[0] != unit[4] and unit[1] == unit[11] and unit[5] == unit[11]:
                #    #duplicate
                #    if G.has_edge(unit[0], unit[4]):
                #         G[unit[0]][unit[4]]['weight'] += 1
                #    else:
                #         #G[unit[0]][unit[4]]['weight']  = 0.5
                #         G.add_edge(unit[0], unit[4], weight=1)
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

    for u,v in G.edges():
        print u, v, G[u][v]['weight']
        if float(G[u][v]['weight']) < 1:
            G.remove_edge(u,v)

    #print G.nodes()
    #print 'number of node: %s' %(len(G.nodes()))
    #print 'number of edges: %s' %(len(G.edges()))   
    #print nx.connected_components(G)
   
    fasta_record = fasta_id(args.fasta)
    ofile = open ('%s.representative.fa' %(os.path.splitext(args.fasta)[0]), 'w') 
    for repeats in nx.connected_components(G):
        if len(repeats) > 1:
            top  = ''
            top1 = 'N'*1000000 #shortest name length
            top2 = top1    #shortest seq length
            for rep in sorted(repeats):
                top1 = rep if len(rep) < len(top1) else top1
                if fasta_record.has_key(top2):
                    top2 = rep if len(str(fasta_record[rep].seq)) < len(str(fasta_record[top2].seq)) else top2
                else:
                    top2 = rep
            
            if top1 == top2:
                top = top1
            else:
                if len(str(fasta_record[top1].seq)) > len(str(fasta_record[top2].seq)):
                    top = top2
                else:
                    top = top1
            SeqIO.write(fasta_record[top], ofile, 'fasta')
            print '%s\tRepresentative\t%s' %(top, ','.join(repeats))
        else:
            SeqIO.write(fasta_record[repeats[0]], ofile, 'fasta')
    ofile.close()

if __name__ == '__main__':
    main()

