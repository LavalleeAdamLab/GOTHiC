import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import csv
import sys
from graph_tool.all import *
import concurrent.futures
import timeit
from tqdm import tqdm
import random
import math
import plotly.graph_objects as go
import time
from goatools import obo_parser
import subprocess
from scipy import stats as sps
import datetime
import gnuplotlib as gp
import os
import importlib
import seaborn as sns
import gcMapExplorer.lib as gmlib
import fanc


if sys.version_info[0] < 3:
    import StringIO
else:
    from io import StringIO


# utility function to quickly get line count of file for tracking progress through file iteration
def count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


# Converts Hi-C SAM files to adjacency list with nodes representing bins of desired size
# to use assemblies other than Hg38.13, define a custom vector of size 24 where entries are chrom lengths in order below
# to use custom bins, define a dict with keys chr1, chr2, ..., and values are a list of bin end positions
def sam_adjacencies(sam1, sam2, m=40000, bins=None, chromlengths=None, verbose=False):  # samt stands for "SAM Table", from the output of samTable function

    # Chromosome lengths based on Hg38.13
    if chromlengths is None:
        cr1 = 248956422
        cr2 = 242193529
        cr3 = 198295559
        cr4 = 190214555
        cr5 = 181538259
        cr6 = 170805979
        cr7 = 159345973
        cr8 = 145138636
        cr9 = 138394717
        cr10 = 133797422
        cr11 = 135086622
        cr12 = 133275309
        cr13 = 114364328
        cr14 = 107043718
        cr15 = 101991189
        cr16 = 90338345
        cr17 = 83257441
        cr18 = 80373285
        cr19 = 58617616
        cr20 = 64444167
        cr21 = 46709983
        cr22 = 50818468
        crX = 156040895
        crY = 57227415

    else:
        cr1 = chromlengths[0]
        cr2 = chromlengths[1]
        cr3 = chromlengths[2]
        cr4 = chromlengths[3]
        cr5 = chromlengths[4]
        cr6 = chromlengths[5]
        cr7 = chromlengths[6]
        cr8 = chromlengths[7]
        cr9 = chromlengths[8]
        cr10 = chromlengths[9]
        cr11 = chromlengths[10]
        cr12 = chromlengths[11]
        cr13 = chromlengths[12]
        cr14 = chromlengths[13]
        cr15 = chromlengths[14]
        cr16 = chromlengths[15]
        cr17 = chromlengths[16]
        cr18 = chromlengths[17]
        cr19 = chromlengths[18]
        cr20 = chromlengths[19]
        cr21 = chromlengths[20]
        cr22 = chromlengths[21]
        crX = chromlengths[22]
        crY = chromlengths[23]

    # chromosome list for checking whether RNAME value is valid or not
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    # initialize empty lists of lists of lists
    # list1: chromosomes (length 24: 22 + 23(X) +24(Y)), list2: bins(length = length(chr)/m),list3:weighted adjacencies
    # list sizes assume Hg38.13
    if bins is None:
        adjlist = [None] * 24
        adjlist[0] = [None] * int((cr1 / m) + 1)
        adjlist[1] = [None] * int((cr2 / m) + 1)
        adjlist[2] = [None] * int((cr3 / m) + 1)
        adjlist[3] = [None] * int((cr4 / m) + 1)
        adjlist[4] = [None] * int((cr5 / m) + 1)
        adjlist[5] = [None] * int((cr6 / m) + 1)
        adjlist[6] = [None] * int((cr7 / m) + 1)
        adjlist[7] = [None] * int((cr8 / m) + 1)
        adjlist[8] = [None] * int((cr9 / m) + 1)
        adjlist[9] = [None] * int((cr10 / m) + 1)
        adjlist[10] = [None] * int((cr11 / m) + 1)
        adjlist[11] = [None] * int((cr12 / m) + 1)
        adjlist[12] = [None] * int((cr13 / m) + 1)
        adjlist[13] = [None] * int((cr14 / m) + 1)
        adjlist[14] = [None] * int((cr15 / m) + 1)
        adjlist[15] = [None] * int((cr16 / m) + 1)
        adjlist[16] = [None] * int((cr17 / m) + 1)
        adjlist[17] = [None] * int((cr18 / m) + 1)
        adjlist[18] = [None] * int((cr19 / m) + 1)
        adjlist[19] = [None] * int((cr20 / m) + 1)
        adjlist[20] = [None] * int((cr21 / m) + 1)
        adjlist[21] = [None] * int((cr22 / m) + 1)
        adjlist[22] = [None] * int((crX / m) + 1)
        adjlist[23] = [None] * int((crY / m) + 1)

    else:
        adjlist = [None] * 24
        adjlist[0] = [None] * len(bins['chr1'])
        adjlist[1] = [None] * len(bins['chr2'])
        adjlist[2] = [None] * len(bins['chr3'])
        adjlist[3] = [None] * len(bins['chr4'])
        adjlist[4] = [None] * len(bins['chr5'])
        adjlist[5] = [None] * len(bins['chr6'])
        adjlist[6] = [None] * len(bins['chr7'])
        adjlist[7] = [None] * len(bins['chr8'])
        adjlist[8] = [None] * len(bins['chr9'])
        adjlist[9] = [None] * len(bins['chr10'])
        adjlist[10] = [None] * len(bins['chr11'])
        adjlist[11] = [None] * len(bins['chr12'])
        adjlist[12] = [None] * len(bins['chr13'])
        adjlist[13] = [None] * len(bins['chr14'])
        adjlist[14] = [None] * len(bins['chr15'])
        adjlist[15] = [None] * len(bins['chr16'])
        adjlist[16] = [None] * len(bins['chr17'])
        adjlist[17] = [None] * len(bins['chr18'])
        adjlist[18] = [None] * len(bins['chr19'])
        adjlist[19] = [None] * len(bins['chr20'])
        adjlist[20] = [None] * len(bins['chr21'])
        adjlist[21] = [None] * len(bins['chr22'])
        adjlist[22] = [None] * len(bins['chrX'])
        adjlist[23] = [None] * len(bins['chrY'])

    if verbose:
        # quickly get file line count for progress bar
        with open(sam1, 'rb') as fp:
            c_generator = count_generator(fp.raw.read)
            # count each \n
            count = sum(buffer.count(b'\n') for buffer in c_generator)
            print('SAM file line count:', count + 1)

    # iterate through SAM rows, and extract RNAME (chromosome) and POS, which is converted to a bin assignment
    with open(sam1) as file1, open(sam2) as file2:

        if verbose:
            pbar = tqdm(total=count+1)  # if verbose option set to True, create a progress bar that tracks parsed lines

        for line1, line2 in zip(file1, file2):
            try:
                if not line1.startswith('@'):  # ignores header lines of SAM file
                    # split input line by whitespace
                    s1 = line1.split()
                    s2 = line2.split()

                    # get info for first contact
                    chrom1 = s1[2]  # extracts chromosome name from 'RNAME' column of SAM
                    pos1 = int(s1[3])  # extracts leftmost mapping position from 'POS' column of SA
                    if chrom1 == '*':  # * means RNAME is unknown, so skip these lines
                        continue

                    bin1 = -1  # initialize bin1 as -1, then reassign according to binning logic
                    if bins is None:  # by default, bin is found by dividing the position by bin size
                        bin1 = int(pos1 / m)
                    else:  # otherwise loop through the bins to find which one it belongs in (by bin index)
                        bincounter = 0
                        for whichbin in bins[chrom1]:
                            binstart, binend = whichbin.split("-")
                            binstart = int(binstart)
                            binend = int(binend)
                            if binstart <= pos1 <= binend:
                                bin1 = bincounter
                                break  # exit loop when correct bin is found
                            bincounter = bincounter + 1
                        if bin1 == -1:  # if bin1 is not reassigned to a bin, exit this loop and go to next line
                            continue

                    # get info for second contact
                    chrom2 = s2[2]  # extracts chromosome name from 'RNAME' column of SAM
                    pos2 = int(s2[3])  # extracts leftmost mapping position from 'POS' column of SAM
                    if chrom2 == '*':  # * means RNAME is unknown, so skip these lines
                        continue

                    bin2 = -1
                    if bins is None:  # by default, bin is found by dividing the position by bin size
                        bin2 = int(pos2 / m)
                    else:  # otherwise loop through the bins to find which one it belongs in (by bin index)
                        bincounter = 0
                        for whichbin in bins[chrom2]:
                            binstart, binend = whichbin.split("-")
                            binstart = int(binstart)
                            binend = int(binend)
                            if binstart <= pos2 <= binend:
                                bin2 = bincounter
                                break  # exit loop when correct bin is found
                            bincounter = bincounter + 1
                        if bin2 == -1:  # if bin2 is not reassigned to a bin, exit this loop and go to next line
                            continue

                    if not (bin1 == bin2 and chrom1 == chrom2):  # filter out reads assigned to same bin
                        if any(item == chrom1 for item in chrlist) and any(item == chrom2 for item in chrlist):
                            # elif block that appends new connection/updates weight for appropriate bins
                            if chrom1 == 'chr1':
                                adjlist[0][bin1] = check_contact(adjlist[0][bin1], chrom2, bin2)
                            elif chrom1 == 'chr2':
                                adjlist[1][bin1] = check_contact(adjlist[1][bin1], chrom2, bin2)
                            elif chrom1 == 'chr3':
                                adjlist[2][bin1] = check_contact(adjlist[2][bin1], chrom2, bin2)
                            elif chrom1 == 'chr4':
                                adjlist[3][bin1] = check_contact(adjlist[3][bin1], chrom2, bin2)
                            elif chrom1 == 'chr5':
                                adjlist[4][bin1] = check_contact(adjlist[4][bin1], chrom2, bin2)
                            elif chrom1 == 'chr6':
                                adjlist[5][bin1] = check_contact(adjlist[5][bin1], chrom2, bin2)
                            elif chrom1 == 'chr7':
                                adjlist[6][bin1] = check_contact(adjlist[6][bin1], chrom2, bin2)
                            elif chrom1 == 'chr8':
                                adjlist[7][bin1] = check_contact(adjlist[7][bin1], chrom2, bin2)
                            elif chrom1 == 'chr9':
                                adjlist[8][bin1] = check_contact(adjlist[8][bin1], chrom2, bin2)
                            elif chrom1 == 'chr10':
                                adjlist[9][bin1] = check_contact(adjlist[9][bin1], chrom2, bin2)
                            elif chrom1 == 'chr11':
                                adjlist[10][bin1] = check_contact(adjlist[10][bin1], chrom2, bin2)
                            elif chrom1 == 'chr12':
                                adjlist[11][bin1] = check_contact(adjlist[11][bin1], chrom2, bin2)
                            elif chrom1 == 'chr13':
                                adjlist[12][bin1] = check_contact(adjlist[12][bin1], chrom2, bin2)
                            elif chrom1 == 'chr14':
                                adjlist[13][bin1] = check_contact(adjlist[13][bin1], chrom2, bin2)
                            elif chrom1 == 'chr15':
                                adjlist[14][bin1] = check_contact(adjlist[14][bin1], chrom2, bin2)
                            elif chrom1 == 'chr16':
                                adjlist[15][bin1] = check_contact(adjlist[15][bin1], chrom2, bin2)
                            elif chrom1 == 'chr17':
                                adjlist[16][bin1] = check_contact(adjlist[16][bin1], chrom2, bin2)
                            elif chrom1 == 'chr18':
                                adjlist[17][bin1] = check_contact(adjlist[17][bin1], chrom2, bin2)
                            elif chrom1 == 'chr19':
                                adjlist[18][bin1] = check_contact(adjlist[18][bin1], chrom2, bin2)
                            elif chrom1 == 'chr20':
                                adjlist[19][bin1] = check_contact(adjlist[19][bin1], chrom2, bin2)
                            elif chrom1 == 'chr21':
                                adjlist[20][bin1] = check_contact(adjlist[20][bin1], chrom2, bin2)
                            elif chrom1 == 'chr22':
                                adjlist[21][bin1] = check_contact(adjlist[21][bin1], chrom2, bin2)
                            elif chrom1 == 'chrX':
                                adjlist[22][bin1] = check_contact(adjlist[22][bin1], chrom2, bin2)
                            elif chrom1 == 'chrY':
                                adjlist[23][bin1] = check_contact(adjlist[23][bin1], chrom2, bin2)
                            else:
                                pass

            except KeyError:  # if mapping is for non-canonical chromosome, skip it
                continue

            if verbose:
                pbar.update(1)  # update pbar (if there is one) after each line

    if verbose:
        pbar.close()  # close pbar if one was created

    return adjlist


# converts adjacency list from sam_adjacencies() to graphml file
def hic_adjlist_to_graphml(adjlist, fileprefix="HicGraph"):
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    mlname = ""
    mlname = mlname.join([fileprefix, ".graphml"])

    g = Graph(directed=False)

    # create adjacency list to then import into graph with add_edgelist()
    contactlist = []  # for storing contacts
    nodelist = []  # for storing  unique node names, to later be  enumerated
    # initialize graph properties
    vname = g.new_vertex_property("string")  # vp map for node names (chromosome-bin pairs)
    weight = g.new_edge_property("double")  # ep map for edge weight

    for i in range(len(adjlist)):  # for each chromosome
        for j in range(len(adjlist[i])):  # for each bin
            # coerce chromosome-index pair to create node string (node u)
            nodeu = ""
            nodeu = nodeu.join([chrlist[i], ":", str(j)])  # get name of first node

            if adjlist[i][j]:
                for k in range(len(adjlist[i][j])):
                    # coerce chromosome-index pair to create node string (node v)
                    nodev = ""
                    nodev = nodev.join([adjlist[i][j][k][0], ":", str(adjlist[i][j][k][1])])  # get name of 2nd node

                    nodelist.append(nodeu)
                    nodelist.append(nodev)

                    # add edge between node string from outer loop and each node from here
                    contactlist.append([nodeu, nodev, weight])

    # once nodes and edges have been established, enumerate list of nodes and create graph using (index1,index2,weight)
    nodelist = list(set(nodelist))  # keep only unique entries in list of nodenames
    nlist = list(enumerate(nodelist))  # enumerate list for passing of index to graph-tool with add_edge_list()
    nodeindexmap = {}  # empty dict for mapping value:index pairs
    g.add_vertex(n=len(nodelist))  # add all vertices to graph first

    for n in nlist:
        nodeindexmap[n[1]] = n[0]  # add node index:name mapping
        vname[n[0]] = n[1]          # add name to vertex properties

    # prune edge list so a->b, b->a duplicates are combined and their weights are summed
    # uses a dictionary for ease of ooking up/ modifying weight value
    nameweightdict = {}
    for i in contactlist:
        n1, n2, w = i  # unpack "tuple" (actually a list rn but)
        # initialize strings for forward and reverse edge direction
        namestring = ""
        namestring = namestring.join([n1, "-", n2])
        rvrstring = ""
        rvrstring = namestring.join([n2, "-", n1])

        # if one name is already in the dictionary, update weight instead of adding
        if namestring in nameweightdict or rvrstring in nameweightdict:
            if namestring in nameweightdict:
                nameweightdict[namestring] = nameweightdict[namestring] + w
            if rvrstring in nameweightdict:
                nameweightdict[rvrstring] = nameweightdict[rvrstring] + w
        # else add as normal
        else:
            nameweightdict[namestring] = w

    # convert [[name1, name2, weight],...] to [(index1, index2, weight),...] so it works with add_edge_list()
    iedgelist = []
    for i in nameweightdict.keys():
        n1, n2 = i.split(sep="-", maxsplit=1)  # splits edge name back into constituent node names
        w = nameweightdict[i]
        i1 = nodeindexmap[n1]  #fetch index mapping by name
        i2 = nodeindexmap[n2]  #fetch index mapping by name

        newtuple = (i1, i2, w)
        iedgelist.append(newtuple)

    g.add_edge_list(iedgelist, eprops=[weight])

    # make properties internal before saving internal
    g.vertex_properties["vname"] = vname
    g.edge_properties["weight"] = weight

    g.save(mlname)


def check_contact(binlist, newchrom, newbin):  # listofchromatin contacts from SAMtoGraph,
    contactflag = 0  # flag raised if adjacency is already in list, so it is not appended

    if not binlist:  # check if list of adjacencies is currently empty
        binlist = [[newchrom, newbin, 1]]
    else:
        for i in binlist:  # iterate through the list of adjacencies, +1 weight if found and flag raised
            if newchrom == i[0] and newbin == i[1]:
                i[2] = i[2] + 1
                contactflag = 1  # set flag so a new contact is not appended
                break  # break out of loop once appropriate bin is found
            else:
                continue

        if contactflag:  # if flag raised, do not append contact to list
            pass
        else:
            binlist.append([newchrom, newbin, 1])  # otherwise, do append

    return binlist


# utility function for writing adjlist to file if not being piped directly into hic_adjlist_to_graphml()
def save_adjlist(adjlist, outfile):
    # chromosome list for mapping adjlist 1st index to chromosome name
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    try:
        with open(outfile, "w") as f:
            for chrom in range(len(adjlist)):
                chromname = chrlist[chrom]
                for thisbin in range(len(adjlist[chrom])):
                    for adjacency in adjlist[chrom][thisbin]:
                        writestring = str(chromname) + ":" + str(thisbin) + "\t" + str(adjacency[0]) + ":" + \
                                      str(adjacency[1]) + "\t" + \
                                      str(adjacency[2]) + "\n"  # writes in format: chr1:bin1 \t chr2:bin2 \t weight \n
                        f.write(writestring)  # write line to file
        print("Adjacency list written...")
    except TypeError:
        print("WARNING: TypeError prevented adjlist from writing to file...")


# takes saved adjlist .tsv file and converts it to
def adjlist_file_to_graphml(adjlistfile, outfile="out.graphml"):
    with open(adjlistfile, "r") as f:

        chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
                   'chr12',
                   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
                   'chrY']

        g = Graph(directed=False)

        # create adjacency list to then import into graph with add_edgelist()
        contactlist = []  # for storing contacts
        nodelist = []  # for storing  unique node names, to later be  enumerated
        # initialize graph properties
        vname = g.new_vertex_property("string")  # vp map for node names (chromosome-bin pairs)
        weight = g.new_edge_property("double")  # ep map for edge weight
        # make properties internal
        g.vertex_properties["vname"] = vname
        g.edge_properties["weight"] = weight

        for line in f:
            nodeu, nodev, strweight = line.split()
            floatweight = float(strweight)
            # append nodes to nodelist and contact info to contactlist
            nodelist.append(nodeu)
            nodelist.append(nodev)
            contactlist.append([nodeu, nodev, floatweight])  # is list of lists now, but will later be list of tuples

        # once nodes and edges have been established, enumerate list of nodes and create graph using (index1,index2,weight)
        nodelist = list(set(nodelist))  # keep only unique entries in list of nodenames
        nlist = list(enumerate(nodelist))  # enumerate list for passing of index to graph-tool with add_edge_list()
        nodeindexmap = {}  # empty dict for mapping value:index pairs
        g.add_vertex(n=len(nodelist))  # add all vertices to graph first

        for n in nlist:
            nodeindexmap[n[1]] = n[0]  # add node index:name mapping
            vname[n[0]] = n[1]  # add name to vertex properties

        # prune edge list so a->b, b->a duplicates are combined and their weights are summed
        # uses a dictionary for ease of ooking up/ modifying weight value
        nameweightdict = {}
        for i in contactlist:
            n1, n2, w = i  # unpack "tuple" (actually a list rn but)
            # initialize strings for forward and reverse edge direction
            namestring = ""
            namestring = namestring.join([n1, "-", n2])
            rvrstring = ""
            rvrstring = rvrstring.join([n2, "-", n1])  # reverse string in case some edges are double encoded

            # if one name is already in the dictionary, update weight instead of adding
            if namestring in nameweightdict or rvrstring in nameweightdict:
                if namestring in nameweightdict:
                    nameweightdict[namestring] = nameweightdict[namestring] + w
                if rvrstring in nameweightdict:
                    nameweightdict[rvrstring] = nameweightdict[rvrstring] + w
            # else add as normal
            else:
                nameweightdict[namestring] = w

        # convert [[name1, name2, weight],...] to [(index1, index2, weight),...] so it works with add_edge_list()
        iedgelist = []
        for i in nameweightdict.keys():
            n1, n2 = i.split(sep="-", maxsplit=1)  # splits edge name back into constituent node names
            w = nameweightdict[i]
            i1 = nodeindexmap[n1]  # fetch index mapping by name
            i2 = nodeindexmap[n2]  # fetch index mapping by name

            newtuple = (i1, i2, w)
            iedgelist.append(newtuple)

        g.add_edge_list(iedgelist, eprops=[weight])

    g.save(outfile)


# function for retrieving all keys with a given value
def get_keys_by_value(adict, value):
    keylist = []
    itemlist = adict.items()
    for item in itemlist:
        for id in item[1]:
            if id == value:
                keylist.append(item[0])
    return keylist


# annotates a given hic network with genes and go terms found in gencodefile and mapfile respectively
def genes_go_annotate(graphmlfile, mapfile="HUMAN_9606_idmapping_selected.tsv", gencodefile="gencode_pcgenes.csv",
                      m=40000, outfile="annotatedgonet.graphml", go_obo="go-basic.obo", binfile=None):

    start_time = time.time()

    go = obo_parser.GODag(go_obo)  # initialize go file

    gencode = pd.read_csv(gencodefile)  # load gencode table
    mapping = pd.read_csv(mapfile, sep='\t')  # load mapping table

    # create attribute dictionaries, dict name will be attribute name in graph
    goterms = {}  # key is node name (chr[#]:[bin]), value is list of goterms
    genes = {}  # key is node name (chr[#]:[bin]), value is list of genes

    try:
        if binfile is None:
            # for each bin, take centroid between start and end pos, then convert pos to bin. annotate that bin
            print("creating gene annotations...")
            print("--- %s seconds since start ---" % (time.time() - start_time))
            for index, row in gencode.iterrows():
                chrom = row['seqid']
                startpos = int(row['start'])
                endpos = int(row['end'])
                centroid = int(startpos + endpos / 2)
                centrebin = int(centroid / m)

                # for each node to be annotated, check if key exists, then create/append new gene as necessary
                nodename = chrom + ':' + str(centrebin)
                if nodename in genes:  # if key already exists, append new gene
                    genes[nodename].append(row['id'])
                else:  # else make new gene list that contains gene
                    genes[nodename] = [row['id']]

        else:  # if binfile is included, use that (should be bed format chr#\wstart\wend)
            # for each bin, take centroid between start and end pos, then convert pos to bin. annotate that bin
            print("creating gene annotations...")
            print("--- %s seconds since start ---" % (time.time() - start_time))

            # make a dict for bin ranges
            bindict = {"chr1": [], "chr2": [], "chr3": [], "chr4": [], "chr5": [], "chr6": [], "chr7": [], "chr8": [],
                       "chr9": [], "chr10": [], "chr11": [], "chr12": [], "chr13": [], "chr14": [], "chr15": [],
                       "chr16": [],
                       "chr17": [], "chr18": [], "chr19": [], "chr20": [], "chr21": [], "chr22": [], "chrX": [],
                       "chrY": [],
                       "chrM": []}
            # populate it with bin start and end positions
            f = open(binfile, "r")
            for line in f:
                chrom, startpos, endpos = line.split()
                binbounds = str(startpos) + "-" + str(endpos)
                bindict[chrom].append(binbounds)
            f.close()

            # iterate through all protein coding genes and add them to appropriate bins
            for index, row in gencode.iterrows():
                chrom = row['seqid']
                startpos = int(row['start'])
                endpos = int(row['end'])
                centroid = int(startpos + endpos / 2)
                genename = row['id']

                indexcounter = 0  # counts iters to keep track of bin index, which is used in vname ("[chr]:[bindex]")
                for whichbin in bindict[chrom]:
                    binstart, binend = whichbin.split("-")  # get bin start and end positions from name
                    if int(binstart) <= centroid <= int(binend):  # if centroid is in bin
                        nodename = chrom + ':' + str(indexcounter)  # nodename set to match vname
                        if nodename in genes:  # if key already exists, append new gene
                            genes[nodename] = genes[nodename] + [genename]
                        else:  # else make new gene list that contains gene
                            genes[nodename] = [genename]
                    indexcounter = indexcounter + 1  # increments index

        # for each node, for each gene, add associated GO terms to new dictionary with node as key, GO term as value
        print("creating GO annotations...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        for ind, ro in tqdm(mapping.iterrows()):  # iterate through mapping table
            ids = str(ro.Ensembl)  # split multiple ids into individual ones
            ids = ids.split("; ")
            for emblid in ids:  # iterate through those ensembl ids
                keys = get_keys_by_value(genes, emblid)  # get all nodes with that gene
                for i in keys:  # iterate through that
                    termstring = str(ro.GO)
                    terms = termstring.split("; ")  # coerce string (list of go terms) to list
                    if i not in goterms:  # if no go terms yet, initialize list
                        goterms[i] = terms
                    else:  # else append to list
                        goterms[i] = goterms[i] + terms

        # for each node, for each GO term add all parent terms to dictionary
        print("adding parent GO terms...")
        print("--- %s seconds since start ---" % (time.time() - start_time))

        for node in tqdm(goterms.keys()):
            try:
                ogterms = goterms[node]  # original list of terms
                if 'nan' in ogterms:
                    ogterms.remove('nan')  # remove NaNs
                allterms = ogterms  # list for all terms to be added to, starting with oglist
                for term in ogterms:
                    if term != 'nan':  # skip NaNs
                        rec = go[term]
                        parents = list(rec.get_all_parents())
                        allterms = allterms + parents
                allterms = list(set(allterms))  # removes duplicates from list
                goterms[node] = allterms
            except KeyError as keyerror:
                print("keyerror: ")
                print(keyerror)  # prints GOTERMS not found in file

            except TypeError:  # ignore nans if they still exist
                pass

        print("adding annotations to graph...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        g = load_graph(graphmlfile)
        # initialize new vertex property objects
        geneprop = g.new_vertex_property("vector<string>")
        goprop = g.new_vertex_property("vector<string>")

        # loop iterates through all gene annotations and adds them to geneprop
        for k in genes.keys():
            try:
                v = g.vertex_index[find_vertex(g, g.vp["vname"], str(k))[0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                geneprop[v] = genes[k]
            except IndexError:
                # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                print("no available vertex for " + str(k) + ". skipping...")

        # loop iterates through all go annotations and addds them to goprop
        for k in goterms.keys():
            try:
                v = g.vertex_index[find_vertex(g, g.vp["vname"], str(k))[0]]  # finds vertex index from vertex name (e.g."chr11:23"->5621)
                goprop[v] = goterms[k]
            except IndexError:
                # TODO it might be bad to ignore empty nodes. Should add missing nodes to the graph? how to connect?
                print("no available vertex for " + str(k) + ". skipping...")

        # make property maps internal
        g.vertex_properties["genes"] = geneprop  # networkx naming convention kept for back compatibility
        g.vertex_properties["goterms"] = goprop  # networkx naming convention kept for back compatibility

        print("writing graph to file...")
        print("--- %s seconds since start ---" % (time.time() - start_time))
        g.save(outfile)

    except MemoryError as memerror:
        print(memerror)

    print("Annotated graph saved as " + outfile)
    print("--- %s seconds since start ---" % (time.time() - start_time))


# use mapfile and genes annotated to nodes to annotate nodes with REACTOME pathways
def reactome_annotate(graph, mapfile="Ensembl2Reactome_All_Levels.txt", outfile=None):
    listoflists = []  # initialize list for coersion to df

    # read lines and append to listoflists
    with open(mapfile, "r") as f:
        for line in f.readlines():
            splitline = line.split(sep="\t")
            splitline[-1] = splitline[-1][:-1]  # remove endline character
            listoflists.append(splitline)

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    # filter table so only human genes included
    df = pd.DataFrame.from_records(listoflists)  # coerce to dataframe
    df = df.rename({0: "EMBLid", 1: "Pathway", 2: "URL", 3: "Notes", 4: "Source", 5: "Species"}, axis=1)
    df = df.loc[df["Species"] == "Homo sapiens"]
    # print(df.loc[df["EMBLid"] == "ENSG00000129038"])
    # print(find_vertex(g, g.vp.genes, 'ENSG00000129038'))

    reactomepathways = g.new_vp("vector<string>")

    # iterate over vertices, extract gene list, use genes to annotate with pathways
    for v in g.vertices():
        genelist = g.vp.genes[v]

        # create list of all reactome paths associated with all genes on a node
        reactomelist = []
        for gene in genelist:
            rows = df.loc[df["EMBLid"] == gene]
            for index, row in rows.iterrows():
                reactomelist.append(row["Pathway"])

        # assign full list to node as vertex property
        reactomepathways[v] = reactomelist

    if outfile == None:
        if type(graph) == str:
            filename = graph[:-3] + "_REACTOMEannotated.gt"
            g.save(filename)
        else:
            filename = "REACTOMEannotatedGraph.gt"
            g.save(filename)
    else:
        g.save(outfile)


# converts gencode gff3 file to parsed, filtered table
def parse_gff3(infile, outfile='./gencode_pcgenes.csv'):
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    mylist = [[]]  # empty list for later coersion to data frame

    with open(infile) as myfile:
        lines = myfile.readlines()
    for line in lines:  # goes line by line and coerces final table to data frame
        if not line.startswith('#'):  # ignores header line
            mylist.append(line.split()[0:9])
    mycolumns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gencode = pd.DataFrame(mylist, columns=mycolumns)

    gencode = gencode.where(gencode.type == "gene")  # only include genes in table
    gencode = gencode.where(gencode.seqid.isin(chrlist))  # only include nuclear chromosomes
    gencode['id'] = gencode.attributes.str.split(';').str.get(0)
    gencode['gene_type'] = gencode.attributes.str.split(';').str.get(2)
    gencode = gencode.dropna()
    gencode = gencode.where(gencode.gene_type == "gene_type=protein_coding")  # only include protein coding genes
    gencode = gencode.dropna()
    gencode['id'] = gencode.id.str.split('=').str.get(1)  # remove 'gene_id=' from id
    gencode['id'] = gencode.id.str.split('.').str.get(0)  # remove decimals from end of gene id

    gencode.to_csv(outfile)


# writes a weighted adjacency list in COO (sparse matrix) format
def output_coo(graph, outfile="graph.coo"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print(
            "bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    f = open(outfile, "w")

    # get edgelist
    # convert to COO and store in file
    for s, t, w in g.iter_edges([g.ep.weight]):  # source, target, weight
        f.write(str(s) + " " + str(t) + " " + str(w) + "\n")
    f.close()


# performs sinkhorn-knopp balancing on a preconstructed chromatin conformation graph
# deprecated
# def skbalance(graphmlfile, outfile="normed.graphml"):
#
#     g = nx.read_graphml(graphmlfile)  # read in graphml
#
#     newg = nx.to_pandas_adjacency(g)  # convert to adjacency matrix
#     sk = skp.SinkhornKnopp()
#     normg = sk.fit(newg)  # run balancing
#     normg = pd.DataFrame(normg, index=newg.index, columns=newg.columns)  # convert back to dataframe
#     normg = nx.from_pandas_adjacency(normg)  # import back to networkx
#     nx.write_graphml(normg, outfile)  # write back to gml


# performs Knight-Ruiz balancing on a preconstructed chromatin conformation graph (faster than skbalance)
# can take either graphml or gt files
def kr_balance(graph, resolution, outfile="normed.graphml"):

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    f = open("tempCOOfile.txt", "w")  # TODO name this properly? or delete after running?

    # get edgelist
    # convert to COO and store in file
    for s, t, w in g.iter_edges([g.ep.weight]):  # source, target, weight
        f.write(str(s) + " " + str(t) + " " + str(w) + "\n")
    f.close()

    # do KR normalization
    del g  # clear up memory for normalization
    cooReader = gmlib.importer.CooMatrixHandler("tempCOOfile.txt", resolution=resolution, coordinate='index')  # import and convert COO edgelist
    cooReader.save_ccmaps('tempccmap.ccmap', xlabels="IndexPosition")  # TODO name this properly? or delete after running?
    del cooReader  # Delete object and generated any temporary files
    gmlib.normalizer.normalizeCCMapByKR('tempccmap.ccmap', outFile='normedccmap.ccmap', memory='RAM')

    # convert back to edgelist?
    ccmap = gmlib.ccmap.load_ccmap('normedccmap.ccmap')
    gmlib.ccmap.export_cmap(ccmap, "normedCOO.txt", doNotWriteZeros=True)
    del ccmap

    # update edge weights with normalized values
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    f = open("normedCOO.txt", "r")  # open normalized COO for reading

    # make new edge weight dict to become new EdgePropertyMap
    eweight_dict = g.new_edge_property("double")
    for line in f:
        splitline = line.split()
        thisedge = g.edge(int(splitline[0]), int(splitline[1]))
        eweight_dict[thisedge] = float(splitline[2])  # add new weight to dict

    # make new edge map internal, replacing old weight property
    g.edge_properties["weight"] = eweight_dict

    # write new graph file
    g.save(outfile)
    # remove temp files
    os.remove("tempCOOfile.txt")
    os.remove("normedCOO.txt")


# make a dictionary of 'GOterm': [nodes with that term]
def make_go_dict(graph, prop="goterms"):
    godict = {}  # initialize dict
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    for v in g.vertices():  # iterate through nodes
        if prop == "goterms":
            terms = g.vp.goterms[v]
        elif prop == "reactomepathways":
            terms = g.vp.reactomepathways[v]
        # clean up string before splitting into list (deprecated from old networkx method)
        # gostring = gostring.replace("\'", "")
        # gostring = gostring.replace(",", "")
        # gostring = gostring.replace("[", "")
        # gostring = gostring.replace("]", "")
        # terms = list(set(gostring.split()))  # TODO why are redundant vertices being added to the list?
        try:
            for term in terms:  # iterate through terms
                if term in godict:
                    godict[term].append(g.vertex_index[v])  # append new node to list if it already exists
                else:
                    godict[term] = [g.vertex_index[v]]  # else make a new list for terms
        except KeyError as keyerror:
            pass
    print(prop + " dictionary created")
    return godict


# make a histogram where y is number of GOTERMs found on x number of nodes
def plot_goterms_per_node_histogram(adict, outfile="GOhistogram.pdf", xlim=250, ylim=250):
    mydict = {}
    dictogram = np.array(list(adict.items()))

    for row in tqdm(dictogram):
        if len(row[1]) in mydict:
            mydict[len(row[1])] += 1
        else:
            mydict[len(row[1])] = 1

    f = plt.bar(list(mydict.keys()), mydict.values())
    plt.xlim(0, xlim)
    plt.ylim(0, ylim)
    plt.show()
    plt.savefig(outfile, bbox_inches='tight')

    for item, value in mydict.items():
        print(item, ":", value)


# creates a nodelist and edgelist for parallelization of dijkstra
# deprecated
# def makenodesedgeslists(graphmlfile, nodesfile="hicGraphNodes.txt", edgesfile="hicGraphEdges.csv"):
#     g = nx.read_graphml(graphmlfile)
#     e = list(g.edges.data('weight', default=0))
#     v = list(g.nodes(data=False))
#
#     with open(nodesfile, 'w+') as f1:  # write nodes to file
#         for item in v:
#             f1.write(str(item) + "\n")
#
#     with open(edgesfile, 'w+') as f2:  # write edges to file
#         writer = csv.writer(f2)
#         writer.writerows(e)


# update all weights in graphs so that lower weight = closer proximity so dijkstra does what we want
def invert_weights(graphfile, outfile="inverted.graphml"):
    g = load_graph(graphfile)
    newweights = g.new_edge_property("double")  # new vp for inverted weights

    for s, t, w in g.iter_edges([g.ep.weight]):  # source, target, weight:
        weight = w
        if 0 > weight > 1:
            print("Please make sure the graph has been normalized before invoking invert_weights()")
            sys.exit()
        else:
            weight = 1 - weight
            newweights[g.edge(s, t)] = weight

    g.ep["weight"] = newweights  # overwrite weight vp with new weights and make internal before saving

    g.save(outfile)


# converts graphml to gt for faster opening in graph tool
def graphml_to_gt(graph, outfile="graph.gt"):
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
    g.save(outfile, fmt="gt")

# converts gt to graphml for utility
def gt_to_graphml(graph, outfile="graph.graphml"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    g.save(outfile, fmt="graphml")


# splits list of goterms (which is a string) to property map for use with find_vertex()
# Is this necessary? Update: No, but it's a good reference for how property maps work in graph-tool
def map_goterms(graph):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    goterm = g.new_vp("vector<string>")
    for v in g.vertices():
        gostring = g.vp.goterms[v]
        # clean up string before splitting into list
        gostring = gostring.replace("\'", "")
        gostring = gostring.replace(",", "")
        gostring = gostring.replace("[", "")
        gostring = gostring.replace("]", "")
        golist = gostring.split()

        goterm[v] = golist

    if type(graph) == str:
        g.save(graph, fmt="gt")
    else:
        return g


# function that obtains sampling vector from network
# each vertex is added a number of times equal to the number of specified annotations it has
# possible annotations are "goterms", "genes", "reactomepathways"
# TODO dont include goterms above threshold that will be counted
def get_sampling_vector(graph, vmax=None, prop="goterms"):
    if prop == "goterms":
        nndict = make_numnodes_dict(graph)
    elif prop == "reactomepathways":
        nndict = make_numnodes_dict(graph, prop="reactomepathways")

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    samplevec = []
    for v in g.vertices():
        golist = g.vp.goterms[v]

        if vmax:  # if vmax is set, don't append node when annotated with goterm with higher number of nodes
            for i in golist:
                if nndict[i] <= vmax:
                    samplevec.append(g.vertex_index[v])

        else:  # otherwise append node for every term
            for i in golist:
                samplevec.append(g.vertex_index[v])

    return samplevec


# pull a weighted sample from the population
# works by pulling samples randomly from the list, coercing to set to remove duplicates, and repeating until full
def get_random_tpd(distarray, samplevec, vlength):
    samplepop = random.sample(samplevec, vlength)  # get a random sample from our weighted vector

    # until there are vlength unique vertices, keep drawing from weighted vector
    while len(set(samplepop)) != vlength:
        samplepop.extend(random.sample(samplevec, vlength - len(set(samplepop))))
    samplepop = list(set(samplepop))  # condense down to appropriately sized list of unique vals

    combs = combinations(samplepop, 2)  # get all possible combinations for sampled nodes
    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory
    row = list(row)
    col = list(col)
    tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
    return tpd


# top percent pairwise distance
# pull a weighted sample from the population
# works by pulling samples randomly from the list, coercing to set to remove duplicates, and repeating until full
def get_random_tptpd(distarray, samplevec, vlength, tpthresh=0.2):
    samplepop = random.sample(samplevec, vlength)  # get a random sample from our weighted vector

    # until there are vlength unique vertices, keep drawing from weighted vector
    while len(set(samplepop)) != vlength:
        samplepop.extend(random.sample(samplevec, vlength - len(set(samplepop))))

    # TODO assess clustering of nodes and cut out bottom% nodes
    # for speed, put all combs into numpy array, vectorize dist lookup?
    # np.fsum(where=) to get separate arrays for each node, then get sum(tpd) for each before sorting and cutting
    samplepop = list(set(samplepop))  # condense down to appropriately sized list of unique vals
    combs = combinations(samplepop, 2)  # get all possible combinations for sampled nodes
    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory

    # get array of combinations with respective shortest path distances
    distlist = distarray[row, col]
    tosumarray = np.array([[[row, col]],
                           [[distlist, distlist]]])  # create a 2x2xncombs np array with axes source x target x dist
    del distlist

    # for each node, get it's individual tpd and add it to a list to be sorted and trimmed
    tpdlist = []
    for i in range(len(samplepop)):
        mask = tosumarray[0, :, :] == i
        tpd = np.sum(tosumarray[1, :, :], where=mask)  # vectorized addition of PD at all combs for TPD calculation
        tpdlist.append(tpd)

    # add nodes with tpds into data frame, sort it, calculate regular tpd for top% of nodes
    tosortframe = pd.DataFrame({"node": samplepop, "tpd": tpdlist})
    sortedframe = tosortframe.sort_values('tpd', ascending=False)

    numtop = int(vlength * tpthresh)
    if numtop <= 3:  # make sure there are always at least 3 nodes in set
        numtop = 3
        # TODO make this so it can handle 2?

    newnodeslist = sortedframe['node'].tolist()
    newnodeslist = newnodeslist[0:numtop]  # cut down list to only top nodes

    #  recalculate tpd for only top nodes
    combs = combinations(newnodeslist, 2)  # get all possible combinations for sampled nodes
    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory
    row = list(row)
    col = list(col)
    tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation

    return tpd


# generate a vector containing the TPD distribution for a given number of vertices (vlength) multiple times
# vlengths can either be an int or a list of ints
def montecarlo_sample_tpds(distmatrix, vlengths, graph, prop="goterms", m=1000000, metric="tpd",
                           ncores=1, outfile=None, tpthresh=0.2, sampling="weighted"):

    start = timeit.default_timer()

    if type(vlengths) == int:  # if only one val given for vlengths, convert to iterable (list)
        vlengths = [vlengths]

    # initialize dict with empty lists for appending of sample tpds
    tpddict = {}
    for i in vlengths:
        tpddict[i] = []

    # get weighted list(vector) of nodes for sampling from
    if sampling == "weighted":  # do weighted sampling
        svec = get_sampling_vector(graph, prop=prop)
    else:  # don't do unweighted sampling
        # load graph
        if type(graph) == str:
            g = load_graph(graph)
        elif type(graph) == Graph:
            g = graph
        else:
            print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")
        svec = list(g.get_vertices())

    # load dist array from file
    distarray = load_dist_matrix(distmatrix)

    print("array created")
    stop = timeit.default_timer()  # Timer Block
    print('Time: ', stop - start)  # Timer Block

    # perform appropriate sampling based on clustering measure selected (tpd or tptpd)
    if metric == "tpd":
        # Perform MC sampling for each vlength. Sampling is done in parallel.
        for numv in vlengths:  # for each vlength, calculate sample dist in parallel
            pbardescr = str("k = " + str(numv))
            pbar = tqdm(total=m, desc=pbardescr)  # set a progress bar to monitor progress
            with concurrent.futures.ThreadPoolExecutor(max_workers=ncores) as executor:
                futuretpd = {executor.submit(get_random_tpd, distarray, svec, numv): i for i in range(m)}
                for future in concurrent.futures.as_completed(futuretpd):
                    tpddict[numv].append(future.result())
                    pbar.update()  # update progress bar for each future

            printstring = ','.join(map(str, tpddict[numv]))
            if outfile is None:
                print("%s,%s\n" % (numv, printstring))
            else:
                with open(outfile, "a+") as mcf:
                    mcf.write("%s,%s\n" % (numv, printstring))
            stop = timeit.default_timer()  # Timer Block
            print(str(numv) + " completed after " + str(stop - start) + " seconds")

    elif metric == "tptpd":
        # Perform MC sampling for each vlength. Sampling is done in parallel.
        for numv in vlengths:  # for each vlength, calculate sample dist in parallel
            pbardescr = str("k = " + str(numv))
            pbar = tqdm(total=m, desc=pbardescr)  # set a progress bar to monitor progress
            with concurrent.futures.ThreadPoolExecutor(max_workers=ncores) as executor:
                futuretpd = {executor.submit(get_random_tptpd, distarray, svec, numv, tpthresh=tpthresh): i for i in
                             range(m)}
                for future in concurrent.futures.as_completed(futuretpd):
                    tpddict[numv].append(future.result())
                    pbar.update()  # update progress bar for each future

            printstring = ','.join(map(str, tpddict[numv]))
            if outfile is None:
                print("%s,%s\n" % (numv, printstring))
            else:
                with open(outfile, "a+") as mcf:
                    mcf.write("%s,%s\n" % (numv, printstring))
            stop = timeit.default_timer()  # Timer Block
            print(str(numv) + " completed after " + str(stop - start) + " seconds")

        stop = timeit.default_timer()  # Timer Block
        print('End time: ', stop - start)  # Timer Block
        print(tpddict)
    return tpddict


# loads distance matrix from file into memory, for use by get_tpd, get_random_tpd, etc.
def load_dist_matrix(distmatrix, indexcol=False):
    start = timeit.default_timer()  # for timing duh

    # open dist matrix
    with open(distmatrix, "r") as f:
        if indexcol:  # old cpp script generated an index column but the new method doesn't, so use this arg if needed
            arraysize = len(f.readline().split()) - 1  # array size is -1 because of index column
        else:
            arraysize = len(f.readline().split())
        distarray = np.zeros((arraysize, arraysize))
        flist = f.readlines()
        del flist[0]
        print("file read")
        stop = timeit.default_timer()  # Timer Block
        print('Time: ', stop - start)  # Timer Block

        # store in memory as an array for speedy access
        for line in range(len(flist)):
            try:
                distarray[line] = np.asarray(flist[0].split())[1:]  # add split line to ndarray, skipping first row/col
                del flist[0]  # clear up memory
            except IndexError as ie:
                print(ie)
                pass
        print("array created")
        stop = timeit.default_timer()  # Timer Block
        print('Time: ', stop - start)  # Timer Block

    return distarray


# utility function to calculcate All Pairs Shortest Paths via graph-tool and save all shortest distances as annotations
def annotate_apsp(graphfile, outfile="distAnnotatedGraph.graphml"):
    g = load_graph(graphfile)  # takes first argument as input graph filename
    dist = shortest_distance(g, weights=g.ep.weight)  # calculate all shortest paths and store as vertex property map
    g.vp["shortest_distances"] = dist  # make property internal
    g.save(outfile)


# creates a distance matrix from the graph-annotated distances calculated with gt.shortest_distance()
# row and column indices are based on vertex indices from the graph
def distmatrix_from_annotations(graphfile, outfile="distmatrix.tsv"):
    g = load_graph(graphfile)  # load graph
    arraylist = []  # initialize empty list for 1d arrays to be appended to

    for v in g.vertices():  # iterate through all vertices
        try:
            distvec = np.array(g.vp.shortest_distances[v])
            arraylist.append(distvec)
        except:
            print("no annotation called 'shortest_distances' found. Please use gt.shortest_distances to find APSP")

    distarray = np.vstack(arraylist)  # create 2D array from the list of 1D arrays
    pd.DataFrame(distarray).to_csv(outfile, sep=" ")
    return(distarray)


# get tpd for a single set
def get_tpd(nodelist, distarray):

    combs = combinations(nodelist, 2)  # get all possible combinations for nodes annotated with term

    col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
    del combs  # delete to clear up memory
    row = list(row)
    col = list(col)
    tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
    print("TPD of set: " + str(tpd))

    return tpd


# get tptpd for a single set
def get_tptpd(nodelist, distmatrix, tpthresh=0.4):
    start = timeit.default_timer()  # Timer Block
    # open dist matrix
    with open(distmatrix, "r") as f:
        arraysize = len(f.readline().split()) - 1  # array size is -1 because of index column
        distarray = np.zeros((arraysize, arraysize))
        flist = f.readlines()
        del flist[0]
        print("file read")
        stop = timeit.default_timer()  # Timer Block
        print('Time: ', stop - start)  # Timer Block

        # store in memory as an array for speedy access
        for line in range(len(flist)):
            try:
                distarray[line] = np.asarray(flist[0].split())[1:]  # add split line to ndarray, skipping first row/col
                del flist[0]  # clear up memory
            except IndexError as ie:
                print(ie)
                pass
        print("array created")

        # do TPTPD calculation
        combs = combinations(nodelist, 2)  # get all possible combinations for nodes annotated with term

        col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
        del combs  # delete to clear up memory

        # get array of combinations with respective shortest path distances
        distlist = distarray[row, col]
        tosumarray = np.array([[[row, col]],
                               [[distlist,
                                 distlist]]])  # create a 2x2xncombs np array with axes source x target x dist
        del distlist

        # for each node, get it's individual tpd and add it to a list to be sorted and trimmed
        tpdlist = []
        for i in range(len(nodelist)):
            mask = tosumarray[0, :, :] == i
            tpd = np.sum(tosumarray[1, :, :],
                         where=mask)  # vectorized addition of PD at all combs for TPD calculation
            tpdlist.append(tpd)

        # add nodes with tpds into data frame, sort it, calculate regular tpd for top% of nodes
        tosortframe = pd.DataFrame({"node": nodelist, "tpd": tpdlist})
        sortedframe = tosortframe.sort_values('tpd', ascending=False)

        numtop = int(len(nodelist) * tpthresh)
        if numtop <= 3:  # make sure there are always at least 3 nodes in set
            numtop = 3
            # TODO make this so it can handle 2?

        newnodeslist = sortedframe['node'].tolist()
        newnodeslist = newnodeslist[0:numtop]  # cut down list to only top nodes

        #  recalculate tpd for only top nodes
        combs = combinations(newnodeslist, 2)  # get all possible combinations for sampled nodes
        col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
        del combs  # delete to clear up memory
        row = list(row)
        col = list(col)
        tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
        print("Top% (" + str(tpthresh * 100) + "%) TPD: " + str(tpd))

    return tpd


# calculate tpd for all goterms in a given godict (created by make_go_dict())
def all_go_tpd(godict, distmatrix, outfile):
    start = timeit.default_timer()  # for timing duh
    keys = list(godict.keys())  # initialize list of keys

    distarray = load_dist_matrix(distmatrix)
    print("Distance matrix loaded...")

    # calculate tpd for all go terms as found in keys
    for term in tqdm(keys):
        combs = combinations(godict[term], 2)  # get all possible combinations for nodes annotated with term

        try:
            col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
            del combs  # delete to clear up memory
            row = list(row)
            col = list(col)
            tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
            print(term + "TPD: " + str(tpd))
            with open(outfile, "a+") as f:
                f.write(term + ", " + str(tpd) + "\n")
                print(term + " TPD written to " + outfile)
        except ValueError as ve:
            print("ERROR ON " + str(term))
            print(ve)
            print("continuing...")
            continue

    stop = timeit.default_timer()  # Timer Block
    print('Total Runtime: ', stop - start)  # Timer Block

    return None


# calculate tpd for all goterms in a given godict (created by make_go_dict())
def all_go_tptpd(godict, distmatrix, outfile, tpthresh=0.2):
    start = timeit.default_timer()  # for timing duh
    keys = list(godict.keys())  # initialize list of keys

    # open dist matrix
    distarray = load_dist_matrix(distmatrix)
    print("Distance matrix loaded...")

    # calculate tpd for all go terms as found in keys
    for term in tqdm(keys):
        combs = combinations(godict[term], 2)  # get all possible combinations for nodes annotated with term

        try:
            col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
            del combs  # delete to clear up memory

            # get array of combinations with respective shortest path distances
            distlist = distarray[row, col]
            tosumarray = np.array([[[row, col]],
                                   [[distlist,
                                     distlist]]])  # create a 2x2xncombs np array with axes source x target x dist
            del distlist

            # for each node, get it's individual tpd and add it to a list to be sorted and trimmed
            tpdlist = []
            for i in range(len(godict[term])):
                mask = tosumarray[0, :, :] == i
                tpd = np.sum(tosumarray[1, :, :],
                             where=mask)  # vectorized addition of PD at all combs for TPD calculation
                tpdlist.append(tpd)

            # add nodes with tpds into data frame, sort it, calculate regular tpd for top% of nodes
            tosortframe = pd.DataFrame({"node": godict[term], "tpd": tpdlist})
            sortedframe = tosortframe.sort_values('tpd', ascending=False)

            numtop = int(len(godict[term]) * tpthresh)
            if numtop <= 3:  # make sure there are always at least 3 nodes in set
                numtop = 3


            newnodeslist = sortedframe['node'].tolist()
            newnodeslist = newnodeslist[0:numtop]  # cut down list to only top nodes

            #  recalculate tpd for only top nodes
            combs = combinations(newnodeslist, 2)  # get all possible combinations for sampled nodes
            col, row = list(zip(*combs))  # unpack combinations into lists so they can be used as indices
            del combs  # delete to clear up memory
            row = list(row)
            col = list(col)
            tpd = math.fsum(distarray[row, col])  # vectorized addition of PD at all combs for TPD calculation
            print(term + "Top% (" + str(tpthresh * 100) + "%) TPD: " + str(tpd))

            with open(outfile, "a+") as f:
                f.write(term + ", " + str(tpd) + "\n")
                print(term + " TPD written to " + outfile)
        except ValueError as ve:
            print("ERROR ON " + str(term))
            print(ve)
            print("continuing...")
            continue

    stop = timeit.default_timer()  # Timer Block
    print('Total Runtime: ', stop - start)  # Timer Block

    return None


# get list of ks for a graph, where k is the number of vertices annotated with a particular GO term
def get_vlengths(graph, prop="goterms"):
    godict = make_go_dict(graph, prop=prop)
    vlengthlist = []
    for key in godict.keys():
        if len(godict[key]) < 2:  # we only care about clusters of size two or more
            continue
        vlengthlist.append(len(godict[key]))
    vlengthlist = sorted(list(set(vlengthlist)))
    return vlengthlist


# pairwise swapping of goterms between random nodes for generation of False Discovery Rate
# 1000000000x(for random pair of nodes, swap one value of 'goterms')
# vars refer to go, as in "golist", but could be reactome or any other annotation set as well
def label_swap(graph, label="goterms", outfile="swappedGraph.gt", m=3000000):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    pbar = tqdm(total=m)  # initialize loading bar

    numnodes = len(g.get_vertices())  # get number of nodes (vertices) in graph
    numswaps = 1  # track number of actual swaps that are made

    while numswaps <= m:  # loop until desired number of swaps is achieved
        nodepair = random.sample(range(numnodes), 2)  # get a random pair of nodes from the graph

        # extract list of goterms/pathways for random node 1
        v1 = g.vertex(nodepair[0])
        if label == "goterms":
            golist1 = list(g.vp.goterms[v1])
        elif label == "reactomepathways":
            golist1 = list(g.vp.reactomepathways[v1])
        elif label == "chromosomes":
            golist1 = list(g.vp.vname[v1])

        # extract list of goterms/pathways for random node 2
        v2 = g.vertex(nodepair[1])
        if label == "goterms":
            golist2 = list(g.vp.goterms[v2])
        elif label == "reactomepathways":
            golist2 = list(g.vp.reactomepathways[v2])
        elif label == "chromosomes":
            golist2 = list(g.vp.vname[v2])

        # check to make sure neither list is empty, if that's true then swap a term in one list for one in the other
        if len(golist1) != 0 and len(golist2) != 0:
            t1 = golist1[random.sample(range(len(golist1)), 1)[0]]
            t2 = golist2[random.sample(range(len(golist2)), 1)[0]]
            ind1 = golist1.index(t1)
            ind2 = golist2.index(t2)
            golist1[ind1] = t2
            golist2[ind2] = t1
            if label == "goterms":
                g.vp.goterms[v1] = golist1
            elif label == "reactomepathways":
                g.vp.reactomepathways[v1] = golist1
            if label == "goterms":
                g.vp.goterms[v2] = golist2
            elif label == "reactomepathways":
                g.vp.reactomepathways[v2] = golist2

            numswaps = numswaps + 1  # count that swap was made
            pbar.update()
        # if either node has no goterms, pass
        else:
            pass

    g.save(outfile, fmt="gt")
    return None


# remove nodes with degree greater than degreecutoff from graph
def trim_nodes_by_degree(graph, outfile="trimmed.gt", degreecutoff=5000):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    del_list = []
    for v in g.vertices():
        if v.out_degree() >= 5000:
            del_list.append(v)
    for v in reversed(sorted(del_list)):
        g.remove_vertex(v)
    g.save(outfile, fmt="graphml")


# generates csv with every GO term and the number of nodes annotated with that term
def make_numnodes_dict(graph, prop="goterms"):
    if prop == "goterms":
        mygodict = make_go_dict(graph)
    elif prop == "reactomepathways":
        mygodict = make_go_dict(graph, prop="reactomepathways")
    nodesdict = {}
    for key in mygodict.keys():
        nodesdict[key] = len(mygodict[key])
    return nodesdict


# get pvals for every go term
def get_go_tpd_pvals(graph, tpdfile, shuftpdfile, distrnfile, approxdist=False, prop="goterms"):
    # read in all files
    if prop == "goterms":
        nndict = make_numnodes_dict(graph)  # nndict is goterm:numnodeswithgoterm
    elif prop == "reactomepathways":
        nndict = make_numnodes_dict(graph, prop="reactomepathways")

    # open Goterm,tpd .csv file as goterm:tpd dict
    tpddict = {}
    with open(tpdfile, "r") as f:
        for line in f:
            splitline = line.split(", ")
            tpddict[splitline[0]] = splitline[1][:-2]

    # open goterm:shuffledtpd csv file as goterm:shuffledtpd dict
    shuftpddict = {}
    with open(shuftpdfile, "r") as f:
        for line in f:
            splitline = line.split(", ")
            shuftpddict[splitline[0]] = splitline[1][:-2]

    # opens knodes:MCdistribution csv as knodes:[MCdistribution] dict
    distdict = {}
    with open(distrnfile, "r") as f:
        for line in f:
            splitline = line.split(",")
            distdict[splitline[0]] = sorted(splitline[2:-3])
            distdict[splitline[0]] = [float(i) for i in distdict[splitline[0]]]  # convert to float

    pvalsdict = {}  # initialize goterm:pval dict
    shufpvalsdict = {}  # initialize goterm:shuffledpval dict

    # for each key (goterm) in our files, look up k, then find that MCdistribution
    # for i in tqdm(tpddict.keys()):
    for i in tpddict.keys():
        if str(nndict[i]) not in distdict:  # skip values of k with no MCdistribution
            print(str(nndict[i]) + " not in MC distribution")
            continue
        if i not in shuftpddict.keys():  # make sure term is in both dicts
            print(str(i) + " not found in both TPD dictionaries")
            continue

        # TODO fix this / make sure it works
        if approxdist:  # approximate the tpd distribution function instead of using the empirical one
            shape, location, scale = sps.gamma.fit(distdict[str(nndict[i])])  # unpack tuple of shape and scale params
            pval = sps.gamma.pdf(np.float64(tpddict[i]), shape, location, scale)  # pulls pval for i from fit distrn
            pvalsdict[i] = pval
            # print("TPD: " + str(tpddict[i]))
            # print("realpval: " + str(pvalsdict[i]))
            # repeat for shuffled set
            pval = sps.gamma.pdf(np.float64(shuftpddict[i]), shape, location, scale)  # pulls pval for i from fit distrn
            shufpvalsdict[i] = pval
            # print("shufTPD: " + str(tpddict[i]))
            # print("shufpval: " + str(pvalsdict[i]))

        else:
            # for every value in distribution, count if it is greater than current goterm's TPD
            # this means a low p-value is associated with our TPD being lower than most vals in MCdist
            counter = 0
            for j in distdict[str(nndict[i])]:
                if float(tpddict[i]) >= float(j):
                    counter += 1
            # print("\n")
            # print(i)
            # print(nndict[i])
            # print("TPD: " + str(tpddict[i]))
            pval = counter / len(distdict[str(nndict[i])])
            pvalsdict[i] = pval
            # print("realpval: " + str(pvalsdict[i]))

            # do the same for the shuffled set
            counter = 0
            for j in distdict[str(nndict[i])]:
                if float(shuftpddict[i]) >= float(j):  # TODO make sure is > not <
                    counter += 1
            # print("shufTPD: " + str(shuftpddict[i]))
            shufpval = counter / len(distdict[str(nndict[i])])
            shufpvalsdict[i] = shufpval
            # print("shufpval: " + str(shufpvalsdict[i]))

        print(str(i) + " pval: " + str(pval) + " / shufpval: " + str(shufpval))
    # coerce to dataframe
    df = pd.DataFrame([nndict, tpddict, pvalsdict, shuftpddict, shufpvalsdict],
                      index=["nnodes", "tpd", "pval", "shuftpd", "shufpval"])
    df = df.dropna(axis=1)

    return df


# get pvals for generic sets from file
# tpdfiles will need to have columns setname, nnodes, clusterscore(tpd)
def get_tpd_pvals(tpdfile, shuftpdfile, distrnfile, approxdist=True):
    # open normal and shuffled clusterscore(tpd) files as setname, nnodes, clusterscore(tpd) df
    df = pd.read_csv(tpdfile, names=["setname", "nnodes", "tpd"], header=0)
    shufdf = pd.read_csv(shuftpdfile, names=["setname", "nnodes", "shuftpd"], header=0)

    # explicitly type all columns so stupid merging will work
    df.setname.astype(str)
    df.nnodes.astype(int)
    df.tpd.astype(float)
    shufdf.setname.astype(str)
    shufdf.nnodes.astype(int)
    shufdf.shuftpd.astype(float)

    df = df.merge(shufdf, on="setname", suffixes=["DROP", None]).filter(regex="^(?!.*DROP)")  # left joins both dfs and drops repeat cols

    # opens knodes:MCdistribution csv as knodes:[MCdistribution] dict
    distdict = {}
    with open(distrnfile, "r") as f:
        for line in f:
            splitline = line.split(",")
            distdict[splitline[0]] = sorted(splitline[2:-3])
            distdict[splitline[0]] = [float(i) for i in distdict[splitline[0]]]  # convert to float

    pvalsdict = {}  # initialize set:pval dict
    shufpvalsdict = {}  # initialize set.:shuffledpval dict

    # port legend
    # i = row[0]
    # tpddict[i] = row[1]
    # nndict[i] = row[2]
    # shuftpddict[i] = row[3]

    # for each set in our files find that MCdistribution and get pval
    for row in df.itertuples(index=False):  # access values like row[0]
        if str(row[2]) not in distdict:  # skip values of k with no MCdistribution
            continue

        if approxdist:  # approximate the tpd distribution function instead of using the empirical one
            shape, location, scale = sps.gamma.fit(
                distdict[str(row[2])])  # unpack tuple of shape and scale params
            pval = sps.gamma.pdf(np.float64(row[1]), shape, location, scale)  # pulls pval for i from fit distrn
            pvalsdict[row[0]] = pval
            # print("TPD: " + str(tpddict[i]))
            # print("realpval: " + str(pvalsdict[i]))
            # repeat for shuffled set
            pval = sps.gamma.pdf(np.float64(row[3]), shape, location, scale)  # pulls pval for i from fit distrn
            shufpvalsdict[row[0]] = pval
            # print("shufTPD: " + str(tpddict[i]))
            # print("shufpval: " + str(pvalsdict[i]))

        else:
            # for every value in distribution, count if it is greater than current goterm's TPD
            # this means a low p-value is associated with our TPD being lower than most vals in MCdist
            counter = 0
            for j in distdict[str(row[2])]:
                if float(row[1]) >= float(j):
                    counter += 1
            # print("\n")
            # print(i)
            # print(nndict[i])
            # print("TPD: " + str(tpddict[i]))
            pval = counter / len(distdict[str(row[2])])
            pvalsdict[row[0]] = pval
            # print("realpval: " + str(pvalsdict[i]))

            # do the same for the shuffled set
            counter = 0
            for j in distdict[str(row[2])]:
                if float(row[3]) >= float(j):
                    counter += 1
            # print("shufTPD: " + str(shuftpddict[i]))
            shufpval = counter / len(distdict[str(row[2])])
            shufpvalsdict[row[0]] = shufpval
            # print("shufpval: " + str(shufpvalsdict[i]))

    # coerce to dataframe
    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
    pvalsdf = pd.DataFrame.from_dict(pvalsdict, orient="index", columns=["pval"])
    pvalsdf.reset_index(inplace=True)
    pvalsdf = pvalsdf.rename(columns={'index': 'setname'})
    pvalsdf['setname'] = pvalsdf['setname'].astype(str)

    shufpvalsdf = pd.DataFrame.from_dict(shufpvalsdict, orient="index", columns=["shufpval"])
    shufpvalsdf.reset_index(inplace=True)
    shufpvalsdf = shufpvalsdf.rename(columns={'index': 'setname'})
    shufpvalsdf['setname'] = shufpvalsdf['setname'].astype(str)

    df['setname'] = df['setname'].astype(str)

    df = df.merge(pvalsdf, on="setname", suffixes=["DROP", None]).filter(regex="^(?!.*DROP)")
    df = df.merge(shufpvalsdf, on="setname", suffixes=["DROP", None]).filter(regex="^(?!.*DROP)")
    df = df.dropna(axis=1)

    return df


# draws histogram of significant pvals at decreasing thresholds to determine FDR cutoff (Panel A)
def plot_fdr_pval_histogram(df, title="FDR Cutoff Plot", stepsize=0.0001):
    pylist = []
    shufpylist = []
    xlist = list(np.arange(0.1, 0, -stepsize))

    # create list of bar
    for i in tqdm(xlist):
        pylist.append(len(list(filter(lambda x: x < i, pd.to_numeric(df["pval"], errors='coerce')))))
        shufpylist.append(len(list(filter(lambda x: x < i, pd.to_numeric(df["shufpval"], errors='coerce')))))

    fig = go.Figure(data=[
        go.Bar(name='Real Graph', x=xlist, y=pylist),
        go.Bar(name='Shuffled Graph', x=xlist, y=shufpylist)
    ], layout=go.Layout(plot_bgcolor='rgba(0,0,0,0)'))
    # Change the bar mode
    fig.update_xaxes(type="log")
    fig.update_layout(barmode='group',
                      title=title,
                      bargap=0,
                      xaxis_title="p-value threshold",
                      yaxis_title="Number of GO terms with p-value lower than threshold")
    return fig


# takes df from gettpdvals() and outputs a comparative distribution plot for TPD pvals (Panel B)
def plot_real_shuffled_tpd_distributions(df, title="Real and Shuffled TPD Distributions"):
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=df["pval"],
        histnorm="probability density",
        name="Real Graph",
        xbins=dict(start=0, end=1, size=0.01),
        marker_color="#2769F4"
    ))

    fig.add_trace(go.Histogram(
        x=df["shufpval"],
        histnorm="",
        name="Shuffled Graph",
        xbins=dict(start=0, end=1, size=0.01),
        marker_color="#E9320B"
    ))

    fig.update_layout(
        title_text=title,  # title of plot
        xaxis_title_text='p-value',  # xaxis label
        yaxis_title_text='count',  # yaxis label
        bargap=0.2,  # gap between bars of adjacent location coordinates
        bargroupgap=0.1  # gap between bars of the same location coordinates
    )

    return fig


# system call to c++ program that calculates distance matrix
def generate_distance_matrix(graphmlfile, outfile="distmatrix.tsv", sofile="gothicAPSP.so", processors=1):
    # check that files have leading "./" so subprocess.run() will recognize them
    if graphmlfile[0] != ".":
        graphmlfile = "./" + graphmlfile
    if outfile[0] != ".":
        outfile = "./" + outfile
    if sofile[0] != ".":
        sofile = "./" + sofile
    # check if file is graphml or gt; if gt convert to graphml
    if graphmlfile.split(".")[-1] == "gt":
        graphmlname = str(".".join(graphmlfile.split(".")[:-1])) + ".graphml"  # so "."s in name don't cause truncated new name
        g = load_graph(graphmlfile)
        g.save(graphmlname)
        graphmlfile = graphmlname  # replaces argument var with new file
    try:
        cpp_process = subprocess.run([sofile, graphmlfile, outfile], capture_output=True, check=True)
        print("stdout:", cpp_process.stdout)
        print("stderr:", cpp_process.stderr)
    except subprocess.CalledProcessError:
        print("ERROR: Djikstra APSP c++ subprocess failed\n")
        print("stdout:", cpp_process.stdout)
        print("stderr:", cpp_process.stderr)

def count_edges(g):
    count = 0
    for edge in g.edges():
        count = count + 1
    return count


def count_vertices(g):
    count = 0
    for v in g.vertices():
        count = count + 1
    return count


# extracts chromosome from graphmlID for each node and applies it as a new vertex property map
def chromosome_annotate(g):
    vprop_chromosome = g.new_vertex_property("string")

    for v in g.vertices():
        node = int(v)
        vprop_chromosome[node] = g.vp.vname[v].split(":")[0]

    g.vp["chromosome"] = vprop_chromosome  # make vprop map internal to g so it will save with it


# annotate whether edges repressent intrachromosomal or interchromosomal contacts
def contact_annotate(g):
    eprop_contactType = g.new_edge_property("string")

    try:
        # annotate graph with chromosomes if not already
        if g.vertex_properties["chromosome"]:
            pass

    except KeyError:
        chromosome_annotate(g)

    # test and annotate edges
    for s, t in g.iter_edges():
        if g.vp.chromosome[s] == g.vp.chromosome[t]:
            eprop_contactType[g.edge(s, t)] = "intrachromosomal"
        else:
            eprop_contactType[g.edge(s, t)] = "interchromosomal"

    # make edge property internal
    g.ep["contactType"] = eprop_contactType  # make eprop map internal so it saves with graph


def create_top_edges_graph(g, threshold=0.05, outfile="top5percentGraph.gt"):
    print("NumEdges: " + str(count_edges(g)))
    # for all edges in graph, add weight to a list, sort the list, then calculate threshold for top%
    eWeightList = []
    for edge in g.edges():
        eWeightList.append(g.ep.weight[edge])

    eWeightList.sort(reverse=True)  # sorts in descending order
    cutoffIndex = int(len(eWeightList) * threshold)
    cutoffVal = eWeightList[cutoffIndex]
    print("sorted")

    # apply a boolean property map where prop=T is edge weight is in the top ?%
    TopPercentEdge = g.new_edge_property("bool")
    g.edge_properties["TopPercentEdge"] = TopPercentEdge
    for e in g.edges():
        if g.ep.weight[e] >= cutoffVal:
            g.ep.TopPercentEdge[e] = True
        else:
            g.ep.TopPercentEdge[e] = False

    # create graph view, and then save the graph
    g.set_edge_filter(TopPercentEdge)
    print("New NumEdges: " + str(count_edges(g)))
    g.save(outfile)


# annotates HiC graph with chromosomes, then creates a subgraph for each (only including intra-chrom contacts)
def split_into_chromosomes(graphfile):
    g = load_graph(graphfile)

    try:
        # annotate graph with chromosomes if not already
        if g.edge_properties["contactType"]:
            pass

    except KeyError:
        contact_annotate(g)

    # chromosome list for checking chroms
    chrlist = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
               'chr13',
               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    # create a subgraph for each chromosome, and save that subgraph
    for chrom in chrlist:
        # create new temporary edge and vertex properties to store whether v or e is correct chromosome
        vprop_isChromosome = g.new_vertex_property("bool")
        eprop_isChromosome = g.new_edge_property("bool")

        # loops to mark each edge and chromosome
        for v in g.vertices():
            if g.vp.chromosome[v] == chrom:
                vprop_isChromosome[v] = True
            else:
                vprop_isChromosome[v] = False

        for e in g.edges():
            s = e.source()
            if g.vp.chromosome[s] == chrom and g.ep.contactType[e] == "intrachromosomal":
                eprop_isChromosome[e] = True
            else:
                eprop_isChromosome[e] = False

        # set filters
        g.set_filters(eprop_isChromosome, vprop_isChromosome)

        # save graph and unset filters
        graphname = graphfile.split(".")[0] + "_" + str(chrom) + "." + graphfile.split(".")[1]
        g.save(graphname)


# converts bin file (.bed format, could be TADs) to dictionary where key is chr# and value is list of bin bounds
def bins_to_dict(binfile):

    bindict = {"chr1":[], "chr2":[], "chr3":[], "chr4":[], "chr5":[], "chr6":[], "chr7":[], "chr8":[], "chr9":[],
               "chr10":[], "chr11":[], "chr12":[], "chr13":[], "chr14":[], "chr15":[], "chr16":[], "chr17":[],
               "chr18":[], "chr19":[], "chr20":[], "chr21":[], "chr22":[], "chrX":[], "chrY":[], "chrM":[]}

    f = open(binfile, "r")
    for line in f:
        chrom, startpos, endpos = line.split()
        binbounds = str(startpos) + "-" + str(endpos)
        bindict[chrom].append(binbounds)

    return bindict


def generate_graph_report(graph, outfile="graphReport.txt"):
    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    chromosome_annotate(g)
    contact_annotate(g)

    # create dict to store all stats to be printed/ used later
    statsdict = {}

    # file size
    if type(graph) == str:
        statsdict["filesize"] = os.path.getsize(graph)

    # number of nodes
    nnodes = count_vertices(g)
    statsdict["nnodes"] = nnodes

    # number of edges
    statsdict["nedges"] = count_edges(g)

    # degree distribution
    ddistCounts, ddistBins = vertex_hist(g, "total")
    statsdict["ddistCounts"] = ddistCounts
    statsdict["ddistBins"] = ddistBins[1:]  # remove first item so this list repressents the upper range of each bin

    # betweenness distribution (too mem/time intensive for a quick report? needs to calculate APSP. Can use pivots to estimate)
    # extract a list of betweenness values from the generated Vertex Property Map so they can be histogrammed later
    sampleSize = int(nnodes / 10)  # set sample size to be 1/10th of the population
    # get list of unmasked vertices (we have to do this as a loop because graph tool is stupid
    realvlist = []
    try:
        for v in g.vertices():
            realvlist.append(v)
        pivotsList = random.sample(realvlist, sampleSize)  # randomly get a sample of nodes to act as pivots
        vBetweenness, eBetweenness = betweenness(g, pivots=pivotsList, weight=g.ep.weight)
        vBetweenList = []
        for key in vBetweenness:
            vBetweenList.append(float(vBetweenness[key]))
        statsdict["EstBetweennessDist"] = np.array(vBetweenList)  # array so it plays nicely with gnuplot (hopefully)
    except IndexError as ie:
        print(ie)
        print("Cannot get betweenness centrality. Continuing...")

    # diameter
    try:
        diameter, dpoints = pseudo_diameter(g, source=realvlist[0], weights=g.ep.weight)
        statsdict["diameter"] = diameter
    except IndexError as ie:
        print(ie)
        print("Cannot get pseudodiameter. Continuing...")

    # global clustering coefficient aka transitivity (weighted)
    try:
        transitivity, transStdDev = global_clustering(g, weight=g.ep.weight)
        statsdict["transitivity"] = transitivity
    except IndexError as ie:
        print(ie)
        print("Cannot get transitivity. Continuing...")

    # number of nodes in largest connected component
    try:
        gv = extract_largest_component(g)  # create a graph view of only the connected component
        lccCount = count_vertices(gv)
        statsdict["nnodesLargestConnectedComponent"] = lccCount
    except ValueError as ve:
        print(ve)
        print("Cannot retrieve largest component. Continuing...")

    # inter vs intra chromosomal contacts ratio
    interCount = 0
    intraCount = 0
    for e in g.edges():
        try:  # catch if graph hasn't been annotated with chromosome / contact info
            if g.ep.contactType[e] == "interchromosomal":
                interCount = interCount + 1
            elif g.ep.contactType[e] == "intrachromosomal":
                intraCount = intraCount + 1
        except:  # if anything goes wrong, try annotating before counting
            chromosome_annotate(g)
            contact_annotate(g)
            if g.ep.contactType[e] == "interchromosomal":
                interCount = interCount + 1
            elif g.ep.contactType[e] == "intrachromosomal":
                intraCount = intraCount + 1

    statsdict["interContactCount"] = interCount
    statsdict["intraContactCount"] = intraCount


    # Generate Report in .txt format for easy readability in terminal
    f = open(outfile, "w")  # open file for writing

    # create list of regular lines to becomethe body of the report
    linelist = []

    # draw ascii logo
    linelist.append("           _/         ,          .                           \n")
    linelist.append("       , -' )         (\-------.,')            (\________________________________________________,\n")
    linelist.append("     , ,-/ |          /\) )     \/            ,' _.--------------------------------------------, /\n")
    linelist.append("   ,',  /, |         /     >--. ,)           / /\\\        _____ _____ _____ _   _ _____ _____  '\n")
    linelist.append("  / ,  //|,'        /'    '\--'\\\)          /,'  \\\      |  __ \  _  |_   _| | | |_   _/  __ \ \n")
    linelist.append(" / ,  // ||       ,'    (.--^( `')         //     \\\     | |  \/ | | | | | | |_| | | | | /  \/\n")
    linelist.append("( ,  //  ||,___,-'    (__\\\  '^^^'        //       \\\    | | __| | | | | | |  _  | | | | |    \n")
    linelist.append(" \  //   ||--.__     (    \`^--)  _____.-'/         \\\   | |_\ \ \_/ / | | | | | |_| |_| \__/\ \n")
    linelist.append("  >'/    ||,        (      \|_(\-'      ,'           \\\   \____/\___/  \_/ \_| |_/\___/ \____/ \n")
    linelist.append(" /,'     ||          \          \      /              \\\         .. _.| \n")
    linelist.append("(/       ||,         \          )  ,'(     `     `    \\\,    /`\n\n")

    # draw real data
    # try/except blocks to ensure the report will run despite errors and will be as complete as possible
    now = datetime.datetime.now()
    try:
        if type(graph) == str:
            linelist.append("Report generated for " + graph + " on " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
        else:
            linelist.append("Report generated on " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    except KeyError:
        print("Did not append datetime information because it is missing")
    try:
        linelist.append("Filesize: " + str(statsdict["filesize"] / 1000000000) + "Gb\n")
    except KeyError:
        print("Did not append filesize information because it is missing")
    try:
        linelist.append("Number of Nodes: " + str(statsdict["nnodes"]) + "\n")
    except KeyError:
        print("Did not append node count information because it is missing")
    try:
        linelist.append("Number of Edges: " + str(statsdict["nedges"]) + "\n")
    except KeyError:
        print("Did not append edge count information because it is missing")
    try:
        linelist.append("Number of Nodes in Largest Connected Component: " + str(statsdict["nnodesLargestConnectedComponent"]) + "\n")
    except KeyError:
        print("Did not append largest component information because it is missing")
    try:
        linelist.append("Number of Intrachromosomal Contacts: " + str(statsdict["intraContactCount"]) + "\n")
    except KeyError:
        print("Did not append intrachromosomal contact information because it is missing")
    try:
        linelist.append("Number of Interchromosomal Contacts: " + str(statsdict["interContactCount"]) + "\n")
    except KeyError:
        print("Did not append interchromosomal contact information because it is missing")
    try:
        linelist.append("Pseudodiameter: " + str(statsdict["diameter"]) + "\n")
    except KeyError:
        print("Did not append pseudodiameter information because it is missing")
    try:
        linelist.append("Global Clustering Coefficient (Weighted Transitivity): " + str(statsdict["transitivity"]) + "\n")
    except KeyError:
        print("Did not append transitivity information because it is missing")
    for line in linelist:  # write lines to file
        f.write(line)

    sys.stdout = f

    # create plots that will be drawn at the bottom of the report
    # plot degree distribution
    try:
        x = statsdict["ddistBins"]
        y = statsdict["ddistCounts"]
        ymax = str(max(y)+5)
        print("\n\n_________________________________ Degree Distribution Plot ________________________________")  # in case title doesn't work
        print("\ny = Vertex Count")  # Since gnuplot doesn't want to print the ylabel in terminal
        gp.plot((x, y, {'with': 'boxes'}),
                     _with    = 'lines',
                     terminal = 'dumb 80,80',
                     unset    = 'grid',
                     title    = 'Degree Distribution',
                     xlabel   = 'Degree (Upper Bin Edge)',
                     ylabel   = 'Count',
                     square = True
                )
    except KeyError:
        print("Could not draw plots due to missing variable")

    # # plot estimated betweenness distribution (old and broken)
    # betwnDat = statsdict["EstBetweennessDist"]
    # betwnCnts, betwnBins = np.histogram(betwnDat, bins=67)  # 67 just because that' how big the above plot is
    # betwnBins = betwnBins[1:]
    # print("\n\n___________________ Estimated Betweenness Distribution Plot (n = 10000) ___________________")  # in case title doesn't work
    # print("\ny = Edge Count")  # Since gnuplot doesn't want to print the ylabel in terminal
    # gp.plot( (betwnBins, betwnCnts, {'with': 'boxes'}),
    #              _with    = 'lines',
    #              terminal = 'dumb 80,40',
    #              unset    = 'grid',
    #              title    = 'Estimated Betweenness Distribution',
    #              xlabel   = 'Betweenness (Upper Bin Edge)',
    #              ylabel   = 'Count'
    #         )
    f.close()
    sys.stdout.close()
    sys.stdout = sys.__stdout__


# function to print list of Top GO terms with FDR and p-val cutoffs for term inclusion
def get_top_goterms(resultscsv, outfile="topterms.csv"):
    # read in and format results csv from getgotpdvals()
    df = pd.read_csv(resultscsv)
    df = df.transpose()
    df.columns = df.iloc[0]
    df["goterm"] = df.index
    df = df[1:]
    df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                    "shufpval": "float", "goterm": "str"})

    sorteddf = df.sort_values("pval")
    shufsorteddf = df.sort_values("shufpval")

    golist = []
    plist = []
    fdrlist = []

    for index, row in sorteddf.iterrows():
        thisgo = row["goterm"]
        thispval = row["pval"]

        # get count in list of pval and shufpval <= pval
        pcount = sum(i <= thispval for i in sorteddf["pval"])
        shufcount = sum(i <= thispval for i in shufsorteddf["shufpval"])  # don't even think I needed to sort for this lol
        denom = pcount + shufcount # denominator for fdr calc TODO make sure this is right (not pcount/shufcount?)
        fdr = pcount / denom

        golist.append(thisgo)
        plist.append(thispval)
        fdrlist.append(fdr)

    # create new df with cols: "goterm, pval, FDR"
    newdf = pd.DataFrame({"goterm": golist, "pval": plist, "FDRatThisCutoff": fdrlist})

    # save df
    print(newdf)
    newdf.to_csv(outfile)


# take results file and return table with cols: pvalThreshold, numRealPassing, numShufPassing, fdr
# applies a monotonic traqnsformation to the FDR such that the FDR never decreases
# ^this is valid because if FDR decreases you can always just pick a higher pval threshold that gave better fdr
# important that startstopstep always decreases for monotonic transformation to work
def get_real_v_null_significants(resultscsv, startstopstep=(0.1, 0, -0.00001),
                                 outfile="SignificantClustersRealVsShuffledTPD.csv", monotonic=True):

    # read in and format results csv from getgotpdvals()
    df = pd.read_csv(resultscsv)
    df = df.transpose()
    df.columns = df.iloc[0]
    df["goterm"] = df.index
    df = df[1:]
    df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                    "shufpval": "float", "goterm": "str"})

    lowestfdr = 1  # initialize lowest fdr as 1, which will be checked against every new fdr value

    pstart, pstop, pstep = startstopstep

    with open(outfile, "w") as f:
        f.write("pvalThreshold, numReal, numShuf, fdr\n")

        for i in np.arange(pstart, pstop, pstep):

            # get count in list of pval and shufpval <= pval
            pcount = sum(j <= i for j in df["pval"])
            shufpcount = sum(j <= i for j in df["shufpval"])
            denom = pcount + shufpcount  # denominator for fdr calc
            fdr = shufpcount / denom

            if monotonic:  # check if monotonic transformation option is selected (True by default)
                if fdr < lowestfdr:  # always ensures fdr that gets written is the lowest possible
                    lowestfdr = fdr

            else:
                lowestfdr = fdr  # if monotonic False, just write new FDR every time

            # write line to file
            f.write(str(i) + ", " + str(pcount) + ", " + str(shufpcount) + ", " + str(lowestfdr) + "\n")


# plots the number of terms passing significance threshold at each FDR value
# input should be csv from get_real_v_null_significants
def plot_passing_terms_v_fdr_fig(inputcsv, outfile="numtermsVfdrFig.png"):

    df = pd.read_csv(inputcsv, sep=", ")  # read input csv into data frame
    df = df.astype({"pvalThreshold": "float", "numReal": "int", "numShuf": "int", "fdr": "float"})
    newdflist = []
    # filter dataframe so there is only one row per FDR value
    currentfdr = 0
    for index, row in df.iterrows():
        if row["fdr"] != currentfdr:  # new fdr at row means we are at the new best threshold (most real hits / fdr)
            newdflist.append({"pvalThreshold": row["pvalThreshold"], "numReal": row["numReal"], "numShuf": row["numShuf"],
                              "fdr": row["fdr"]})
            currentfdr = row["fdr"]  # update to new fdr for row

    newdf = pd.DataFrame.from_records(newdflist)  # convert from list of dicts to df

    # make plot using new filtered df
    plt.figure(figsize=(8, 6))  # Optional: Adjust the figure size
    plt.scatter(newdf["fdr"], newdf["numReal"], marker='o', color='b', alpha=0.5)
    #plt.title("Number of real significant GO terms versus False Discovery Rate")
    plt.xlabel("False Discovery Rate")
    plt.xlim(0, 0.25)
    plt.ylim(0, 150)
    plt.ylabel("Number of statistically significantly clustered GO terms")
    plt.savefig(outfile)


# takes a graph and outputs a list of node ids annotated with a given go term
def get_nodes_with_term(graphfile, term):
    g = load_graph(graphfile)

    nodelist = []
    for v in g.vertices():
        for t in g.vp.goterms[v]:
            if t == term:
                nodelist.append("n" + str(int(v)))

    return nodelist


# calculates A/B compartments from .hic file and annotates graph
def annotate_ab_compartments(graphfile, hicfile, outfile, genomefile):
    hic = fanc.load(hicfile)  # load hic file into mem
    ab = fanc.ABCompartmentMatrix.from_hic(hic)  # get
    ev = ab.eigenvector()
    domains = ab.domains()

    g = load_graph(graphfile)

    for region in domains.regions:
        # extract region information
        chr = region.chromosome
        start = region.start
        end = region.end
        compartment = region.name

        print(str(chr) + ":" + str(start) + "-" + str(end) + " -> " + str(compartment))


# opens tad file (.bed format) and creates a {[chrom]:[list of region bounds], ...} dictionary
def tads2dict(tadfile):

    taddict = {"chr1":[], "chr2":[], "chr3":[], "chr4":[], "chr5":[], "chr6":[], "chr7":[], "chr8":[], "chr9":[],
               "chr10":[], "chr11":[], "chr12":[], "chr13":[], "chr14":[], "chr15":[], "chr16":[], "chr17":[],
               "chr18":[], "chr19":[], "chr20":[], "chr21":[], "chr22":[], "chrX":[], "chrY":[], "chrM":[]}

    f = open(tadfile, "r")
    for line in f:
        chrom, startpos, endpos = line.split()
        binbounds = str(startpos) + "-" + str(endpos)
        taddict[chrom].append(binbounds)

    return taddict

# plots correlation of TPD vs # of interchromosomal contacts for each go term or other annotation
def plot_interchromosomal_contact_tpd_correlation(resultscsv, graph, outcsv=None, outplot=None, prop="goterms"):
    # load results df to get tpds
    df = pd.read_csv(resultscsv)
    df = df.transpose()
    df.columns = df.iloc[0]
    df["goterm"] = df.index
    df = df[1:]
    df = df.astype({"nnodes": "int", "tpd": "float", "pval": "float", "shuftpd": "float",
                    "shufpval": "float", "goterm": "str"})

    # make godict to get list of nodes for each go term
    gd = make_go_dict(graph, prop=prop)  # TODO put after loading graph once I fix graph loading in all functions

    # load graph
    if type(graph) == str:
        g = load_graph(graph)
    elif type(graph) == Graph:
        g = graph
    else:
        print("bad argument: Graph should be a file name (<class 'str'>) or graph-tool graph object (<class ‘graph_tool.Graph‘>)")

    # annotate edges with contact info
    contact_annotate(g)

    # initialize data structures
    contactsdfcol = []
    intercontactsdfcol = []
    contactratiodfcol = []


    # for each GO term, for each node, count number of interchromosomal contacts (as a function of all contacts?)
    # also fetch TPD from results df for go term while in this loop
    for index, row in df.iterrows():
        goterm = row["goterm"]
        print(goterm)
        tpd = row["tpd"]
        intercontactcounter = 0
        contactcounter = 0

        # for each node in list from godict, count intercontact edges and total edges
        for node in gd[goterm]:  # iterate through list of nodes with goterm
            v = g.vertex(node)
            for e in v.out_edges():
                if g.ep.contactType[e] == "interchromosomal":
                    intercontactcounter = intercontactcounter + 1
                contactcounter = contactcounter + 1
            contactratio = intercontactcounter / contactcounter
            print("contacts:" + str(contactcounter) + " intercontacts:" +
                  str(intercontactcounter) + " ratio:" + str(contactratio))

        # add stats to respective lists to later append to df
        contactsdfcol.append(contactcounter)
        intercontactsdfcol.append(intercontactcounter)
        contactratiodfcol.append(contactratio)

    # add to df
    df["totalContacts"] = contactsdfcol
    df["interchromosomalContacts"] = intercontactsdfcol
    df["contactRatio"] = contactratiodfcol

    # write to csv
    if outcsv is not None:
        df.to_csv(outcsv)

    # plot
    if outplot is not None:
        sns.set_theme(style="ticks")
        corplot = sns.lmplot(data=df, x="tpd", y="contactRatio",
        palette="muted", ci=None,
        height=4, scatter_kws={"s": 5, "alpha": 1})

    corplot.savefig(outplot)


# full run of the pipeline, from sam files to final results dataframe plus figures
# files that need to be in filedir:
# - map(HUMAN_9606_idmapping_selected.tsv)
# - gencode(gencode_pcgenes.csv)
# - go obo(go-basic.obo)
# - GothicAPSP.so
# - Ensembl2Reactome_All_Levels.txt
def gothic_full_run(runname, sam1=None, sam2=None, filedir=".", binsize=80000, vmax=100, degreecutoff=5000, ncores=1,
                    numshuffles=3000000, step=0, endstep=12, binfile=None, saveadjlist=False):
    start = timeit.default_timer()

    if binfile:
        bindict = bins_to_dict(binfile)
        print("Bin dict created from file...")
    else:
        bindict = None

    # create graph
    if step <= 1 <= endstep:
        print("Creating Graph...")
        adjlist = sam_adjacencies(sam1, sam2, m=binsize, bins=bindict, verbose=True)
        print("Adjacency list created...")
        if saveadjlist:  # write adjacency list to file if saveadjlist option is true
            print("Writing adjacency list to file...")
            adjlistoutfile = runname + "_ADJLIST.tsv"
            save_adjlist(adjlist, adjlistoutfile)

    # convert adjlist to graphml
    if step <= 2 <= endstep:
        print("Converting adjlist to .graphml...")
        graphname = runname + "_noannotation.graphml"
        if saveadjlist:
            adjlist = []  # lil hack to overwrite the adjlist or else create one if there is none (ie start w this step)
            del adjlist  # clear up memory from filled adjlist dictionary object
            adjlistoutfile = runname + "_ADJLIST.tsv"
            adjlist_file_to_graphml(adjlistoutfile, outfile=graphname)
        else:
            hic_adjlist_to_graphml(adjlist, fileprefix=graphname)
            del adjlist  # clear up memory from filled adjlist dictionary object
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(GraphStep): ', stop - start)  # Timer Block

    # annotate graph
    if step <= 3 <= endstep:
        print("Annotating graph...")
        graphname = runname + "_noannotation.graphml"
        mymap = filedir + "/HUMAN_9606_idmapping_selected.tsv"
        gencode = filedir + "/gencode_pcgenes.csv"
        goobo = filedir + "/go-basic.obo"
        reactomemap = filedir + "/Ensembl2Reactome_All_Levels.txt"
        newgraphname = runname + "_incompleteannotated.graphml"
        genes_go_annotate(graphname, mapfile=mymap, gencodefile=gencode, m=binsize,
                          binfile=binfile, outfile=newgraphname, go_obo=goobo)
        oldgraphname = newgraphname
        newgraphname = runname + "_annotated.graphml"
        reactome_annotate(oldgraphname, mapfile=reactomemap, outfile=newgraphname)

        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(AnnotateStep): ', stop - start)  # Timer Block

    # normalize, invert, trim by degree
    if step <= 4 <= endstep:
        print("Normalizing...")
        graphname = runname + "_incompleteannotated.graphml"
        newgraphname = runname + "_normalized.graphml"
        # skbalance(graphname, outfile=newgraphname)  # deprecated
        krres = str(binsize/1000) + "kb"  # resolution string for gcmapexplorer
        kr_balance(graphname, krres, outfile=newgraphname)
        print("Inverting...")
        graphname = newgraphname
        newgraphname = runname + "_inverted.graphml"
        invert_weights(graphname, outfile=newgraphname)
        # TODO reintroduce this block if we decide we still want to trim high degree nodes (edge trimmming should be enough)
        # print("Trimming...")
        # graphname = newgraphname
        # newgraphname = runname + "_trimmed.graphml"
        # trim_nodes_by_degree(graphname, newgraphname, degreecutoff=degreecutoff)
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(NormalizeStep): ', stop - start)  # Timer Block

    # calculate All Pairs Shortest Paths (APSP) and create distance matrix
    if step <= 5 <= endstep:
        newgraphname = runname + "_inverted.graphml"
        distgraphname = runname + "_wdistances.graphml"
        annotate_apsp(newgraphname, distgraphname)  # calculate all pairs shortest paths as graph annotations
        matfilename = runname + "_distmatrix.tsv"
        distmatrix_from_annotations(distgraphname, outfile=matfilename)  # convert/format distances to distance matrix
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(DijkstraStep): ', stop - start)  # Timer Block

    # convert from .graphml to .gt
    if step <= 6 <= endstep:
        print("Converting to .gt...")
        graphname = runname + "_inverted.graphml"
        newgraphname = runname + ".gt"  # naming convention assumes basic graph name is the complete one TODO is that ok?
        graphml_to_gt(graphname, newgraphname)

    # generate graph report
    if step <= 7 <= endstep:
        graphname = runname + ".gt"
        print("Generating graph report...")
        generate_graph_report(graphname, outfile=runname + "graphReport.txt")

    # get sampling vector and perform MC sampling
    if step <= 8 <= endstep:
        graphname = runname + ".gt"
        print("Starting MC sampling...")
        matfilename = runname + "_distmatrix.tsv"
        vlens = get_vlengths(graphname)
        distrnfile = runname + "_distributions.csv"
        montecarlo_sample_tpds(matfilename, vlens, graphname, ncores=ncores, outfile=distrnfile)
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(MCsamplingStep): ', stop - start)  # Timer Block

    # create shuffled graph
    if step <= 9 <= endstep:
        graphname = runname + ".gt"
        newgraphname = runname + "_goshuffled.gt"
        label_swap(graphname, outfile=newgraphname, m=numshuffles)
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(ShuffleStep): ', stop - start)  # Timer Block

    # get real and shuffled TPDs for all go terms
    if step <= 10 <= endstep:
        graphname = runname + ".gt"
        newgraphname = runname + "_goshuffled.gt"
        matfilename = runname + "_distmatrix.tsv"
        gd = make_go_dict(graphname)  # need to make godict before doing all_go_tpd()
        shufgd = make_go_dict(newgraphname)  # need to make godict before doing all_go_tpd()
        all_go_tpd(gd, matfilename, runname + "_realTPDs.csv")
        all_go_tpd(shufgd, matfilename, runname + "_shuffledTPDs.csv")
        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(TPDstep): ', stop - start)  # Timer Block

    # get final dataframe with "nnodes", "tpd", "pval", "shuftpd", "shufpval" for all go terms
    if step <= 11 <= endstep:
        graphname = runname + ".gt"
        distrnfile = runname + "_distributions.csv"
        resultscsvname = runname + "_resultsDF.csv"
        topgolistname = runname + "_topGOterms.csv"

        df = get_go_tpd_pvals(graphname, runname + "_realTPDs.csv", runname + "_shuffledTPDs.csv", distrnfile)
        df.to_csv(resultscsvname)

        # also get list of top go terms with FDR
        get_top_goterms(resultscsvname, outfile=topgolistname)

        # plots
        pvalhistofig = plot_real_shuffled_tpd_distributions(df)
        pvalhistofig.write_image(runname + "_pvalHistogram.png")
        fdrfig = plot_fdr_pval_histogram(df)
        fdrfig.write_image(runname + "_FDRpvalComparePlot.png")

        stop = timeit.default_timer()  # Timer Block
        print('Global Runtime(EndStep): ', stop - start)  # Timer Block
