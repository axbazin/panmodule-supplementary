#!/usr/bin/env python3
#coding: utf-8

# Author: Adelme Bazin

#default libraries
from collections import defaultdict, Counter
import argparse
from itertools import combinations
from math import sqrt
import os
import glob

#installed libraries
from tqdm import tqdm

class Module:
    def __init__(self, modname, start, stop, category=None, subclass=None):
        self.name = modname
        self.start = int(start)
        self.stop = int(stop)
        self.category = category
        self.subclass=subclass

class gene:
    def __init__(self, genome, start, stop, mods, nb=0):
        self.genome = genome
        self.start = int(start)
        self.stop = int(stop)
        self.modules = set(mods.split(',')) if mods != "" else set()
        self.number_of_genes = int(nb)

def writeMods_CGVIEW(genome2mod):
    for genome, mods in genome2mod.items():
        f = open("CGVIEW/" + genome + "_cgview.tab","w")
        f.write("name\ttype\tstart\tstop\tstrand\n")
        for mod in mods:
            f.write('\t'.join([mod.category, mod.name,mod.start, mod.stop, "+"]) + "\n")
        f.close()

def genomeAndid():
    ID_matchings = ["GCF_000013305.1\t536\tNC_008253.1","GCF_000014845.1\tAPEC\tNC_008563.1","GCF_000007445.1\tCFT073\tNC_004431.1","GCF_000026305.1\tED1a\tNC_011745.1","GCF_000017765.1\tHS\tNC_009800.1","GCF_000026265.1\tIAI1\tNC_011741.1","GCF_000026345.1\tIAI39\tNC_011750.1","GCF_000005845.2\tK12\tNC_000913.3","GCF_000026285.1\tS88\tNC_011742.1","GCF_000010385.1\tSE11\tNC_011415.1","GCF_000026325.1\tUMN026\tNC_011751.1","GCF_000013265.1\tUTI89\tNC_007946.1"]
    genome2id = {}
    id2genome = {}
    for line in ID_matchings:
        line = line.strip().split("\t")
        genome2id[line[1]] = line[2]
        id2genome[line[2]] = line[1]
    return genome2id, id2genome

def genome2mods(references, unknown = True ):

    f = open(references,"r")
    f.readline()
    genome2mod=defaultdict(set)
    for line in f:
        line = line.split("\t")
        if line[4] == "unknown" and not unknown:
            continue#if unknown is not kept, we do not store the module
        genome2mod[line[0]].add(Module(modname=line[1], start=line[2], stop=line[3], category=line[4], subclass=line[5]))
    f.close()
    return genome2mod

def readGenomicIslands(genomeids, gifile):
    genome2gi = defaultdict(set)
    infile = open(gifile, "r")
    for line in infile:
        line = line.strip().split("\t")
        if line[2] in genomeids:
            genome2gi[line[2]].add(range(int(line[3]), int(line[4])))#start, stop
    if len(genome2gi) != len(genomeids):
        raise Exception(f"Found  {len(genome2gi)} reference genomes only, while we expected {len(genomeids)}")
    return genome2gi

def readModules(genomeids, modfile):
    genome2genes = defaultdict(set)
    infile = open(modfile, "r")
    for line in infile:
        line = line.split("\t")
        if line[1] in genomeids:
            genome2genes[line[1]].add(gene(genome=line[1], start=line[3], stop=line[4], mods=line[5].strip()))
    if len(genome2genes) != len(genomeids):
        raise Exception(f"Found  {len(genome2genes)} reference genomes only, while we expected {len(genomeids)}")
    #print(f"found {len(genome2genes)} reference genomes among {len(genomeids)} for the modules files")
    return genome2genes

def filterRefMods(ref, gis):
    kept_refs = defaultdict(dict)
    
    for genome, refmodules in ref.items():
        kept_nb = 0
        for refmod in refmodules:
            for gi in gis[genome]:
                if refmod.start in gi or refmod.stop in gi:
                    kept_nb+=1
                    try:
                        kept_refs[genome][gi].add(refmod)
                    except KeyError:
                        kept_refs[genome][gi] = set([refmod])
        #print(f"Kept {kept_nb} reference modules among {len(ref[genome])} for {genome}.")
    return kept_refs

def getGIgenes(gis, modules):
    kept_genes = defaultdict(dict)
    
    for genome, genelist in modules.items():
        kept_nb=0
        tot_genes = set()
        for gene in genelist:
            for gi in gis[genome]:
                if gene.start in gi and gene.stop in gi:
                    kept_nb+=1
                    tot_genes.add(gene)
                    try:
                        kept_genes[genome][gi].add(gene)
                    except KeyError:
                        kept_genes[genome][gi] = set([gene])

    return kept_genes

def write_metrics_occurrences(met_occ):
    output = open("metrics_occurrences.tsv","w")
    for metric, counter in met_occ.items():
        for nb_genes, value in counter.items():
            output.write(f"{metric}\t{nb_genes}\t{value}\n")
    output.close()

def calculateMetrics(refs, gi_genes):

    FN_desc = Counter()

    metrics_occurrences = {"TP":Counter(), "FN":Counter(), "FP": Counter(), "TN":Counter()}

    TP = 0
    FN = 0
    FP = 0
    TN = 0
    for genome, gis in gi_genes.items():
        for gi, genes in gis.items():
            try:
                curr_refs = refs[genome][gi]
            except KeyError:
                continue
            gene2ref = defaultdict(set)
            curr_regions_mods = {}
            for ref in curr_refs:
                for gene in genes:
                    ## associate genes to a reference module, if they are in a reference's positions
                    if gene.start in range(ref.start, ref.stop) and gene.stop in range(ref.start, ref.stop):
                        gene2ref[gene].add(ref)
                    elif gene.start in range(ref.start, ref.stop):
                        if (gene.stop - gene.start)*0.8 <= ref.stop - gene.start:#if more than 80% of the gene is in the ref, add it to ref
                            gene2ref[gene].add(ref)
                    elif gene.stop in range(ref.start, ref.stop):
                        if (gene.stop - gene.start)*0.8 <=  gene.stop - ref.start:#if more than 80% of the gene is in the ref, add it to ref
                            gene2ref[gene].add(ref)

                    ## get predicted modules position in the current GI
                    for mod in gene.modules:
                        if mod not in curr_regions_mods:
                            curr_regions_mods[mod] = [gene.start, gene.stop]
                        else:
                            if curr_regions_mods[mod][0] > gene.start:
                                curr_regions_mods[mod][0] = gene.start
                            if curr_regions_mods[mod][1] < gene.stop:
                                curr_regions_mods[mod][1] = gene.stop

            #all genes in refs
            #compute metrics for all pairs of gene in the GI
            
            for gene1, gene2 in combinations(gene2ref.keys(), 2):
                if len(gene2ref[gene1] & gene2ref[gene2]) > 0:
                    if len(gene1.modules & gene2.modules) > 0:
                        TP+=1
                        metrics_occurrences["TP"][min([gene1.number_of_genes, gene2.number_of_genes])] += 1
                    else:
                        FN+=1
                        metrics_occurrences["FN"][min([gene1.number_of_genes, gene2.number_of_genes])] += 1
                        if len(gene1.modules) == 0 and len(gene2.modules) == 0:
                            FN_desc["both_not_in_mods"] +=1
                        elif len(gene1.modules) == 0 or len(gene2.modules) == 0 :
                            FN_desc["one_not_in_mod"] +=1

                            curr_inside_mod_pos = False
                            if len(gene1.modules) == 0:
                                #gene2 has module(s)
                                checked_gene = gene1
                                checked_mods = gene2.modules
                            else:
                                checked_gene = gene2
                                checked_mods = gene1.modules
                            
                            for mod in checked_mods:
                                if curr_regions_mods[mod][0] < checked_gene.start and curr_regions_mods[mod][1] > checked_gene.stop:
                                    curr_inside_mod_pos = True
                                    break

                            if curr_inside_mod_pos:
                                FN_desc["inside_mod_position"]+=1
                        elif len(gene1.modules) > 0 and len(gene2.modules) > 0:
                            FN_desc["in_different_mod"] +=1
                            intersection = False
                            for mod1 in gene1.modules:
                                if intersection:
                                    break
                                for mod2 in gene2.modules:
                                    if (curr_regions_mods[mod1][0] > curr_regions_mods[mod2][0] and curr_regions_mods[mod1][0] < curr_regions_mods[mod2][1]) or (curr_regions_mods[mod2][0] > curr_regions_mods[mod1][0] and curr_regions_mods[mod2][0]<curr_regions_mods[mod1][1]):
                                        intersection=True
                                        break
                            if intersection:
                                FN_desc["mods_overlap"]+=1
                else:
                    if len(gene1.modules & gene2.modules) > 0:
                        FP+=1
                        metrics_occurrences["FP"][min([gene1.number_of_genes, gene2.number_of_genes])] += 1
                    else:
                        TN+=1
                        metrics_occurrences["TN"][min([gene1.number_of_genes, gene2.number_of_genes])] += 1

    #write_metrics_occurrences(metrics_occurrences)

    return TP, FN, FP, TN, FN_desc


def computeBench(ref, gis, modules):
    ## remove reference modules that are not in gis
    reference = filterRefMods(ref, gis)

    gi_genes = getGIgenes(gis, modules)

    TP, FN, FP, TN, FN_desc = calculateMetrics(reference, gi_genes)

    recall = TP / (TP + FN)
    precision = TP / (TP + FP) if (TP+FP) > 0 else 0
    accuracy = (TP + TN) / (TP + FN + FP + TN)
    F1scores =  2*TP / (2*TP + FN + FP)
    curr_MCC = (( TP * TN ) - (FP * FN )) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0 else 0

    results = {"TP":str(TP), "FN":str(FN), "FP":str(FP), "TN":str(TN), "MCC":str(curr_MCC), "recall":str(recall), "precision":str(precision), "accuracy":str(accuracy), "F1":str(F1scores)}
    results["both_not_in_mods"] = str(FN_desc["both_not_in_mods"])
    results["one_not_in_mod"] = str(FN_desc["one_not_in_mod"])
    results["inside_mod_position"] = str(FN_desc["inside_mod_position"])
    results["mods_overlap"] = str(FN_desc["mods_overlap"])
    #print(TP, FN, FP, TN, FN_desc)
    return results

def cmdLine():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--genome_islands', type=str, required=True, help="tsv file describing genomic island positions. Expected format is the one used for 'plastic_regions.tsv' file in PPanGGOLiN. The used columns (in 'normal' indenting) are 2nd (organism genbank id), 4th (GI start position) and 5th (GI stop position).")
    parser.add_argument("-r","--references",type=str, required=True, help = "The reference to use to compute the benchmark")
    parser.add_argument("-m", '--modules', type=str, required=True, help="directory containing predicted modules tsv files. Expected file naming is 'softname_parameter1-val1[_parameterN-valN_][...].tsv. One line per genes, all genes of all genomes are expected. Expected format is, in tsv:\n genbank id, gene id, start, stop, comma-separated modules list if any.")
    parser.add_argument("-o", "--output", type=str, required=True, help="output file name for result files")
    parser.add_argument("-u","--unknown", action='store_false', help="Use this option to remove the 'unknown' functional category from the bench")
    return parser.parse_args()


def main():
    args = cmdLine()
    genome2mod = genome2mods(args.references, args.unknown)
    genome2id , id2genome = genomeAndid()
    id2mod = {}
    for genome, mods in genome2mod.items():
        id2mod[genome2id[genome]] = mods
    # get the genomic island positions
    gis = readGenomicIslands(id2genome.keys(), args.genome_islands)
    # get the modules for each genome.

    fout = open(args.output,"w")
    c = 0
    ##glob files for each module file:
    print("done loading data. Computing bench for all files...")
    for fname in tqdm(glob.glob(args.modules + '/*'), unit="file"):
        modules = readModules(id2genome.keys(), fname)
        # compare for each genome one by one.

        results = computeBench(id2mod, gis, modules)

        file_dependent_stuff = '.'.join(os.path.basename(fname).split('.')[:-1]).split("_")

        for element in file_dependent_stuff[1:]:
            key, val = element.split("-")
            results[key] = val

        if c == 0:
            fout.write("software\t" + "\t".join(results.keys())+"\n")
        fout.write(file_dependent_stuff[0] + "\t" + "\t".join(results.values())+"\n")
        fout.flush()
        c+=1
    fout.close()

if __name__ == "__main__":
    main()