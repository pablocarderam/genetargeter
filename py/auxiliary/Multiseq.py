from __future__ import print_function
# Used to process batch download from plasmoDB to be input into GeneTargeter

from builtins import str
from py.utils.GenBankToolbox import *
from py.utils.BioUtils import *
from copy import deepcopy

# use to process plasmoDB fasta download (no introns, use gff_to_genbank now
# instead!)
def processFastas(fastaFilepath):
    fastas = loadFastas(fastaFilepath);
    genes = [];
    for f in fastas:
        info = f.split("|");
        geneName = (info[0].strip()+"_"+info[2].strip()).replace("/",",").replace(" ","_");
        genIndexes = [1000,len(fastas[f])-1000]
        geneSeq = fastas[f][genIndexes[0]:genIndexes[1]];
        geneAnn = GenBankAnn(geneName, "gene", geneSeq, False, genIndexes);
        gene = GenBank();
        gene.name = geneName;
        gene.origin = fastas[f];
        gene.features.append(geneAnn);

        gene.save("gb/"+geneName+".gb",True);

        if "WARNING:" in f:
            print("WARNING in seq " + geneName);


# Used to process gff_to_genbank output
def processGenBank(gbFilePath):
    txt = open(gbFilePath); # Access given file
    d = txt.read(); # Read file

    gbTxts = d.split("//"); # Split each sequence (only first occurrence of \n) into lines: first line contains sequence name. Returns d as a list of lists (sequences) containing strings (seq name, seq).
    count = 0;
    total = d.count("gene");

    for gbTxt in gbTxts:
        if len(gbTxt.strip()) > 0:
            gb = GenBank();
            gb.load(gbTxt+"//",loadFromFile=False);

            for annOriginal in gb.features:
                annGB = deepcopy(gb);
                annsLabel = annGB.findAnnsLabel(annOriginal.label)
                ann = GenBankAnn()
                for a in annsLabel:
                    if a.label == annOriginal.label:
                        ann = a
                        break

                if len(ann.label) == 0:
                    print("ERROR: no annotation found")

                if ann.type == "gene":
                    annGB.removeSeq([0,max(ann.index[0]-1000,0)]);
                    annGB.removeSeq([min(ann.index[1]+1000,len(annGB.origin)-1),len(annGB.origin)]);
                    annGB.name = ann.label;
                    annGB.save("geneFiles/"+ann.label+".gb",saveToFile=True);
                    count+=1;
                    print([str(count)+"/"+str(total)]);
