# Used to process batch download from plasmoDB to be input into GeneTargeter

from py.GenBankToolbox import *

def processFastas(fastas):
    genes = [];
    for f in fastas:
        info = f.split("|");
        geneName = (info[0].strip()+"_"+info[2].strip()).replace("/",",").replace(" ","_");
        genIndexes = [1000,len(fastas[f])-1000]
        geneSeq = fastas[f][genIndexes[0]:genIndexes[1]];
        geneAnn = GenBankAnn(geneName, "misc", geneSeq, False, genIndexes);
        gene = GenBank();
        gene.name = geneName;
        gene.origin = fastas[f];
        gene.features.append(geneAnn);

        gene.save("gb/"+geneName+".gb",True);
