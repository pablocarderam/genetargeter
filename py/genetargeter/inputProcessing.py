
import os; # needed for file handling
import copy; # Import object copying methods for deep copies

from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants

"""
Accepts GenBank string possibly with multiple mRNA splice forms.
Returns dictionary of GenBank objects with only one mRNA and its exons/introns
annotated. Keys are transcript names.
"""
def preprocessInputFile(geneName, geneFileStr, useFileStrs=False):
    gbDict = {}; # will store GenBank objects to be returned

    geneGB = GenBank(); # initializes variable to hold gene GenBank object
    if useFileStrs: # if we are using file strings,
        geneGB.load(geneFileStr, loadFromFile=False); # load gene from file string
        path = ""; # sets an empty path, needed for save functions of class GenBank
    else: # If we're using files,
        geneGB.load(geneFileStr, loadFromFile=True); # load gene file
        path = "output/" + geneName; # string with path to this gene's directory in output folder
        if not os.path.exists(path): # if a directory for this gene does not exist already,
            os.mkdir(path); # makes new directory


    geneGB.setAllColors(annColors['otherAnnColor']); # set all annotation colors to grey

    geneList = geneGB.findAnnsLabel(geneName); # search for annotations with gene name
    if len(geneList) > 0: # if genes found,
        gene = geneList[0]; # stores gene annotation
        geneAnns = geneGB.findAnnsLabel(geneName); # stores all gene annotations with gene name in them
        for a in geneAnns: # loops through found annotations
            if a.type == "gene": # if this annotation's type is gene,
                gene = a; # keep it as the gene annotation
                break; # stop loop

        gene.color = annColors['targetGeneColor'] # set gene to orange, using Okabe and Ito's colorblind palette from 2002
        mRNAs = geneGB.findAnnsType("mRNA",True)+geneGB.findAnnsType("tRNA",True)+geneGB.findAnnsType("rRNA",True); # list of all mRNAs in GB file. EDIT 11 oct 2017: Include rRNA and tRNA, not just mRNAs!
        for mRNA in mRNAs: # loop over all mRNAs
            if gene.index[0] <= mRNA.index[0] < mRNA.index[1] <= gene.index[1]: # if mRNA is inside gene,
                newGB = copy.deepcopy(geneGB); # copy original GB object,
                newMRNAs = newGB.findAnnsType("mRNA")+newGB.findAnnsType("tRNA",True)+newGB.findAnnsType("rRNA",True); # list of all mRNAs in new GB file. EDIT 11 oct 2017: Include rRNA and tRNA, not just mRNAs!
                for newMRNA in newMRNAs: # loop over mRNAs of new GB file
                    if gene.index[0] <= newMRNA.index[0] < newMRNA.index[1] <= gene.index[1] and mRNA.label != newMRNA.label: # if inside gene and different from current mRNA,
                        newGB.features.remove(newMRNA); # remove newMRNA from new GB file


                newGB.name = mRNA.label; # set new GB object's name to mRNA name
                gbDict[mRNA.label] = newGB; # saves new GB object to output dictionary



    if len(gbDict) == 0: # if no output saved until now,
        gbDict[geneName] = geneGB; # save original gb file as output

    return gbDict;

'''
Check if gene's name is in list of predicted signal peptides
'''
def chkSignalPeptide5Prime(geneName):
    sigPep = False; # TODO
    return sigPep;
