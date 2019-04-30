# -*- coding: utf-8 -*-
"""
Created on June 7 2016
GeneTargeterMethods
Methods for creating constructs aimed at inserting regulatory elements into
chromosomal DNA.
For now, creates a construct based on pSN054 (updated) V5 or pSN150 (5'
targeting) developed by the Niles Lab at MIT tagged designed to deliver the
regulatory element 3' UTR payload to a specific gene, given as a parameter.
@author: Pablo the awesome molecular jedi
"""

#### Libraries ####
import os; # needed for file handling
import copy; # needed for deepcopy
from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants
from py.genetargeter.gRNASelection import *; # Imports methods
from py.genetargeter.HRSelection import *; # Imports methods
from py.genetargeter.recodeSelection import *; # Imports methods
from py.genetargeter.plasmidConstruction import *; # Imports methods
from py.genetargeter.primerSelection import *; # Imports methods

#### Methods ####
"""
Inserts targeting elements into plasmid in order to target it to the given gene
through CRISPR homologous recombination. Input: Gene and gRNAs
Returns dictionary with a bunch of stuff
[0] New annotated plasmid targeted at given gene, containing primers
necessary for Gibson assembly, necessary Klenow oligos and gBlock fragments, if
necessary (annotates a Klenow reaction or gBlock fragment in plasmid given if
codon recoded region is more than minGBlockSize bp long).
[1]
Designed for homologous recombination in Plasmodium falciparum D7.

USAGE: download GenBank sequence file from Benchling with at least one
annotation of the gene to be targeted in 5' to 3' sense, and at least annotation
of the guide RNA chosen for CRISPR recombination. The gene annotation label must
be unique within the sequence. The gRNA must have "gRNA 1" in its label. If more
than one gRNA is to be tested using the same vector (just by replacing the gRNA
in the Cas9 RNA sequence, thus making the recoded region necessarily longer),
annotate all gRNAs in order of preferrence ("gRNA 1", "gRNA 2", etc.).

geneName is a string with the exact label of the gene to be annotated.
TODO: update geneFileName is the name of the .gb file (including ".gb") containing the gene
    to be targeted. The file must be within the input folder. Alternatively, it
    can contain the raw data of an already read file, if useFileStrs is True.
gibsonHomRange is a list of three parameters [min, preferred, max] defining the
    range and preferrence of size for homologous regions used in Gibson Assembly
    (used in vector construction, not homologous regions used for homologous
    recombination in final target). Usually 30, but due to P. falciparum's AT-
    rich genome, we amp it up to 40 bp.
minGBlockSize defines the minimum size a gBlock can have (125 for IDT, june
    2016).
HRannotated is true if the GenBank file given as input includes manual LHR and
    RHR annotations, false otherwise (in which case the function determines
    them).
filterCutSites is a list of strings containing cut sequences to be filtered if
    user provides LHR and RHR.
"""
def targetGene(geneName, geneGB, codonOptimize="T. gondii", HRannotated=False, lengthLHR=[450,500,650], lengthRHR=[450,500,750], gibsonHomRange=[30,40,50], optimRangeLHR=[-20,10], optimRangeRHR=[-20,20], endSizeLHR=40, endSizeRHR=40, endTempLHR=55, endTempRHR=59, gibTemp=65, gibTDif=5, maxDistLHR=500, maxDistRHR=500, minGBlockSize=10, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI], useFileStrs=False, codonSampling=False, minGRNAGCContent=0.3, onTargetMethod="azimuth", minOnTargetScore=30, offTargetMethod="cfd", offTargetThreshold=0.5, maxOffTargetHitScore=35, enzyme="Cas9", PAM="NGG", gBlockDefault=True, plasmidType="pSN054", haTag=True, sigPep=False): # cfd, 0.5, 35; hsu, 75, 5

    outputDicHA = {"geneName":geneName, "newGene":GenBank(), "editedLocus":GenBank(), "newPlasmid":GenBank(), "geneFileStr":"", "plasmidFileStr":"", "oligoFileStr":"", "logFileStr":"", "editedLocusFileStr":"", "gRNATable":""}; # dictionary containing keys to all values being returned for HA-containing design
    outputDic = {"geneName":geneName, "newGene":GenBank(), "editedLocus":GenBank(), "newPlasmid":GenBank(), "geneFileStr":"", "plasmidFileStr":"", "oligoFileStr":"", "logFileStr":"", "editedLocusFileStr":"", "gRNATable":"", "outputHA":outputDicHA}; # dictionary containing keys to all values being returned
    outputDic["logFileStr"] = outputDic["logFileStr"] + " **** Message log for " + geneName + "-targeting construct based on plasmid " + plasmidType + "_" + enzyme + " **** \n\n"; # starts message log to file

    if sigPep: # if gene in signal peptide list,
        outputDic["logFileStr"] = outputDic["logFileStr"] + "Gene contains putative signal peptide according to SignalP information downloaded from PlasmoDB.\n\n"; # say so in output log

    if useFileStrs: # if we are using file strings,
        path = ""; # sets an empty path, needed for save functions of class GenBank
    else: # If we're using files,
        path = "output/" + geneName; # string with path to this gene's directory in output folder
        if not os.path.exists(path): # if a directory for this gene does not exist already,
            os.mkdir(path); # makes new directory

    geneOrientationNegative = False; # used to flip GenBank objects at end if originally on comp strand
    geneList = geneGB.findAnnsLabel(geneName); # search for annotations with gene name
    if len(geneList) > 0: # if genes found,
        gene = geneList[0]; # stores gene annotation
        geneAnns = geneGB.findAnnsLabel(geneName); # stores all gene annotations with gene name in them
        for a in geneAnns: # loops through found annotations
            if a.label == geneName and (a.type == "gene" or a.type == "mRNA" or a.type == "tRNA" or a.type == "rRNA"): # if this annotation's type is gene or mRNA, rRNA, or tRNA,
                gene = a; # keep it as the gene annotation
                if a.comp: # if gene is on complementary strand,
                    geneOrientationNegative = True; # used to flip GenBank objects at end
                    geneGB = geneGB.revComp(); # flip whole GenBank object for processing
                    newAnns = geneGB.findAnnsLabel(a.label); # search new GenBank object for gene annotations
                    for newA in newAnns: # loop through search products
                        if a.label == newA.label and (newA.type == "gene" or newA.type == "mRNA" or a.type == "tRNA" or a.type == "rRNA"): # if this annotation has the exact same name as the original and is a gene or mRNA, rRNA, or tRNA,
                            gene = newA; # save it as the gene annotation object

                break; # stop loop



        closestGene = 0; # initialize variable storing index inside closest relevant gene (depends on targeting position)
        target3Prime = True; # initialize position being targeted, default 3' end of gene
        plasmid = GenBank(); # stores plasmid map variable
        if plasmidType == "pSN054": # if using 3' plasmid
            if enzyme == "Cas9": # if Cas9,
                plasmid = pSN054_V5_Cas9; # set plasmid
            elif enzyme == "Cas12": # if Cas12,
                plasmid = pSN054_V5_Cas12; # set plasmid

            closestGene = len(geneGB.origin); # by default, assume next gene downstream is after end of file
            target3Prime = True; # this plasmid targets 3' end
            for ann in geneGB.features: # for every annotation,
                if ann.index[0] > gene.index[1] and ann.index[0] < closestGene and not "LHR" in ann.label.upper() and not "RHR" in ann.label.upper() and not "GRNA" in ann.label.upper() and not "PRIMER" in ann.label.upper() and not "OLIGO" in ann.label.upper(): # if this annotation is downstream of the gene but upstream of the previous annotation found as possibly immediately downstream and the annotation is not a HR, a gRNA, or an oligo,
                    closestGene = ann.index[0]; # save this annotation's start as the start of the next gene



        elif plasmidType == "pSN150": # if using 5' plasmid
            if enzyme == "Cas9": # if Cas9,
                plasmid = pSN150_Cas9; # set plasmid
            elif enzyme == "Cas12": # if Cas12,
                plasmid = pSN150_Cas12; # set plasmid

            closestGene = 0; # by default, assume next gene downstream is before start of file
            target3Prime = False; # this plasmid targets 3' end
            for ann in geneGB.features: # for every annotation,
                if ann.index[1] > gene.index[0] and ann.index[1] > closestGene and not "LHR" in ann.label.upper() and not "RHR" in ann.label.upper() and not "GRNA" in ann.label.upper() and not "PRIMER" in ann.label.upper() and not "OLIGO" in ann.label.upper(): # if this annotation is downstream of the gene but upstream of the previous annotation found as possibly immediately downstream and the annotation is not a HR, a gRNA, or an oligo,
                    closestGene = ann.index[1]; # save this annotation's start as the start of the next gene



        codingGene = True; # used further on to decide whether to take stop codons into account
        if not (compareSeq(gene.seq[len(gene.seq)-3:len(gene.seq)],"TAA") or compareSeq(gene.seq[len(gene.seq)-3:len(gene.seq)],"TAG") or compareSeq(gene.seq[len(gene.seq)-3:len(gene.seq)],"TGA") ) or ( "RNA" in (gene.type + "_" + gene.label).upper() and not "MRNA" in (gene.type + "_" + gene.label).upper() ): # if gene ends on a stop codon or has the word RNA in its name or type but isn't an mRNA,
            outputDic["logFileStr"] = outputDic["logFileStr"] + "Note: this gene appears not to be a protein-coding sequence (it does not end on a stop codon).\n"; # starts message log to file
            codingGene = False;

        if HRannotated: # if using manual annotations
            gRNA = findGRNA(geneGB, gene); # finds gRNA most upstream annotated manually.
            if len(gRNA["out"].label) == 0: # if no manual annotation found,
                gRNA = chooseGRNA(geneGB, gene, PAM=PAM, minGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, onTargetMethod=onTargetMethod, minOffTargetScore=offTargetThreshold, offTargetMethod=offTargetMethod, maxOffTargetHitScore=maxOffTargetHitScore, gBlockOverlapSize=gibsonHomRange[1], codingGene=codingGene, enzyme=enzyme, closestGene=closestGene, target3Prime=target3Prime, filterCutSites=filterCutSites); # chooses gRNA.

        else: # if not,
            gRNA = chooseGRNA(geneGB, gene, PAM=PAM, minGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, onTargetMethod=onTargetMethod, minOffTargetScore=offTargetThreshold, offTargetMethod=offTargetMethod, maxOffTargetHitScore=maxOffTargetHitScore, gBlockOverlapSize=gibsonHomRange[1], codingGene=codingGene, enzyme=enzyme, closestGene=closestGene, target3Prime=target3Prime, filterCutSites=filterCutSites); # chooses gRNA.

        outputDic["logFileStr"] = outputDic["logFileStr"] + gRNA["log"]; # add logs
        outputDic["gRNATable"] = gRNA["gRNATable"]; # saves gRNA output values
        gRNA = gRNA["out"]; # saves actual data

        if len(gRNA.label) > 0: # if gRNA found,
            LHR = GenBankAnn(); # init LHR var
            RHR = GenBankAnn(); # init RHR var

            # pick HRs first
            if target3Prime: # if going for 3' payload,
                LHR = chooseHR(geneGB, gene, doingHR='LHR', targetExtreme='end', lengthHR=lengthLHR, minTmEnds=endTempLHR, endsLength=endSizeLHR, gBlockDefault=gBlockDefault, minGBlockSize=minGBlockSize, codingGene=codingGene, filterCutSites=filterCutSites); # chooses an LHR
                RHR = chooseHR(geneGB, gene, doingHR='RHR', targetExtreme='end', lengthHR=lengthRHR, minTmEnds=endTempRHR, endsLength=endSizeRHR, filterCutSites=filterCutSites); # chooses RHR
            else: # if going for 5' payload, TODO: switch params?
                LHR = chooseHR(geneGB, gene, doingHR='LHR', targetExtreme='start', lengthHR=lengthLHR, minTmEnds=endTempLHR, endsLength=endSizeLHR, filterCutSites=filterCutSites); # chooses LHR
                RHR = chooseHR(geneGB, gene, doingHR='RHR', targetExtreme='start', lengthHR=lengthRHR, minTmEnds=endTempRHR, endsLength=endSizeRHR, gBlockDefault=gBlockDefault, minGBlockSize=minGBlockSize, codingGene=codingGene, filterCutSites=filterCutSites); # chooses an RHR

            if HRannotated: # if LHR and RHR are already annotated,
                LHRlist = geneGB.findAnnsLabel("LHR"); # overwrite LHR annotations
                if len(LHRlist) > 0: # if LHR found,
                    LHR = LHRlist[0]; # save first LHR annotation as LHR
                    outputDic["logFileStr"] = outputDic["logFileStr"] + "\nFound user LHR annotation, replaced automatic annotation with it." + "\n"; # add warning to log
                else: # if no LHR found,
                    outputDic["logFileStr"] = outputDic["logFileStr"] + "\nWarning: Did not find user LHR annotation, used automatic annotation instead." + "\n"; # add warning to log
                    outputDic["logFileStr"] = outputDic["logFileStr"] + LHR["log"]; # add logs
                    LHR = LHR["out"]; # saves actual data
                    geneGB.features.append(LHR); # adds LHR to gene annotations

                RHRlist = geneGB.findAnnsLabel("RHR"); # saves RHR annotation
                if len(RHRlist) > 0: # if RHR found,
                    RHR = RHRlist[0]; # save first RHR annotation as RHR
                    outputDic["logFileStr"] = outputDic["logFileStr"] + "\nFound user RHR annotation, replaced automatic annotation with it." + "\n"; # add warning to log
                else: # if no RHR found,
                    outputDic["logFileStr"] = outputDic["logFileStr"] + "\nWarning: Did not find user RHR annotation, used automatic annotation instead." + "\n"; # add warning to log
                    outputDic["logFileStr"] = outputDic["logFileStr"] + RHR["log"]; # add logs
                    RHR = RHR["out"]; # saves actual data
                    geneGB.features.append(RHR); # adds RHR to gene annotations

                for site in filterCutSites: # for every cut site being filtered
                    if findFirst(site,LHR.seq) > -1 or findFirst(revComp(site),LHR.seq) > -1: # if cut site found,
                        outputDic["logFileStr"] = outputDic["logFileStr"] + "\nWarning: LHR sequence for gene " + gene.label + ": \n" + LHR + "\ncontains restriction site " + site + "\n"; # add warning to log
                    if findFirst(site,RHR.seq) > -1 or findFirst(revComp(site),RHR.seq) > -1: # if cut site found,
                        outputDic["logFileStr"] = outputDic["logFileStr"] + "\nWarning: RHR sequence for gene " + gene.label + ": \n" + RHR + "\ncontains restriction site " + site + "\n"; # add warning to log


            else: # if HRs not annotated,
                outputDic["logFileStr"] = outputDic["logFileStr"] + LHR["log"]; # add logs
                LHR = LHR["out"]; # saves actual data
                outputDic["logFileStr"] = outputDic["logFileStr"] + RHR["log"]; # add logs
                RHR = RHR["out"]; # saves actual data
                geneGB.features.append(LHR); # adds LHR to gene annotations
                geneGB.features.append(RHR); # adds RHR to gene annotations



            recoded = chooseRecodeRegion(geneGB, gene, offTargetMethod, pamType=PAM, orgCodonTable=codonUsageTables[codonOptimize],codonSampling=codonSampling, gRNATableString=outputDic["gRNATable"], target3Prime=target3Prime, filterCutSites=filterCutSites); # defines region to be recoded, returns recoded sequence
            recodedHA = {}; # will contain recoded region with HA tag
            if haTag: # if using HA tags,
                recodedHA = chooseRecodeRegion(geneGB, gene, offTargetMethod, pamType=PAM, orgCodonTable=codonUsageTables[codonOptimize],codonSampling=codonSampling, gRNATableString=outputDic["gRNATable"], target3Prime=target3Prime, haTag=True, filterCutSites=filterCutSites); # defines region to be recoded with HA tag, returns recoded sequence
                recodedHA = recodedHA["out"]; # saves actual data

            outputDic["logFileStr"] = outputDic["logFileStr"] + recoded["log"]; # add logs
            outputDic["gRNATable"] = recoded["gRNATable"]; # saves gRNA output values
            recoded = recoded["out"]; # saves actual data

            plasmidArmed = insertTargetingElements(plasmid, gene.label, gRNA.seq, LHR.seq, recoded.seq, RHR.seq, plasmidType=plasmidType, haTag=False); # inserts targeting elements
            plasmidArmedHA = ""; # will contain plasmid with HA tags
            outputDicHA = copy.deepcopy(outputDic); # will store outputs
            outputDicHA["logFileStr"] = outputDicHA["logFileStr"].replace(" **** \n\n", "_HA_Tag **** \n\n")
            outputDic = postProcessPlasmid(geneName, geneGB, gene, plasmidArmed, recoded, outputDic, path, useFileStrs, geneOrientationNegative=geneOrientationNegative, plasmidType=plasmidType, enzyme=enzyme, gibsonHomRange=gibsonHomRange, gibTemp=gibTemp, gibTDif=gibTDif, minGBlockSize=minGBlockSize, haTag=False); # generate and annotate assembly info
            if haTag: # if using HA tags,
                plasmidArmedHA = insertTargetingElements(plasmid, gene.label, gRNA.seq, LHR.seq, recodedHA.seq, RHR.seq, plasmidType=plasmidType, haTag=True); # inserts targeting elements
                outputDicHA = postProcessPlasmid(geneName, geneGB, gene, plasmidArmedHA, recodedHA, outputDicHA, path, useFileStrs, geneOrientationNegative=geneOrientationNegative, plasmidType=plasmidType, enzyme=enzyme, gibsonHomRange=gibsonHomRange, gibTemp=gibTemp, gibTDif=gibTDif, minGBlockSize=minGBlockSize, haTag=True); # generate and annotate assembly info
                outputDic["outputHA"] = outputDicHA; # if using HA tags, save design inside output dictionary
                if not useFileStrs: # if saving to files,
                    output(outputDicHA["oligoFileStr"], path + "/" + geneName + plasmidType + "_" + enzyme + "HA_Tags_Oligos.csv",wipe=True); # saves oligos to file
                    output(outputDicHA["logFileStr"], path + "/" + geneName + plasmidType + "_" + enzyme + "HA_Tags_Message_Log.txt",wipe=True); # saves message log to file

            if not useFileStrs: # if saving to files,
                output(outputDic["oligoFileStr"], path + "/" + geneName + plasmidType + "_" + enzyme + "_Oligos.csv",wipe=True); # saves oligos to file
                output(outputDic["logFileStr"], path + "/" + geneName + plasmidType + "_" + enzyme + "_Message_Log.txt",wipe=True); # saves message log to file



    else: # if no gene found,
        outputDic["logFileStr"] = outputDic["logFileStr"] + "\nERROR: No gene annotations found in this file with name " + geneName + "\nProcess terminated.\n";

    return outputDic; # returns output dictionary

def postProcessPlasmid(geneName, geneGB, gene, plasmidArmed, recoded, outputDic, path, useFileStrs, geneOrientationNegative, plasmidType="pSN054", enzyme="Cas9", gibsonHomRange=[30,40,50], gibTemp=65, gibTDif=5, minGBlockSize=10, haTag=False):
    gRNAOnPlasmid = plasmidArmed.findAnnsLabel(gene.label + " gRNA")[0]; # saves gRNA annotation actually on plasmid
    if len(recoded.seq) > 0: # if there actually is a recoded region,
        recodedOnPlasmid = plasmidArmed.findAnnsLabel(gene.label + " Recoded")[0]; # saves recoded annotation actually on plasmid

    LHROnPlasmid = plasmidArmed.findAnnsLabel(gene.label + " LHR")[0]; # saves LHR annotation actually on plasmid
    RHROnPlasmid = plasmidArmed.findAnnsLabel(gene.label + " RHR")[0]; # saves RHR annotation actually on plasmid

    primerString = "Primer name,Sequence"; # "OLIGOS for construct targeting gene " + geneName + "\n\n"; # String will save all primer information to be written to file

    #TODO: Check what changes with pSN150 from here to end of method (check createPrimers, createGBlock, createKlenowOligos, createGibsonPrimers, editLocus)
    if len(recoded.seq) > 0 and len(recoded.seq) + gibsonHomRange[1]*2 >= minGBlockSize: # if there is a recoded region and length of recoded region plus homology regions necessary for Gibson Assembly is greater or equal to minimum gBlock size,
        gBlock = createGBlock(plasmidArmed,recodedOnPlasmid,gibsonHomRange[1]); # annotates gBlock on plasmid
        outputDic["logFileStr"] = outputDic["logFileStr"] + gBlock["log"]; # add logs
        gBlock = gBlock["out"]; # saves actual data
        plasmidArmed.features.append(gBlock); # add to plasmid annotations

        primGBlock = createPrimers(plasmidArmed, gBlock); # creates gBlock primers
        outputDic["logFileStr"] = outputDic["logFileStr"] + primGBlock["log"]; # add logs
        primGBlock = primGBlock["out"]; # saves actual data
        plasmidArmed.features.append(primGBlock[0]); # add fwd primer to plasmid annotations
        plasmidArmed.features.append(primGBlock[1]); # add rev primer to plasmid annotations

        primerString = primerString + "\n" + geneName + " gBlock primer (fwd)," + primGBlock[0].seq + "\n" + geneName + " gBlock primer (rev)," + primGBlock[1].seq; # write primers to output string
    elif len(recoded.seq) >= gibsonHomRange[0]/2: # if length of recoded region greater or equal to half of minimum size of homology region,
        outputDic["logFileStr"] = outputDic["logFileStr"] + "gBlock deemed not feasible for recoded region of construct targeting gene " + geneName + ", used Klenow instead.\n"; # say so
        klenowRecoded = createKlenowOligos(plasmidArmed, recodedOnPlasmid, gibsonHomRange[1]); # creates Klenow oligos
        outputDic["logFileStr"] = outputDic["logFileStr"] + klenowRecoded["log"]; # add logs
        klenowRecoded = klenowRecoded["out"]; # saves actual data
        plasmidArmed.features.append(klenowRecoded[0]); # add fwd primer to plasmid annotations
        plasmidArmed.features.append(klenowRecoded[1]); # add rev primer to plasmid annotations
        primerString = primerString + "\n" + geneName + " Recoded region Klenow oligo (fwd)," + klenowRecoded[0].seq + "\n" + geneName + " Recoded region Klenow oligo (rev)," + klenowRecoded[1].seq; # write oligos to output string
    else: # if Klenow unnecesary too,
        outputDic["logFileStr"] = outputDic["logFileStr"] + "\ngBlock and Klenow for recoded region deemed unnecesary for construct targeting gene " + geneName +".\n\n"; # say so

    primLHR = createGibsonPrimers(plasmidArmed, LHROnPlasmid, rangeHom=gibsonHomRange, minMeltTemp=gibTemp, maxTempDif=gibTDif); # creates LHR Gibson primers
    outputDic["logFileStr"] = outputDic["logFileStr"] + primLHR["log"]; # add logs
    primLHR = primLHR["out"]; # saves actual data
    plasmidArmed.features.append(primLHR[0]); # add fwd primer to plasmid annotations
    plasmidArmed.features.append(primLHR[1]); # add rev primer to plasmid annotations
    primerString = primerString + "\n" + geneName + " LHR primer (fwd)," + primLHR[0].seq + "\n" + geneName + " LHR primer (rev)," + primLHR[1].seq; # write primers to output string

    primRHR = createGibsonPrimers(plasmidArmed, RHROnPlasmid, rangeHom=gibsonHomRange, minMeltTemp=gibTemp, maxTempDif=gibTDif); # creates RHR Gibson primers
    outputDic["logFileStr"] = outputDic["logFileStr"] + primRHR["log"]; # add logs
    primRHR = primRHR["out"]; # saves actual data
    plasmidArmed.features.append(primRHR[0]); # add fwd primer to plasmid annotations
    plasmidArmed.features.append(primRHR[1]); # add rev primer to plasmid annotations
    primerString = primerString  + "\n" + geneName + " RHR primer (fwd)," + primRHR[0].seq + "\n" + geneName + " RHR primers (rev)," + primRHR[1].seq; # write primers to output string

    klenow = createKlenowOligos(plasmidArmed, gRNAOnPlasmid, gibsonHomRange[1]); # creates Klenow oligos
    outputDic["logFileStr"] = outputDic["logFileStr"] + klenow["log"]; # add logs
    klenow = klenow["out"]; # saves actual data
    plasmidArmed.features.append(klenow[0]); # add fwd primer to plasmid annotations
    plasmidArmed.features.append(klenow[1]); # add rev primer to plasmid annotations
    primerString = primerString + "\n" + geneName + " gRNA Klenow oligo (fwd)," + klenow[0].seq + "\n" + geneName + " gRNA Klenow oligo (rev)," + klenow[1].seq; # write oligos to output string

    gRNACassetteStart = -1 # will be filled in below
    gRNACassetteEnd = -1

    if plasmidType == "pSN054": # if using pSN150 instead of pSN054,
        gRNACassetteStart = plasmidArmed.findAnnsLabel("Lox")[0].index[0]; # gBlock starts at first Lox
        gRNACassetteEnd = plasmidArmed.findAnnsLabel("RHR_vector overlap_lef")[0].index[1]; # gBlock ends at RHR_vector overlap_lef
    if plasmidType == "pSN150": # if using pSN150 instead of pSN054,
        gRNACassetteStart = findFirst(plasmidArmed.origin, cut_BsiWI) + 1; # gBlock starts at BsiWI cut
        gRNACassetteEnd = findFirst(plasmidArmed.origin, cut_AsiSI) + 4; # gBlock ends at AsiSI cut

    gRNACassette = GenBankAnn("sgRNA cassette", "misc", plasmidArmed.origin[gRNACassetteStart:gRNACassetteEnd], False, [gRNACassetteStart,gRNACassetteEnd], annColors["otherAnnColor"]); # create cassette annotation
    gRNAGBlock = createGBlock(plasmidArmed,gRNACassette,gibsonHomRange[1]); # annotates gBlock on plasmid
    outputDic["logFileStr"] = outputDic["logFileStr"] + gRNAGBlock["log"]; # add logs
    gRNAGBlock = gRNAGBlock["out"]; # saves actual data
    plasmidArmed.features.append(gRNAGBlock); # add to plasmid annotations

    primGRNAGBlock = createPrimers(plasmidArmed, gRNAGBlock); # creates gBlock primers
    outputDic["logFileStr"] = outputDic["logFileStr"] + primGRNAGBlock["log"]; # add logs
    primGRNAGBlock = primGRNAGBlock["out"]; # saves actual data
    plasmidArmed.features.append(primGRNAGBlock[0]); # add fwd primer to plasmid annotations
    plasmidArmed.features.append(primGRNAGBlock[1]); # add rev primer to plasmid annotations

    primerString = primerString + "\n" + geneName + "sgRNA cassette gBlock primer (fwd)," + primGRNAGBlock[0].seq + "\n" + geneName + "sgRNA cassette gBlock primer (rev)," + primGRNAGBlock[1].seq; # write primers to output string

    primerString = shortenOligoNames(primerString) + "\n"; # abbreviates primer names to fit on commercial tube labels

    editedLocus = editLocus(geneName, geneGB, plasmidArmed); # inserts construct into genomic context
    outputDic["logFileStr"] += editedLocus["log"]; # add logs
    editedLocus = editedLocus["out"]; # saves actual data

    outputDic["logFileStr"] = outputDic["logFileStr"] + "\nVector constructed and written to file. End of process.\n"; # saves message log to file
    plasmidArmed.definition = (plasmidArmed.definition + "  " + outputDic["logFileStr"]).replace("\n","   "); # save logs to file definition to be viewed in benchling

    if geneOrientationNegative: # if gene was originally on comp strand,
        geneGB = geneGB.revComp(); # flip pre-editing locus
        editedLocus = editedLocus.revComp(); # flip post-editing locus

    haName = ""; # default no HA tags notated in name
    if haTag: # if using HA tags,
        haName = "HA_tags"; # include in name

    geneGB.name = geneName + "_" + plasmidType + "_" + enzyme + "_" + haName + "_Locus_Pre-editing"; # save file type in genbank locus
    plasmidArmed.name = geneName + "_" + plasmidType + "_" + enzyme + "_" + haName; # save file type in genbank locus
    editedLocus.name = geneName + "_" + plasmidType + "_" + enzyme + "_" + haName + "_Locus_Post-editing"; # save file type in genbank locus

    outputDic["geneFileStr"] = geneGB.save(path + "/" + geneName + plasmidType + "_" + enzyme + "_Locus_Pre-editing.gb", saveToFile=(not useFileStrs)); # saves annotated gene
    outputDic["plasmidFileStr"] = plasmidArmed.save(path + "/" + plasmidType + "_" + enzyme + "_" + geneName + haName, saveToFile=(not useFileStrs)); # saves plasmid
    outputDic["editedLocusFileStr"] = editedLocus.save(path + "/" + geneName + plasmidType + "_" + enzyme + haName + "_Locus_Post-editing.gb", saveToFile=(not useFileStrs)); # saves edited locus
    outputDic["oligoFileStr"] = primerString; # saves primers to file
    outputDic["newPlasmid"] = plasmidArmed; # saves new plasmid to output dictionary
    outputDic["newGene"] = geneGB; # saves new plasmid to output dictionary
    outputDic["editedLocus"] = editedLocus; # saves edited locus to output dictionary
    outputDic["geneName"] = geneName; # saves gene name to output

    return outputDic;