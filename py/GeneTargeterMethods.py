# -*- coding: utf-8 -*-
"""
Created on June 7 2016
GeneTargeterMethods
Methods for creating constructs aimed at inserting regulatory elements into
chromosomal DNA.
For now, creates a construct based on pSN054 (updated) V5 developed by the Niles
Lab at MIT tagged designed to deliver the regulatory element 3' UTR payload to a
specific gene, given as a parameter.
@author: Pablo the awesome molecular jedi
"""

#### Libraries ####
from copy import deepcopy; # Import object copying methods for deep copies
from BioUtils import *; # Imports utils
from GenBankToolbox import *; # Imports utils
import os; # needed for file handling
import operator; # Needed for sorting
from gRNAScores.gRNAScoring import * # import methods for gRNA scoring
from gRNAScores.gRNAScoreDB import * # import methods for gRNA scoring from precalculated DB

#### Constants ####
# Restriction enzyme cut sequences used
cut_FseI = "ggccggcc"; # FseI cut site
cut_AsiSI = "gcgatcgc"; # AsiSI cut site
cut_IPpoI = "ctctcttaaggtagc"; # I-PpoI cut site
cut_ISceI = "attaccctgttatcccta"; # I-SceI cut site
cut_AflII = "cttaag"; # AflII cut site
# Plasmids
pSN054_V5_Cas9 = GenBank();
pSN054_V5_Cas9.load("input/plasmids/psn054-updated_v5-tagged-jn_final.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN054_V5_Cpf1 = GenBank();
pSN054_V5_Cpf1.load("input/plasmids/psn054_v5ha-tags_lbcpf1.gb",loadFromFile=True); # load Cpf1 plasmid sequence from GenBank format
# Codon usage tables
codonUsageTables = {
    'P. falciparum 3D7': codonUsage('input/codonUsageTables/Pfalciparum3D7.txt'),
    'P. vivax': codonUsage('input/codonUsageTables/Pvivax.txt'),
    'T. gondii': codonUsage('input/codonUsageTables/Tgondii.txt'),
    'E. coli K12': codonUsage('input/codonUsageTables/EcoliK12.txt'),
    'S. cerevisiae': codonUsage('input/codonUsageTables/Scerevisiae.txt'),
    'H. sapiens': codonUsage('input/codonUsageTables/Hsapiens.txt'),
    'R. norvegicus': codonUsage('input/codonUsageTables/Rnorvegicus.txt'),
    'scramble': codonUsage()
}


#### Methods ####

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


    geneList = geneGB.findAnnsLabel(geneName); # search for annotations with gene name
    if len(geneList) > 0: # if genes found,
        gene = geneList[0]; # stores gene annotation
        geneAnns = geneGB.findAnnsLabel(geneName); # stores all gene annotations with gene name in them
        for a in geneAnns: # loops through found annotations
            if a.type == "gene": # if this annotation's type is gene,
                gene = a; # keep it as the gene annotation
                break; # stop loop

        mRNAs = geneGB.findAnnsType("mRNA")+geneGB.findAnnsType("tRNA")+geneGB.findAnnsType("rRNA"); # list of all mRNAs in GB file. EDIT 11 oct 2017: Include rRNA and tRNA, not just mRNAs!
        for mRNA in mRNAs: # loop over all mRNAs
            if gene.index[0] <= mRNA.index[0] < mRNA.index[1] <= gene.index[1]: # if mRNA is inside gene,
                newGB = deepcopy(geneGB); # copy original GB object,
                newMRNAs = newGB.findAnnsType("mRNA")+newGB.findAnnsType("tRNA")+newGB.findAnnsType("rRNA"); # list of all mRNAs in new GB file. EDIT 11 oct 2017: Include rRNA and tRNA, not just mRNAs!
                for newMRNA in newMRNAs: # loop over mRNAs of new GB file
                    if gene.index[0] <= newMRNA.index[0] < newMRNA.index[1] <= gene.index[1] and mRNA.label != newMRNA.label: # if inside gene and different from current mRNA,
                        newGB.features.remove(newMRNA); # remove newMRNA from new GB file


                newGB.name = mRNA.label; # set new GB object's name to mRNA name
                gbDict[mRNA.label] = newGB; # saves new GB object to output dictionary



    if len(gbDict) == 0: # if no output saved until now,
        gbDict[geneName] = geneGB; # save original gb file as output

    return gbDict;


"""
Inserts targeting elements into pSN054 in order to target it to the given gene
through CRISPR homologous recombination. Input: Gene and gRNAs
Returns dictionary with
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
def pSN054TargetGene(geneName, geneGB, codonOptimize="T. gondii", HRannotated=False, lengthLHR=[450,500,650], lengthRHR=[450,500,750], gibsonHomRange=[30,40,50], optimRangeLHR=[-20,10], optimRangeRHR=[-20,20], endSizeLHR=40, endSizeRHR=40, endTempLHR=55, endTempRHR=59, gibTemp=65, gibTDif=5, maxDistLHR=500, maxDistRHR=500, minGBlockSize=10, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AflII], useFileStrs=False, codonSampling=False, minGRNAGCContent=0.3, onTargetMethod="azimuth", minOnTargetScore=30, offTargetMethod="cfd", offTargetThreshold=0.5, maxOffTargetHitScore=35, enzyme="Cas9", PAM="NGG", gBlockDefault=True): # cfd, 0.5, 35; hsu, 75, 5
    outputDic = {"geneName":geneName, "newGene":GenBank(), "editedLocus":GenBank(), "newPlasmid":GenBank(), "geneFileStr":"", "plasmidFileStr":"", "oligoFileStr":"", "logFileStr":"", "editedLocusFileStr":"", "gRNATable":""}; # dictionary containing keys to all values being returned
    outputDic["logFileStr"] = outputDic["logFileStr"] + " **** Message log for " + geneName + "-targeting construct based on plasmid pSN054_V5_"+ enzyme +" **** \n\n"; # starts message log to file

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



        if not (compareSeq(gene.seq[len(gene.seq)-3:len(gene.seq)],"TAA") or compareSeq(gene.seq[len(gene.seq)-3:len(gene.seq)],"TAG") or compareSeq(gene.seq[len(gene.seq)-3:len(gene.seq)],"TGA")): # if gene ends on a stop codon,
            outputDic["logFileStr"] = outputDic["logFileStr"] + "Note: this gene appears not to be a protein-coding sequence (it does not end on a stop codon).\n"; # starts message log to file

        if HRannotated: # if using manual annotations
            gRNA = findGRNA(geneGB, gene); # finds gRNA most upstream annotated manually.
            if len(gRNA["out"].label) == 0: # if no manual annotation found,
                gRNA = chooseGRNA(geneGB, gene, PAM=PAM, minGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, onTargetMethod=onTargetMethod, minOffTargetScore=offTargetThreshold, offTargetMethod=offTargetMethod, maxOffTargetHitScore=maxOffTargetHitScore, enzyme=enzyme); # chooses gRNA.

        else: # if not,
            gRNA = chooseGRNA(geneGB, gene, PAM=PAM, minGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, onTargetMethod=onTargetMethod, minOffTargetScore=offTargetThreshold, offTargetMethod=offTargetMethod, maxOffTargetHitScore=maxOffTargetHitScore, enzyme=enzyme); # chooses gRNA.

        outputDic["logFileStr"] = outputDic["logFileStr"] + gRNA["log"]; # add logs
        outputDic["gRNATable"] = gRNA["gRNATable"]; # saves gRNA output values
        gRNA = gRNA["out"]; # saves actual data

        if len(gRNA.label) > 0: # if gRNA found,
            LHR = GenBankAnn(); # init LHR var
            RHR = GenBankAnn(); # init RHR var

            # pick HRs first
            LHR = chooseLHR(geneGB, gene, lengthLHR=lengthLHR, minTmEnds=endTempLHR, endsLength=endSizeLHR, optimizeRange=optimRangeLHR, maxDistanceFromGRNA=maxDistLHR, gBlockDefault=gBlockDefault, minGBlockSize=minGBlockSize); # chooses an LHR
            RHR = chooseRHR(geneGB, gene, lengthRHR=lengthRHR, minTmEnds=endTempRHR, endsLength=endSizeRHR, optimizeRange=optimRangeRHR, maxDistanceFromGene=maxDistRHR); # chooses RHR

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
                if len(RHRlist) > 0: # if LHR found,
                    RHR = RHRlist[0]; # save first RHR annotation as RHR
                    outputDic["logFileStr"] = outputDic["logFileStr"] + "\nFound user RHR annotation, replaced automatic annotation with it." + "\n"; # add warning to log
                else: # if no LHR found,
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



            recoded = chooseRecodeRegion(geneGB, gene, offTargetMethod, pamType=PAM, orgCodonTable=codonUsageTables[codonOptimize],codonSampling=codonSampling, gRNATableString=outputDic["gRNATable"]); # defines region to be recoded, returns recoded sequence
            outputDic["logFileStr"] = outputDic["logFileStr"] + recoded["log"]; # add logs
            outputDic["gRNATable"] = recoded["gRNATable"]; # saves gRNA output values
            recoded = recoded["out"]; # saves actual data

            plasmid = pSN054_V5_Cas9; # stores plasmid map variable
            if enzyme == "Cas9": # if Cas9,
                plasmid = pSN054_V5_Cas9; # set plasmid
            elif enzyme == "Cpf1": # if Cpf1,
                plasmid = pSN054_V5_Cpf1; # set plasmid

            pSN054_ARMED = insertTargetingElementsPSN054(pSN054_V5_Cas9, gene.label, gRNA.seq, LHR.seq, recoded.seq, RHR.seq); # inserts targeting elements

            gRNAOnPlasmid = pSN054_ARMED.findAnnsLabel(gene.label + " gRNA")[0]; # saves gRNA annotation actually on plasmid
            if len(recoded.seq) > 0: # if there actually is a recoded region,
                recodedOnPlasmid = pSN054_ARMED.findAnnsLabel(gene.label + " Recoded")[0]; # saves recoded annotation actually on plasmid

            LHROnPlasmid = pSN054_ARMED.findAnnsLabel(gene.label + " LHR")[0]; # saves LHR annotation actually on plasmid
            RHROnPlasmid = pSN054_ARMED.findAnnsLabel(gene.label + " RHR")[0]; # saves RHR annotation actually on plasmid

            primerString = "Primer name,Sequence"; # "OLIGOS for construct targeting gene " + geneName + "\n\n"; # String will save all primer information to be written to file

            if len(recoded.seq) > 0 and len(recoded.seq) + gibsonHomRange[1]*2 >= minGBlockSize: # if there is a recoded region and length of recoded region plus homology regions necessary for Gibson Assembly is greater or equal to minimum gBlock size,
                gBlock = createGBlock(pSN054_ARMED,recodedOnPlasmid); # annotates gBlock on plasmid
                outputDic["logFileStr"] = outputDic["logFileStr"] + gBlock["log"]; # add logs
                gBlock = gBlock["out"]; # saves actual data
                pSN054_ARMED.features.append(gBlock); # add to plasmid annotations

                primGBlock = createPrimers(pSN054_ARMED, gBlock); # creates gBlock primers
                outputDic["logFileStr"] = outputDic["logFileStr"] + primGBlock["log"]; # add logs
                primGBlock = primGBlock["out"]; # saves actual data
                pSN054_ARMED.features.append(primGBlock[0]); # add fwd primer to plasmid annotations
                pSN054_ARMED.features.append(primGBlock[1]); # add rev primer to plasmid annotations

                primerString = primerString + "\n" + geneName + " gBlock primer (fwd)," + primGBlock[0].seq + "\n" + geneName + " gBlock primer (rev)," + primGBlock[1].seq; # write primers to output string
            elif len(recoded.seq) >= gibsonHomRange[0]/2: # if length of recoded region greater or equal to half of minimum size of homology region,
                outputDic["logFileStr"] = outputDic["logFileStr"] + "gBlock deemed not feasible for recoded region of construct targeting gene " + geneName + ", used Klenow instead.\n"; # say so
                klenowRecoded = createKlenowOligos(pSN054_ARMED, recodedOnPlasmid, gibsonHomRange[1]); # creates Klenow oligos
                outputDic["logFileStr"] = outputDic["logFileStr"] + klenowRecoded["log"]; # add logs
                klenowRecoded = klenowRecoded["out"]; # saves actual data
                pSN054_ARMED.features.append(klenowRecoded[0]); # add fwd primer to plasmid annotations
                pSN054_ARMED.features.append(klenowRecoded[1]); # add rev primer to plasmid annotations
                primerString = primerString + "\n" + geneName + " Recoded region Klenow oligo (fwd)," + klenowRecoded[0].seq + "\n" + geneName + " Recoded region Klenow oligo (rev)," + klenowRecoded[1].seq; # write oligos to output string
            else: # if Klenow unnecesary too,
                outputDic["logFileStr"] = outputDic["logFileStr"] + "\ngBlock and Klenow for recoded region deemed unnecesary for construct targeting gene " + geneName +".\n\n"; # say so

            primLHR = createGibsonPrimers(pSN054_ARMED, LHROnPlasmid, rangeHom=gibsonHomRange, minMeltTemp=gibTemp, maxTempDif=gibTDif); # creates LHR Gibson primers
            outputDic["logFileStr"] = outputDic["logFileStr"] + primLHR["log"]; # add logs
            primLHR = primLHR["out"]; # saves actual data
            pSN054_ARMED.features.append(primLHR[0]); # add fwd primer to plasmid annotations
            pSN054_ARMED.features.append(primLHR[1]); # add rev primer to plasmid annotations
            primerString = primerString + "\n" + geneName + " LHR primer (fwd)," + primLHR[0].seq + "\n" + geneName + " LHR primer (rev)," + primLHR[1].seq; # write primers to output string

            primRHR = createGibsonPrimers(pSN054_ARMED, RHROnPlasmid, rangeHom=gibsonHomRange, minMeltTemp=gibTemp, maxTempDif=gibTDif); # creates RHR Gibson primers
            outputDic["logFileStr"] = outputDic["logFileStr"] + primRHR["log"]; # add logs
            primRHR = primRHR["out"]; # saves actual data
            pSN054_ARMED.features.append(primRHR[0]); # add fwd primer to plasmid annotations
            pSN054_ARMED.features.append(primRHR[1]); # add rev primer to plasmid annotations
            primerString = primerString  + "\n" + geneName + " RHR primer (fwd)," + primRHR[0].seq + "\n" + geneName + " RHR primers (rev)," + primRHR[1].seq; # write primers to output string

            klenow = createKlenowOligos(pSN054_ARMED, gRNAOnPlasmid, gibsonHomRange[1]); # creates Klenow oligos
            outputDic["logFileStr"] = outputDic["logFileStr"] + klenow["log"]; # add logs
            klenow = klenow["out"]; # saves actual data
            pSN054_ARMED.features.append(klenow[0]); # add fwd primer to plasmid annotations
            pSN054_ARMED.features.append(klenow[1]); # add rev primer to plasmid annotations
            primerString = primerString + "\n" + geneName + " gRNA Klenow oligo (fwd)," + klenow[0].seq + "\n" + geneName + " gRNA Klenow oligo (rev)," + klenow[1].seq; # write oligos to output string

            primerString = shortenOligoNames(primerString) + "\n"; # abbreviates primer names to fit on commercial tube labels

            editedLocus = editLocus(geneName, geneGB, pSN054_ARMED); # inserts construct into genomic context
            outputDic["logFileStr"] += editedLocus["log"]; # add logs
            editedLocus = editedLocus["out"]; # saves actual data

            outputDic["logFileStr"] = outputDic["logFileStr"] + "\nVector constructed and written to file. End of process.\n"; # saves message log to file
            pSN054_ARMED.definition = (pSN054_ARMED.definition + "  " + outputDic["logFileStr"]).replace("\n","   "); # save logs to file definition to be viewed in benchling

            if geneOrientationNegative: # if gene was originally on comp strand,
                geneGB = geneGB.revComp(); # flip pre-editing locus
                editedLocus = editedLocus.revComp(); # flip post-editing locus

            outputDic["geneFileStr"] = geneGB.save(path + "/" + geneName + "_Locus_Pre-editing.gb", saveToFile=(not useFileStrs)); # saves annotated gene
            outputDic["plasmidFileStr"] = pSN054_ARMED.save(path + "/" +  "pSN054_V5_Cas9_targeting" + geneName, saveToFile=(not useFileStrs)); # saves plasmid
            outputDic["editedLocusFileStr"] = editedLocus.save(path + "/" + geneName+"_Locus_Post-editing.gb", saveToFile=(not useFileStrs)); # saves edited locus
            outputDic["oligoFileStr"] = primerString; # saves primers to file
            outputDic["newPlasmid"] = pSN054_ARMED; # saves new plasmid to output dictionary
            outputDic["newGene"] = geneGB; # saves new plasmid to output dictionary
            outputDic["editedLocus"] = editedLocus; # saves edited locus to output dictionary
            outputDic["geneName"] = geneName; # saves gene name to output

            if not useFileStrs: # if saving to files,
                output(outputDic["oligoFileStr"], path + "/" + geneName + "_Oligos.csv",wipe=True); # saves oligos to file
                output(outputDic["logFileStr"], path + "/" + geneName + "_Message_Log.txt",wipe=True); # saves message log to file



    else: # if no gene found,
        outputDic["logFileStr"] = outputDic["logFileStr"] + "\nERROR: No gene annotations found in this file with name " + geneName + "\nProcess terminated.\n";

    return outputDic; # returns output dictionary


"""
Inserts targeting elements given as arguments into pSN054_V5_Cas9 at predetermined
sites. Elements given as arguments must be strings, not GenBankAnn objects.
The gRNA, LHR and RHR given must be in 5' to 3' sense and in the
positive strand and must not contain RE cut sites cut_FseI, cut_AsiSI,
cut_IPpoI, cut_ISceI, or cut_AflII. Returns new GenBank object with targeting
elements.
Specifically:
Inserts gRNA used by Cas9 in I-PpoI cut site. Removes recognition sequence. Adds
    GG at 5’-end of the gRNA sequence. The GG is required for the T7 RNA
    polymerase to efficiently transcribe the gRNA to high levels.
Inserts the LHR between FseI and AsiSI cut sites, but leaves the sites intact
    with their recognition sequences.
Inserts the recoded region after the LHR, before the AsiSI cut site (leaves the
    site intact).
Inserts RHR after I-SceI cut site, leaves the site intact with its recognition
    sequence.
"""
def insertTargetingElementsPSN054(plasmid, geneName, gRNA, LHR, recodedRegion, RHR):
    plas = deepcopy(plasmid); # makes a copy of the plasmid object to modify without altering the original

    startLHR = findFirst(plas.origin, cut_FseI) + len(cut_FseI); # index of LHR start site (at end of FseI cut sequence)
    endLHR = findFirst(plas.origin, cut_AsiSI); # index of LHR end site (at start of AsiSI cut sequence)
    plas.removeSeq([startLHR, endLHR]); # removes sequence that LHR will replace
    plas.insertSeq(LHR, startLHR); # inserts LHR sequence
    annLHR = GenBankAnn(geneName+" LHR", "misc_feature", LHR, False, [startLHR,startLHR+len(LHR)]); # annotation object
    plas.features.append(annLHR); # adds annotation

    if len(recodedRegion) > 0: # if there is a recoded region,
        inRecode = annLHR.index[1]; # index of recoded region start site (at end of LHR)
        plas.insertSeq(recodedRegion, inRecode); # inserts recoded region sequence
        annRecoded = GenBankAnn(geneName+" Recoded Region", "misc_feature", recodedRegion, False, [inRecode,inRecode+len(recodedRegion)]); # annotation object
        plas.features.append(annRecoded); # adds annotation

    startGRNA = findFirst(plas.origin, cut_IPpoI); # index of gRNA start site (at start of I-PpoI cut sequence)
    endGRNA = startGRNA + len(cut_IPpoI); # index of gRNA end site (at end of I-PpoI cut sequence)
    plas.removeSeq([startGRNA, endGRNA]); # removes sequence that gRNA will replace
    plas.insertSeq("gg" + gRNA, startGRNA); # inserts gRNA sequence with gg sequence used by T7 polymerase
    annGRNA = GenBankAnn(geneName+" gRNA", "misc_feature", gRNA, False, [startGRNA+2,startGRNA+2+len(gRNA)]); # annotation object. Note that gRNA starts after "gg" added for T7 polymerase
    plas.features.append(annGRNA); # adds annotation

    inRHR = findFirst(plas.origin, cut_ISceI); # index of RHR insertion site (at start of I-SceI cut sequence)
    plas.insertSeq(RHR, inRHR); # inserts RHR sequence
    annRHR = GenBankAnn(geneName+" RHR", "misc_feature", RHR, False, [inRHR,inRHR+len(RHR)]); # annotation object
    plas.features.append(annRHR); # adds annotation

    return plas; # returns modified plasmid


"""
Selects an appropriate gRNA for the given gene. GenBank gene sequence given as
argument must have an annotation with "gRNA 1" as a part of its label, and the
gene must be at least searchRange[0] bp long. The file must include at least
searchRange[1] bp of 3' UTR. searchRange is a list indicating search start and
end indexes, counted with the last bp in the gene's stop codon as index 0.
side3Prime is false if PAM sequence is at the 5' end of gRNA, true if at 3'.
gRNA: guide RNA used by CRISPR enzyme.
"""
def chooseGRNA(geneGB, gene, searchRange=[-500,125], PAM="NGG", side3Prime=True, minGCContent=0.3, minOnTargetScore=25, minOffTargetScore=75, maxOffTargetHitScore=35, onTargetMethod="azimuth", offTargetMethod="hsu", gLength=20, maxDistanceBetweenGRNAS=50, enzyme="Cas9", gBlockDefault=True, maxTier1GBlockSize="250", filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII]): # could've been useful at some point: http://grna.ctegd.uga.edu/ http://www.broadinstitute.org/rnai/public/software/sgrna-scoring-help http://crispr.mit.edu/about
    log = "Choosing gRNA with PAM sequence " + PAM + " for use with enzyme " + enzyme; # init log
    gRNATable = []; # will store information on each gRNA evaluated. Format: Label, Status, Enzyme, Position, Strand, GC_content, On-target_score, On-target_method, Aggregated_off-target_score, Max_pairwise_off-target_score, Off-target_method, >9_consecutive_A/T, 4-Homopolymer, Triple_T, Sequence, Recoded_sequence, Recoded_sequence_pairwise_off-target_score
    geneList = geneGB.findAnnsLabel(gene.label); # stores gene annotation
    gene = geneList[0]; # stores gene annotation
    for g in geneList: # search for right gene
        if g.type == "gene": # if found,
            gene = g; # save it


    for g in geneList: # search for right gene
        if g.label == gene.label and g.type == "gene": # if found,
            gene = g; # save it


    gRNAs = []; # list of GenBankAnn objects for gRNAs.
    backupGRNAs = []; # stores gRNAs that fail off-target score as possible backups
    gRNAUpstream = GenBankAnn(); # will store gRNA most upstream

    searchRange[0] = max(searchRange[0],-1*len(geneGB.origin)); # adjusts search range in case it goes beyond the total length of the gene

    pamSeqs = ambiguousSeqs(PAM); # store actual PAM sequences
    extGRNASeqIndexes = []; # stores indexes of extended gRNA sequence (NNNN-gRNA 20mer-PAM 3mer-NNN for Cas9, PAM 4mer-gRNA 23mer-NNN for Cpf1) relative to PAM start point.
    realGRNASeqIndexes = []; # stores indexes of actual gRNA sequence (gRNA 20mer for Cas9, gRNA 23mer for Cpf1) relative to PAM start point.
    pamIndexes = []; # stores indexes of PAM sequences within extended gRNA

    if enzyme == "Cas9": # if enzyme is Cas9,
        extGRNASeqIndexes = [-24,6,-3,27]; # Start and end indexes for NNNN-gRNA 20mer-PAM 3mer-NNN extended sequence relative to PAM start position. Second two numbers are for rev comp strand.
        realGRNASeqIndexes = [4,24]; # Start and end indexes for gRNA 20mer sequence within extended sequence.
        pamIndexes = [24,24+len(PAM)]; # Stores PAM seq indexes, including Ns and Vs
    elif enzyme == "Cpf1": # if enzyme is Cpf1,
        extGRNASeqIndexes = [0,30,-26,4]; # Start and end indexes for PAM 4mer-gRNA 23mer extended sequence relative to PAM start position. Second two numbers are for rev comp strand.
        realGRNASeqIndexes = [len(PAM),23+len(PAM)]; # Start and end indexes for PAM gRNA 23mer sequence within extended sequence.
        pamIndexes = [0,len(PAM)]; # Stores PAM seq indexes, including Ns and Vs
        # assign values

    if len(gene.label) > 0: # if gene found,
        searchStart = gene.index[side3Prime]+searchRange[0]; # Start searching for gRNAs in a position relative to gene start or end point
        searchEnd = gene.index[side3Prime]+searchRange[1]; # Finish searching for gRNAs in a position relative to gene start or end point
        searchSeq = geneGB.origin[searchStart:searchEnd].upper(); # get sequence where gRNAs will be searched for, centered around start or end of gene
        i = len(searchSeq)-extGRNASeqIndexes[1]+1; # start searching for gRNAs by searching for PAM in searchSeq. Start at downstream end (allow five bases downstream to accomodate on-target 30-mer gRNA sequence), go upstream
        while i > len(PAM) and len(gRNAs) < 3: # iterate through all searchSeq until done with searchSeq (allow enough bases upstream of search end point to accomodate on-target 30-mer minus-strand gRNA sequence) or three candidates found.
            comp = False; # set sense to plus strand
            extGRNASeq = ""; # stores extended sequence sequence (NNNN-gRNA 20mer-PAM 3mer-NNN for Cas9, PAM 4mer-gRNA 23mer for Cpf1)
            gRNAIndexes = []; # store gRNA indexes
            if searchSeq[i:i+len(PAM)] in pamSeqs: # if PAM found on plus strand,
                extGRNASeq = searchSeq[i+extGRNASeqIndexes[0]:i+extGRNASeqIndexes[1]]; # store extended sequence
                gRNAIndexes = [searchStart+i+extGRNASeqIndexes[0]+realGRNASeqIndexes[0],searchStart+i+extGRNASeqIndexes[0]+realGRNASeqIndexes[1]]; # store gRNA indexes
            elif searchSeq[i:i+len(PAM)] in [revComp(p) for p in pamSeqs]: # if PAM found on minus strand, +extGRNASeqIndexes[0]+realGRNASeqIndexes[0]
                extGRNASeq = revComp(searchSeq[i+extGRNASeqIndexes[2]:i+extGRNASeqIndexes[3]]); # store extended sequence on comp strand
                gRNAIndexes = [searchStart+i+extGRNASeqIndexes[3]-realGRNASeqIndexes[1],searchStart+i+extGRNASeqIndexes[3]-realGRNASeqIndexes[0]]; # store gRNA indexes on comp strand
                comp = True; # set sense to complementary strand

            if len(extGRNASeq) == extGRNASeqIndexes[1]-extGRNASeqIndexes[0]: # if extended gRNA is right size (doesn't overstep boundaries)
                pamSeq = extGRNASeq[pamIndexes[0]:pamIndexes[1]]; # stores PAM sequence
                gRNASeq = extGRNASeq[realGRNASeqIndexes[0]:realGRNASeqIndexes[1]]; # store actual gRNA seq without PAM
                gc = gcContent(gRNASeq); # store gc content
                strandString = '+'; # stores string denoting strand orientation
                if comp: # if on opposite strand,
                    strandString = '-'; # save that info

                gRNATable.append(['Unlabeled','Rejected (low GC content)',enzyme,str(gRNAIndexes).replace(',',' to'),strandString,str(gc),'Not evaluated',onTargetMethod,'Not evaluated','Not evaluated',offTargetMethod,str(findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") > -1),'Not evaluated','Not evaluated',gRNASeq,'Not recoded','-']); # starts storing info
                if gc >= minGCContent and findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") < 0: # if gc content is acceptable and does not contain 10 or more consecutive As or Ts,
                    gRNAInfo = getGRNAInfoFromDB(extGRNASeq,enzyme); # access gRNA scores from DB
                    onTarget = 0;
                    if len(gRNAInfo) == 0: # if not found in DB
                        onTarget = onTargetScore(extGRNASeq,onTargetMethod); # store on-target score
                    else: # if found in DB,
                        if onTargetMethod == "azimuth": # if using Azimuth (Doench et al., 2016)
                            onTarget = gRNAInfo[1]; # use it
                        elif onTargetMethod == "ruleset2": # if using Rule Set 2 (Doench et al., 2016)
                            onTarget = onTargetScore(extGRNASeq,onTargetMethod); # store on-target score
                        elif onTargetMethod == "cindel": # if using CINDEL (Kim et al., 2017)
                            onTarget = gRNAInfo[2]; # use it


                    gRNATable[len(gRNATable)-1][1] = 'Rejected (low on-target score)'; # Edit this gRNA's status
                    gRNATable[len(gRNATable)-1][6] = str(onTarget); # Edit this gRNA's on-target score
                    if onTarget > minOnTargetScore: # if on-target score is passable,
                        offTargetScores = [];
                        if len(gRNAInfo) == 0: # if not found in DB
                            offTargetScores = offTargetScore(gRNASeq,offTargetMethod,enzyme,pamSeq,PAM); # store off-target score
                        else: # if found in DB,
                            if offTargetMethod == "hsu": # if using Hsu
                                offTargetScores = gRNAInfo[3:6]; # use it
                            elif offTargetMethod == "cfd": # if using CFD
                                offTargetScores = gRNAInfo[9:12]; # use it


                        newGRNA = GenBankAnn(gene.label + " " + enzyme + " gRNA ","misc",gRNASeq,comp,gRNAIndexes); # create annotation object with gRNA information
                        newGRNA.onTarget = onTarget; # Add on-target score as attribute
                        newGRNA.offTarget = offTargetScores; # Add off-target scores as attribute
                        newGRNA.gc = gc; # Add gc content as attribute
                        newGRNA.homopolymer = (findFirst(gRNASeq,"AAAA") > 0 or findFirst(gRNASeq,"TTTT") > 0 or findFirst(gRNASeq,"CCCC") > 0 or findFirst(gRNASeq,"GGGG") > 0 ); # 4-homopolymers are bad https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1581-4
                        newGRNA.tripleT = (findFirst(gRNASeq,"TTT") > 0); # Triple Ts (triple Us) are bad because they're an RNApol stop codon https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1581-4

                        gRNATable[len(gRNATable)-1][1] = 'Possible backup (low aggregate off-target)'; # Edit this gRNA's status
                        gRNATable[len(gRNATable)-1][8] = str(offTargetScores[0]); # Edit this gRNA's aggregate off-target score
                        gRNATable[len(gRNATable)-1][9] = str(offTargetScores[1]); # Edit this gRNA's max pairwise off-target score
                        gRNATable[len(gRNATable)-1][12] = str(newGRNA.homopolymer); # Edit this gRNA's homopolymer
                        gRNATable[len(gRNATable)-1][13] = str(newGRNA.tripleT); # Edit this gRNA's tripleT
                        gRNATable[len(gRNATable)-1][15] = 'No recoded region necessary'; # Edit this gRNA's recoded status

                        if offTargetScores[0] >= minOffTargetScore and offTargetScores[1] <= maxOffTargetHitScore and not newGRNA.homopolymer and not newGRNA.tripleT: # if total off-target score and max hit score are passable, and if no homopolymer or triple T
                            gRNAs.append(newGRNA); # add this gRNA's information to list.
                            gRNATable[len(gRNATable)-1][1] = 'Valid'; # Edit this gRNA's status
                        else: # if failed off-target score culling.
                            newGRNA.label = newGRNA.label + "(backup)"; # notify backup on label
                            backupGRNAs.append(newGRNA); # add this gRNA's information to backup list.
                            if offTargetScores[1] > maxOffTargetHitScore: # if failed max pairwise off-target check,
                                gRNATable[len(gRNATable)-1][1] = 'Possible backup (high max pairwise off-target)'; # Edit this gRNA's status
                            elif newGRNA.homopolymer: # if failed max pairwise off-target check,
                                gRNATable[len(gRNATable)-1][1] = 'Possible backup (4-homopolymer)'; # Edit this gRNA's status
                            elif newGRNA.tripleT: # if failed max pairwise off-target check,
                                gRNATable[len(gRNATable)-1][1] = 'Possible backup (triple T)'; # Edit this gRNA's status



                elif gc > minGCContent : # if rejected due to low gc,
                        gRNATable.append(['Unlabeled','Rejected (>9 consecutive A/Ts)',enzyme,str(gRNAIndexes).replace(',',' to'),strandString,str(gc),'Not evaluated',onTargetMethod,'Not evaluated','Not evaluated',offTargetMethod,str(findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") > -1),'Not evaluated','Not evaluated',gRNASeq,'Not recoded','-']); # starts storing info






            i -= 1; # advance indexer

        inUTR = False; # will be used to explain status in table

        if len(gRNAs) == 0: # if no gRNAs found,
            log = log + "\nWarning: no acceptable gRNA with at least " + str(minGCContent*100) + "% GC content, " + str(minOnTargetScore) + " on-target score, " + str(minOffTargetScore) + " off-target total score, no 4-homopolymer sequences, and no TTT sequences found on gene " + gene.label + ".\n" + "Will use backup gRNA with highest GC content, if there are any backups."; # add warning to log
        else: # if gRNAs found,
            newList = []; # new list will only contain gRNAs within acceptable range
            if gRNAs[0].index[0] > gene.index[1]-3: # if most downstream gRNA starts after start of stop codon,
                inUTR = True; # note downstream gRNA is in UTR
                for g in gRNAs: # loop through gRNAs
                    if g.index[0] > gene.index[1]-3: # if this gRNA is still in the 3' UTR or stop codon,
                        newList.append(g); # add it to the new list



            else: # if most downstream gRNA is upstream of stop codon,
                for g in gRNAs: # loop through gRNAs
                    if (gRNAs[0].index[0] - g.index[0] <= maxDistanceBetweenGRNAS) or (gBlockDefault and gene.index[1]-3 - g.index[0] < maxTier1GBlockSize): # if this gRNA is within the max distance from the most downstream gRNA or gBlocks are the default and the gRNA is within the cheapest gBlock range from the end of the gene
                        newList.append(g); # add it to the new list




            gRNAs = newList; # keep the new list of gRNAs as the main list
            gRNAs.sort(reverse=True,key=lambda g: g.gc); # sorts gRNA list according to custom function (GC content)
            count = 1; # counter for numbering candidates according to their quality
            for g in gRNAs: # loop thorugh gRNAs ordered by GC content
                g.label = g.label + str(count); # add number to gRNA label
                count += 1; # advance counter

            bestGRNA = gRNAs[0]; # will store best gRNA
            gRNAUpstream = gRNAs[0]; # will find most upstream gRNA
            for g in gRNAs: # loop through gRNAs
                if g.index[0] < gRNAUpstream.index[0]: # if more upstream than previous most upstream
                    gRNAUpstream = g; # set this gRNA as most upstream


            log = log + "\n" + str(len(gRNAs)) + " acceptable gRNAs were selected automatically on gene " + gene.label + ". \ngRNA 1 has GC content of " + str(bestGRNA.gc*100) + "%, on-target score of " + str(bestGRNA.onTarget)  + ", \nand aggregated off-target score of " + str(bestGRNA.offTarget[0]) + " (Method: " + offTargetMethod + ", Max. Hit Score: " + str(bestGRNA.offTarget[1]) + ", Num. hits: " + str(bestGRNA.offTarget[2]) + ")."; # add warning to log


    geneGB.features = geneGB.features + gRNAs; # add gRNAs to gene GenBank object features list
    countBackups = 0; # counts how many backup gRNAs are included
    allGRNAS = gRNAs; # will store both valid and backup gRNAs
    if len(gRNAs) > 0: # if there is at least one gRNA
        for g in backupGRNAs: # for every possible backupGRNAs
            if g.index[0] > gRNAUpstream.index[0]: # if this backup is downstream of most upstream gRNA,
                geneGB.features.append(g); # add to features list
                allGRNAS.append(g); # add to full gRNA list
                countBackups +=1; # advances counter


    else: # if there are no gRNAs,
        maxGC = 0; # will track max gc content of gRNA
        for g in backupGRNAs: # for every possible backupGRNAs
            geneGB.features.append(g); # add to features list
            allGRNAS.append(g); # add to full gRNA list
            countBackups +=1; # advances counter
            if g.gc >= maxGC: # if this gRNA has a greater or equal GC content,
                gRNAUpstream = g; # set as most upstream gRNA



    for gArr in gRNATable: # loop over all gRNAs in table
        found = False; # determine whether the gRNA is in labeled list
        for g in allGRNAS: # loop over list of labeled gRNAs
            if gArr[14] == g.seq: # if not included in annotations
                found = True; # set found to True
                break; # break inner loop


        if not found: # if not included in annotations
            gArr[15] = 'Not recoded'; # Edit this gRNA's recoded status
            if gArr[1][0:8] == "Rejected": # if status is rejected,
                gArr[1] = gArr[1] + " unlabeled (due to rejection)"; # Add unlabeled explanation to status
            else: # if not rejected,
                if inUTR: # if downstream gRNA in UTR
                    gArr[1] = gArr[1] + " unlabeled (lead gRNA in UTR and this one in gene)"; # Add unlabeled explanation to status
                else: # if downstream gRNA in gene,
                    gArr[1] = gArr[1] + " unlabeled (more than " + str(maxDistanceBetweenGRNAS) + " bp upstream of lead gRNA)"; # Add unlabeled explanation to status




    for g in allGRNAS: # for every gRNA annotated
        for gArr in gRNATable: # find this gRNA in table
            if g.seq == gArr[14]: # if the sequence is the same,
                gArr[0] = g.label.replace(',',' '); # add this label to this row




    gRNATableNew = []; # will store new ordered list
    for g in gRNATable: # loop over list,
        if g[1][0:5] == "Valid": # if status is valid
            gRNATableNew.append(g); # add to new list


    gRNATableNew = sorted(gRNATableNew,key=operator.itemgetter(0)); # sort according to name (gRNA priority)
    for g in gRNATable: # loop over list,
        if g[1][0:14] == "Kept as backup": # if status is Kept as backup
            gRNATableNew.append(g); # add to new list


    for g in gRNATable: # loop over list,
        if g[1][0:15] == "Possible backup": # if status is Possible backup
            gRNATableNew.append(g); # add to new list


    for g in gRNATable: # loop over list,
        if g[1][0:23] == "Rejected (low on-target": # if status is Rejected due to on-target
            gRNATableNew.append(g); # add to new list


    for g in gRNATable: # loop over list,
        if g[1][0:12] == "Rejected (>9": # if status is Rejected due to >9 consecutive A/Ts
            gRNATableNew.append(g); # add to new list


    for g in gRNATable: # loop over list,
        if g[1][0:16] == "Rejected (low GC": # if status is Rejected due to gc content
            gRNATableNew.append(g); # add to new list


    gRNATableString = "\n".join([",".join(g) for g in gRNATableNew]); # join array into csv string
    gRNATableString = "Values for rejected gRNAs, Rejected, "+enzyme+", "+str(searchRange[0]+gene.index[side3Prime])+" to "+str(searchRange[1]+gene.index[side3Prime])+", +/-, <" + str(minGCContent) + ", <"+str(minOnTargetScore)+", "+onTargetMethod+", Not evaluated, Not evaluated, -, True, Not evaluated, Not evaluated, -, Not recoded, -\n" + gRNATableString; # Add rejected threshold
    gRNATableString = "Values for backup gRNAs, Backup, "+enzyme+", "+str(searchRange[0]+gene.index[side3Prime])+" to "+str(searchRange[1]+gene.index[side3Prime])+", +/-, >=" + str(minGCContent) + ", >="+str(minOnTargetScore)+", "+onTargetMethod+", <"+str(minOffTargetScore)+", <"+str(maxOffTargetHitScore)+", "+offTargetMethod+", False, True, True, -, Recoded if upstream of stop codon, >=threshold\n" + gRNATableString; # Add backup threshold
    gRNATableString = "Values for valid gRNAs, Valid, "+enzyme+", "+str(searchRange[0]+gene.index[side3Prime])+" to "+str(searchRange[1]+gene.index[side3Prime])+", +/-, >=" + str(minGCContent) + ", >="+str(minOnTargetScore)+", "+onTargetMethod+", >="+str(minOffTargetScore)+", >="+str(maxOffTargetHitScore)+", "+offTargetMethod+", False, False, False, -, Recoded if upstream of stop codon, >=threshold\n" + gRNATableString; # Add valid threshold
    gRNATableString = "Label, Status, Enzyme, Position, Strand, GC_content, On-target_score, On-target_method, Aggregated_off-target_score, Max_pairwise_off-target_score, Off-target_method, >9_consecutive_A/T, 4-Homopolymer, Triple_T, Sequence, Recoded_sequence, Recoded_sequence_pairwise_off-target_score\n" + gRNATableString; # Add column heads

    log = log + "\n" + str(countBackups) + " backup gRNAs with possible off-target effects annotated.\n\n"
    if len(backupGRNAs) + len(gRNAs) == 0: # If there were absolutely no gRNAs under these settings,
        log = log + "\n" + "ERROR: no gRNAs found. Please modify your criteria or select and annotate one manually.\n\n"; # say so

    return {"out":gRNAUpstream, "log":log, "gRNATable":gRNATableString}; # returns gRNA and log


"""
Finds gRNA already annotated on gene. Filters for given restriction enzyme cut
sites. If no gRNA found on gene, nothing will happen (logs this).
"""
def findGRNA(geneGB, gene, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII]):
    log = ""; # init log
    gRNAs = geneGB.findAnnsLabel("gRNA 1", True); # List of all gRNAs
    gRNAUpstream = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA1", True); # List of all gRNAs
        if len(gRNAs) == 0: # if no gRNAs found
            gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
            if len(gRNAs) == 0: # if no gRNAs found
                log = log + "\nWarning: no gRNA found on gene " + gene.label + ", selecting one automatically.\n"; # add warning to log
            else: # if gRNAs found,
                gRNAUpstream = gRNAs[0]; # will store gRNA most upstream
                for g in gRNAs: # loops across gRNAs
                    if g.index[0] < gRNAUpstream.index[0]: # if more upstream
                        gRNAUpstream = g; # replace as most upstream





    if len(gRNAs) > 0: # if gRNAs found
        gRNAUpstream = gRNAs[0]; # will store gRNA most upstream
        gRNAUpstream = deepcopy(gRNAUpstream); # fixes referencing issue. We want this to be a genuinenly new annotation
        for site in filterCutSites: # for every cut site being filtered
            if findFirst(site,gRNAUpstream.seq) > -1 or findFirst(revComp(site),gRNAUpstream.seq) > -1: # if cut site found,
                log = log + "\nWarning: gRNA sequence for gene " + gene.label + ": \n" + gRNAUpstream + "\ncontains restriction site " + site + "\n"; # add warning to log


        gRNAUpstream.label = gene.label + " gRNA"; # renames gRNA according to this program's convention

        log = log + "gRNA for gene " + gene.label + " found on gene.\n\n"; # logs this process finished
        gRNATableString = "gRNAs not evaluated if they are user-defined.\nIf you want to check their scores, run the gene in automatic mode!\n"; # add disclaimer

    return {"out":gRNAUpstream, "log":log, "gRNATable":gRNATableString}; # returns gRNA and log



"""
Selects an appropriate LHR for the given gene. GenBank sequence object given as
argument (geneGB) must have an annotation with label=gene, at least one
annotation with gRNA in label (or some similar indicator), and the length of the
gene must be at least lengthLHR[2]. lengthLHR is [min, preferred, max] length in
bp of LHR. minTmEnds is min melting temp in extremes of LHR. endsLength is
length of extremes of the LHR contemplated in analysis. optimizeRange is a list
with the index positions relative to the start of the LHR over which the start
of the LHR is optimized once it is chosen. maxDistanceFromGRNA is maximum
distance in bp between end of LHR and start of gRNA.
LHR: Left Homologous Region used for chromosomal integration by homologous
recombination during repair.
"""
def chooseLHR(geneGB, gene, lengthLHR=[450,500,650], minTmEnds=55, endsLength=40, optimizeRange=[-20,10], maxDistanceFromGRNA=500, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII]):
    #TODO: debug cases in which LHR has to be shifted, do case in which LHR is in intron?
    log = ""; # init log
    gRNAs = []; # List of all gRNAs
    gRNAUpstream = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
        gRNAUpstream = gRNAs[0]; # will store gRNA most upstream
        for g in gRNAs: # loops across gRNAs
            if g.index[0] < gRNAUpstream.index[0]: # if more upstream
                gRNAUpstream = g; # replace as most upstream


    else: # if gRNAs found,
        gRNAUpstream = gRNAs[0]; # will store gRNA most upstream

    endLHR = min(gRNAUpstream.index[0],gene.index[1]-3); # saves end index of LHR as whatever is more upstream between the start of the gRNA or the end of the gene (minus the stop codon) (in Python indexing, i.e. not included in LHR).
    if gBlockDefault and endLHR < gene.index[1]-3 and endLHR > minGBlockSize: # if defaulting to a gBlock for any recoded region,
        endLHR = gene.index[0]-3 - minGBlockSize; # extend recoded region to minimum gBlock size

    while not geneGB.checkInExon(endLHR) and endLHR > lengthLHR[2]: # Loop as long as the end of LHR is not in an exon and the end of the LHR is inside the max length
        endLHR -= 1; # shift LHR end upstream one bp

    startLHR = max(endLHR - lengthLHR[1],0); # stores LHR start index according to preferred length, or 0 if not enough space before gene
    failSearchEnd = False; # true if no suitable LHR end region is found within given parameters.

    # First search for end region
    while meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchEnd and geneGB.checkInExon(endLHR-1): # while no suitable end region found, still within max distance from gRNA, the search for a suitable end region hasn't failed yet, and still within exon,
        endLHR -= 1; # shift endLHR upstream

    if meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) < minTmEnds: # if no suitable end region found,
        endLHR = min(gRNAUpstream.index[0],gene.index[1]-3); # saves end index of LHR by default as whatever is more upstream between the start of the gRNA or the end of the gene (minus the stop codon)
        while not geneGB.checkInExon(endLHR) and endLHR > lengthLHR[2]: # Loop as long as the end of LHR is not in an exon and the end of the LHR is inside the max length
            endLHR -= 1; # shift LHR end upstream one bp

        failSearchEnd = True; # changes failSearchEnd status
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of LHR and gRNA. \nDefaulted to ending right before start of gRNA most upstream." + "\n"; # give a warning

    # Then search for start region; modify end region if necessary
    while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and (gRNAUpstream.index[0]-endLHR <= maxDistanceFromGRNA and geneGB.checkInExon(endLHR-1)): # if no start region has been found, and still within max distance from gRNA and inside exon,
        # find new end region if necessary
        while meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchEnd and geneGB.checkInExon(endLHR-1): # while no suitable end region found, still within max distance from gRNA, the search for a suitable end region hasn't failed yet, and still within exon,
            endLHR -= 1; # shift endLHR upstream

        # search for starts upstream
        while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR <= lengthLHR[2] and startLHR > 0: # while no suitable start region found and still within max length of LHR and inside gene file,
            startLHR -= 1; # shift startLHR upstream

        # if not found upstream
        if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR > lengthLHR[2]: # if no start region found this way,
            startLHR = endLHR - lengthLHR[1]; # return to preferred position
            # search downstream
            while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR >= lengthLHR[0]: # while no suitable start region found and still within min length of LHR,
                startLHR += 1; # shift startLHR downstream

        if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and geneGB.checkInExon(endLHR-1): # if still not found and still inside exon,
            endLHR -= 1; # shifts end of LHR upstream



    if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds: # if no suitable start region found,
        endLHR = gRNAUpstream.index[0]; # resets endLHR
        while not geneGB.checkInExon(endLHR) and endLHR > lengthLHR[2]: # Loop as long as the end of LHR is not in an exon and the end of the LHR is inside the max length
            endLHR -= 1; # shift LHR end upstream one bp

        while meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchEnd and geneGB.checkInExon(endLHR-1): # while no suitable end region found, still within max distance from gRNA, the search for a suitable end region hasn't failed yet, and still within exon,
            endLHR -= 1; # shift endLHR upstream

        startLHR = endLHR - lengthLHR[1]; # stores LHR start index according to preferred length by default
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of LHR and gRNA. \nDefaulted to starting so that the LHR size is equal to preferred size of " + str(lengthLHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting LHR by adjusting start sequence around the given range
    searchIndexes = [int(startLHR+optimizeRange[0]), int(startLHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexes[0], searchIndexes[1]): # iterates across optimization range
        if i >= endsLength: # if inside gene file,
            if meltingTemp(geneGB.origin[i:(i+endsLength)]) > meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) and lengthLHR[0] >= endLHR-i >= lengthLHR[2]: # if this start point has a better Tm and is still within bounds,
                startLHR = i; # make this the starting position



    searchIndexesEnd = [int(endLHR+optimizeRange[0]), int(endLHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better end point.
    for i in range(searchIndexesEnd[0], searchIndexesEnd[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[(i-endsLength):i]) > meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) and lengthLHR[2] >= i-startLHR >= lengthLHR[0] and i < gRNAUpstream.index[0] and geneGB.checkInExon(i): # if this start point has a better Tm, and is still within bounds, before gRNA and within exon,
            endLHR = i; # make this the ending position




    for site in filterCutSites: # for every cut site being filtered
        if findFirst(site,geneGB.origin[startLHR:endLHR]) > -1 or findFirst(revComp(site),geneGB.origin[startLHR:endLHR]) > -1: # if cut site found,
            log = log + "\nWarning: LHR sequence for gene " + gene.label + ": \n" + geneGB.origin[startLHR:endLHR] + "\ncontains restriction site " + site + "\n"; # add warning to log


    # creates var to store finished LHR as annotation
    LHR = GenBankAnn(gene.label + " LHR", "misc_feature", geneGB.origin[startLHR:endLHR], False, [startLHR,endLHR]); # creates GenBankAnn object to hold LHR

    log = log + "LHR for gene " + gene.label + " selected.\n\n"; # logs this process finished
    return {"out":LHR, "log":log}; # returns LHR GenBankAnn object


"""
Selects an appropriate RHR for the given gene. GenBank sequence object given as
argument (geneGB) must have an annotation with label=gene, at least one
annotation with "gRNA" in label (or some similar indicator), and at least
lengthRHR[2] + maxDistanceFromGene - optimizeRange[0] + optimizeRange[1]
bp of 3' UTR (about 1000 bp is fine). lengthRHR is [min, preferred, max]
length in bp of RHR. minTmEnds is min melting temp in extremes of RHR.
endsLength is length of extremes of the RHR contemplated in analysis.
optimizeRange is a list with the index positions relative to the start of the
RHR over which the start of the LHR is optimized once it is chosen.
maxDistanceFromGRNA is maximum distance in bp between end of RHR and start of
gRNA.
RHR: Right Homologous Region used for chromosomal integration by homologous
recombination during repair.
"""
def chooseRHR(geneGB, gene, lengthRHR=[450,500,750], minTmEnds=59, endsLength=40, optimizeRange=[-20,20], maxDistanceFromGene=500, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII]):
    #TODO: debug cases in which RHR has to be shifted
    log = ""; # init log
    gRNAs = []; # List of all gRNAs
    gRNADownstream = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
        gRNADownstream = gRNAs[0]; # will store gRNA most downstream
        for g in gRNAs: # loops across gRNAs
            if g.index[0] > gRNADownstream.index[0]: # if more downstream
                gRNADownstream = g; # replace as most downstream


    else: # if gRNAs found,
        gRNADownstream = gRNAs[0]; # will store gRNA most downstream

    genes = geneGB.findAnnsType("gene"); # list of all genes, used to verify RHR doesn't start inside any genes (truncating them)
    def checkInGene(pIndex): # checks whether the index given is inside a gene.
        inGene = False; # default is false
        for g in genes: # loop through all genes
            inGene = g.index[0] <= pIndex < g.index[1]; # if inside gene, set to True.

        return inGene; # return

    startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or after most downstream gRNA.
    endRHR = min(startRHR + lengthRHR[1],len(geneGB.origin)-1); # stores RHR start index according to preferred length, or the last base in file if not enough space after gene
    failSearchStart = False; # true if no suitable RHR start region is found within given parameters.

    # First search for start region
    while meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) < minTmEnds and startRHR - gene.index[1] <= maxDistanceFromGene and startRHR + lengthRHR[0] < gene.index[1] and not checkInGene(startRHR): # while no suitable start region found, still within max distance from gene, not beyond the border of the chromosome, and the search for a suitable start region hasn't failed yet,
        startRHR += 1; # shift startRHR downstream

    if meltingTemp(geneGB.origin[startRHR:(min(startRHR+endsLength,len(geneGB.origin))+endsLength)]) < minTmEnds and ( startRHR - gene.index[1] > maxDistanceFromGene or checkInGene(startRHR) ): # if no suitable start region found,
        startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or 1 bp after most downstream gRNA by default
        failSearchStart = True; # changes failSearchStart status
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of gene and start of RHR. \nDefaulted to starting right after end of gene." + "\n"; # give a warning

    # search for end downstream
    while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2] and endRHR < len(geneGB.origin)-1: # while no suitable end region found, still within max length of RHR, and still inside gene file,
        endRHR += 1; # shift endRHR downstream

    # if not found downstream
    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
        endRHR = startRHR + lengthRHR[1]; # return to preferred position
        # search upstream
        while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR >= lengthRHR[0]: # while no suitable start region found and still within min length of RHR,
            endRHR -= 1; # shift endRHR downstream


    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and checkInGene(startRHR): # if still not found and within bounds,
        # Then modify start region, search for end region;
        while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and not checkInGene(startRHR): # if no end region has been found,
            startRHR += 1; # shifts start of RHR downstream

            # find new start region if necessary
            while meltingTemp(geneGB.origin[startRHR:(startRHR+endsLength)]) < minTmEnds and startRHR - gene.index[1] <= maxDistanceFromGene and checkInGene(startRHR) and not failSearchStart: # while no suitable start region found, still within max distance from gene, and the search for a suitable start region hasn't failed yet,
                startRHR += 1; # shift startRHR downstream

            # search for end downstream
            while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2]: # while no suitable end region found and still within max length of RHR,
                endRHR += 1; # shift endRHR downstream

            # if not found downstream
            if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
                endRHR = startRHR + lengthRHR[1]; # return to preferred position
                # search upstream
                while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR >= lengthRHR[0]: # while no suitable start region found and still within min length of RHR,
                    endRHR -= 1; # shift endRHR downstream



    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and ( startRHR-gene.index[1] > maxDistanceFromGene or checkInGene(startRHR) ): # if no suitable end region found,
        startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or 1 bp after most downstream gRNA by default
        while meltingTemp(geneGB.origin[startRHR:(startRHR+endsLength)]) < minTmEnds and startRHR-gene.index[0] <= maxDistanceFromGene and not failSearchStart: # while no suitable start region found, still within max distance from gene, and the search for a suitable start region hasn't failed yet,
            startRHR += 1; # shift startRHR downstream

        endRHR = min([startRHR + lengthRHR[1],len(geneGB.origin)]); # stores RHR end index according to preferred length by default, or max length available if there is not enough sequence space on the chromosome after the gene end
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of RHR and end of gene. \nDefaulted to ending at preferred size of " + str(lengthRHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting RHR by adjusting start and stop sequences around the given range
    searchIndexesEnd = [int(endRHR+optimizeRange[0]), int(endRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better end point.
    for i in range(searchIndexesEnd[1], searchIndexesEnd[0], -1): # iterates across optimization range in reverse
        if meltingTemp(geneGB.origin[ min([i-endsLength,len(geneGB.origin)-endsLength]):min([i,len(geneGB.origin)]) ]) > meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) and lengthRHR[2] >= i-startRHR >= lengthRHR[0]: # if this end point has a better Tm and is still within bounds,
            if i < len(geneGB.origin): # if inside gene file,
                endRHR = i; # make this the ending position




    searchIndexesStart = [int(startRHR+optimizeRange[0]), int(startRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexesStart[0], searchIndexesStart[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[i:(i+endsLength)]) > meltingTemp(geneGB.origin[startRHR:(startRHR+endsLength)]) and lengthRHR[2] >= endRHR-i >= lengthRHR[0] and i >= gRNADownstream.index[1] and not checkInGene(i): # if this start point has a better Tm and is still within bounds and after gRNA and not inside a gene,
            startRHR = i; # make this the starting position


    if endRHR-startRHR < lengthRHR[0]+searchIndexesEnd[0]-searchIndexesEnd[1]: # if RHR is smaller than it should be (gene too close to the end of chromosome, for example)
        log = log + "\ERROR: RHR found for gene " + geneGB.name + " \nis smaller than allowed by design, with a length of " + str(endRHR-startRHR) + " bp. Maybe gene was too close to start/end of chromosome?\n"; # give a warning

    for site in filterCutSites: # for every cut site being filtered
        if findFirst(site,geneGB.origin[startRHR:endRHR]) > -1 or findFirst(revComp(site),geneGB.origin[startRHR:endRHR]) > -1: # if cut site found,
            log = log + "\nWarning: RHR sequence for gene " + gene.label + ": \n" + geneGB.origin[startRHR:endRHR] + "\ncontains restriction site " + site + "\n"; # add warning to log

    # creates var to store finished RHR as annotation
    RHR = GenBankAnn(gene.label + " RHR", "misc_feature", geneGB.origin[startRHR:endRHR], False, [startRHR,endRHR]); # creates GenBankAnn object to hold RHR

    log = log + "RHR for gene " + gene.label + " selected.\n\n"; # logs this process finished

    return {"out":RHR, "log":log}; # returns RHR GenBankAnn object


"""
Chooses the region to be recoded to avoid gRNA targeting in already transfected
regions. Returns GenBankAnn object with recoded sequence and indexes between
which it should go. GenBank object given as argument should contain one gene
with geneName included in its label, and at least one annotation with "LHR" in
its label. Also needs all gRNAs to be annotated in the file. Returns empty
region if LHR end is at or downstream of gene stop codon. Checks against
restriction sites given as parameters. Checks that gRNA recoded sequence has a
pairwise off-target score lower than the given threshold with respect to the
original gRNA.
"""
def chooseRecodeRegion(geneGB, gene, offTargetMethod="cfd", pamType="NGG", orgCodonTable=codonUsage(), filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII], codonSampling=False, offScoreThreshold=10, minGCEnd5Prime=0.375, gRNATableString=""):
    #TODO: debug #TODO: Recoded if upstream of stop codon add recode values to table
    gRNAs = geneGB.findAnnsLabel("gRNA", True); # List of all gRNAs
    gRNATable = gRNATableString.split('\n'); # split string into lines
    gRNATable = [g.split(',') for g in gRNATable]; # split each line into values

    if offTargetMethod == "hsu": # if off-target scoring with Hsu scores
        offScoreThreshold = 1; # set threshold to 1%

    log = ""; # init log
    LHR = geneGB.findAnnsLabel("LHR")[0]; # LHR annotation object

    annRecoded = GenBankAnn(); # creates GenBankAnn object to hold recoded region
    if LHR.index[1] < gene.index[1]: # if end of LHR is inside gene
        startRecode = LHR.index[1]; # start of recode region (start of gRNA most upstream)
        while not geneGB.checkInExon(startRecode): # while recode region start is in intron,
            startRecode += 1; # shift downstream

        intronStartIndices = []; # stores start indexes of introns starting after recode sequence start
        intronEndIndices = []; # stores end indexes of introns starting after recode sequence start
        for ann in geneGB.features: # loop through all annotations
            if ann.type == "exon": # if annotation is exon
                if gene.index[1] > ann.index[1] > startRecode: # if annotation is an exon ending before gene end and after recode start,
                    intronStartIndices.append(ann.index[1]); # add this intron start index
                if gene.index[1] > ann.index[0] > startRecode: # if annotation is an exon starting after recode start,
                    intronEndIndices.append(ann.index[0]); # add this intron end index


            elif ann.type == "intron": # if annotation is intron,
                if ann.index[0] > startRecode: # if annotation is an intron starting after recode start,
                    intronStartIndices.append(ann.index[0]); # add this intron start index
                    intronEndIndices.append(ann.index[1]); # add this intron start index



        intronIndices = []; # will contain final indexes of introns downstream of recode start (introns to be removed from recoded region)
        intronStartIndices = sorted(intronStartIndices); # sort
        intronEndIndices = sorted(intronEndIndices); # sort
        for i in range(len(intronEndIndices)): # for every intron end,
            if intronEndIndices[i] >= LHR.index[1]: # if after LHR end,
                start = 0; # will store corresponding intron start (largest underneath intron end)
                for startIndex in intronStartIndices: # loop through starts
                    if startIndex > intronEndIndices[i]: # if start surpasses this end,
                        break; # stop loop
                    else: # if not,
                        start = startIndex; # set as start


                intronIndices.append([start,intronEndIndices[i]]); # add these coordinates to intron splice list
            #for g in gRNAs: # for every gRNA,
                #if not (intronStartIndices[i] >= g.index[1] or intronEndIndices[i] <= g.index[0]): # if intron and gRNA overlap,
                    #intronIndices.append([intronStartIndices[i],intronEndIndices[i]]); # add these coordinates to intron splice list



        endRecode = gene.index[1] - 3; # end of recode region (end of gene, exclude stop codon)
        recodeSeq = geneGB.origin[startRecode:endRecode]; # will contain sequence to be recorded
        if len(intronIndices) > 0: # if there are introns,
            recodeSeq = geneGB.origin[startRecode:intronIndices[0][0]]; # get recode sequence until first intron
            for i in range(len(intronIndices)-1): # for every intron except last one,
                recodeSeq = recodeSeq + geneGB.origin[intronIndices[i][1]:intronIndices[i+1][0]]; # add next exon to recode seq

            recodeSeq = recodeSeq + geneGB.origin[intronIndices[len(intronIndices)-1][1]:endRecode]; # get rest of recode sequence until endRecode

        frame = len(recodeSeq) % 3; # stores reading frame, index from start of sequence to be recoded
        startRecode += frame; # modify recode start site according to reading frame
        nonRecodedStart = recodeSeq[0:frame]; # stores 0, 1 or 2 nucleotides not recoded due to reading frame
        recodeSeq = recodeSeq[frame:len(recodeSeq)]; # adjust recode region

        cutSeqs = filterCutSites + [g.seq for g in gRNAs]; # list of all cut seqs. all gRNAs in gene are to be included as cut sequences
        cutCheck = 0; # variable used to check if cut sequences are present. Initially greater than -1*len(cutSeqs) since all gRNAs are present.
        offScore = 100; # stores off-target score. Default is 100% due to the fact that gRNA sequence is the same.
        count = 0; # iteration counter
        recodedSeq = recodeSeq; # assign recoded sequence to same as original
        bestRecodedSeq = recodedSeq; # will store best candidate sequence
        if len(recodeSeq) > 2: # if recodeSeq contains at least one codon,
            tricky = False; # True if suspected to be hard to synthesize
            badStart = False; # True if first bases have low melting temp (important for Gibson assembly)
            candidateFound = False; # signal possible candidate found
            bestRecodedSeq = recodedSeq; # will store best candidate sequence
            while cutCheck > -2*len(cutSeqs) or offScore > offScoreThreshold or tricky or badStart: # while cutCheck is greater than what you would expect for no hits in all cut sequences plus the gRNAs on both positive and comp strands, or while the pairwise off-target score is over the threshold, or while there are difficult-to-synthesize structures in the recoded region, or while the first 40 bp have a bad gc content
                if count > 0: # if recoded region has failed checks once,
                    codonSampling = True; # forces codonSampling to true if so

                cutCheck = 0; # reset cutCheck
                offScore = 0; # reset offScore
                tricky = False; # reset tricky Boolean
                badStart = False; # reset badStart Boolean
                recodedSeq = optimizeCodons(recodeSeq,orgCodonTable,codonSampling=codonSampling); # optimize codons.
                for g in gRNAs: # for every gRNA candidate within recoded region,
                    if g.index[0] >= startRecode-frame and g.index[1] <= endRecode: # if grna is inside recoded region
                        gOnSeq = g.seq; # get original gRNA sequence
                        wholeRecSeq = nonRecodedStart + recodedSeq; # add initial bases
                        gOffSeq = "";
                        anchor = -1; # will store index of gRNA bp most to the left (whichever strand). Default to -1 to indicate excision
                        if geneGB.checkInExon(g.index[0]) or geneGB.checkInExon(g.index[1]): # if the gRNA hasn't been completely excised,
                            if pamType == "NGG" and g.comp or pamType == "TTTV" and not g.comp: # if PAM is to the left of the rest of the gRNA sequence (on whichever strand),
                                anchor = g.index[0]-startRecode+frame; # stores index of gRNA bp most to the left (whichever strand)
                                for intron in intronIndices: # for every intron,
                                    if g.index[0] > intron[1]: # if anchor after end of intron,
                                        anchor -= intron[1]-intron[0]; # substract intron length from anchor index
                                    elif intron[0] >= g.index[0] >= intron[1]: # if anchor inside intron,
                                        anchor -= g.index[0] - intron[0]; # substract distance between intron start and anchor from anchor


                                gOffSeq = wholeRecSeq[anchor:anchor+len(g.seq)]; # get recoded sequence that used to be gRNA
                                if g.comp: # if on comp strand
                                    gOffSeq = revComp(gOffSeq); # save as reverse complement

                            else: # if PAM is to the right,
                                anchor = g.index[1]-startRecode+frame; # stores index of gRNA bp most to the right (whichever strand)
                                for intron in intronIndices: # for every intron,
                                    if g.index[1] > intron[1]: # if anchor after end of intron,
                                        anchor -= intron[1]-intron[0]; # substract intron length from anchor index
                                    elif intron[0] >= g.index[1] >= intron[1]: # if anchor inside intron,
                                        anchor -= g.index[1] - intron[0]; # substract distance between intron start and anchor from anchor


                                gOffSeq = wholeRecSeq[anchor-len(g.seq):anchor]; # get recoded sequence that used to be gRNA
                                if g.comp: # if on comp strand
                                    gOffSeq = revComp(gOffSeq); # save as reverse complement



                        gNewPAM = ""; # will store new PAM sequence
                        if pamType == "NGG" and anchor > -1: # if using NGG PAM and gRNA not completely excised,
                            if  (g.index[1]+3 >= endRecode and not g.comp) or (g.index[0]-3 >= startRecode and g.comp): # if PAM is within recoded region,
                                if not g.comp: # if on positive strand,
                                    gNewPAM = wholeRecSeq[anchor+len(g.seq):anchor+len(g.seq)+3]; # retrieve PAM downstream of gRNA sequence
                                else: # if on negative strand,
                                    gNewPAM = revComp(wholeRecSeq[anchor+len(g.seq)-3:anchor+len(g.seq)]); # retrieve PAM upstream of gRNA sequence, on comp strand

                            else: # if outside recoded region,
                                if g.comp: # if on comp strand,
                                    gNewPAM = geneGB.origin[g.index[1]:g.index[1]+3]; # will store new PAM sequence
                                else: # if on positive strand,
                                    gNewPAM = revComp(geneGB.origin[g.index[0]-3:g.index[0]]); # will store new PAM sequence


                        elif pamType == "TTTV" and anchor > -1: # if using TTTV PAM and gRNA not completely excised,
                            if (g.index[1]+4 >= endRecode and g.comp) or (g.index[0]-4 >= startRecode and not g.comp): # if PAM is inside recoded region,
                                if not g.comp: # if on positive strand,
                                    gNewPAM = wholeRecSeq[anchor+len(g.seq)-4:anchor+len(g.seq)]; # retrieve PAM upstream of gRNA sequence
                                else: # if on negative strand,
                                    gNewPAM = revComp(wholeRecSeq[anchor+len(g.seq):anchor+len(g.seq)+4]); # retrieve PAM downstream of gRNA sequence, on comp strand

                            else: # if outside recoded region,
                                if g.comp: # if on comp strand,
                                    gNewPAM = geneGB.origin[g.index[1]:g.index[1]+4]; # will store new PAM sequence
                                else: # if on positive strand,
                                    gNewPAM = revComp(geneGB.origin[g.index[0]-4:g.index[0]]); # will store new PAM sequence


                        newOffScore = 0; # Assume gRNA was excised
                        if offTargetMethod == "cfd" and len(gOffSeq) > 22: # if using cfd and gRNA not completely excised,
                            newOffScore = pairScoreCFD(gOnSeq,gOffSeq,gNewPAM,pamType); # calculate pairwise off-target score
                        elif offTargetMethod == "hsu" and len(gOffSeq) > 22: # if using hsu and gRNA not completely excised,
                            newOffScore = pairScoreHsu(gOnSeq,gOffSeq,gNewPAM,pamType); # calculate pairwise off-target score

                        offScore = max(offScore,newOffScore); # set offscore for next iteration

                        for g in gRNATable: # find this gRNA in table
                            if g[14] == gOnSeq: # if found,
                                g[15] = gOffSeq; # store recoded sequence
                                g[16] = str(newOffScore); # store recoded sequence's pair score



                    else: # if gRNA is not entirely contained,
                        offScore = max(offScore,0); # assume recoded


                for site in cutSeqs: # for every cut site being filtered,
                    cutCheck += findFirst(site,recodedSeq); # Find cut site, register in cutCheck
                    cutCheck += findFirst(revComp(site),recodedSeq); # Find cut site in comp strand, register in cutCheck

                if findFirst(recodedSeq,"TATATATATATATATATATA") > -1: # if 10 TA repeats found,
                    tricky = True; # it's tricky
                elif findFirst(recodedSeq,"GCGCGCGCGCGCGC") > -1: # if 7 GC repeats found,
                    tricky = True; # it's tricky
                elif findFirst(recodedSeq,"AAAAAAAAAAAAA") > -1: # if 13 A repeats found,
                    tricky = True; # it's tricky
                elif findFirst(recodedSeq,"TTTTTTTTTTTTT") > -1: # if 13 T repeats found,
                    tricky = True; # it's tricky
                elif findFirst(recodedSeq,"GGGGGGGGG") > -1: # if 9 G repeats found,
                    tricky = True; # it's tricky
                elif findFirst(recodedSeq,"CCCCCCCCC") > -1: # if 9 C repeats found,
                    tricky = True; # it's tricky

                if gcContent(recodedSeq[0:40]) < minGCEnd5Prime: # if the first bases don't have enough gc content
                    badStart = True;

                if not tricky and offScore <= offScoreThreshold and cutCheck <= -2*len(cutSeqs): # if parameters other than badStart are ok and this sequence has better start than previous best,
                    if not candidateFound: # if no candidate found until now,
                        bestRecodedSeq = recodedSeq; # make this new best
                    elif gcContent(recodedSeq[0:40]) > gcContent(bestRecodedSeq[0:40]):
                        bestRecodedSeq = recodedSeq; # make this new best

                    candidateFound = True; # signal possible candidate found

                count += 1; # advances iteration counter
                if count > 6000: # if out of iteration limit,
                    if not candidateFound: # if no candidate without cut sequences found,
                        log = log + "\nWarning: Recoded region for gene " + gene.label + " could not reshuffle gRNA cut sites enough to fulfill the maximum off-target score threshold, or contains at least one of the following cut sequences: \n" + str(cutSeqs) + "\n\n"; # log warning

                    break; # escape loop


                #print [gOnSeq+"NGG",gOffSeq+gNewPAM,pairScoreCFD(gOnSeq,gOffSeq,gNewPAM,pamType),pairScoreHsu(gOnSeq,gOffSeq,gNewPAM,pamType)]

        recodedSeq = nonRecodedStart + bestRecodedSeq; # adds initial bases from reading frame adjustment to best candidate
        annRecoded = GenBankAnn(gene.label + " Recoded", "misc_feature", recodedSeq, False, [startRecode,endRecode]); # creates var to store finished recodedSeq as annotation
        log = log + "Recoded region with size " + str(len(recodedSeq)) + " for gene " + gene.label + " selected.\n\n"; # logs this process finished

    else: # if no recoded region necessary,
        log = log + "Recoded region not deemed necessary for gene " + gene.label + ".\n\n"; # logs this process finished

    gRNATableString = "\n".join([",".join(g) for g in gRNATable]); # Creates string from grna array
    gRNATableString = gRNATableString.replace(">=threshold",">="+str(offScoreThreshold)); # adds pairwise recoded threshold values
    return {"out":annRecoded, "log":log, "gRNATable":gRNATableString}; # returns recoded region GenBankAnn object


"""
Creates list with GenBankAnn objects for forward and reverse primers for
obtaining part given. Poor design, user should check with other tools
afterwards.
"""
def createPrimers(plasmid, part, rangeSize=[18,22,50], rangeMeltTemp=[55,62,65], maxTempDif=3): #TODO: this code is a little crap.
    log = ""; # init log
    startPF = part.index[0]; # Rev primer preferred start position
    endPF = part.index[0] + rangeSize[1]; # Rev primer preferred end position
    primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

    if meltingTemp(primFwdSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq) > rangeMeltTemp[2] or not primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C": # if out of Tm range or no GC clamp
        endPF = part.index[0] + rangeSize[0]; # Smallest fwd primer end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
        maxIndexes = [startPF, endPF]; # store start and end positions of best primer in search range
        while (meltingTemp(primFwdSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq) > rangeMeltTemp[2] or not primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C") and rangeSize[0] <= len(primFwdSeq) <= rangeSize[2]: # while still no suitable Tm found or no gc clamp and still within length parameters,
            endPF = endPF + 1; # shift primer start position upstream
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            if meltingTemp(primFwdSeq)-rangeMeltTemp[1] > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]])-rangeMeltTemp[1] and primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C": # if this primer has Tm closer to the preferred Tm, and has gc clamp
                maxIndexes = [startPF, endPF]; # store start and end positions of this primer as best

        if meltingTemp(primFwdSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq) > rangeMeltTemp[2]: # if still no use
            startPF = maxIndexes[0]; # Fwd primer default start position
            endPF = maxIndexes[1]; # Fwd primer default end position
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            log = log + "\nWarning: Best fwd primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primFwdSeq)) + "\n"; # give warning


    startPR = part.index[1] - rangeSize[1]; # Rev primer start position
    endPR = part.index[1]; # Rev primer end position
    primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

    if meltingTemp(primRevSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif or not primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C": # if out of Tm range or no gc clamp
        startPR = part.index[1] - rangeSize[0]; # Smallest rev primer end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
        maxIndexes = [startPR, endPR]; # store start and end positions of best primer in search range
        while (meltingTemp(primRevSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif or not primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C") and rangeSize[0] <= len(primRevSeq) <= rangeSize[2]: # while still no suitable Tm found o no gc clamp, and still within length parameters,
            startPR = startPR - 1; # shift primer start position upstream
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            if meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) < meltingTemp(primFwdSeq)-meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]]) and primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C": # if this primer has Tm closer to the fwd primer's,
                maxIndexes = [startPR, endPR]; # store start and end positions of this primer as best

        if meltingTemp(primRevSeq) < rangeMeltTemp[0] or meltingTemp(primRevSeq) > rangeMeltTemp[2]: # if still no use and it's the reverse primer's fault,
            startPR = maxIndexes[0]; # Rev primer default start position
            endPR = maxIndexes[1]; # Rev primer default end position
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            log = log + "\nWarning: Best rev primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primRevSeq)) + "\n"; # give warning
        if meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif: # if temp difference exceeds specs
            startPR = maxIndexes[0]; # Rev primer default start position
            endPR = maxIndexes[1]; # Rev primer default end position
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            lastBase = len(primFwdSeq); # stores possible end points of fwd primer
            primFwdSeq = primFwdSeq.upper();
            while meltingTemp(primFwdSeq[0:lastBase])-meltingTemp(primRevSeq) > maxTempDif and meltingTemp(primFwdSeq[0:lastBase]) > minMeltTemp and lastBase > rangeHom[0]*2: # while T diff is still out of bounds and still within bounds of Fwd primer,
                lastBase -= 1;
                while not primFwdSeq[lastBase-1].upper().replace("G","C") == "C" and lastBase > rangeHom[0]: # find next G or C upstream
                    lastBase -= 1;


            if meltingTemp(primFwdSeq[0:lastBase])-meltingTemp(primRevSeq) < maxTempDif and meltingTemp(primFwdSeq[0:lastBase]) > rangeMeltTemp[0]: # while T diff is still out of bounds and still within bounds of Fwd primer,
                endPF = endPF - (len(primFwdSeq)-lastBase);
                primFwdSeq = plasmid.origin[startPF:endPF];
            else: # if temp difference still exceeds specs
                log = log + "\nWarning: Primers for sequence " + part.label + " under given constraints have a Tm difference of " + str(meltingTemp(primFwdSeq)-meltingTemp(primRevSeq)) + ", above the given threshold of " + str(maxTempDif)  + "\n"; # give warning

    annPrimFwd = GenBankAnn(part.label + " Primer (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF]); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Primer (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR]); # creates GenBankAnn object to hold rev primer

    log = log + "Primers for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":[annPrimFwd, annPrimRev], "log":log}; # return primer annotation objects


"""
Creates list with GenBankAnn objects for forward and reverse primers for
obtaining part for Gibson Assembly in plasmid given. The plasmid must have the
part as an annotation. rangeHom gives the length of homology on each side of the
primer [min,preferred,max].
"""
def createGibsonPrimers(plasmid, part, rangeHom=[30,40,50], minMeltTemp=68, maxTempDif=5): #TODO: debug.
    log = ""; # init log
    startPF = part.index[0] - rangeHom[1]; # Fwd primer preferred start position
    endPF = part.index[0] + rangeHom[1]; # Fwd primer preferred end position
    primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

    if meltingTemp(primFwdSeq) < minMeltTemp or not primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C": # if still no use
        startPF = part.index[0] - rangeHom[0]; # Smallest fwd primer start position
        endPF = part.index[0] + rangeHom[0]; # Smallest fwd primer end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

        maxIndexes = [startPF, endPF]; # store start and end positions of best primer in search range
        while (meltingTemp(primFwdSeq) < minMeltTemp or not primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C") and len(primFwdSeq)/2 <= rangeHom[2]: # while still no suitable Tm found and still within length parameters,
            startPF = startPF - 1; # shift primer start position upstream
            endPF = endPF + 1; # shift primer start position upstream
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            if (meltingTemp(primFwdSeq) > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]]) or not plasmid.origin[maxIndexes[1]-1].upper().replace("G","C") == "C") and primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C": # if this primer has higher Tm than the max or the current max has no gc clamp, and this one does have a gc clamp,
                maxIndexes = [startPF, endPF]; # store start and end positions of this primer

        startPF = maxIndexes[0]; # Fwd primer default start position
        endPF = maxIndexes[1]; # Fwd primer default end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
        if meltingTemp(primFwdSeq) < minMeltTemp: # if still no use
            log = log + "\nWarning: Best Gibson fwd primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primFwdSeq)) + ", below the given threshold of " + str(minMeltTemp) + "\n"; # give warning


    startPR = part.index[1] - rangeHom[1]; # Rev primer start position
    endPR = part.index[1] + rangeHom[1]; # Rev primer end position
    primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

    if meltingTemp(primRevSeq) < minMeltTemp or meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif or not primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C": # if still no use
        startPR = part.index[1] - rangeHom[0]; # Smallest fwd primer start position
        endPR = part.index[1] + rangeHom[0]; # Smallest fwd primer end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

        maxIndexes = [startPR, endPR]; # store start and end positions of best primer in search range
        while (meltingTemp(primRevSeq) < minMeltTemp or meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif or not primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C") and rangeHom[0] <= len(primRevSeq)/2 <= rangeHom[2]: # while still no suitable Tm found and still within length parameters,
            startPR = startPR - 1; # shift primer start position upstream
            endPR = endPR + 1; # shift primer start position upstream
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            if meltingTemp(primRevSeq) > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]]) and primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C": # if this primer has higher Tm than the max and has gc clamp,
                maxIndexes = [startPR, endPR]; # store start and end positions of this primer

        startPR = maxIndexes[0]; # Rev primer default start position
        endPR = maxIndexes[1]; # Rev primer default end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
        if meltingTemp(primRevSeq) < minMeltTemp: # if still no use
            log = log + "\nWarning: Best Gibson rev primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primRevSeq)) + ", below the given threshold of " + str(minMeltTemp) + "\n"; # give warning
        elif meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif: # if temp difference exceeds specs
            startPR = maxIndexes[0]; # Rev primer default start position
            endPR = maxIndexes[1]; # Rev primer default end position
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            primFwdSeq = primFwdSeq.upper(); # to uppercase
            lastBase = len(primFwdSeq)-2; # stores possible end points of fwd primer
            while not primFwdSeq[lastBase-1].upper().replace("G","C") == "C" and lastBase > rangeHom[0]: # find next G or C upstream
                lastBase -= 1;

            while meltingTemp(primFwdSeq[0:lastBase])-meltingTemp(primRevSeq) > maxTempDif and meltingTemp(primFwdSeq[0:lastBase]) > minMeltTemp and lastBase > rangeHom[0]*2: # while T diff is still out of bounds and still within bounds of Fwd primer,
                lastBase -= 1;
                while not primFwdSeq[lastBase-1].upper().replace("G","C") == "C" and lastBase > rangeHom[0]: # find next G or C upstream
                    lastBase -= 1;


            if meltingTemp(primFwdSeq[0:lastBase])-meltingTemp(primRevSeq) < maxTempDif and meltingTemp(primFwdSeq[0:lastBase]) > minMeltTemp: # while T diff is still out of bounds and still within bounds of Fwd primer,
                endPF = endPF - (len(primFwdSeq)-lastBase);
                primFwdSeq = plasmid.origin[startPF:endPF];
            else: # if temp difference still exceeds specs
                log = log + "\nWarning: Gibson primers for sequence " + part.label + " under given constraints have a Tm difference of " + str(meltingTemp(primFwdSeq)-meltingTemp(primRevSeq)) + ", above the given threshold of " + str(maxTempDif) + "\n"; # give warning

    annPrimFwd = GenBankAnn(part.label + " Gibson Primer (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF]); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Gibson Primer (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR]); # creates GenBankAnn object to hold rev primer

    log = log + "Gibson primers for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":[annPrimFwd, annPrimRev],"log":log}; # return list of primers


"""
Returns gBlock GenBankAnn object in plasmid given for part (a GenBankAnn object)
to be synthesized and inserted via Gibson in plasmid. Will give a warning if it
suspects the gBlock won't be able to be synthesized by IDT, reverse-engineering
from the IDT gene synthesis webpage.
"""
def createGBlock(plasmid, part):
    log = ""; # init log
    startGBlock = part.index[0] - 40; # gBlock start position
    endGBlock = part.index[1] + 40; # gBlock end position
    gBlockSeq = plasmid.origin[startGBlock:endGBlock]; # gBlock sequence
    tricky = False; # True if suspected to be hard to synthesize
    '''
    for i in range(0,len(gBlockSeq)-20):
        if gcContent(gBlockSeq[i:i+20]) < 0.05: # If gc content of 20 bp gBlock window is too low or tricky sequences are present
            tricky = True; # might be tricky'''

    if findFirst(gBlockSeq,"TATATATATATATATATATA") > -1: # if 10 TA repeats found,
        tricky = True; # it's tricky
    elif findFirst(gBlockSeq,"GCGCGCGCGCGCGC") > -1: # if 7 GC repeats found,
        tricky = True; # it's tricky
    elif findFirst(gBlockSeq,"AAAAAAAAAAAAA") > -1: # if 13 A repeats found,
        tricky = True; # it's tricky
    elif findFirst(gBlockSeq,"TTTTTTTTTTTTT") > -1: # if 13 T repeats found,
        tricky = True; # it's tricky
    elif findFirst(gBlockSeq,"GGGGGGGGG") > -1: # if 9 G repeats found,
        tricky = True; # it's tricky
    elif findFirst(gBlockSeq,"CCCCCCCCC") > -1: # if 9 C repeats found,
        tricky = True; # it's tricky

    if tricky: # if sequence seems tricky
        log = log + "\nWarning: I suspect IDT might reject the gBlock created: \n" + gBlockSeq + "\n"; # warn user

    annGBlock = GenBankAnn(part.label + " gBlock", "misc_feature", gBlockSeq, False, [startGBlock,endGBlock]); # creates GenBankAnn object to hold gBlock

    log = log + "gBlock for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":annGBlock, "log":log}; # returns annotation


"""
Creates list with GenBankAnn objects for forward and reverse oligos for
obtaining part for Gibson assembly in plasmid given via Klenow reaction. The
plasmid must have part as an annotation. lengthHom gives
the length of homology on each side of the part sequence.
"""
def createKlenowOligos(plasmid, part, lengthHom=40): #TODO: debug.
    log = ""; # init log
    startPF = part.index[0] - lengthHom; # Fwd primer preferred start position
    endPF = part.index[1]; # Fwd primer preferred end position
    primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

    startPR = part.index[0]; # Rev primer start position
    endPR = part.index[1] + lengthHom; # Rev primer end position
    primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

    annPrimFwd = GenBankAnn(part.label + " Klenow oligo (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF]); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Klenow oligo (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR]); # creates GenBankAnn object to hold rev primer

    log = log + "Klenow oligos for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":[annPrimFwd, annPrimRev], "log":log}; # return list of primers

"""
Returns dictionary with ["log":messageLog, "out":method out]
Out: GenBank object with the annotated, edited locus;
with all the inserted sequences in their genomic context.
"""
def editLocus(geneName, gene, construct):
    log = ""; # message log
    LHRlistAll = construct.findAnnsLabel(geneName + " LHR"); # find LHR annotations
    LHRlist = []; # contain only exact matches
    LHR = ""; # init LHR var
    for ann in LHRlistAll: # will keep only those exact matches (exclude primers, for example)
        if ann.label == geneName + " LHR":
            LHRlist.append(ann);

    if len(LHRlist) == 1: # if LHR found,
        LHR = LHRlist[0]; # save first LHR annotation as LHR
        log = log + "Found LHR annotation for locus editing." + "\n"; # add warning to log
    elif len(LHRlist) == 0: # if no LHR found,
        log = log + "ERROR: Did not find LHR annotations, genomic context file could not be built." + "\n"; # add warning to log
    else:
        LHR = LHRlist[0]; # save first LHR annotation as RHR
        log = log + "Warning: Found multiple LHR annotations, genomic context file built with the first one found." + "\n"; # add warning to log

    RHRlistAll = construct.findAnnsLabel("RHR"); # saves RHR annotation
    RHRlist = []; # contain only exact matches
    RHR = ""; # init LHR var
    for ann in RHRlistAll: # will keep only those exact matches (exclude primers, for example)
        if ann.label == geneName + " RHR":
            RHRlist.append(ann);

    if len(RHRlist) == 1: # if LHR found,
        RHR = RHRlist[0]; # save first RHR annotation as RHR
        log = log + "Found RHR annotation for locus editing." + "\n"; # add warning to log
    elif len(RHRlist) == 0: # if no LHR found,
        log = log + "ERROR: Did not find RHR annotations, genomic context file could not be built." + "\n"; # add warning to log
    else:
        RHR = RHRlist[0]; # save first RHR annotation as RHR
        log = log + "Warning: Found multiple RHR annotations, genomic context file built with the first one found." + "\n"; # add warning to log

    startInsertChrom = findFirst(gene.origin, LHR.seq) + len(LHR.seq); # find index of insertion start in the chromosomal DNA given
    endInsertChrom = findFirst(gene.origin, RHR.seq); # find index of insertion end in the chromosomal DNA given

    startInsertConst = findFirst(construct.origin, LHR.seq) + len(LHR.seq); # find index of insertion start in the chromosomal DNA given
    endInsertConst = findFirst(construct.origin, RHR.seq); # find index of insertion end in the chromosomal DNA given
    insert = construct.origin[startInsertConst:endInsertConst]; # saves string containing insert

    editedLocus = deepcopy(gene); # var will contain the edited locus
    editedLocus.removeSeq([startInsertChrom,endInsertChrom]); # remove deleted bases between homologous regions
    editedLocus.insertSeq(insert, startInsertChrom); # inserts construct into edited locus

    # Now we add all the plasmid's annotations to the edited locus
    for ann in construct.features: # loop through plasmid annotations
        if ann.index[0] >= startInsertConst and ann.index[1] <= endInsertConst: # if annotation is completely inside integrating region
            editedAnn = deepcopy(ann); # annotation on edited chromosome
            editedAnn.index[0] = editedAnn.index[0] + startInsertChrom - startInsertConst; # shift index according to new context
            editedAnn.index[1] = editedAnn.index[1] + startInsertChrom - startInsertConst; # shift index according to new context
            editedLocus.features.append(editedAnn); # adds annotation to edited locus


    log += "Edited locus (genomic context) file built.\n"; # add log entry

    return {"out":editedLocus, "log":log}; # return dictionary

"""
Abbreviates primer names to fit on commercial tube labels with the format:
Seven Digit Gene Identifier_Oligo type_Orientation

The seven digit gene identifier code follows "PF3D7_".
Oligo types include:
LHR (Gibson overhang PCR primers)
RHR (Gibson overhang PCR primers)
gRNA (Klenow oligos for gRNA sequence)
gBlock (gBlock sequencing primer)
RecKlen (Klenow oligos for recoded region, if the region is small enough)

Orientation refers to forward (F) and reverse (R) primers.
"""
def shortenOligoNames(primerString):
    mat = primerString.split("\n"); # split string into lines
    mat = [l.split(",") for l in mat]; # split lines into cells. mat is now 2D array

    for primer in mat: # iterates across all primers (rows in array)
        name = primer[0]; # gets primer name
        if name[0:6] == "PF3D7_": # if primer format start is correct
            newName = name[6:13] + "_"; # Adds numerical identifier to new name
            # Add oligo type to new name:
            if name.find("gBlock") > -1:
                newName = newName + "gBlock" + "_";
            elif name.find("LHR") > -1:
                newName = newName + "LHR" + "_";
            elif name.find("RHR") > -1:
                newName = newName + "RHR" + "_";
            elif name.find("gRNA") > -1:
                newName = newName + "gRNA" + "_";
            elif name.find("Recoded region Klenow") > -1:
                newName = newName + "RecKlen" + "_";

            # Add oligo orientation
            if name.find("fwd") > -1:
                newName = newName + "F";
            elif name.find("rev") > -1:
                newName = newName + "R";

            primer[0] = newName; # replace name


    mat = [",".join(l) for l in mat]; # join cells into lines
    outStr = "\n".join(mat); # join lines into string

    return outStr;
