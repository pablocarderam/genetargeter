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
from gRNAScores.gRNAScoring import * # import methods for gRNA scoring

#### Constants ####
# Restriction enzyme cut sequences used
cut_FseI = "ggccggcc"; # FseI cut site
cut_AsiSI = "gcgatcgc"; # AsiSI cut site
cut_IPpoI = "ctctcttaaggtagc"; # I-PpoI cut site
cut_ISceI = "attaccctgttatcccta"; # I-SceI cut site
# Plasmids
pSN054_V5 = GenBank();
pSN054_V5.load("input/plasmids/psn054-updated_v5-tagged-jn_final.gb",loadFromFile=True); # load plasmid sequence from GenBank format
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
geneFileName is the name of the .gb file (including ".gb") containing the gene
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
def pSN054TargetGene(geneName, geneFileName, codonOptimize="T. gondii", HRannotated=False, lengthLHR=[450,500,650], lengthRHR=[450,500,750], gibsonHomRange=[30,40,50], optimRangeLHR=[-20,10], optimRangeRHR=[-20,20], endSizeLHR=40, endSizeRHR=40, endTempLHR=55, endTempRHR=59, gibTemp=65, gibTDif=5, maxDistLHR=500, maxDistRHR=500, minGBlockSize=125, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI], useFileStrs=False, codonSampling=False, minGRNAGCContent=0.3, minOnTargetScore=30, offTargetMethod="cfd", offTargetThreshold=0.5, maxOffTargetHitScore=35): # cfd, 0.5, 35; hsu, 75, 5
    outputDic = {"geneName":geneName, "newGene":GenBank(), "editedLocus":GenBank(), "newPlasmid":GenBank(), "geneFileStr":"", "plasmidFileStr":"", "oligoFileStr":"", "logFileStr":"", "editedLocusFileStr":""}; # dictionary containing keys to all values being returned
    outputDic["logFileStr"] = outputDic["logFileStr"] + " **** Message log for " + geneName + "-targeting construct based on plasmid pSN054_V5 **** \n\n"; # starts message log to file

    geneGB = GenBank(); # initializes variable to hold gene GenBank object
    if useFileStrs: # if we are using file strings,
        geneGB.load(geneFileName, loadFromFile=False); # load gene from file string
        path = ""; # sets an empty path, needed for save functions of class GenBank
    else: # If we're using files,
        geneGB.load(geneFileName, loadFromFile=True); # load gene file
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



        if HRannotated: # if not using manual annotations
            gRNA = findGRNA(geneGB, gene); # finds gRNA most upstream annotated manually.
            if len(gRNA["out"].label) == 0: # if no manual annotation found,
                gRNA = chooseGRNA(geneGB, gene, minGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, minOffTargetScore=offTargetThreshold, offTargetMethod=offTargetMethod, maxOffTargetHitScore=maxOffTargetHitScore); # chooses gRNA.

        else: # if not,
            gRNA = chooseGRNA(geneGB, gene, minGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, minOffTargetScore=offTargetThreshold, offTargetMethod=offTargetMethod, maxOffTargetHitScore=maxOffTargetHitScore); # chooses gRNA.

        outputDic["logFileStr"] = outputDic["logFileStr"] + gRNA["log"]; # add logs
        gRNA = gRNA["out"]; # saves actual data

        if len(gRNA.label) > 0: # if gRNA found,
            LHR = GenBankAnn(); # init LHR var
            RHR = GenBankAnn(); # init RHR var

            # pick HRs first
            LHR = chooseLHR(geneGB, gene, lengthLHR=lengthLHR, minTmEnds=endTempLHR, endsLength=endSizeLHR, optimizeRange=optimRangeLHR, maxDistanceFromGRNA=maxDistLHR); # chooses an LHR
            RHR = chooseRHR(geneGB, gene, lengthRHR=lengthRHR, minTmEnds=endTempLHR, endsLength=endSizeRHR, optimizeRange=optimRangeRHR, maxDistanceFromGene=maxDistRHR); # chooses RHR

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



            recoded = chooseRecodeRegion(geneGB, gene, orgCodonTable=codonUsageTables[codonOptimize],codonSampling=codonSampling); # defines region to be recoded, returns recoded sequence
            outputDic["logFileStr"] = outputDic["logFileStr"] + recoded["log"]; # add logs
            recoded = recoded["out"]; # saves actual data

            pSN054_ARMED = insertTargetingElementsPSN054(pSN054_V5, gene.label, gRNA.seq, LHR.seq, recoded.seq, RHR.seq); # inserts targeting elements

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
                outputDic["logFileStr"] = outputDic["logFileStr"] + "\ngBlock deemed not feasible for recoded region of construct targeting gene " + geneName + ", used Klenow instead.\n"; # say so
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

            outputDic["logFileStr"] = outputDic["logFileStr"] + "\n\nVector constructed and written to file. End of process.\n"; # saves message log to file
            pSN054_ARMED.definition = (pSN054_ARMED.definition + "  " + outputDic["logFileStr"]).replace("\n","   "); # save logs to file definition to be viewed in benchling

            outputDic["geneFileStr"] = geneGB.save(path + "/" + geneName + "_Locus_Pre-editing.gb", saveToFile=(not useFileStrs)); # saves annotated gene
            outputDic["plasmidFileStr"] = pSN054_ARMED.save(path + "/" +  "pSN054_V5_targeting" + geneName, saveToFile=(not useFileStrs)); # saves plasmid
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
        outputDic["logFileStr"] = outputDic["logFileStr"] + "\nERROR: No gene annotations found in this file with name " + geneName + "\nProcess terminated.";

    return outputDic; # returns output dictionary


"""
Inserts targeting elements given as arguments into pSN054_V5 at predetermined
sites. Elements given as arguments must be strings, not GenBankAnn objects.
The gRNA, LHR and RHR given must be in 5' to 3' sense and in the
positive strand and must not contain RE cut sites cut_FseI, cut_AsiSI,
cut_IPpoI or cut_ISceI. Returns new GenBank object with targeting elements.
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
def chooseGRNA(geneGB, gene, searchRange=[-500,125], PAM="NGG", side3Prime=True, minGCContent=0.3, minOnTargetScore=25, minOffTargetScore=75, maxOffTargetHitScore=35, offTargetMethod="hsu", gLength=20, maxDistanceBetweenGRNAS=50, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI]): # could've been useful at some point: http://grna.ctegd.uga.edu/ http://www.broadinstitute.org/rnai/public/software/sgrna-scoring-help http://crispr.mit.edu/about
    log = ""; # init log
    pamSeq = "GG"; # should be derived from PAM, but replacing N would involve regex and I'm lazy. The Doench et al. (2016) gRNA scoring technologies only work with NGG PAMs anyway.
    gene = geneGB.findAnnsLabel(gene.label)[0]; # stores gene annotation
    gRNAs = []; # list of GenBankAnn objects for gRNAs.
    backupGRNAs = []; # stores gRNAs that fail off-target score as possible backups
    gRNAUpstream = GenBankAnn(); # will store gRNA most upstream

    if len(gene.label) > 0: # if gene found,
        searchStart = gene.index[side3Prime]+searchRange[0]; # Start searching for gRNAs in a position relative to gene start or end point
        searchEnd = gene.index[side3Prime]+searchRange[1]; # Finish searching for gRNAs in a position relative to gene start or end point
        searchSeq = geneGB.origin[searchStart:searchEnd].upper(); # get sequence where gRNAs will be searched for, centered around start or end of gene
        i = len(searchSeq)-6; # start searching for gRNAs by searching for PAM in searchSeq. Start at downstream end (allow five bases downstream to accomodate on-target 30-mer gRNA sequence), go upstream
        while i > 3 and len(gRNAs) < 3: # iterate through all searchSeq until done with searchSeq (allow 3 bases upstream of search end point to accomodate on-target 30-mer minus-strand gRNA sequence) or three candidates found.
            comp = False; # set sense to plus strand
            extGRNASeq = ""; # stores extended sequence sequence NNNN-gRNA-PAM-NNN
            gRNAIndexes = []; # store gRNA indexes
            if searchSeq[i:i+2] == pamSeq: # if PAM found on plus strand,
                extGRNASeq = searchSeq[i-25:i+5]; # store extended sequence NNNN-gRNA-PAM-NNN
                gRNAIndexes = [searchStart+i-21,searchStart+i-1]; # store gRNA indexes
            elif searchSeq[i:i+2] == revComp(pamSeq): # if PAM found on minus strand,
                extGRNASeq = revComp(searchSeq[i-3:i+27]); # store extended sequence NNNN-gRNA-PAM-NNN on comp strand
                gRNAIndexes = [searchStart+i+3,searchStart+i+23]; # store gRNA indexes on comp strand
                comp = True; # set sense to complementary strand

            if len(extGRNASeq) == 30: # if extended gRNA is right size
                gRNASeqPAM = extGRNASeq[4:27]; # store actual gRNA seq with PAM
                gRNASeq = extGRNASeq[4:24]; # store actual gRNA seq without PAM
                gc = gcContent(gRNASeq); # store gc content
                if gc >= minGCContent and findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") < 0: # if gc content is acceptable and does not contain 10 or more consecutive As or Ts,
                    onTarget = onTargetScore(extGRNASeq); # store on-target score
                    if onTarget > minOnTargetScore: # if on-target score is passable,
                        offTargetScores = offTargetScore(gRNASeqPAM,offTargetMethod); # store off-target score
                        newGRNA = GenBankAnn(gene.label+" gRNA ","misc",gRNASeq,comp,gRNAIndexes); # create annotation object with gRNA information
                        newGRNA.onTarget = onTarget; # Add on-target score as attribute
                        newGRNA.offTarget = offTargetScores; # Add off-target scores as attribute
                        newGRNA.gc = gc; # Add gc content as attribute
                        if offTargetScores[0] > minOffTargetScore and offTargetScores[1] < maxOffTargetHitScore: # if total off-target score and max hit score are passable
                            gRNAs.append(newGRNA); # add this gRNA's information to list.
                        else: # if failed off-target score culling.
                            newGRNA.label = newGRNA.label + "(backup)"; # notify backup on label
                            backupGRNAs.append(newGRNA); # add this gRNA's information to backup list.





            i -= 1; # advance indexer

        if len(gRNAs) == 0: # if no gRNAs found,
            log = log + "\nERROR: no acceptable gRNA with at least " + str(minGCContent) + " GC content, " + str(minOnTargetScore) + " on-target score and " + str(minOffTargetScore) + " off-target total score found on gene " + gene.label + ", please modify your criteria or select and annotate one manually.\n\n"; # add warning to log
        else: #if gRNAs found,
            newList = []; # new list will only contain gRNAs within acceptable range
            if gRNAs[0].index[0] > gene.index[1]-3: # if most downstream gRNA starts after start of stop codon,
                for g in gRNAs: # loop through gRNAs
                    if g.index[0] > gene.index[1]-3 and gRNAs[0].index[0] - g.index[0] < maxDistanceBetweenGRNAS: # if this gRNA is still in the 3' UTR or stop codon and is within the max distance from the most downstream gRNA,
                        newList.append(g); # add it to the new list



            else: # if most downstream gRNA is upstream of stop codon,
                for g in gRNAs: # loop through gRNAs
                    if gRNAs[0].index[0] - g.index[0] < maxDistanceBetweenGRNAS: # if this gRNA is within the max distance from the most downstream gRNA
                        newList.append(g); # add it to the new list




            gRNAs = newList; # keep the new list of gRNAs as the main list
            gRNAs.sort(key=lambda g: g.gc); # sorts gRNA list according to custom function (GC content)
            count = 1; # counter for numbering candidates according to their quality
            for g in gRNAs: # loop thorugh gRNAs ordered by GC content
                g.label = g.label + str(count); # add number to gRNA label
                count += 1; # advance counter

            bestGRNA = gRNAs[0]; # will store best gRNA
            gRNAUpstream = gRNAs[0]; # will find most upstream gRNA
            for g in gRNAs: # loop through gRNAs
                if g.index[0] < gRNAUpstream.index[0]: # if more upstream than previous most upstream
                    gRNAUpstream = g; # set this gRNA as most upstream


            log = log + "\n" + str(len(gRNAs)) + " acceptable gRNAs were selected automatically on gene " + gene.label + ". \ngRNA 1 has GC content of " + str(bestGRNA.gc) + ", on-target score of " + str(bestGRNA.onTarget)  + ", \nand aggregated off-target score of " + str(bestGRNA.offTarget[0]) + " (Method: " + offTargetMethod + ", Max. Hit Score: " + str(bestGRNA.offTarget[1]) + ", Num. hits: " + str(bestGRNA.offTarget[2]) + "."; # add warning to log


    geneGB.features = geneGB.features + gRNAs; # add gRNAs to gene GenBank object features list
    countBackups = 0; # counts how many backup gRNAs are included
    for g in backupGRNAs: # for every possible backupGRNAs
        if g.index[0] > gRNAUpstream.index[0]: # if this backup is downstream of most upstream gRNA,
            geneGB.features.append(g); # add to features list
            countBackups +=1; # advances counter

    log = log + "\n" + str(len(backupGRNAs)) + " backup gRNAs with possible off-target effects annotated.\n\n"

    return {"out":gRNAUpstream, "log":log}; # returns gRNA and log


"""
Finds gRNA already annotated on gene. Filters for given restriction enzyme cut
sites. If no gRNA found on gene, nothing will happen (logs this).
"""
def findGRNA(geneGB, gene, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI]):
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

    return {"out":gRNAUpstream, "log":log}; # returns gRNA and log



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
def chooseLHR(geneGB, gene, lengthLHR=[450,500,650], minTmEnds=55, endsLength=40, optimizeRange=[-20,10], maxDistanceFromGRNA=500, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI]):
    #TODO: debug cases in which LHR has to be shifted, do case in which LHR is in intron?
    log = ""; # init log
    gRNAs = geneGB.findAnnsLabel("gRNA 1", True); # List of all gRNAs
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
    startLHR = endLHR - lengthLHR[1]; # stores LHR start index according to preferred length
    failSearchEnd = False; # true if no suitable LHR end region is found within given parameters.

    # First search for end region
    while meltingTemp(geneGB.origin[(endLHR-40):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchEnd: # while no suitable end region found, still within max distance from gRNA, and the search for a suitable end region hasn't failed yet,
        endLHR -= 1; # shift endLHR upstream

    if meltingTemp(geneGB.origin[(endLHR-40):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR >= maxDistanceFromGRNA: # if no suitable end region found,
        endLHR = min(gRNAUpstream.index[0],gene.index[1]-3); # saves end index of LHR by default as whatever is more upstream between the start of the gRNA or the end end of the gene (minus the stop codon)
        failSearchEnd = True; # changes failSearchEnd status
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of LHR and gRNA. \nDefaulted to ending right before start of gRNA most upstream." + "\n"; # give a warning

    # Then search for start region; modify end region if necessary
    while meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) < minTmEnds and gRNAUpstream.index[0]-endLHR <= maxDistanceFromGRNA: # if no start region has been found and still within max distance from gRNA,
        # find new end region if necessary
        while meltingTemp(geneGB.origin[(endLHR-40):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR <= maxDistanceFromGRNA and not failSearchEnd: # while no suitable end region found, still within max distance from gRNA, and the search for a suitable end region hasn't failed yet,
            endLHR -= 1; # shift endLHR upstream

        # search for starts upstream
        while meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) < minTmEnds and endLHR-startLHR <= lengthLHR[2]: # while no suitable start region found and still within max length of LHR,
            startLHR -= 1; # shift startLHR upstream

        # if not found upstream
        if meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) < minTmEnds and endLHR-startLHR > lengthLHR[2]: # if no start region found this way,
            startLHR = endLHR - lengthLHR[1]; # return to preferred position
            # search downstream
            while meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) < minTmEnds and endLHR-startLHR >= lengthLHR[0]: # while no suitable start region found and still within min length of LHR,
                startLHR += 1; # shift startLHR downstream

        if meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) < minTmEnds: # if still not found
            endLHR -= 1; # shifts end of LHR upstream


    if meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) < minTmEnds and gRNAUpstream.index[0]-endLHR > maxDistanceFromGRNA: # if no suitable start region found,
        endLHR = gRNAUpstream.index[0]; # resets endLHR
        while meltingTemp(geneGB.origin[(endLHR-40):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR <= maxDistanceFromGRNA and not failSearchEnd: # while no suitable end region found, still within max distance from gRNA, and the search for a suitable end region hasn't failed yet,
            endLHR -= 1; # shift endLHR upstream

        startLHR = endLHR - lengthLHR[1]; # stores LHR start index according to preferred length by default
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of LHR and gRNA. \nDefaulted to starting so that the LHR size is equal to preferred size of " + str(lengthLHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting LHR by adjusting start sequence around the given range
    searchIndexes = [int(startLHR+optimizeRange[0]), int(startLHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexes[0], searchIndexes[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[i:(i+40)]) > meltingTemp(geneGB.origin[startLHR:(startLHR+40)]) and lengthLHR[0] >= endLHR-i >= lengthLHR[2]: # if this start point has a better Tm and is still within bounds,
            startLHR = i; # make this the starting position

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
def chooseRHR(geneGB, gene, lengthRHR=[450,500,750], minTmEnds=59, endsLength=40, optimizeRange=[-20,20], maxDistanceFromGene=500, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI]):
    #TODO: debug cases in which RHR has to be shifted
    log = ""; # init log
    gRNAs = geneGB.findAnnsLabel("gRNA 1", True); # List of all gRNAs
    gRNADownstream = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
        gRNADownstream = gRNAs[0]; # will store gRNA most upstream
        for g in gRNAs: # loops across gRNAs
            if g.index[0] > gRNADownstream.index[0]: # if more upstream
                gRNADownstream = g; # replace as most upstream


    else: # if gRNAs found,
        gRNADownstream = gRNAs[0]; # will store gRNA most upstream

    genes = geneGB.findAnnsType("gene"); # list of all genes, used to verify RHR doesn't start inside any genes (truncating them)
    def checkInGene(pIndex): # checks whether the index given is inside a gene.
        inGene = False; # default is false
        for g in genes: # loop through all genes
            inGene = g.index[0] <= pIndex < g.index[1]; # if inside gene, set to True.

        return inGene; # return

    startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or after most downstream gRNA.
    print startRHR

    endRHR = startRHR + lengthRHR[1]; # stores RHR start index according to preferred length

    failSearchStart = False; # true if no suitable RHR start region is found within given parameters.

    # First search for start region
    while meltingTemp(geneGB.origin[startRHR:(startRHR+40)]) < minTmEnds and startRHR - gene.index[1] <= maxDistanceFromGene and not checkInGene(startRHR): # while no suitable start region found, still within max distance from gene, and the search for a suitable start region hasn't failed yet,
        startRHR += 1; # shift startRHR downstream

    if meltingTemp(geneGB.origin[startRHR:(startRHR+40)]) < minTmEnds and ( startRHR - gene.index[1] > maxDistanceFromGene or checkInGene(startRHR) ): # if no suitable start region found,
        startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or 1 bp after most downstream gRNA by default
        failSearchStart = True; # changes failSearchStart status
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of gene and start of RHR. \nDefaulted to starting right after end of gene." + "\n"; # give a warning

    # search for end downstream
    while meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2]: # while no suitable end region found and still within max length of RHR,
        endRHR += 1; # shift endRHR downstream

    # if not found downstream
    if meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
        endRHR = startRHR + lengthRHR[1]; # return to preferred position
        # search upstream
        while meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and endRHR-startRHR >= lengthRHR[0]: # while no suitable start region found and still within min length of RHR,
            endRHR -= 1; # shift endRHR downstream


    if meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and checkInGene(startRHR): # if still not found and within bounds,
        # Then modify start region, search for end region;
        while meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and not checkInGene(startRHR): # if no end region has been found,
            startRHR += 1; # shifts start of RHR downstream

            # find new start region if necessary
            while meltingTemp(geneGB.origin[startRHR:(startRHR+40)]) < minTmEnds and startRHR - gene.index[1] <= maxDistanceFromGene and checkInGene(startRHR) and not failSearchStart: # while no suitable start region found, still within max distance from gene, and the search for a suitable start region hasn't failed yet,
                startRHR += 1; # shift startRHR downstream

            # search for end downstream
            while meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2]: # while no suitable end region found and still within max length of RHR,
                endRHR += 1; # shift endRHR downstream

            # if not found downstream
            if meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
                endRHR = startRHR + lengthRHR[1]; # return to preferred position
                # search upstream
                while meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and endRHR-startRHR >= lengthRHR[0]: # while no suitable start region found and still within min length of RHR,
                    endRHR -= 1; # shift endRHR downstream



    if meltingTemp(geneGB.origin[(endRHR-40):endRHR]) < minTmEnds and ( startRHR-gene.index[1] > maxDistanceFromGene or checkInGene(startRHR) ): # if no suitable end region found,
        startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or 1 bp after most downstream gRNA by default
        while meltingTemp(geneGB.origin[startRHR:(startRHR+40)]) < minTmEnds and startRHR-gene.index[0] <= maxDistanceFromGene and not failSearchStart: # while no suitable start region found, still within max distance from gene, and the search for a suitable start region hasn't failed yet,
            startRHR += 1; # shift startRHR downstream

        endRHR = startRHR + lengthRHR[1]; # stores RHR start index according to preferred length by default
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of RHR and end of gene. \nDefaulted to ending at preferred size of " + str(lengthRHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting RHR by adjusting start and stop sequences around the given range
    searchIndexesEnd = [int(endRHR+optimizeRange[0]), int(endRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better end point.
    for i in range(searchIndexesEnd[1], searchIndexesEnd[0], -1): # iterates across optimization range in reverse
        if meltingTemp(geneGB.origin[(i-40):i]) > meltingTemp(geneGB.origin[(endRHR-40):endRHR]) and lengthRHR[2] >= i-startRHR >= lengthRHR[0]: # if this end point has a better Tm and is still within bounds,
            if not checkInGene(i): # if not inside a gene,
                endRHR = i; # make this the ending position



    searchIndexesStart = [int(startRHR+optimizeRange[0]), int(startRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexesStart[0], searchIndexesStart[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[i:(i+40)]) > meltingTemp(geneGB.origin[startRHR:(startRHR+40)]) and lengthRHR[2] >= endRHR-i >= lengthRHR[0] and i > gRNADownstream.index[1]: # if this start point has a better Tm and is still within bounds and after gRNA,
            if not checkInGene(i): # if not inside a gene,
                startRHR = i; # make this the starting position



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
restriction sites given as parameters.
"""
def chooseRecodeRegion(geneGB, gene, orgCodonTable=codonUsage(), filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI], codonSampling=False):
    #TODO: debug
    log = ""; # init log
    LHR = geneGB.findAnnsLabel("LHR")[0]; # LHR annotation object

    annRecoded = GenBankAnn(); # creates GenBankAnn object to hold recoded region
    if LHR.index[1] < gene.index[1]: # if end of LHR is inside gene
        startRecode = LHR.index[1]; # start of recode region (start of gRNA most upstream)
        endRecode = gene.index[1] - 3; # end of recode region (end of gene, exclude stop codon)
        frame = (endRecode-startRecode) % 3; # stores reading frame, index from start of sequence to be recoded
        startRecode += frame; # modify recode start site according to reading frame
        recodeSeq = geneGB.origin[startRecode:endRecode]; # get sequence to be recoded

        gRNAs = geneGB.findAnnsLabel("gRNA 1", True); # List of all gRNAs
        gRNAUpstream = GenBankAnn(); # init var to hold gRNA
        if len(gRNAs) == 0: # if no gRNAs found
            gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
            gRNAUpstream = gRNAs[0]; # will store gRNA most upstream
            for g in gRNAs: # loops across gRNAs
                if g.index[0] < gRNAUpstream.index[0]: # if more upstream
                    gRNAUpstream = g; # replace as most upstream


        else: # if gRNAs found,
            gRNAUpstream = gRNAs[0]; # will store gRNA most upstream

        cutSeqs = filterCutSites + [g.seq for g in gRNAs]; # list of all cut seqs. all gRNAs in gene are to be included as cut sequences
        cutCheck = 0; # variable used to check if cut sequences are present. Initially greater than -1*len(cutSeqs) since all gRNAs are present.
        count = 0; # iteration counter
        recodedSeq = recodeSeq; # assign recoded sequence to same as original
        if len(recodeSeq) > 2: # if recodeSeq contains at least one codon,
            while cutCheck > -2*len(cutSeqs): # while cutCheck is greater than what you would expect for no hits in all cut sequences plus the gRNAs on both positive and comp strands,
                cutCheck = 0; # reset cutCheck
                recodedSeq = optimizeCodons(recodeSeq,orgCodonTable,codonSampling=codonSampling); # optimize codons.
                for site in cutSeqs: # for every cut site being filtered,
                    cutCheck += findFirst(site,recodedSeq); # Find cut site, register in cutCheck
                    cutCheck += findFirst(revComp(site),recodedSeq); # Find cut site in comp strand, register in cutCheck

                count += 1; # advances iteration counter
                if count > 10000: # if out of iteration limit,
                    log = log + "\nWarning: Recoded region for gene " + gene.label + " contains at least one of the following cut sequences: \n" + str(cutSeqs) + "\n\n"; # log warning
                    break; # escape loop



        recodedSeq = geneGB.origin[LHR.index[1]:LHR.index[1]+frame] + recodedSeq; # adds initial bases from reading frame adjustment
        # creates var to store finished recodedSeq as annotation
        annRecoded = GenBankAnn(gene.label + " Recoded", "misc_feature", recodedSeq, False, [startRecode,endRecode]); # creates GenBankAnn object to hold RHR
        log = log + "\nRecoded region for gene " + gene.label + " selected.\n\n"; # logs this process finished

    else: # if no recoded region necessary,
        log = log + "\nRecoded region for gene " + gene.label + " not deemed necessary.\n\n"; # logs this process finished

    return {"out":annRecoded, "log":log}; # returns recoded region GenBankAnn object


"""
Creates list with GenBankAnn objects for forward and reverse primers for
obtaining part given. Poor design, user should check with other tools
afterwards.
"""
def createPrimers(plasmid, part, rangeSize=[18,20,50], rangeMeltTemp=[59,60,65], maxTempDif=3): #TODO: this code is crap.
    log = ""; # init log
    startPF = part.index[0]; # Rev primer preferred start position
    endPF = part.index[0] + rangeSize[1]; # Rev primer preferred end position
    primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

    if meltingTemp(primFwdSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq) > rangeMeltTemp[2]: # if out of Tm range
        endPF = part.index[0] + rangeSize[0]; # Smallest fwd primer end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
        maxIndexes = [startPF, endPF]; # store start and end positions of best primer in search range
        while (meltingTemp(primFwdSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq) > rangeMeltTemp[2]) and rangeSize[0] <= len(primFwdSeq) <= rangeSize[2]: # while still no suitable Tm found and still within length parameters,
            endPF = endPF + 1; # shift primer start position upstream
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            if meltingTemp(primFwdSeq)-rangeMeltTemp[1] > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]])-rangeMeltTemp[1]: # if this primer has Tm closer to the preferred Tm,
                maxIndexes = [startPF, endPF]; # store start and end positions of this primer as best

        if meltingTemp(primFwdSeq) < rangeMeltTemp[0] or meltingTemp(primFwdSeq) > rangeMeltTemp[2]: # if still no use
            startPF = maxIndexes[0]; # Fwd primer default start position
            endPF = maxIndexes[1]; # Fwd primer default end position
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            log = log + "\nWarning: Best fwd primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primFwdSeq)) + "\n"; # give warning


    startPR = part.index[1] - rangeSize[1]; # Rev primer start position
    endPR = part.index[1]; # Rev primer end position
    primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

    if meltingTemp(primRevSeq) < rangeMeltTemp[0] or meltingTemp(primRevSeq) > rangeMeltTemp[2] or meltingTemp(primRevSeq)-meltingTemp(primFwdSeq) > maxTempDif: # if out of Tm range
        startPR = part.index[1] - rangeSize[0]; # Smallest rev primer end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
        maxIndexes = [startPR, endPR]; # store start and end positions of best primer in search range
        while (meltingTemp(primRevSeq) < rangeMeltTemp[0] or meltingTemp(primRevSeq) > rangeMeltTemp[2] or meltingTemp(primRevSeq)-meltingTemp(primFwdSeq) > maxTempDif) and len(primRevSeq) <= rangeSize[2] <= len(primRevSeq): # while still no suitable Tm found and still within length parameters,
            startPR = startPR - 1; # shift primer start position upstream
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            if meltingTemp(primRevSeq)-meltingTemp(primFwdSeq) > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]])-meltingTemp(primFwdSeq): # if this primer has Tm closer to the fwd primer's,
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

    if meltingTemp(primFwdSeq) < minMeltTemp: # if still no use
        startPF = part.index[0] - rangeHom[0]; # Smallest fwd primer start position
        endPF = part.index[0] + rangeHom[0]; # Smallest fwd primer end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

        maxIndexes = [startPF, endPF]; # store start and end positions of best primer in search range
        while meltingTemp(primFwdSeq) < minMeltTemp and rangeHom[0] <= len(primFwdSeq) <= rangeHom[2]: # while still no suitable Tm found and still within length parameters,
            startPF = startPF - 1; # shift primer start position upstream
            endPF = endPF + 1; # shift primer start position upstream
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            if meltingTemp(primFwdSeq) > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]]): # if this primer has higher Tm than the max,
                maxIndexes = [startPF, endPF]; # store start and end positions of this primer

        startPF = maxIndexes[0]; # Fwd primer default start position
        endPF = maxIndexes[1]; # Fwd primer default end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
        if meltingTemp(primFwdSeq) < minMeltTemp: # if still no use
            log = log + "\nWarning: Best Gibson fwd primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primFwdSeq)) + ", below the given threshold of " + str(minMeltTemp) + "\n"; # give warning


    startPR = part.index[1] - rangeHom[1]; # Rev primer start position
    endPR = part.index[1] + rangeHom[1]; # Rev primer end position
    primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

    if meltingTemp(primRevSeq) < minMeltTemp or meltingTemp(primRevSeq)-meltingTemp(primFwdSeq) > maxTempDif: # if still no use
        startPR = part.index[1] - rangeHom[0]; # Smallest fwd primer start position
        endPR = part.index[1] + rangeHom[0]; # Smallest fwd primer end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

        maxIndexes = [startPR, endPR]; # store start and end positions of best primer in search range
        while (meltingTemp(primRevSeq) < minMeltTemp or meltingTemp(primRevSeq)-meltingTemp(primFwdSeq) > maxTempDif) and rangeHom[0] <= len(primRevSeq) <= rangeHom[2]: # while still no suitable Tm found and still within length parameters,
            startPR = startPR - 1; # shift primer start position upstream
            endPR = endPR + 1; # shift primer start position upstream
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            if meltingTemp(primRevSeq) > meltingTemp(plasmid.origin[maxIndexes[0]:maxIndexes[1]]): # if this primer has higher Tm than the max,
                maxIndexes = [startPR, endPR]; # store start and end positions of this primer

        startPR = maxIndexes[0]; # Rev primer default start position
        endPR = maxIndexes[1]; # Rev primer default end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
        if meltingTemp(primRevSeq) < minMeltTemp: # if still no use
            log = log + "\nWarning: Best Gibson rev primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(primRevSeq)) + ", below the given threshold of " + str(minMeltTemp) + "\n"; # give warning
        elif meltingTemp(primFwdSeq)-meltingTemp(primRevSeq) > maxTempDif: # if temp difference exceeds specs
            log = log + "\nWarning: Gibson primers for sequence " + part.label + " under given constraints have a Tm difference of " + str(meltingTemp(primFwdSeq)-meltingTemp(primRevSeq)) + ", above the given threshold of " + str(maxTempDif) + "\n"; # give warning

    annPrimFwd = GenBankAnn(part.label + " Gibson Primer (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF]); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Gibson Primer (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR]); # creates GenBankAnn object to hold rev primer

    log = log + "Gibson primers for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":[annPrimFwd, annPrimRev],"log":log}; # return list of primers


"""
Returns gBlock GenBankAnn object in plasmid given for part (a GenBankAnn object)
to be synthesized and inserted via Gibson in plasmid. Will give a warning if it
suspects the gBlock won't be able to be synthesized by IDT.
"""
def createGBlock(plasmid, part):
    log = ""; # init log
    startGBlock = part.index[0] - 40; # gBlock start position
    endGBlock = part.index[1] + 40; # gBlock end position
    gBlockSeq = plasmid.origin[startGBlock:endGBlock]; # gBlock sequence
    tricky = False; # True if suspected to be hard to synthesize
    for i in range(0,len(gBlockSeq)-20):
        if gcContent(gBlockSeq[i:i+20]) < 0.05: # If gc content of 20 bp gBlock window is too low or tricky sequences are present
            tricky = True; # might be tricky

    if findFirst("TATATATATATATATATA", gBlockSeq) > -1: # if 9 TA repeats found,
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
    startPF = part.index[0]; # Fwd primer preferred start position
    endPF = part.index[1] + lengthHom; # Fwd primer preferred end position
    primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

    startPR = part.index[0] - lengthHom; # Rev primer start position
    endPR = part.index[1]; # Rev primer end position
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
        log = log + "\nFound LHR annotation for locus editing." + "\n"; # add warning to log
    elif len(LHRlist) == 0: # if no LHR found,
        log = log + "\nERROR: Did not find LHR annotations, genomic context file could not be built." + "\n"; # add warning to log
    else:
        LHR = LHRlist[0]; # save first LHR annotation as RHR
        log = log + "\nWarning: Found multiple LHR annotations, genomic context file built with the first one found." + "\n"; # add warning to log

    RHRlistAll = construct.findAnnsLabel("RHR"); # saves RHR annotation
    RHRlist = []; # contain only exact matches
    RHR = ""; # init LHR var
    for ann in RHRlistAll: # will keep only those exact matches (exclude primers, for example)
        if ann.label == geneName + " RHR":
            RHRlist.append(ann);

    if len(RHRlist) == 1: # if LHR found,
        RHR = RHRlist[0]; # save first RHR annotation as RHR
        log = log + "\nFound RHR annotation for locus editing." + "\n"; # add warning to log
    elif len(RHRlist) == 0: # if no LHR found,
        log = log + "\nERROR: Did not find RHR annotations, genomic context file could not be built." + "\n"; # add warning to log
    else:
        RHR = RHRlist[0]; # save first RHR annotation as RHR
        log = log + "\nWarning: Found multiple RHR annotations, genomic context file built with the first one found." + "\n"; # add warning to log

    startInsertChrom = findFirst(gene.origin, LHR.seq); # find index of insertion start in the chromosomal DNA given
    endInsertChrom = findFirst(gene.origin, RHR.seq) + len(RHR.seq); # find index of insertion end in the chromosomal DNA given

    startInsertConst = findFirst(construct.origin, LHR.seq); # find index of insertion start in the chromosomal DNA given
    endInsertConst = findFirst(construct.origin, RHR.seq) + len(RHR.seq); # find index of insertion end in the chromosomal DNA given
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


    log += "\nEdited locus (genomic context) file built."; # add log entry

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
