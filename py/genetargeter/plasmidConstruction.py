
import copy; # Import object copying methods for deep copies

from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants

"""
Inserts targeting elements given as arguments into pSN054_V5_Cas9 at predetermined
sites. Elements given as arguments must be strings, not GenBankAnn objects.
The gRNA, LHR and RHR given must be in 5' to 3' sense and in the
positive strand and must not contain RE cut sites cut_FseI, cut_AsiSI,
cut_IPpoI, cut_ISceI, or cut_AflII. Returns new GenBank object with targeting
elements.
Specifically:
Inserts gRNA used by Cas9 in I-PpoI cut site. Removes recognition sequence. Adds
    GG at 5'-end of the gRNA sequence. The GG is required for the T7 RNA
    polymerase to efficiently transcribe the gRNA to high levels.
Inserts the LHR between FseI and AsiSI cut sites, but leaves the sites intact
    with their recognition sequences.
Inserts the recoded region after the LHR, before the AsiSI cut site (leaves the
    site intact).
Inserts RHR after I-SceI cut site, leaves the site intact with its recognition
    sequence.
"""
def insertTargetingElementsPSN054(plasmid, geneName, gRNA, LHR, recodedRegion, RHR):
    plas = copy.deepcopy(plasmid); # makes a copy of the plasmid object to modify without altering the original

    startLHR = findFirst(plas.origin, cut_FseI) + len(cut_FseI); # index of LHR start site (at end of FseI cut sequence)
    endLHR = findFirst(plas.origin, cut_AsiSI); # index of LHR end site (at start of AsiSI cut sequence)
    plas.removeSeq([startLHR, endLHR]); # removes sequence that LHR will replace
    plas.insertSeq(LHR.upper(), startLHR); # inserts LHR sequence
    annLHR = GenBankAnn(geneName+" LHR", "misc_feature", LHR, False, [startLHR,startLHR+len(LHR)], annColors['LHRColor']); # annotation object
    plas.features.append(annLHR); # adds annotation

    if len(recodedRegion) > 0: # if there is a recoded region,
        inRecode = annLHR.index[1]; # index of recoded region start site (at end of LHR)
        plas.insertSeq(recodedRegion.upper(), inRecode); # inserts recoded region sequence
        annRecoded = GenBankAnn(geneName+" Recoded Region", "misc_feature", recodedRegion, False, [inRecode,inRecode+len(recodedRegion)], annColors['recodedRegionColor']); # annotation object
        plas.features.append(annRecoded); # adds annotation

    startGRNA = findFirst(plas.origin, cut_IPpoI); # index of gRNA start site (at start of I-PpoI cut sequence)
    endGRNA = startGRNA + len(cut_IPpoI); # index of gRNA end site (at end of I-PpoI cut sequence)
    plas.removeSeq([startGRNA, endGRNA]); # removes sequence that gRNA will replace
    plas.insertSeq("gg" + gRNA.upper(), startGRNA); # inserts gRNA sequence with gg sequence used by T7 polymerase
    annGRNA = GenBankAnn(geneName+" gRNA", "misc_feature", gRNA, False, [startGRNA+2,startGRNA+2+len(gRNA)], annColors['gRNAColor']); # annotation object. Note that gRNA starts after "gg" added for T7 polymerase
    plas.features.append(annGRNA); # adds annotation

    inRHR = findFirst(plas.origin, cut_ISceI); # index of RHR insertion site (at start of I-SceI cut sequence)
    plas.insertSeq(RHR.upper(), inRHR); # inserts RHR sequence
    annRHR = GenBankAnn(geneName+" RHR", "misc_feature", RHR, False, [inRHR,inRHR+len(RHR)], annColors['RHRColor']); # annotation object
    plas.features.append(annRHR); # adds annotation

    return plas; # returns modified plasmid



"""
Inserts targeting elements given as arguments into pSN150 at predetermined
sites. Elements given as arguments must be strings, not GenBankAnn objects.
The gRNA, LHR and RHR given must be in 5' to 3' sense and in the
positive strand and must not contain RE cut sites cut_FseI, cut_AsiSI,
cut_IPpoI, cut_ISceI, or cut_AflII. Returns new GenBank object with targeting
elements.
Specifically:
Inserts gRNA used by Cas9 in I-PpoI cut site. Removes recognition sequence. Adds
    GG at 5'-end of the gRNA sequence. The GG is required for the T7 RNA
    polymerase to efficiently transcribe the gRNA to high levels.
Inserts the LHR between FseI and AsiSI cut sites, but leaves the sites intact
    with their recognition sequences.
Inserts the recoded region after the LHR, before the AsiSI cut site (leaves the
    site intact).
Inserts RHR after I-SceI cut site, leaves the site intact with its recognition
    sequence.
"""
def insertTargetingElementsPSN150(plasmid, geneName, gRNA, LHR, recodedRegion, RHR, haTag=True, KO=False):
    plas = copy.deepcopy(plasmid); # makes a copy of the plasmid object to modify without altering the original

    inLHR = findFirst(plas.origin, cut_FseI) + len(cut_FseI); # index of LHR start site (at end of FseI cut sequence)
    plas.insertSeq(LHR.upper() + cut_FseI, inLHR); # inserts LHR sequence
    annLHR = GenBankAnn(geneName+" LHR", "misc_feature", LHR, False, [inLHR,inLHR+len(LHR)], annColors['LHRColor']); # annotation object
    plas.features.append(annLHR); # adds annotation

    if KO: # if knocking out,
        endRHR = findFirst(plas.origin, cut_XmaI); # index of RHR end site is start of XmaI cut sequence
        startRHR = findFirst(plas.origin, cut_AscI) + len(cut_AscI); # index of RHR end site is end of AscI cut sequence
        plas.removeSeq([startRHR, endRHR]); # removes sequence that LHR will replace
    else:
        endRHR = findFirst(plas.origin, cut_AhdI) + 6; # index of RHR end site is middle of AhdI cut sequence, conserves frame
        startRHR = endRHR; # assume keeping HA tag
        if not haTag: # if deleting HA tag,
            startRHR = findFirst(plas.origin, cut_NheI) + len(cut_NheI); # index of RHR start site (end of NheI cut sequence)
            plas.insertSeq(cut_NheI, startRHR); # inserts other NheI cut sequence downstream

    if len(recodedRegion) > 0: # if there is a recoded region,
        inRecode = startRHR; # index of recoded region start site (middle of AhdI cut sequence)
        plas.insertSeq(recodedRegion.upper(), inRecode); # inserts recoded region sequence
        annRecoded = GenBankAnn(geneName+" Recoded Region", "misc_feature", recodedRegion, False, [inRecode,inRecode+len(recodedRegion)], annColors['recodedRegionColor']); # annotation object
        if haTag: # if recoded region contains HA tag,
            annHATag = GenBankAnn("HA tag (recoded)", "misc_feature", recodedRegion[3:len(ha_tag)+3], False, [inRecode+3,inRecode+len(ha_tag)+3], annColors['otherAnnColor']); # annotation object for HA tag
            plas.features.append(annHATag); # adds annotation

        plas.features.append(annRecoded); # adds annotation
        startRHR = inRecode + len(recodedRegion); # shift RHR start site downstream of recoded region

    plas.insertSeq(RHR.upper(), startRHR); # inserts RHR sequence
    annRHR = GenBankAnn(geneName+" RHR", "misc_feature", RHR, False, [startRHR,startRHR+len(RHR)], annColors['RHRColor']); # annotation object
    plas.features.append(annRHR); # adds annotation

    startGRNA = findFirst(plas.origin, cut_IPpoI); # index of gRNA start site (at start of I-PpoI cut sequence)
    endGRNA = startGRNA + len(cut_IPpoI); # index of gRNA end site (at end of I-PpoI cut sequence)
    plas.removeSeq([startGRNA, endGRNA]); # removes sequence that gRNA will replace
    plas.insertSeq("gg" + gRNA.upper(), startGRNA); # inserts gRNA sequence with gg sequence used by T7 polymerase
    annGRNA = GenBankAnn(geneName+" gRNA", "misc_feature", gRNA, False, [startGRNA+2,startGRNA+2+len(gRNA)], annColors['gRNAColor']); # annotation object. Note that gRNA starts after "gg" added for T7 polymerase
    plas.features.append(annGRNA); # adds annotation

    return plas; # returns modified plasmid


"""
Inserts targeting elements given as arguments into a custom plasmid.
Elements given as arguments must be strings, not GenBankAnn objects.
Returns new GenBank object with targeting elements.
"""
def insertTargetingElementsCustom(plasmid, geneName, gRNA, LHR, recodedRegion, RHR):
    plas = copy.deepcopy(plasmid); # makes a copy of the plasmid object to modify without altering the original

    LHRcomp = False;
    RHRcomp = False;
    gRNAcomp = False;
    RRcomp = False;

    if len(LHR) == 0 and len(plas.findAnnsLabel("LHR")) > 0:
        ann = plas.findAnnsLabel("LHR")[0];
        plas.removeSeq([ann.index[0], ann.index[1]]);
    elif len(plas.findAnnsLabel("LHR")) > 0:
        if plas.findAnnsLabel("LHR")[0].comp:
            LHR = revComp(LHR)
            LHRcomp = True;

        inLHR = plas.findAnnsLabel("LHR")[0].index[0]; # index of LHR start site
        plas.removeSeq([inLHR, plas.findAnnsLabel("LHR")[0].index[1]]); # removes sequence that LHR will replace
        plas.insertSeq(LHR.upper(), inLHR); # inserts LHR sequence

        if LHRcomp:
            LHR = revComp(LHR)

        annLHR = GenBankAnn(geneName+" LHR", "misc_feature", LHR, LHRcomp, [inLHR,inLHR+len(LHR)], annColors['LHRColor']); # annotation object
        plas.features.append(annLHR); # adds annotation

    if len(RHR) == 0 and len(plas.findAnnsLabel("RHR")) > 0:
        ann = plas.findAnnsLabel("RHR")[0];
        plas.removeSeq([ann.index[0], ann.index[1]]);
    elif len(plas.findAnnsLabel("RHR")) > 0:
        if plas.findAnnsLabel("RHR")[0].comp:
            RHR = revComp(RHR)
            RHRcomp = True;

        inRHR = plas.findAnnsLabel("RHR")[0].index[0]; # index of LHR start site
        plas.removeSeq([inRHR, plas.findAnnsLabel("RHR")[0].index[1]]); # removes sequence that RHR will replace
        plas.insertSeq(RHR.upper(), inRHR); # inserts LHR sequence

        if RHRcomp:
            RHR = revComp(RHR)

        annRHR = GenBankAnn(geneName+" RHR", "misc_feature", RHR, RHRcomp, [inRHR,inRHR+len(RHR)], annColors['RHRColor']); # annotation object
        plas.features.append(annRHR); # adds annotation

    if len(gRNA) == 0 and len(plas.findAnnsLabel("sgRNA Sequence")) > 0:
        ann = plas.findAnnsLabel("sgRNA Sequence")[0];
        plas.removeSeq([ann.index[0], ann.index[1]]);
    elif len(plas.findAnnsLabel("sgRNA Sequence")) > 0:
        if plas.findAnnsLabel("sgRNA Sequence")[0].comp:
            gRNA = revComp(gRNA)
            gRNAcomp = True;

        ingRNA = plas.findAnnsLabel("sgRNA Sequence")[0].index[0]; # index of LHR start site
        plas.removeSeq([ingRNA, plas.findAnnsLabel("sgRNA Sequence")[0].index[1]]); # removes sequence that RHR will replace
        plas.insertSeq(gRNA.upper(), ingRNA); # inserts LHR sequence

        if gRNAcomp:
            gRNA = revComp(gRNA)

        annGRNA = GenBankAnn(geneName+" gRNA", "misc_feature", gRNA, gRNAcomp, [ingRNA,ingRNA+len(gRNA)], annColors['gRNAColor']); # annotation object. Note that gRNA starts after "gg" added for T7 polymerase
        plas.features.append(annGRNA); # adds annotation

    if len(recodedRegion) == 0 and len(plas.findAnnsLabel("Recoded Region")) > 0:
        ann = plas.findAnnsLabel("Recoded Region")[0];
        plas.removeSeq([ann.index[0], ann.index[1]]);
    elif len(plas.findAnnsLabel("Recoded Region")) > 0:
        if plas.findAnnsLabel("Recoded Region")[0].comp:
            recodedRegion = revComp(recodedRegion)
            RRcomp = True;

        inRR = plas.findAnnsLabel("Recoded Region")[0].index[0]; # index of RR start site
        plas.removeSeq([inRR, plas.findAnnsLabel("Recoded Region")[0].index[1]]); # removes sequence that RR will replace
        plas.insertSeq(recodedRegion.upper(), inRR); # inserts RR sequence

        if RRcomp:
            recodedRegion = revComp(recodedRegion)

        annRecoded = GenBankAnn(geneName+" Recoded Region", "misc_feature", recodedRegion, RRcomp, [inRR,inRR+len(recodedRegion)], annColors['recodedRegionColor']); # annotation object
        plas.features.append(annRecoded); # adds annotation

    return plas; # returns modified plasmid


"""
Inserts targeting elements given as arguments into plasmid at predetermined
sites. Elements given as arguments must be strings, not GenBankAnn objects.
The gRNA, LHR and RHR given must be in 5' to 3' sense and in the
positive strand and must not contain RE cut sites cut_FseI, cut_AsiSI,
cut_IPpoI, cut_ISceI, or cut_AflII. Returns new GenBank object with targeting
elements.
"""
def insertTargetingElements(plasmid, geneName, gRNA, LHR, recodedRegion, RHR, plasmidType, haTag = False):
    out = {}; # will contain method output
    if plasmidType == 'pSN054': # if using pSN054,
        out = insertTargetingElementsPSN054(plasmid, geneName, gRNA, LHR, recodedRegion, RHR); # use this method
    elif plasmidType == 'pSN150': # if using pSN150,
        out = insertTargetingElementsPSN150(plasmid, geneName, gRNA, LHR, recodedRegion, RHR, haTag); # use other method
    elif plasmidType == 'pSN150-KO': # if using pSN150-KO,
        out = insertTargetingElementsPSN150(plasmid, geneName, gRNA, LHR, recodedRegion, RHR, haTag, KO=True); # use other method
    elif plasmidType == 'custom':
        out = insertTargetingElementsCustom(plasmid, geneName, gRNA, LHR, recodedRegion, RHR); # use other method

    return out;


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
        log = log + "Found LHR annotation for locus editing." + "\n\n"; # add warning to log
    elif len(LHRlist) == 0: # if no LHR found,
        log = log + "ERROR: Did not find LHR annotations, genomic context file could not be built." + "\n\n"; # add warning to log
    else:
        LHR = LHRlist[0]; # save first LHR annotation as RHR
        log = log + "Warning: Found multiple LHR annotations, genomic context file built with the first one found." + "\n\n"; # add warning to log

    RHRlistAll = construct.findAnnsLabel("RHR"); # saves RHR annotation
    RHRlist = []; # contain only exact matches
    RHR = ""; # init LHR var
    for ann in RHRlistAll: # will keep only those exact matches (exclude primers, for example)
        if ann.label == geneName + " RHR":
            RHRlist.append(ann);

    if len(RHRlist) == 1: # if LHR found,
        RHR = RHRlist[0]; # save first RHR annotation as RHR
        log = log + "Found RHR annotation for locus editing." + "\n\n"; # add warning to log
    elif len(RHRlist) == 0: # if no LHR found,
        log = log + "ERROR: Did not find RHR annotations, genomic context file could not be built." + "\n\n"; # add warning to log
    else:
        RHR = RHRlist[0]; # save first RHR annotation as RHR
        log = log + "Warning: Found multiple RHR annotations, genomic context file built with the first one found." + "\n\n"; # add warning to log

    startInsertChrom = findFirst(gene.origin, LHR.seq) + len(LHR.seq); # find index of insertion start in the chromosomal DNA given
    endInsertChrom = findFirst(gene.origin, RHR.seq); # find index of insertion end in the chromosomal DNA given

    startInsertConst = findFirst(construct.origin, LHR.seq) + len(LHR.seq); # find index of insertion start in the chromosomal DNA given
    endInsertConst = findFirst(construct.origin, RHR.seq); # find index of insertion end in the chromosomal DNA given
    insert = construct.origin[startInsertConst:endInsertConst]; # saves string containing insert

    editedLocus = copy.deepcopy(gene); # var will contain the edited locus
    editedLocus.removeSeq([startInsertChrom,endInsertChrom]); # remove deleted bases between homologous regions
    editedLocus.insertSeq(insert, startInsertChrom); # inserts construct into edited locus

    # Now we add all the plasmid's annotations to the edited locus
    for ann in construct.features: # loop through plasmid annotations
        if ann.index[0] >= startInsertConst and ann.index[1] <= endInsertConst: # if annotation is completely inside integrating region
            editedAnn = copy.deepcopy(ann); # annotation on edited chromosome
            editedAnn.index[0] = editedAnn.index[0] + startInsertChrom - startInsertConst; # shift index according to new context
            editedAnn.index[1] = editedAnn.index[1] + startInsertChrom - startInsertConst; # shift index according to new context
            editedLocus.features.append(editedAnn); # adds annotation to edited locus


    log += "Edited locus (genomic context) file built.\n\n"; # add log entry

    return {"out":editedLocus, "log":log}; # return dictionary
