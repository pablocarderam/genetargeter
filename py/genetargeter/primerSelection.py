
from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants

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

    annPrimFwd = GenBankAnn(part.label + " Primer (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF], annColors['primerColor']); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Primer (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR], annColors['primerColor']); # creates GenBankAnn object to hold rev primer

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

    if meltingTemp(plasmid.origin[part.index[0]:endPF]) < minMeltTemp or not primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C": # if still no use
        startPF = part.index[0] - rangeHom[0]; # Smallest fwd primer start position
        endPF = part.index[0] + rangeHom[0]; # Smallest fwd primer end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence

        maxIndexes = [startPF, endPF]; # store start and end positions of best primer in search range
        while (meltingTemp(plasmid.origin[part.index[0]:endPF]) < minMeltTemp or not primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C") and len(primFwdSeq) <= rangeHom[2]*2: # while still no suitable Tm found and still within length parameters,
            startPF = startPF - 1; # shift primer start position upstream
            endPF = endPF + 1; # shift primer start position upstream
            primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
            if (meltingTemp(plasmid.origin[part.index[0]:endPF]) > meltingTemp(plasmid.origin[part.index[0]:maxIndexes[1]]) or not plasmid.origin[maxIndexes[1]-1].upper().replace("G","C") == "C") and primFwdSeq[len(primFwdSeq)-1].upper().replace("G","C") == "C": # if this primer has higher Tm than the max or the current max has no gc clamp, and this one does have a gc clamp,
                maxIndexes = [startPF, endPF]; # store start and end positions of this primer

        startPF = maxIndexes[0]; # Fwd primer default start position
        endPF = maxIndexes[1]; # Fwd primer default end position
        primFwdSeq = plasmid.origin[startPF:endPF]; # Fwd primer sequence
        if meltingTemp(plasmid.origin[part.index[0]:endPF]) < minMeltTemp: # if still no use
            log = log + "\nWarning: Best Gibson fwd primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(plasmid.origin[part.index[0]:endPF])) + ", below the given threshold of " + str(minMeltTemp) + "\n"; # give warning


    startPR = part.index[1] - rangeHom[1]; # Rev primer start position
    endPR = part.index[1] + rangeHom[1]; # Rev primer end position
    primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

    if meltingTemp(plasmid.origin[startPR:part.index[1]]) < minMeltTemp or meltingTemp(plasmid.origin[startPR:part.index[1]])-meltingTemp(plasmid.origin[startPR:part.index[1]]) > maxTempDif or not primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C": # if still no use
        startPR = part.index[1] - rangeHom[0]; # Smallest fwd primer start position
        endPR = part.index[1] + rangeHom[0]; # Smallest fwd primer end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence

        maxIndexes = [startPR, endPR]; # store start and end positions of best primer in search range
        while (meltingTemp(plasmid.origin[startPR:part.index[1]]) < minMeltTemp or meltingTemp(plasmid.origin[startPR:part.index[1]])-meltingTemp(plasmid.origin[startPR:part.index[1]]) > maxTempDif or not primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C") and rangeHom[0]*2 <= len(primRevSeq) <= rangeHom[2]*2: # while still no suitable Tm found and still within length parameters,
            startPR = startPR - 1; # shift primer start position upstream
            endPR = endPR + 1; # shift primer start position upstream
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            if meltingTemp(plasmid.origin[startPR:part.index[1]]) > meltingTemp(plasmid.origin[maxIndexes[0]:part.index[1]]) and primRevSeq[len(primRevSeq)-1].upper().replace("G","C") == "C": # if this primer has higher Tm than the max and has gc clamp,
                maxIndexes = [startPR, endPR]; # store start and end positions of this primer

        startPR = maxIndexes[0]; # Rev primer default start position
        endPR = maxIndexes[1]; # Rev primer default end position
        primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
        if meltingTemp(plasmid.origin[startPR:part.index[1]]) < minMeltTemp: # if still no use
            log = log + "\nWarning: Best Gibson rev primer for sequence " + part.label + " under given constraints has a Tm of " + str(meltingTemp(plasmid.origin[startPR:part.index[1]])) + ", below the given threshold of " + str(minMeltTemp) + "\n"; # give warning
        elif meltingTemp(plasmid.origin[part.index[0]:endPF])-meltingTemp(plasmid.origin[startPR:part.index[1]]) > maxTempDif: # if temp difference exceeds specs
            startPR = maxIndexes[0]; # Rev primer default start position
            endPR = maxIndexes[1]; # Rev primer default end position
            primRevSeq = revComp(plasmid.origin[startPR:endPR]); # Rev primer sequence
            primFwdSeq = primFwdSeq.upper(); # to uppercase
            lastBase = endPF-2; # stores possible end points of fwd primer
            while not plasmid.origin[lastBase-1].upper().replace("G","C") == "C" and lastBase-part.index[0] > rangeHom[0]: # find next G or C upstream
                lastBase -= 1;

            while meltingTemp(plasmid.origin[part.index[0]:lastBase])-meltingTemp(plasmid.origin[startPR:part.index[1]]) > maxTempDif and meltingTemp(plasmid.origin[part.index[0]:lastBase]) > minMeltTemp and lastBase-part.index[0] > rangeHom[0]: # while T diff is still out of bounds and still within bounds of Fwd primer,
                lastBase -= 1;
                while not plasmid.origin[lastBase-1].upper().replace("G","C") == "C" and lastBase-part.index[0] > rangeHom[0]: # find next G or C upstream
                    lastBase -= 1;


            if meltingTemp(plasmid.origin[part.index[0]:lastBase])-meltingTemp(plasmid.origin[startPR:part.index[1]]) < maxTempDif and meltingTemp(plasmid.origin[0:lastBase]) > minMeltTemp: # while T diff is still out of bounds and still within bounds of Fwd primer,
                endPF = lastBase;
                primFwdSeq = plasmid.origin[startPF:endPF];
            else: # if temp difference still exceeds specs
                log = log + "\nWarning: Gibson primers for sequence " + part.label + " under given constraints have a Tm difference of " + str(meltingTemp(plasmid.origin[part.index[0]:endPF])-meltingTemp(plasmid.origin[startPR:part.index[1]])) + ", above the given threshold of " + str(maxTempDif) + "\n"; # give warning

    annPrimFwd = GenBankAnn(part.label + " Gibson Primer (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF], annColors['primerColor']); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Gibson Primer (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR], annColors['primerColor']); # creates GenBankAnn object to hold rev primer

    log = log + "Gibson primers for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":[annPrimFwd, annPrimRev],"log":log}; # return list of primers


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

    annPrimFwd = GenBankAnn(part.label + " Klenow oligo (Fwd)", "misc_feature", primFwdSeq, False, [startPF,endPF], annColors['primerColor']); # creates GenBankAnn object to hold fwd primer
    annPrimRev = GenBankAnn(part.label + " Klenow oligo (Rev)", "misc_feature", primRevSeq, True, [startPR,endPR], annColors['primerColor']); # creates GenBankAnn object to hold rev primer

    log = log + "Klenow oligos for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":[annPrimFwd, annPrimRev], "log":log}; # return list of primers


"""
Returns gBlock GenBankAnn object in plasmid given for part (a GenBankAnn object)
to be synthesized and inserted via Gibson in plasmid. Will give a warning if it
suspects the gBlock won't be able to be synthesized by IDT, reverse-engineering
from the IDT gene synthesis webpage.
"""
def createGBlock(plasmid, part, overlapSize):
    log = ""; # init log
    startGBlock = part.index[0] - overlapSize; # gBlock start position
    endGBlock = part.index[1] + overlapSize; # gBlock end position
    gBlockSeq = plasmid.origin[startGBlock:endGBlock]; # gBlock sequence
    tricky = False; # True if suspected to be hard to synthesize
    '''
    for i in range(0,len(gBlockSeq)-20):
        if gcContent(gBlockSeq[i:i+20]) < 0.05: # If gc content of 20 bp gBlock window is too low or tricky sequences are present
            tricky = True; # might be tricky'''

    if findFirst(gBlockSeq.replace("T", "A"),"AAAAAAAAAAAAAAAAAAAA") > -1: # if 20 continuous T/A nucleotides found,
        tricky = True; # it's tricky
    elif findFirst(gBlockSeq.replace("G", "C"),"CCCCCCCCCCCCCC") > -1: # if 14 continuous G/C nucleotides found,
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

    annGBlock = GenBankAnn(part.label + " gBlock", "misc_feature", gBlockSeq, False, [startGBlock,endGBlock], annColors['gBlockColor']); # creates GenBankAnn object to hold gBlock

    log = log + "gBlock for part " + part.label + " selected.\n\n"; # logs this process finished
    return {"out":annGBlock, "log":log}; # returns annotation


"""
Abbreviates primer names to fit on commercial tube labels with the format:
Seven Digit Gene Identifier_Oligo type_Orientation

The seven digit gene identifier code follows "PF3D7_".
Oligo types include:
LHR (Gibson overhang PCR primers)
RHR (Gibson overhang PCR primers)
gRNA (Klenow oligos for gRNA sequence)
gBlock (gBlock sequencing primer)
RecAE (Anneal-extension oligos for recoded region, if the region is small enough)

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
                newName = newName + "RecAE" + "_";

            # Add oligo orientation
            if name.find("fwd") > -1:
                newName = newName + "F";
            elif name.find("rev") > -1:
                newName = newName + "R";

            primer[0] = newName; # replace name


    mat = [",".join(l) for l in mat]; # join cells into lines
    outStr = "\n".join(mat); # join lines into string

    return outStr;
