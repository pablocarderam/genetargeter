
from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants

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
def chooseLHR3Prime(geneGB, gene, lengthLHR=[450,500,650], minTmEnds=55, endsLength=40, optimizeRange=[-20,10], maxDistanceFromGRNA=500, gBlockDefault=True, minGBlockSize=125, codingGene=True, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI]):
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



    endLHR = min(gRNAUpstream.index[0],gene.index[1]-(3*codingGene)); # saves end index of LHR as whatever is more upstream between the start of the gRNA or the end of the gene (minus the stop codon) (in Python indexing, i.e. not included in LHR).
    if gBlockDefault and ((endLHR < gene.index[1]-(3*codingGene) and endLHR > gene.index[1]-(3*codingGene) - minGBlockSize) or (meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds and gene.index[1]-(3*codingGene)-endLHR < minGBlockSize)): # if defaulting to a gBlock for any recoded region and a) there is a recoded region, and it is under the minimum gBlock size, or b) if the LHR is supposed to end at the end of the gene, but this region does not make for a good LHR ending point,
        endLHR = gene.index[1]-(3*codingGene) - minGBlockSize; # extend recoded region to minimum gBlock size

    while (not geneGB.checkInExon(endLHR)) and endLHR > lengthLHR[2]: # Loop as long as the end of LHR is not in an exon and the end of the LHR is inside the max length
        endLHR -= 1; # shift LHR end upstream one bp

    startLHR = max(endLHR - lengthLHR[1],0); # stores LHR start index according to preferred length, or 0 if not enough space before gene
    failSearchStart = False; # true if no suitable LHR end region is found within given parameters.

    # First search for end region
    while meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchStart and geneGB.checkInExon(endLHR-1): # while no suitable end region found, still within max distance from gRNA, the search for a suitable end region hasn't failed yet, and still within exon,
        endLHR -= 1; # shift endLHR upstream

    if meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds: # if no suitable end region found,
        endLHR = min(gRNAUpstream.index[0],gene.index[1]-(3*codingGene)); # saves end index of LHR by default as whatever is more upstream between the start of the gRNA or the end of the gene (minus the stop codon)
        while not geneGB.checkInExon(endLHR) and endLHR > lengthLHR[2]: # Loop as long as the end of LHR is not in an exon and the end of the LHR is inside the max length
            endLHR -= 1; # shift LHR end upstream one bp

        failSearchStart = True; # changes failSearchStart status
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of LHR and gRNA. \nDefaulted to ending right before start of gRNA most upstream." + "\n"; # give a warning

    # Then search for start region; modify end region if necessary
    while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and gRNAUpstream.index[0]-endLHR <= maxDistanceFromGRNA and geneGB.checkInExon(endLHR-1): # if no start region has been found, and still within max distance from gRNA and inside exon,
        # find new end region if necessary
        while meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchStart and geneGB.checkInExon(endLHR-1): # while no suitable end region found, still within max distance from gRNA, the search for a suitable end region hasn't failed yet, and still within exon,
            endLHR -= 1; # shift endLHR upstream

        if gBlockDefault and gene.index[1]-(3*codingGene)-endLHR < minGBlockSize: # if defaulting to a gBlock for any recoded region and the current recoded region is smaller than the minimum gBlock size,
            endLHR = gene.index[1]-(3*codingGene) - minGBlockSize; # extend recoded region to minimum gBlock size

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
        # resets endLHR
        endLHR = min(gRNAUpstream.index[0],gene.index[1]-(3*codingGene)); # saves end index of LHR as whatever is more upstream between the start of the gRNA or the end of the gene (minus the stop codon) (in Python indexing, i.e. not included in LHR).
        if gBlockDefault and endLHR < gene.index[1]-(3*codingGene) and endLHR > gene.index[1]-(3*codingGene) - minGBlockSize: # if defaulting to a gBlock for any recoded region, there is a recoded region, and it is under the minimum gBlock size,
            endLHR = gene.index[1]-(3*codingGene) - minGBlockSize; # extend recoded region to minimum gBlock size

        while meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds and gRNAUpstream.index[0]-endLHR < maxDistanceFromGRNA and not failSearchStart and geneGB.checkInExon(endLHR-1) and endLHR > lengthLHR[2]: # while no suitable end region found, still within max distance from gRNA, the search for a suitable end region hasn't failed yet, still within exon, and the end of the LHR is inside the max length
            endLHR -= 1; # shift endLHR upstream

        startLHR = endLHR - lengthLHR[1]; # stores LHR start index according to preferred length by default
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of LHR and gRNA. \nDefaulted to starting so that the LHR size is equal to preferred size of " + str(lengthLHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting LHR by adjusting start sequence around the given range
    searchIndexes = [int(startLHR+optimizeRange[0]), int(startLHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexes[0], searchIndexes[1]): # iterates across optimization range
        if i >= endsLength: # if inside gene file,
            if meltingTemp(geneGB.origin[i:(i+endsLength)]) > meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) and lengthLHR[0] >= endLHR-i >= lengthLHR[2]: # if this start point has a better Tm and is still within bounds,
                startLHR = i; # make this the starting position



    searchIndexesEnd = [int(endLHR+optimizeRange[0]), min(int(endLHR+optimizeRange[1]),gene.index[1]-(3*codingGene))]; # list with indexes across which we are going to search for a better end point, can't be after end of gene (or start of stop codon, to be precise)
    for i in range(searchIndexesEnd[0], searchIndexesEnd[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[max(i-endsLength,0):i]) > meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) and lengthLHR[2] >= i-startLHR >= lengthLHR[0] and i < gRNAUpstream.index[0] and geneGB.checkInExon(i) and i <= gene.index[1]-(3*codingGene) and (not gBlockDefault or (i <= gene.index[1]-(3*codingGene) - minGBlockSize) or (i == gene.index[1]-(3*codingGene)) ): # if this end point has a better Tm, and is still within bounds, before gRNA and within exon, and if it allows for a gBlock of the minimum size or no gBlock at all if gBlocks are being used as default, and not in stop codon,
            endLHR = i; # make this the ending position




    for site in filterCutSites: # for every cut site being filtered
        if findFirst(site,geneGB.origin[startLHR:endLHR]) > -1 or findFirst(revComp(site),geneGB.origin[startLHR:endLHR]) > -1: # if cut site found,
            log = log + "\nWarning: LHR sequence for gene " + gene.label + ": \n" + geneGB.origin[startLHR:endLHR] + "\ncontains restriction site " + site + "\n"; # add warning to log


    # creates var to store finished LHR as annotation
    LHR = GenBankAnn(gene.label + " LHR", "misc_feature", geneGB.origin[startLHR:endLHR], False, [startLHR,endLHR], annColors['LHRColor']); # creates GenBankAnn object to hold LHR

    log = log + "LHR for gene " + gene.label + " selected.\n\n"; # logs this process finished
    return {"out":LHR, "log":log}; # returns LHR GenBankAnn object


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
def chooseLHR5Prime(geneGB, gene, lengthLHR=[450,500,750], minTmEnds=59, endsLength=40, optimizeRange=[-20,20], maxDistanceFromGene=500, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI]):
    #TODO: debug cases in which LHR has to be shifted
    log = ""; # init log
    gRNAs = []; # List of all gRNAs
    gRNAUpstream = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
        gRNAUpstream = gRNAs[0]; # will store gRNA most upstream
        for g in gRNAs: # loops across gRNAs
            if g.index[0] < gRNAUpstream.index[0]: # if more upstream
                gRNAUpstream = g; # replace as most upstream



    genes = geneGB.findAnnsType("gene"); # list of all genes, used to verify LHR doesn't start inside any genes (truncating them)

    endLHR = min(gene.index[0], gRNAUpstream.index[0]); # saves end index of LHR as one bp before gene or before most upstream gRNA.
    failSearchEnd = False; # true if no suitable LHR start region is found within given parameters.

    # First search for end region
    while meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds and gene.index[0] - endLHR <= maxDistanceFromGene and endLHR - lengthLHR[0] > 0 and not geneGB.checkInExon(endLHR-1): # while no suitable start region found, still within max distance from gene, not beyond the border of the chromosome, the search for a suitable start region hasn't failed yet, and start of LHR is not downstream of the start of a gene,
        endLHR -= 1; # shift endLHR upstream

    if meltingTemp(geneGB.origin[max(endLHR-endsLength,0):endLHR]) < minTmEnds and ( gene.index[0] - endLHR > maxDistanceFromGene or geneGB.checkInExon(endLHR) or endLHR - lengthLHR[0] <= 0 ): # if no suitable start region found,
        endLHR = min(gene.index[0], gRNAUpstream.index[0]); # saves start index of LHR as one bp before gene or 1 bp before most downstream gRNA by default
        failSearchEnd = True; # changes failSearchEnd status
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of LHR and start of gene. \nDefaulted to ending right before start of gene or most upstream gRNA." + "\n"; # give a warning

    startLHR = max(endLHR - lengthLHR[1],0); # stores LHR start index according to preferred length, or the first base in file if not enough space before gene
    # search for start upstream
    while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR <= lengthLHR[2] and startLHR > 0: # while no suitable end region found, still within max length of LHR, and still inside gene file,
        startLHR -= 1; # shift startLHR upstream

    # if not found upstream
    if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR > lengthLHR[2]: # if no end region found this way,
        startLHR = endLHR - lengthLHR[1]; # return to preferred position
        # search upstream
        while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR >= lengthLHR[0]: # while no suitable start region found and still within min length of LHR,
            startLHR += 1; # shift startLHR downstream


    if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and gene.index[0]-endLHR <= maxDistanceFromGene and not geneGB.checkInExon(endLHR) and endLHR-lengthLHR[0] > 0: # if still not found and within bounds and start of LHR is not upstream of the start of a gene and not too close to start of gene file,
        # Then modify end region, search for start region;
        while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and gene.index[0]-endLHR <= maxDistanceFromGene and not geneGB.checkInExon(endLHR-1): # if no end region has been found and start of LHR is not downstream of the start of a gene,
            endLHR -= 1; # shifts end of LHR upstream

            # find new end region if necessary
            while meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) < minTmEnds and gene.index[0] - endLHR <= maxDistanceFromGene and not geneGB.checkInExon(endLHR-1) and not failSearchEnd: # while no suitable start region found, still within max distance from gene, the search for a suitable start region hasn't failed yet, and start of LHR is not downstream of the start of a gene,
                endLHR -= 1; # shift endLHR upstream

            # search for start upstream
            while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR <= lengthLHR[2] and startLHR > 0: # while no suitable end region found and still within max length of LHR and still inside gene file,
                startLHR -= 1; # shift startLHR upstream

            # if not found upstream
            if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR > lengthLHR[2]: # if no end region found this way,
                startLHR = endLHR - lengthLHR[1]; # return to preferred position
                # search downstream
                while meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and endLHR-startLHR >= lengthLHR[0]: # while no suitable start region found and still within min length of LHR,
                    startLHR += 1; # shift startLHR downstream



    if meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) < minTmEnds and ( gene.index[0] - endLHR > maxDistanceFromGene or geneGB.checkInExon(endLHR) ): # if no suitable start region found,
        endLHR = min(gene.index[0], gRNAUpstream.index[0]); # saves end index of LHR as one bp before gene or 1 bp before most upstream gRNA by default
        while meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) < minTmEnds and gene.index[0]-endLHR <= maxDistanceFromGene and not failSearchEnd and not geneGB.checkInExon(endLHR-1): # while no suitable start region found, still within max distance from gene, the search for a suitable start region hasn't failed yet, and start of LHR is not downstream of the start of a gene,
            endLHR -= 1; # shift endLHR upstream

        startLHR = max([endLHR - lengthLHR[1],0]); # stores LHR end index according to preferred length by default, or max length available if there is not enough sequence space on the chromosome before the gene start
        log = log + "\nWarning: No LHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of LHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of LHR and start of gene. \nDefaulted to starting at preferred size of " + str(lengthLHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting LHR by adjusting start and stop sequences around the given range
    searchIndexesEnd = [int(startLHR+optimizeRange[0]), int(startLHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexesEnd[1], searchIndexesEnd[0]): # iterates across optimization range
        if meltingTemp(geneGB.origin[ max([i,0]):min([i+endsLength,endsLength]) ]) > meltingTemp(geneGB.origin[startLHR:(startLHR+endsLength)]) and lengthLHR[2] >= endLHR-i >= lengthLHR[0]: # if this end point has a better Tm and is still within bounds,
            if i > 0: # if inside gene file,
                startLHR = i; # make this the starting position




    searchIndexesStart = [int(endLHR+optimizeRange[0]), int(endLHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better end point.
    for i in range(searchIndexesStart[0], searchIndexesStart[1], -1): # iterates across optimization range in reverse
        if meltingTemp(geneGB.origin[(i-endsLength):i]) > meltingTemp(geneGB.origin[(endLHR-endsLength):endLHR]) and lengthLHR[2] >= i-startLHR >= lengthLHR[0] and i < gRNAUpstream.index[0] and not geneGB.checkInExon(i): # if this start point has a better Tm and is still within bounds and after gRNA and not after a start of a gene,
            endLHR = i; # make this the ending position


    if endLHR-startLHR < lengthLHR[0]+searchIndexesEnd[0]-searchIndexesEnd[1]: # if LHR is smaller than it should be (gene too close to the end of chromosome, for example)
        log = log + "\ERROR: LHR found for gene " + geneGB.name + " \nis smaller than allowed by design, with a length of " + str(endLHR-startLHR) + " bp. Maybe gene was too close to start/end of chromosome?\n"; # give a warning

    for site in filterCutSites: # for every cut site being filtered
        if findFirst(site,geneGB.origin[startLHR:endLHR]) > -1 or findFirst(revComp(site),geneGB.origin[startLHR:endLHR]) > -1: # if cut site found,
            log = log + "\nWarning: LHR sequence for gene " + gene.label + ": \n" + geneGB.origin[startLHR:endLHR] + "\ncontains restriction site " + site + "\n"; # add warning to log

    # creates var to store finished LHR as annotation
    LHR = GenBankAnn(gene.label + " LHR", "misc_feature", geneGB.origin[startLHR:endLHR], False, [startLHR,endLHR], annColors['LHRColor']); # creates GenBankAnn object to hold LHR

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
def chooseRHR3Prime(geneGB, gene, lengthRHR=[450,500,750], minTmEnds=59, endsLength=40, optimizeRange=[-20,20], maxDistanceFromGene=500, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI]):
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



    genes = geneGB.findAnnsType("gene"); # list of all genes, used to verify RHR doesn't start inside any genes (truncating them)

    startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or after most downstream gRNA.
    failSearchStart = False; # true if no suitable RHR start region is found within given parameters.

    # First search for start region
    while meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) < minTmEnds and startRHR - gene.index[1] <= maxDistanceFromGene and startRHR + lengthRHR[0] < len(geneGB.origin) and not geneGB.checkInExon(startRHR): # while no suitable start region found, still within max distance from gene, not beyond the border of the chromosome, the search for a suitable start region hasn't failed yet, and start of RHR is not downstream of the start of a gene,
        startRHR += 1; # shift startRHR downstream

    if meltingTemp(geneGB.origin[startRHR:(min(startRHR+endsLength,len(geneGB.origin))+endsLength)]) < minTmEnds and ( startRHR - gene.index[1] > maxDistanceFromGene or geneGB.checkInExon(startRHR-1) ): # if no suitable start region found,
        startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or 1 bp after most downstream gRNA by default
        failSearchStart = True; # changes failSearchStart status
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between end of gene and start of RHR. \nDefaulted to starting right after end of gene or most downstream gRNA." + "\n"; # give a warning

    endRHR = min(startRHR + lengthRHR[1],len(geneGB.origin)-1); # stores RHR start index according to preferred length, or the last base in file if not enough space after gene
    # search for end downstream
    while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2] and endRHR < len(geneGB.origin)-1: # while no suitable end region found, still within max length of RHR, and still inside gene file,
        endRHR += 1; # shift endRHR downstream

    # if not found downstream
    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
        endRHR = startRHR + lengthRHR[1]; # return to preferred position
        # search upstream
        while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR >= lengthRHR[0]: # while no suitable start region found and still within min length of RHR,
            endRHR -= 1; # shift endRHR upstream


    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and not geneGB.checkInExon(startRHR-1) and startRHR+lengthRHR[0] < len(geneGB.origin)-1: # if still not found and within bounds and start of RHR is not downstream of the start of a gene and not too close to end of gene file,
        # Then modify start region, search for end region;
        while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and not geneGB.checkInExon(startRHR): # if no end region has been found and start of RHR is not downstream of the start of a gene,
            startRHR += 1; # shifts start of RHR downstream

            # find new start region if necessary
            while meltingTemp(geneGB.origin[startRHR:(startRHR+endsLength)]) < minTmEnds and startRHR - gene.index[1] <= maxDistanceFromGene and not geneGB.checkInExon(startRHR) and not failSearchStart: # while no suitable start region found, still within max distance from gene, the search for a suitable start region hasn't failed yet, and start of RHR is not downstream of the start of a gene,
                startRHR += 1; # shift startRHR downstream

            # search for end downstream
            while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2] and endRHR < len(geneGB.origin)-1: # while no suitable end region found and still within max length of RHR and still inside gene file,
                endRHR += 1; # shift endRHR downstream

            # if not found downstream
            if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
                endRHR = startRHR + lengthRHR[1]; # return to preferred position
                # search upstream
                while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR >= lengthRHR[0]: # while no suitable start region found and still within min length of RHR,
                    endRHR -= 1; # shift endRHR upstream



    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and ( startRHR-gene.index[1] > maxDistanceFromGene or geneGB.checkInExon(startRHR-1) ): # if no suitable end region found,
        startRHR = max(gene.index[1], gRNADownstream.index[1]); # saves start index of RHR as one bp after gene or 1 bp after most downstream gRNA by default
        while meltingTemp(geneGB.origin[startRHR:(startRHR+endsLength)]) < minTmEnds and startRHR-gene.index[1] <= maxDistanceFromGene and not failSearchStart and not geneGB.checkInExon(startRHR): # while no suitable start region found, still within max distance from gene, the search for a suitable start region hasn't failed yet, and start of RHR is not downstream of the start of a gene,
            startRHR += 1; # shift startRHR downstream

        endRHR = min([startRHR + lengthRHR[1],len(geneGB.origin)]); # stores RHR end index according to preferred length by default, or max length available if there is not enough sequence space on the chromosome after the gene end
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGene) + " bp between start of RHR and end of gene. \nDefaulted to ending at preferred size of " + str(lengthRHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting RHR by adjusting start and stop sequences around the given range
    searchIndexesEnd = [int(endRHR+optimizeRange[0]), int(endRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better end point.
    for i in range(searchIndexesEnd[1], searchIndexesEnd[0], -1): # iterates across optimization range in reverse
        if meltingTemp(geneGB.origin[ min([i-endsLength,len(geneGB.origin)-endsLength]):min([i,len(geneGB.origin)]) ]) > meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) and lengthRHR[2] >= i-startRHR >= lengthRHR[0]: # if this end point has a better Tm and is still within bounds,
            if i < len(geneGB.origin): # if inside gene file,
                endRHR = i; # make this the ending position




    searchIndexesStart = [int(startRHR+optimizeRange[0]), int(startRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point.
    for i in range(searchIndexesStart[0], searchIndexesStart[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[i:(i+endsLength)]) > meltingTemp(geneGB.origin[startRHR:(startRHR+endsLength)]) and lengthRHR[2] >= endRHR-i >= lengthRHR[0] and i >= gRNADownstream.index[1] and not geneGB.checkInExon(i-1): # if this start point has a better Tm and is still within bounds and after gRNA and not after a start of a gene,
            startRHR = i; # make this the starting position


    if endRHR-startRHR < lengthRHR[0]+searchIndexesEnd[0]-searchIndexesEnd[1]: # if RHR is smaller than it should be (gene too close to the end of chromosome, for example)
        log = log + "\ERROR: RHR found for gene " + geneGB.name + " \nis smaller than allowed by design, with a length of " + str(endRHR-startRHR) + " bp. Maybe gene was too close to start/end of chromosome?\n"; # give a warning

    for site in filterCutSites: # for every cut site being filtered
        if findFirst(site,geneGB.origin[startRHR:endRHR]) > -1 or findFirst(revComp(site),geneGB.origin[startRHR:endRHR]) > -1: # if cut site found,
            log = log + "\nWarning: RHR sequence for gene " + gene.label + ": \n" + geneGB.origin[startRHR:endRHR] + "\ncontains restriction site " + site + "\n"; # add warning to log

    # creates var to store finished RHR as annotation
    RHR = GenBankAnn(gene.label + " RHR", "misc_feature", geneGB.origin[startRHR:endRHR], False, [startRHR,endRHR], annColors['RHRColor']); # creates GenBankAnn object to hold RHR

    log = log + "RHR for gene " + gene.label + " selected.\n\n"; # logs this process finished

    return {"out":RHR, "log":log}; # returns RHR GenBankAnn object


"""
Selects an appropriate RHR for the given gene. GenBank sequence object given as
argument (geneGB) must have an annotation with label=gene, at least one
annotation with "gRNA" in label (or some similar indicator), and at least
lengthRHR[2] + maxDistanceFromGene - optimizeRange[0] + optimizeRange[1]
bp of 3' UTR (about 1000 bp is fine). lengthRHR is [min, preferred, max]
length in bp of RHR. minTmEnds is min melting temp in extremes of RHR.
endsLength is length of extremes of the RHR contemplated in analysis.
optimizeRange is a list with the index positions relative to the start of the
RHR over which the start of the RHR is optimized once it is chosen.
maxDistanceFromGRNA is maximum distance in bp between end of RHR and start of
gRNA.
RHR: Right Homologous Region used for chromosomal integration by homologous
recombination during repair.
"""
def chooseRHR5Prime(geneGB, gene, lengthRHR=[450,500,650], minTmEnds=55, endsLength=40, optimizeRange=[-20,10], maxDistanceFromGRNA=500, gBlockDefault=True, minGBlockSize=125, codingGene=True, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI]):
    #TODO: debug cases in which RHR has to be shifted, do case in which RHR is in intron?
    log = ""; # init log
    gRNAs = []; # List of all gRNAs
    gRNADownstream = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
        gRNADownstream = gRNAs[0]; # will store gRNA most downstream
        for g in gRNAs: # loops across gRNAs
            if g.index[0] > gRNADownstream.index[0]: # if more downstream
                gRNADownstream = g; # replace as most downstream



    startRHR = max(gRNADownstream.index[1],gene.index[0]); # saves start index of RHR as whatever is more downstream between the end of the gRNA or the start of the gene (in Python indexing, i.e. not included in RHR).
    if gBlockDefault and ((startRHR > gene.index[0] and startRHR < gene.index[0] - minGBlockSize) or (meltingTemp(geneGB.origin[startRHR:startRHR+endsLength]) < minTmEnds and startRHR-gene.index[0] < minGBlockSize)): # if defaulting to a gBlock for any recoded region and a) there is a recoded region, and it is under the minimum gBlock size, or b) if the RHR is supposed to start at the start of the gene, but this region does not make for a good RHR starting point,
        startRHR = gene.index[0] + minGBlockSize; # extend recoded region to minimum gBlock size

    while (not geneGB.checkInExon(startRHR)) and startRHR < gene.index[1] and startRHR < len(geneGB.origin)-lengthRHR[0]: # Loop as long as the start of RHR is not in an exon and the start of the RHR is inside the gene and there is enough room for RHR in file
        startRHR += 1; # shift RHR start downstream one bp

    endRHR = min(startRHR + lengthRHR[1],len(geneGB.origin)); # stores RHR start index according to preferred length, or end of file if not enough space after gene
    failSearchEnd = False; # true if no suitable RHR end region is found within given parameters.

    # First search for start region
    while meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) < minTmEnds and startRHR-gRNADownstream.index[1] < maxDistanceFromGRNA and not failSearchEnd and geneGB.checkInExon(startRHR) and startRHR < gene.index[1]: # while no suitable start region found, still within max distance from gRNA, the search for a suitable start region hasn't failed yet, and still within exon,
        startRHR += 1; # shift startRHR downstream

    if meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) < minTmEnds: # if no suitable start region found,
        startRHR = max(gRNADownstream.index[1],gene.index[0]); # saves end index of RHR by default as whatever is more downstream between the start of the gRNA or the start of the gene
        while not geneGB.checkInExon(startRHR) and startRHR < gene.index[1] and startRHR < len(geneGB.origin)-lengthRHR[2]: # Loop as long as the end of RHR is not in an exon and the end of the RHR is inside the max length
            startRHR += 1; # shift RHR end downstream one bp

        failSearchEnd = True; # changes failSearchEnd status
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the first " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between end of gRNA and RHR. \nDefaulted to starting right after end of gRNA most upstream." + "\n"; # give a warning

    # Then search for start region; modify end region if necessary
    while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and (startRHR-gRNADownstream.index[1] <= maxDistanceFromGRNA and geneGB.checkInExon(startRHR)): # if no end region has been found, and still within max distance from gRNA and inside exon,
        # find new start region if necessary
        while meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) < minTmEnds and startRHR-gRNADownstream.index[1] < maxDistanceFromGRNA and not failSearchEnd and geneGB.checkInExon(startRHR): # while no suitable start region found, still within max distance from gRNA, the search for a suitable start region hasn't failed yet, and still within exon,
            startRHR += 1; # shift startRHR downstream

        if gBlockDefault and startRHR-gene.index[0] < minGBlockSize: # if defaulting to a gBlock for any recoded region and the current recoded region is smaller than the minimum gBlock size,
            startRHR = gene.index[0] + minGBlockSize; # extend recoded region to minimum gBlock size

        # search for ends upstream
        while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR <= lengthRHR[2] and endRHR < len(geneGB.origin): # while no suitable end region found and still within max length of RHR and inside gene file,
            endRHR += 1; # shift endRHR downstream

        # if not found downstream
        if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[2]: # if no end region found this way,
            endRHR = startRHR - lengthRHR[1]; # return to preferred position
            # search downstream
            while meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and endRHR-startRHR > lengthRHR[0]: # while no suitable end region found and still within min length of RHR,
                endRHR -= 1; # shift endRHR downstream

        if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds and geneGB.checkInExon(startRHR): # if still not found and still inside exon,
            startRHR += 1; # shifts start of RHR downstream



    if meltingTemp(geneGB.origin[(endRHR-endsLength):endRHR]) < minTmEnds: # if no suitable end region found,
        # resets startRHR
        startRHR = max(gRNADownstream.index[1],gene.index[0]); # saves end index of RHR as whatever is more downstream between the end of the gRNA or the start of the gene (in Python indexing, i.e. not included in RHR).
        if gBlockDefault and startRHR > gene.index[0] and startRHR < gene.index[0] + minGBlockSize: # if defaulting to a gBlock for any recoded region, there is a recoded region, and it is under the minimum gBlock size,
            startRHR = gene.index[0] + minGBlockSize; # extend recoded region to minimum gBlock size

        while meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) < minTmEnds and startRHR-gRNADownstream.index[1] < maxDistanceFromGRNA and not failSearchEnd and geneGB.checkInExon(startRHR) and startRHR < len(geneGB.origin) - lengthRHR[2]: # while no suitable start region found, still within max distance from gRNA, the search for a suitable start region hasn't failed yet, still within exon, and the end of the RHR is inside the max length
            startRHR += 1; # shift startRHR downstream

        endRHR = startRHR + lengthRHR[1]; # stores RHR start index according to preferred length by default
        log = log + "\nWarning: No RHR found for gene " + geneGB.name + " \nwith more than " + str(minTmEnds) + " C melting temperature in the last " + str(endsLength) + " bp of RHR, \nwith a max distance of " + str(maxDistanceFromGRNA) + " bp between start of RHR and gRNA. \nDefaulted to ending so that the RHR size is equal to preferred size of " + str(lengthRHR[1]) + "\n"; # give a warning

    # Now we optimize the resulting RHR by adjusting end sequence around the given range
    searchIndexes = [int(endRHR+optimizeRange[0]), int(min(endRHR+optimizeRange[1],gene.index[1]))]; # list with indexes across which we are going to search for a better end point.
    for i in range(searchIndexes[0], searchIndexes[1]): # iterates across optimization range
        if i >= endsLength: # if inside gene file,
            if meltingTemp(geneGB.origin[(i-endsLength):i]) > meltingTemp(geneGB.origin[(endRHR-endsLength):]) and lengthRHR[0] >= i-startRHR >= lengthRHR[2]: # if this end point has a better Tm and is still within bounds,
                endRHR = i; # make this the starting position



    searchIndexesEnd = [max(int(startRHR+optimizeRange[0]),gene.index[0]), int(startRHR+optimizeRange[1])]; # list with indexes across which we are going to search for a better start point, can't be before start of gene
    for i in range(searchIndexesEnd[0], searchIndexesEnd[1]): # iterates across optimization range
        if meltingTemp(geneGB.origin[i:min(i+endsLength,len(geneGB.origin))]) > meltingTemp(geneGB.origin[startRHR:min(startRHR+endsLength,len(geneGB.origin))]) and lengthRHR[2] >= endRHR-i >= lengthRHR[0] and i > gRNADownstream.index[1] and geneGB.checkInExon(i) and i >= gene.index[0] and (not gBlockDefault or (i >= gene.index[0] + minGBlockSize) or (i == gene.index[0]) ): # if this end point has a better Tm, and is still within bounds, after gRNA and within exon, and if it allows for a gBlock of the minimum size or no gBlock at all if gBlocks are being used as default
            startRHR = i; # make this the starting position




    for site in filterCutSites: # for every cut site being filtered
        if findFirst(site,geneGB.origin[startRHR:endRHR]) > -1 or findFirst(revComp(site),geneGB.origin[startRHR:endRHR]) > -1: # if cut site found,
            log = log + "\nWarning: RHR sequence for gene " + gene.label + ": \n" + geneGB.origin[startRHR:endRHR] + "\ncontains restriction site " + site + "\n"; # add warning to log


    # creates var to store finished RHR as annotation
    RHR = GenBankAnn(gene.label + " RHR", "misc_feature", geneGB.origin[startRHR:endRHR], False, [startRHR,endRHR], annColors['RHRColor']); # creates GenBankAnn object to hold RHR

    log = log + "RHR for gene " + gene.label + " selected.\n\n"; # logs this process finished
    return {"out":RHR, "log":log}; # returns RHR GenBankAnn object
