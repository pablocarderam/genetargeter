
from builtins import str
from builtins import range
from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants
from gRNAScores.gRNAScoring import *; # Import scoring metrics

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
def chooseRecodeRegion3Prime(geneGB, gene, offTargetMethod="cfd", pamType="NGG", orgCodonTable=codonUsage(), avoidPeptidesAtTermini=[[],[]], filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI], codonSampling=False, offScoreThreshold=10, minGCEnd=0.375, gRNATableString=""):
    #TODO: debug #TODO: Recoded if upstream of stop codon add recode values to table
    gRNAs = geneGB.findAnnsLabel("gRNA", True); # List of all gRNAs
    gRNATable = gRNATableString.split('\n'); # split string into lines
    gRNATable = [g.split(',') for g in gRNATable]; # split each line into values

    if offTargetMethod == "hsu": # if off-target scoring with Hsu scores
        offScoreThreshold = 1; # set threshold to 1%

    log = ""; # init log
    LHR = geneGB.findAnnsLabel("LHR")[0]; # LHR annotation object

    annRecoded = GenBankAnn(); # creates GenBankAnn object to hold recoded region
    if LHR.index[1] < gene.index[1]: # if end of LHR is inside gene (or before)
        startRecode = max(LHR.index[1], gene.index[0]); # start of recode region (end of LHR or start of gene, most downstream)
        while not geneGB.checkInExon(startRecode) and startRecode <= len(geneGB.origin): # while recode region start is in intron,
            startRecode += 1; # shift downstream

        intronStartIndices = []; # stores start indexes of introns starting after recode sequence start
        intronEndIndices = []; # stores end indexes of introns starting after recode sequence start
        for ann in geneGB.findAnnsLabel(gene.label): # loop through annotations associated with transcript
            if ann.type == "CDS": # if annotation is cds
                if gene.index[1] > ann.index[1] > startRecode: # if annotation is an exon ending before gene end and after recode start,
                    intronStartIndices.append(ann.index[1]); # add this intron start index
                if gene.index[1] > ann.index[0] > startRecode: # if annotation is an exon starting after recode start,
                    intronEndIndices.append(ann.index[0]); # add this intron end index




        # if len(intronStartIndices) == 0: # if no CDS exons found in this way,
        #     for ann in geneGB.findAnnsLabel(gene.label.split('.')[0]): # loop through annotations associated with gene
        #         if ann.type == "exon": # if annotation is exon
        #             if gene.index[1] > ann.index[1] > startRecode: # if annotation is an exon ending before gene end and after recode start,
        #                 intronStartIndices.append(ann.index[1]); # add this intron start index
        #             if gene.index[1] > ann.index[0] > startRecode: # if annotation is an exon starting after recode start,
        #                 intronEndIndices.append(ann.index[0]); # add this intron end index
        #
        #
        #         elif ann.type == "intron": # if annotation is intron,
        #             if ann.index[0] > startRecode: # if annotation is an intron starting after recode start,
        #                 intronStartIndices.append(ann.index[0]); # add this intron start index
        #                 intronEndIndices.append(ann.index[1]); # add this intron start index



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
        cutCheck = True; # variable used to check if no cut sequences are present.
        offScore = 100; # stores off-target score. Default is 100% due to the fact that gRNA sequence is the same.
        count = 0; # iteration counter
        recodedSeq = recodeSeq; # assign recoded sequence to same as original
        bestRecodedSeq = recodedSeq; # will store best candidate sequence
        if len(recodeSeq) > 2: # if recodeSeq contains at least one codon,
            tricky = False; # True if suspected to be hard to synthesize
            badStart = False; # True if first bases have low melting temp (important for Gibson assembly)
            candidateFound = False; # signal possible candidate found
            bestRecodedSeq = recodedSeq; # will store best candidate sequence
            moCloProblem = translate(nonRecodedStart+recodedSeq[0:(6-frame)]) in avoidPeptidesAtTermini[0] or translate(recodedSeq[-6:]) in avoidPeptidesAtTermini[1]; # True if terminal dipeptides are in corresponding lists of incompatible dipeptides
            while not cutCheck or offScore > offScoreThreshold or tricky or badStart or moCloProblem: # while cutCheck shows hits in a cut sequences, or while the pairwise off-target score is over the threshold, or while there are difficult-to-synthesize structures in the recoded region, or while the first 40 bp have a bad gc content, or while the terminal dipeptides are incompatible with MoClo,
                if count == 1: # if recoded region has failed checks once,
                    log = log + "Defaulted recoded region recodonization to codon sampling due to possible difficulties in synthesis or enzyme cut sites.\n\n"; # log warning
                    codonSampling = True; # forces codonSampling to true
                    if count == 10:
                        log = log + "Defaulted recoded region recodonization to random codon sampling due to possible difficulties in synthesis or enzyme cut sites.\n\n"; # log warning
                        orgCodonTable = codonUsage(); # forces random codon selection to true

                cutCheck = True; # reset cutCheck
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
                    cutCheck = cutCheck * ( findFirst(recodedSeq,site) < 0 ); # Find cut site, register in cutCheck
                    cutCheck = cutCheck * ( findFirst(recodedSeq,revComp(site)) < 0 ); # Find cut site in comp strand, register in cutCheck

                tricky = isTricky(recodedSeq); # check if tricky to synthesize
                trickyCount = 0
                while tricky and tricky < len(recodedSeq)-9 and trickyCount < 10000: # targeted recoding of problematic fragments
                    recodedSeq = recodedSeq[0:tricky-tricky%3] + optimizeCodons(recodedSeq[tricky-tricky%3:tricky-tricky%3+9]) + recodedSeq[tricky-tricky%3+9:]; # optimize codons.
                    tricky = isTricky(recodedSeq); # check if tricky to synthesize
                    trickyCount += 1

                if gcContent(recodedSeq[0:40]) < minGCEnd: # if the first bases don't have enough gc content
                    badStart = True;

                moCloProblem = translate(nonRecodedStart+recodedSeq[0:(6-frame)]) in avoidPeptidesAtTermini[0] or translate(recodedSeq[-6:]) in avoidPeptidesAtTermini[1]; # True if terminal dipeptides are in corresponding lists of incompatible dipeptides

                if offScore <= offScoreThreshold and cutCheck and not moCloProblem: # if parameters other than badStart are ok and this sequence has better start than previous best,
                    if not candidateFound or isTricky(bestRecodedSeq): # if no candidate found until now or current best is already tricky,
                        bestRecodedSeq = recodedSeq; # make this new best
                    elif not tricky and gcContent(recodedSeq[0:40]) > gcContent(bestRecodedSeq[0:40]):
                        bestRecodedSeq = recodedSeq; # make this new best

                    if not tricky:
                        candidateFound = True; # signal possible candidate found

                count += 1; # advances iteration counter
                if count > 5000: # if out of iteration limit,
                    if not candidateFound: # if no candidate without cut sequences found,
                        if tricky:
                            log = log + "Warning: Recoded region for gene " + gene.label + " could not reshuffle enough to avoid repeated sequences or low-complexity regions.\n\n"; # log warning
                        else:
                            log = log + "Warning: Recoded region for gene " + gene.label + " could not reshuffle enough to fulfill the maximum off-target sgRNA score threshold, or avoid all the following cut sequences: \n" + str(cutSeqs) + "\n\n"; # log warning

                    break; # escape loop


                #print [gOnSeq+"NGG",gOffSeq+gNewPAM,pairScoreCFD(gOnSeq,gOffSeq,gNewPAM,pamType),pairScoreHsu(gOnSeq,gOffSeq,gNewPAM,pamType)]

        recodedSeq = nonRecodedStart + bestRecodedSeq; # adds initial bases from reading frame adjustment to best candidate
        annRecoded = GenBankAnn(gene.label + " Recoded", "misc_feature", recodedSeq, False, [startRecode,endRecode], annColors['recodedRegionColor']); # creates var to store finished recodedSeq as annotation
        log = log + "Recoded region with size " + str(len(recodedSeq)) + " for gene " + gene.label + " selected.\n\n"; # logs this process finished

    else: # if no recoded region necessary,
        log = log + "Recoded region not deemed necessary for gene " + gene.label + ".\n\n"; # logs this process finished

    gRNATableString = "\n".join([",".join(g) for g in gRNATable]); # Creates string from grna array
    gRNATableString = gRNATableString.replace(">=threshold",">="+str(offScoreThreshold)); # adds pairwise recoded threshold values
    return {"out":annRecoded, "log":log, "gRNATable":gRNATableString}; # returns recoded region GenBankAnn object



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
def chooseRecodeRegion5Prime(geneGB, gene, offTargetMethod="cfd", pamType="NGG", orgCodonTable=codonUsage(), avoidPeptidesAtTermini=[[],[]], filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI], codonSampling=False, offScoreThreshold=10, minGCEnd=0.375, gRNATableString="", haTag=True):
    #TODO: debug #TODO: Recoded if upstream of stop codon add recode values to table
    gRNAs = geneGB.findAnnsLabel("gRNA", True); # List of all gRNAs
    gRNATable = gRNATableString.split('\n'); # split string into lines
    gRNATable = [g.split(',') for g in gRNATable]; # split each line into values

    if offTargetMethod == "hsu": # if off-target scoring with Hsu scores
        offScoreThreshold = 1; # set threshold to 1%

    log = ""; # init log
    RHR = geneGB.findAnnsLabel("RHR")[0]; # RHR annotation object

    annRecoded = GenBankAnn(); # creates GenBankAnn object to hold recoded region
    if RHR.index[0] > gene.index[0]: # if end of RHR is inside gene
        endRecode = min(RHR.index[0],gene.index[1]); # end of recode region (start of RHR or end of gene, most upstream)
        while not geneGB.checkInExon(endRecode): # while recode region end is in intron,
            endRecode -= 1; # shift upstream

        intronStartIndices = []; # stores start indexes of introns starting after recode sequence start
        intronEndIndices = []; # stores end indexes of introns starting after recode sequence start
        for ann in geneGB.findAnnsLabel(gene.label): # loop through annotations associated with transcript
            if ann.type == "CDS": # if annotation is cds
                if gene.index[0] < ann.index[1] < endRecode: # if annotation is an exon ending after gene start and before recode end,
                    intronStartIndices.append(ann.index[1]); # add this intron start index
                if gene.index[0] < ann.index[0] < endRecode: # if annotation is an exon starting before recode end,
                    intronEndIndices.append(ann.index[0]); # add this intron end index




        # if len(intronStartIndices) == 0: # if no CDS exons found in this way,
        #     for ann in geneGB.findAnnsLabel(gene.label.split('.')[0]): # loop through annotations associated with gene
        #         if ann.type == "exon": # if annotation is exon
        #             if gene.index[0] < ann.index[1] < endRecode: # if annotation is an exon ending after gene start and before recode end,
        #                 intronStartIndices.append(ann.index[1]); # add this intron start index
        #             if gene.index[0] < ann.index[0] < endRecode: # if annotation is an exon starting before recode end,
        #                 intronEndIndices.append(ann.index[0]); # add this intron end index
        #
        #
        #         elif ann.type == "intron": # if annotation is intron,
        #             if ann.index[1] < endRecode: # if annotation is an intron ending before recode end,
        #                 intronStartIndices.append(ann.index[0]); # add this intron start index
        #                 intronEndIndices.append(ann.index[1]); # add this intron start index



        intronIndices = []; # will contain final indexes of introns downstream of recode start (introns to be removed from recoded region)
        intronStartIndices = sorted(intronStartIndices); # sort
        intronEndIndices = sorted(intronEndIndices); # sort
        for i in range(len(intronEndIndices)): # for every intron end,
            if intronEndIndices[i] <= RHR.index[0]: # if before RHR start,
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



        startRecode = gene.index[0]; # end of recode region (end of gene, exclude stop codon)
        recodeSeq = geneGB.origin[startRecode:endRecode]; # will contain sequence to be recorded
        if len(intronIndices) > 0: # if there are introns,
            recodeSeq = geneGB.origin[startRecode:intronIndices[0][0]]; # get recode sequence until first intron
            for i in range(len(intronIndices)-1): # for every intron except last one,
                recodeSeq = recodeSeq + geneGB.origin[intronIndices[i][1]:intronIndices[i+1][0]]; # add next exon to recode seq

            recodeSeq = recodeSeq + geneGB.origin[intronIndices[len(intronIndices)-1][1]:endRecode]; # get rest of recode sequence until endRecode

        frame = len(recodeSeq) % 3; # stores reading frame, index from start of sequence to be recoded
        endRecode -= frame; # modify recode end site according to reading frame
        nonRecodedEnd = ""; # stores 0, 1 or 2 nucleotides not recoded due to reading frame
        if frame != 0: # if frame shift is not zero, to avoid listing all of recode region with recodeSeq[-0:],
            nonRecodedEnd = recodeSeq[-frame:]; # stores 0, 1 or 2 nucleotides not recoded due to reading frame

        recodeSeq = recodeSeq[0:len(recodeSeq)-frame]; # adjust recode region
        if haTag: # if adding an HA tag,
            recodeSeq = ha_tag + recodeSeq; # add HA tag to start of recoded region

        cutSeqs = filterCutSites + [g.seq for g in gRNAs]; # list of all cut seqs. all gRNAs in gene are to be included as cut sequences
        cutCheck = True; # variable used to check if cut sequences are present. Initially false since all gRNAs are present.
        offScore = 100; # stores off-target score. Default is 100% due to the fact that gRNA sequence is the same.
        count = 0; # iteration counter
        recodedSeq = recodeSeq; # assign recoded sequence to same as original
        bestRecodedSeq = recodedSeq; # will store best candidate sequence
        if len(recodeSeq) > 2: # if recodeSeq contains at least one codon,
            tricky = False; # True if suspected to be hard to synthesize
            badStart = False; # True if first bases have low melting temp (important for Gibson assembly)
            candidateFound = False; # signal possible candidate found
            bestRecodedSeq = recodedSeq; # will store best candidate sequence
            moCloProblem = translate(nonRecodedStart+recodedSeq[0:(6-frame)]) in avoidPeptidesAtTermini[0] or translate(recodedSeq[-6:]) in avoidPeptidesAtTermini[1]; # True if terminal dipeptides are in corresponding lists of incompatible dipeptides
            while not cutCheck or offScore > offScoreThreshold or tricky or badStart or moCloProblem: # while cutCheck is greater than what you would expect for no hits in all cut sequences plus the gRNAs on both positive and comp strands, or while the pairwise off-target score is over the threshold, or while there are difficult-to-synthesize structures in the recoded region, or while the first 40 bp have a bad gc content
                if count == 1: # if recoded region has failed checks once,
                    log = log + "Defaulted recoded region recodonization to codon sampling due to possible difficulties in synthesis or enzyme cut sites.\n\n"; # log warning
                    codonSampling = True; # forces codonSampling to true
                    if count == 10:
                        log = log + "Defaulted recoded region recodonization to random codon sampling due to possible difficulties in synthesis or enzyme cut sites.\n\n"; # log warning
                        orgCodonTable = codonUsage(); # forces random codon selection to true

                cutCheck = True; # reset cutCheck
                offScore = 0; # reset offScore
                tricky = False; # reset tricky Boolean
                badStart = False; # reset badStart Boolean
                recodedSeq = optimizeCodons(recodeSeq,orgCodonTable,codonSampling=codonSampling); # optimize codons.
                for g in gRNAs: # for every gRNA candidate within recoded region,
                    if g.index[0] >= startRecode and g.index[1] <= endRecode-frame: # if grna is inside recoded region
                        gOnSeq = g.seq; # get original gRNA sequence
                        wholeRecSeq = recodedSeq + nonRecodedEnd; # add initial bases
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
                    cutCheck = cutCheck * ( findFirst(recodedSeq,site) < 0 ); # Find cut site, register in cutCheck
                    cutCheck = cutCheck * ( findFirst(recodedSeq,revComp(site)) < 0 ); # Find cut site in comp strand, register in cutCheck

                tricky = isTricky(recodedSeq); # check if tricky to synthesize
                trickyCount = 0
                while tricky and tricky < len(recodedSeq)-9 and trickyCount < 10000: # targeted recoding of problematic fragments
                    recodedSeq = recodedSeq[0:tricky-tricky%3] + optimizeCodons(recodedSeq[tricky-tricky%3:tricky-tricky%3+9]) + recodedSeq[tricky-tricky%3+9:]; # optimize codons.
                    tricky = isTricky(recodedSeq); # check if tricky to synthesize
                    trickyCount += 1

                if gcContent(recodedSeq[-40:]) < minGCEnd: # if the last bases don't have enough gc content
                    badStart = True;

                moCloProblem = translate(nonRecodedStart+recodedSeq[0:(6-frame)]) in avoidPeptidesAtTermini[0] or translate(recodedSeq[-6:]) in avoidPeptidesAtTermini[1]; # True if terminal dipeptides are in corresponding lists of incompatible dipeptides

                if offScore <= offScoreThreshold and cutCheck and not moCloProblem: # if parameters other than badStart are ok and this sequence has better start than previous best,
                    if not candidateFound or isTricky(bestRecodedSeq): # if no candidate found until now or current best is already tricky,
                        bestRecodedSeq = recodedSeq; # make this new best
                    elif not tricky and gcContent(recodedSeq[-40:]) > gcContent(bestRecodedSeq[-40:]):
                        bestRecodedSeq = recodedSeq; # make this new best

                    if not tricky:
                        candidateFound = True; # signal possible candidate found

                count += 1; # advances iteration counter
                if count > 5000: # if out of iteration limit,
                    if not candidateFound: # if no candidate without cut sequences found,
                        if tricky:
                            log = log + "Warning: Recoded region for gene " + gene.label + " could not reshuffle enough to avoid repeated sequences or low-complexity regions.\n\n"; # log warning
                        else:
                            log = log + "Warning: Recoded region for gene " + gene.label + " could not reshuffle enough to fulfill the maximum off-target sgRNA score threshold, or avoid all the following cut sequences: \n" + str(cutSeqs) + "\n\n"; # log warning

                    break; # escape loop


                #print [gOnSeq+"NGG",gOffSeq+gNewPAM,pairScoreCFD(gOnSeq,gOffSeq,gNewPAM,pamType),pairScoreHsu(gOnSeq,gOffSeq,gNewPAM,pamType)]

        recodedSeq = bestRecodedSeq + nonRecodedEnd; # adds end bases from reading frame adjustment to best candidate
        annRecoded = GenBankAnn(gene.label + " Recoded", "misc_feature", recodedSeq, False, [startRecode,endRecode], annColors['recodedRegionColor']); # creates var to store finished recodedSeq as annotation
        haTagMsg = ""; # used to output message
        if haTag: # if using an HA tag,
            haTagMsg = " with a recoded HA tag"; # msg modifier

        log = log + "Recoded region with size " + str(len(recodedSeq)) + " for gene " + gene.label + haTagMsg + " selected.\n\n"; # logs this process finished

    else: # if no recoded region necessary,
        log = log + "Recoded region not deemed necessary for gene " + gene.label + ".\n\n"; # logs this process finished

    gRNATableString = "\n".join([",".join(g) for g in gRNATable]); # Creates string from grna array
    gRNATableString = gRNATableString.replace(">=threshold",">="+str(offScoreThreshold)); # adds pairwise recoded threshold values

    return {"out":annRecoded, "log":log, "gRNATable":gRNATableString}; # returns recoded region GenBankAnn object


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
#TODO: different cut sites for different plasmids
def chooseRecodeRegion(geneGB, gene, offTargetMethod="cfd", pamType="NGG", orgCodonTable=codonUsage(), avoidPeptidesAtTermini=[[],[]], filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI], codonSampling=False, offScoreThreshold=10, minGCEnd=0.375, gRNATableString="", target3Prime=True, haTag=False):
    out = {}; # will contain method output
    if target3Prime: # if targeting 3'
        out = chooseRecodeRegion3Prime(geneGB, gene, offTargetMethod, pamType=pamType, orgCodonTable=orgCodonTable,codonSampling=codonSampling, gRNATableString=gRNATableString); # defines region to be recoded, returns recoded sequence
    else: # if using pSN150,
        out = chooseRecodeRegion5Prime(geneGB, gene, offTargetMethod, pamType=pamType, orgCodonTable=orgCodonTable,codonSampling=codonSampling, gRNATableString=gRNATableString, haTag=haTag); # defines region to be recoded, returns recoded sequence

    return out;
